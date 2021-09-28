package Bio::Adventure::Map;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use File::Basename;
use File::Spec;
use File::Which qw"which";

=head1 NAME

Bio::Adventure::Map - Perform highthroughput sequence alignments with tools like bowtie/tophat/etc

=head1 SYNOPSIS

use Bio::Adventure;
my $hpgl = new Bio::Adventure;
$hpgl->Bowtie();

=head1 METHODS

=head2 C<Bowtie>

Perform a bowtie alignment.  Unless instructed otherwise, it will do so with 0
mismatches and with multi-matches randomly placed 1 time among the
possibilities. (options -v 0 -M 1)

It checks to see if a bowtie1 compatible index is in
$libdir/$libtype/indexes/$species, if not it attempts to create
them.

It will continue on to convert the bowtie sam output to a compressed, sorted,
indexed bam file, and pass that to htseq-count using a gff file of the same
species.

=cut
sub Bowtie {
    my ($class, %args) = @_;
    my $check = which('bowtie-build');
    die("Could not find bowtie in your PATH.") unless($check);
    print "Recall that you can change the bowtie arguments via 'bt_type'.\n";
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        bt_type => 'v0M1',
        count => 1,
        libtype => 'genome',
        htseq_type => 'gene',
        htseq_id => 'ID',
        jprefix => '10',
        modules => ['bowtie1', 'samtools', 'htseq'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $start_species = $options->{species};
    my $species = $start_species;
    if ($species =~ /\:/) {
        my @species_lst = split(/:/, $species);
        my @result_lst = ();
        foreach my $sp (@species_lst) {
            print "Invoking bowtie on ${sp}\n";
            $class->{species} = $sp;
            my $result = Bio::Adventure::Map::Bowtie($class, %{$options});
            push (@result_lst, $result);
        }
        $options->{species} = $start_species;
        return(@result_lst);
    }

    my $bt_type = $options->{bt_type};
    my $bt_args = $options->{bt_args}->{$bt_type};
    $bt_args = ' --best -v 0 -M 1 ' if (!defined($bt_args));

    my $sleep_time = 3;
    my $bt_input = $options->{input};

    my $test_file = "";
    if ($bt_input =~ /\:|\;|\,|\s+/) {
        my @pair_listing = split(/\:|\;|\,|\s+/, $bt_input);
        $pair_listing[0] = File::Spec->rel2abs($pair_listing[0]);
        $pair_listing[1] = File::Spec->rel2abs($pair_listing[1]);
        $bt_input = qq" <(less $pair_listing[0]) <(less $pair_listing[1]) ";
        $test_file = $pair_listing[0];
    } else {
        $test_file = File::Spec->rel2abs($bt_input);
        $bt_input = qq" <(less ${test_file}) ";
    }

    my $jname = qq"bt${bt_type}_${species}";
    ##$jname = $options->{jname} if ($options->{jname});
    my $libtype = $options->{libtype};
    my $count = $options->{count};
    my $bt_dir = qq"outputs/bowtie_${species}";
    $bt_dir = $options->{bt_dir} if ($options->{bt_dir});

    ## Check that the indexes exist
    my $bt_reflib = "$options->{libdir}/${libtype}/indexes/${species}";
    my $bt_reftest = qq"${bt_reflib}.1.ebwt";
    my $current_prefix = 10;
    if (defined($options->{jprefix})) {
        $current_prefix = $options->{jprefix};
    }
    my $index_job;
    if (!-r $bt_reftest && !$options->{bt1_indexjobs}) {
        $options = $class->Set_Vars(bt1_indexjobs => 1);
        $index_job = $class->Bio::Adventure::Map::BT1_Index(
            jdepends => $options->{jdepends},
            jprefix => $current_prefix,
            libtype => $libtype,);
        $options->{jdepends} = $index_job->{jobid};
    }

    my $bowtie_input_flag = "-q"; ## fastq by default
    $bowtie_input_flag = "-f" if ($options->{input} =~ /\.fasta/);

    my $cpus = $options->{cpus};
    my $error_file = qq"${bt_dir}/$options->{jbasename}-${bt_type}.err";
    my $comment = qq!## This is a bowtie1 alignment of ${bt_input} against
## ${bt_reflib} using arguments: ${bt_args}.
!;
    my $aligned_filename = qq"${bt_dir}/$options->{jbasename}-${bt_type}_aligned_${species}.fastq";
    my $unaligned_filename = qq"${bt_dir}/$options->{jbasename}-${bt_type}_unaligned_${species}.fastq";
    my $sam_filename = qq"${bt_dir}/$options->{jbasename}-${bt_type}.sam";
    my $jstring = qq!mkdir -p ${bt_dir} && sleep ${sleep_time} && bowtie \\
  ${bt_reflib} \\
  ${bt_args} \\
  -p ${cpus} \\
  ${bowtie_input_flag} ${bt_input} \\
  --un ${unaligned_filename} \\
  --al ${aligned_filename} \\
  -S ${sam_filename} \\
  2>${error_file} \\
  1>${bt_dir}/$options->{jbasename}-${bt_type}.out
!;

    my $bt_job = $class->Submit(
        aligned => $aligned_filename,
        comment => $comment,
        input => $bt_input,
        jdepends => $options->{jdepends},
        jname => $jname,
        jprefix => $options->{jprefix},
        jstring => $jstring,
        output => $sam_filename,
        modules => $options->{modules},
        postscript => $options->{postscript},
        prescript => $options->{prescript},
        jqueue => 'workstation',
        unaligned => $unaligned_filename,);
    if (defined($index_job)) {
        $bt_job->{index} = $index_job;
    }

    my $un_comp = $class->Bio::Adventure::Compress::Recompress(
        comment => qq"## Compressing the sequences which failed to align against ${bt_reflib} using options ${bt_args}.\n",
        jdepends => $bt_job->{job_id},
        jname => "xzun",
        jprefix => $options->{jprefix} + 1,
        xz_input => "${bt_dir}/$options->{jbasename}-${bt_type}_unaligned_${species}.fastq",);
    $bt_job->{unaligned_compression} = $un_comp;

    my $al_comp = $class->Bio::Adventure::Compress::Recompress(
        comment => qq"## Compressing the sequences which successfully aligned against ${bt_reflib} using options ${bt_args}.",
        jdepends => $bt_job->{job_id},
        jname => "xzal",
        jprefix => $options->{jprefix} + 2,
        xz_input => "${bt_dir}/$options->{jbasename}-${bt_type}_aligned_${species}.fastq",);
    $bt_job->{aligned_compression} = $al_comp;

    ## BT1_Stats also reads the trimomatic output, which perhaps it should not.
    my $trim_output_file = qq"outputs/trimomatic_stats.csv";

    my $sam_job = $class->Bio::Adventure::Convert::Samtools(
        input => $sam_filename,
        jdepends => $bt_job->{job_id},
        jname => "s2b_${jname}",
        jprefix => $options->{jprefix} + 3,);
    $bt_job->{samtools} = $sam_job;
    my $htmulti;
    if ($count) {
        if ($libtype eq 'rRNA') {
            $htmulti = $class->Bio::Adventure::Count::HTSeq(
                htseq_id => $options->{htseq_id},
                htseq_input => $sam_job->{output},
                htseq_type => $options->{htseq_type},
                jdepends => $sam_job->{job_id},
                jname => "ht_${jname}",
                jprefix => $options->{jprefix} + 4,
                libtype => $libtype,
                jqueue => 'workstation',
                suffix => $bt_type,
                mapper => 'bowtie1',);
            $bt_job->{rRNA_count} = $htmulti;
        } else {
            $htmulti = $class->Bio::Adventure::Count::HT_Multi(
                htseq_id => $options->{htseq_id},
                htseq_input => $sam_job->{output},
                htseq_type => $options->{htseq_type},
                jdepends => $sam_job->{job_id},
                jname => "ht_${jname}",
                jprefix => $options->{jprefix} + 4,
                libtype => $libtype,
                jqueue => 'workstation',
                suffix => $bt_type,
                mapper => 'bowtie1',);
            $bt_job->{htseq} = $htmulti;
        }
    }  ## End if ($count)

    my $stats = $class->Bio::Adventure::Map::BT1_Stats(
        %args,
        bt_input => $error_file,
        bt_type => $bt_type,
        count_table => qq"$options->{jbasename}-${bt_type}.count.xz",
        jdepends => $bt_job->{job_id},
        jname => "${jname}_stats",
        jprefix => $options->{jprefix} + 5,
        trim_input => ${trim_output_file},);
    $bt_job->{stats} = $stats;

    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    
    return($bt_job);
}

=head2 C<Bowtie2>

Perform a bowtie2 alignment.  Unless instructed otherwise, it will do so with 0
mismatches and with multi-matches randomly placed 1 time among the
possibilities. (options -v 0 -M 1)

It checks to see if a bowtie2 compatible index is in
$libdir/$libtype/indexes/$species, if not it attempts to create
them.

It will continue on to convert the bowtie sam output to a
compressed, sorted, indexed bam file, and pass that to htseq-count
using a gff file of the same species.

=cut
sub Bowtie2 {
    my ($class, %args) = @_;
    my $check = which('bowtie2-build');
    die("Could not find bowtie2 in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input',],
        do_htseq => 1,
        htseq_type => 'gene',
        htseq_id => 'ID',
        jmem => 28,
        jprefix => 20,
        modules => ['bowtie']);

    if ($options->{species} =~ /\:/) {
        my @species_lst = split(/:/, $options->{species});
        my @result_lst = ();
        my $start_species = $options->{species};
        foreach my $sp (@species_lst) {
            print "Invoking bowtie2 on ${sp}\n";
            $class->{species} = $sp;
            my $result = Bio::Adventure::Map::Bowtie2($class %{$options});
            push (@result_lst, $result);
        }
        $options->{species} = $start_species;
        return(@result_lst);
    }

    my $ready = $class->Check_Input(
        files => $options->{input},
    );
    my $sleep_time = 3;
    my %bt_jobs = ();
    my $libtype = 'genome';
    my $bt2_args = $options->{bt2_args};

    my $prefix_name = qq"bt2";
    my $bt2_name = qq"${prefix_name}_$options->{species}";
    my $suffix_name = $prefix_name;
    if ($options->{jname}) {
        $bt2_name .= qq"_$options->{jname}";
        $suffix_name .= qq"_$options->{jname}";
    }

    my $bt_dir = qq"outputs/bowtie2_$options->{species}";
    if ($args{bt_dir}) {
        $bt_dir = $args{bt_dir};
    }
    my $bt_input = $options->{input};

    my $test_file = "";
    if ($bt_input =~ /\:|\;|\,|\s+/) {
        my @pair_listing = split(/\:|\;|\,|\s+/, $bt_input);
        $pair_listing[0] = File::Spec->rel2abs($pair_listing[0]);
        $pair_listing[1] = File::Spec->rel2abs($pair_listing[1]);
        $bt_input = qq" -1 <(less $pair_listing[0]) -2 <(less $pair_listing[1]) ";
        $test_file = $pair_listing[0];
    } else {
        $test_file = File::Spec->rel2abs($bt_input);
        $bt_input = qq" <(less ${test_file}) ";
    }

    ## Check that the indexes exist
    my $bt_reflib = "$options->{libdir}/$options->{libtype}/indexes/$options->{species}";
    my $bt_reftest = qq"${bt_reflib}.1.bt2";
    my $bt_reftest_large = qq"${bt_reflib}.1.bt2l";
    my $index_job;
    if (!-r $bt_reftest && !-r $bt_reftest_large) {
        print "Hey! The Indexes do not appear to exist, check this out: ${bt_reftest}\n";
        sleep(20);
        $index_job = Bio::Adventure::Map::BT2_Index(
            $class,
            jdepends => $options->{jdepends},
            jprefix => $options->{jprefix} - 1,
            libtype => $libtype,);
        $options->{jdepends} = $index_job->{jobid};
    }
    my $bowtie_input_flag = "-q "; ## fastq by default
    $bowtie_input_flag = "-f " if (${bt_input} =~ /\.fasta$/);

    my $cpus = $options->{cpus};
    my $error_file = qq"${bt_dir}/$options->{jbasename}.err";
    my $comment = qq!## This is a bowtie2 alignment of ${bt_input} against
## ${bt_reflib} using arguments: ${bt2_args}.
## This jobs depended on: $options->{jdepends}
!;
    my $aligned_filename = qq"${bt_dir}/$options->{jbasename}_aligned_$options->{species}.fastq";
    my $unaligned_filename = qq"${bt_dir}/$options->{jbasename}_unaligned_$options->{species}.fastq";
    my $sam_filename = qq"${bt_dir}/$options->{jbasename}.sam";
    my $jstring = qq!mkdir -p ${bt_dir} && \\
  sleep ${sleep_time} && \\
  bowtie2 -x ${bt_reflib} ${bt2_args} \\
    -p ${cpus} \\
    ${bowtie_input_flag} ${bt_input} \\
    --un ${unaligned_filename} \\
    --al ${aligned_filename} \\
    -S ${sam_filename} \\
    2>${error_file} \\
    1>${bt_dir}/$options->{jbasename}.out
!;

    my $bt2_job = $class->Submit(
        aligned => $aligned_filename,
        comment => $comment,
        input => $bt_input,
        jname => $bt2_name,
        jdepends => $options->{jdepends},
        jstring => $jstring,
        jprefix => $options->{jprefix},
        jmem => $options->{jmem},
        modules => $options->{modules},
        output => $sam_filename,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        unaligned => $unaligned_filename,);

    my $un_comp = $class->Bio::Adventure::Compress::Recompress(
        comment => qq"## Compressing the sequences which failed to align against ${bt_reflib} using options ${bt2_args}\n",
        xz_input => "${bt_dir}/$options->{jbasename}_unaligned_$options->{species}.fastq",
        jdepends => $bt2_job->{job_id},
        jname => "xzun_${suffix_name}",
        jprefix => $options->{jprefix} + 1,);
    $bt2_job->{unaligned_compression} = $un_comp;

    my $al_comp = $class->Bio::Adventure::Compress::Recompress(
        comment => qq"## Compressing the sequences which successfully aligned against ${bt_reflib} using options ${bt2_args}",
        xz_input => "${bt_dir}/$options->{jbasename}_aligned_$options->{species}.fastq",
        jname => "xzal_${suffix_name}",
        jprefix => $options->{jprefix} + 2,
        depends => $bt2_job->{job_id},);
    $bt2_job->{aligned_compression} = $al_comp;

    ## BT1_Stats also reads the trimomatic output, which perhaps it should not.
    ## my $trim_output_file = qq"outputs/$options->{jbasename}-trimomatic.out";
    my $stats = $class->Bio::Adventure::Map::BT2_Stats(
        bt_input => $error_file,
        count_table => qq"$options->{jbasename}.count.xz",
        jdepends => $bt2_job->{job_id},
        jname => "bt2st_${suffix_name}",
        jprefix => $options->{jprefix} + 3,
        ## trim_input => ${trim_output_file},
        );
    $bt2_job->{stats} = $stats;
    my $sam_job = $class->Bio::Adventure::Convert::Samtools(
        input => $sam_filename,
        jdepends => $bt2_job->{job_id},
        jname => "s2b_${suffix_name}",
        jprefix => $options->{jprefix} + 4,);
    $bt2_job->{samtools} = $sam_job;
    my $htseq_input = $sam_job->{output};
    my $htmulti;
    if ($options->{do_htseq}) {
        if ($libtype eq 'rRNA') {
            $htmulti = $class->Bio::Adventure::Count::HTSeq(
                htseq_input => $sam_job->{output},
                htseq_type => $options->{htseq_type},
                htseq_id => $options->{htseq_id},
                jdepends => $sam_job->{job_id},
                jname => $suffix_name,
                jprefix => $options->{jprefix} + 5,
                libtype => $libtype,
                mapper => 'bowtie2',);
        } else {
            $htmulti = $class->Bio::Adventure::Count::HT_Multi(
                htseq_input => $sam_job->{output},
                htseq_type => $options->{htseq_type},
                htseq_id => $options->{htseq_id},
                jdepends => $sam_job->{job_id},
                jname => $suffix_name,
                jprefix => $options->{jprefix} + 6,
                libtype => $libtype,
                mapper => 'bowtie2',);
            $bt2_job->{htseq} = $htmulti;
        }
    }
    return($bt2_job);
}

=head2 C<Bowtie_RRNA>

Perform an alignment against a home-curated set of ribosomal RNA/tRNA sequences.
The alignment requires a fastq input library and fasta library found in
'libraries/rRNA/$class->{species}.fasta'

Example:
  my $rrna = $hpgl->Bowtie_RRNA();
  ## If you want to exclude the rRNA sequences from future alignments:
  my $rrna = $hpgl->Bowtie_RRNA(exclude => 1);

=cut
sub Bowtie_RRNA {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ["species"],);
    my $job_name = qq"rRNA_$options->{jbasename}";
    my $exclude = 0;
    $exclude = $options->{exclude} if ($options->{exclude});
    my $species = $options->{species};
    my ${bt_dir} = qq"outputs/bowtie_$options->{species}";
    my $job = $class->Bio::Adventure::Map::Bowtie(
        jdepends => $options->{jdepends},
        jname => $job_name,
        libtype => 'rRNA',
        prescript => $args{prescript},
        postscript => $args{postscript},);
    ## Return the basename back to normal so that future tasks don't
    ## get confuseled.
    return($job);
}

=head2 C<BT_Multi>

Attempts to run multiple bowtie1 runs for a given species.  One run is performed
for each of a few parameters which are kept in the variable '$hpgl->{bt_types}'
and generally include: 0 mismatch, 1 mismatch, 2 mismatches, 1 randomly placed
hit, 0 randomly placed hits, or default options.

=cut
sub BT_Multi {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ["species", "input", "htseq_type"],);
    my $bt_input = $options->{input};
    my $species = $options->{species};
    my %bt_types = %{$options->{bt_args}};
    my $count = 0;
    my $job;
    foreach my $type (keys %bt_types) {
        my $jname = qq"bt${type}_${species}";
        my $bt_job = $class->Bio::Adventure::Map::Bowtie(
            input => $bt_input,
            bt_type => $type,
            jdepends => $options->{jdepends},
            jname => $jname,
            prescript => $args{prescript},
            postscript => $args{postscript},);
        if ($count == 0) {
            $job = $bt_job;
        } else {
            $job->{$count} = $bt_job;
        }
        $count++;
    }
    return($job);
}

=head2 C<BT1_Index>

Create a bowtie1 index using $hpgl->{species}.fasta and leaves it in the
indexes/ directory.

=cut
sub BT1_Index {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ["species"],);

    my $jstring = qq!bowtie-build $options->{libdir}/$options->{libtype}/$options->{species}.fasta \\
  $options->{libdir}/$options->{libtype}/indexes/$options->{species}
!;
    my $comment = qq!## Generating bowtie1 indexes for species: $options->{species} in $options->{libdir}/$options->{libtype}/indexes!;
    my $bt1_index = $class->Submit(
        comment => $comment,
        jname => "bt1idx_$options->{species}",
        jdepends => $options->{jdepends},
        jstring => $jstring,
        jprefix => "10",
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($bt1_index);
}

=head2 C<BT2_Index>

Create a bowtie2 index using ${species}.fasta and leaves it in the indexes/
directory.

=cut
sub BT2_Index {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ["species"],
        modules => ['bowtie'],);
    my $libtype = $options->{libtype};
    my $libdir = File::Spec->rel2abs($options->{libdir});
    my $jstring = qq!
if test \! -e "${libdir}/genome/$options->{species}.fa"; then
  ln -sf ${libdir}/genome/$options->{species}.fasta ${libdir}/genome/$options->{species}.fa
fi
if test \! -e "${libdir}/genome/indexes/$options->{species}.fa"; then
  ln -sf ${libdir}/genome/$options->{species}.fasta ${libdir}/genome/indexes/$options->{species}.fa
fi

bowtie2-build $options->{libdir}/genome/$options->{species}.fasta \\
  $options->{libdir}/${libtype}/indexes/$options->{species}
!;
    my $comment = qq!## Generating bowtie2 indexes for species: $options->{species} in $options->{libdir}/${libtype}/indexes!;
    my $indexer = $class->Submit(
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => "bt2idx_$options->{species}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        modules => $options->{modules},
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($indexer);
}

=head2 C<BT1_Stats>

Collect some alignment statistics from bowtie1.

=cut
sub BT1_Stats {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args);
    my $bt_input = $options->{bt_input};
    my $bt_type = "";
    $bt_type = $options->{bt_type} if ($options->{bt_type});
    my $jname = "stats";
    $jname = $options->{jname} if ($options->{jname});
    my $jobid = qq"$options->{jbasename}_stats";
    my $count_table = "";
    $count_table = $options->{count_table} if ($options->{count_table});
    my $stat_output = qq"outputs/bowtie_stats.csv";
    my $comment = qq!## This is a stupidly simple job to collect alignment statistics.!;
    my $jstring = qq!
if [ \! -r "${stat_output}" ]; then
  echo "name,type,original_reads,reads,one_hits,failed,samples,rpm,count_table" > ${stat_output}
fi
original_reads=0
if [ -r "outputs/trimomatic_stats.csv" ]; then
  original_reads_tmp=\$(tail -n 1 outputs/trimomatic_stats.csv | awk -F, '{print \$2}')
  original_reads=\${original_reads_tmp:-0}
fi
reads_tmp=\$(grep "^# reads processed" ${bt_input} | awk -F: '{print \$2}' | sed 's/ //g')
reads=\${reads_tmp:-0}
one_align_tmp=\$(grep "^# reads with at least one reported" ${bt_input} | awk -F": " '{print \$2}' | sed 's/ .*//g')
one_align=\${one_align_tmp:-0}
failed_tmp=\$(grep "^# reads that failed to align" ${bt_input} | awk -F": " '{print \$2}' | sed 's/ .*//g')
failed=\${failed_tmp:-0}
sampled_tmp=\$(grep "^# reads with alignments sampled" ${bt_input} | awk -F": " '{print \$2}' | sed 's/ .*//g')
sampled=\${sampled_tmp:-0}
rpm_tmp=\$(perl -e "printf(1000000 / \${one_align})" 2>/dev/null)
rpm=\${rpm_tmp:-0}
stat_string=\$(printf "$options->{jbasename},${bt_type},%s,%s,%s,%s,%s,%s,${count_table}" "\${original_reads}" "\${reads}" "\${one_align}" "\${failed}" "\${sampled}" "\$rpm")
echo "\$stat_string" >> ${stat_output}!;
    my $stats = $class->Submit(
        comment => $comment,
        input => $bt_input,
        jdepends => $options->{jdepends},
        jname => $jname,
        jprefix => $options->{jprefix},
        jstring => $jstring,
        cpus => 1,
        jmem => 1,
        jqueue => 'throughput',);
    return($stats);
}

=head2 C<BT2_Stats>

Collects alignment statistics from bowtie 2.

=cut
sub BT2_Stats {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $bt_input = $options->{bt_input};
    my $bt_type = $options->{bt_type};
    my $jname = "bt2_stats";
    $jname = $options->{jname} if ($options->{jname});
    my $jobid = qq"$options->{jbasename}_stats";
    my $count_table = "";
    $count_table = $options->{count_table} if ($options->{count_table});
    my $comment = qq!## This is a stupidly simple job to collect alignment statistics.!;
    my $output = "outputs/bowtie2_stats.csv";
    my $jstring = qq!
if [ \! -e "${output}" ]; then
    echo "original reads, single hits, failed reads, multi-hits, rpm" > ${output}
fi
original_reads_tmp=\$(grep " reads; of these" "${bt_input}" 2>/dev/null | awk '{print \$1}' | sed 's/ //g')
original_reads=\${original_reads_tmp:-0}
one_align_tmp=\$(grep " aligned exactly 1 time" "${bt_input}" | awk '{print \$1}' | sed 's/ .*//g')
one_align=\${one_align_tmp:-0}
failed_tmp=\$(grep " aligned 0 times" "${bt_input}" | awk '{print \$1}' | sed 's/ .*//g')
failed=\${failed_tmp:-0}
sampled_tmp=\$(grep " aligned >1 times" "${bt_input}" | awk '{print \$1}' | sed 's/ .*//g')
sampled=\${sampled_tmp:-0}
rpm_tmp=\$(perl -e "printf(1000000 / \$(( \${one_align} + \${sampled} )) ) " 2>/dev/null)
rpm=\${rpm_tmp:-0}
stat_string=\$(printf "$options->{jbasename},${bt_type},%s,%s,%s,%s,%s" "\${original_reads}" "\${one_align}" "\${failed}" "\${sampled}" "\${rpm}")
echo "\$stat_string" >> ${output}!;
    my $stats = $class->Submit(
        comment => $comment,
        input => $bt_input,
        jname => $jname,
        jdepends => $options->{jdepends},
        jprefix => $options->{jprefix},
        jstring => $jstring,
        cpus => 1,
        jmem => 1,
        jqueue => 'throughput',);
    return($stats);
}

=head2 C<BWA>

Perform a bwa alignment using both the sam(s|p)e and aln algorithms.  It then
converts the output (when appropriate) to sorted/indexed bam and passes them to
htseq.

=cut
sub BWA {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        species => 'lmajor',
        libtype => 'genome',
        jprefix => 30,
        modules => ['bwa'],);
    my $check = which('bwa');
    die("Could not find bwa in your PATH.") unless($check);

    if ($options->{species} =~ /\:/) {
        my @species_lst = split(/:/, $options->{species});
        my @result_lst = ();
        my $start_species = $options->{species};
        foreach my $sp (@species_lst) {
            print "Invoking bwa on ${sp}\n";
            $options->{species} = $sp;
            my $result = Bio::Adventure::Map::BWA($class, %{$options});
            push (@result_lst, $result);
        }
        $options->{species} = $start_species;
        return(@result_lst);
    }

    my $sleep_time = 3;
    my $bwa_input = $options->{input};

    my $test_file = "";
    my $forward_reads = "";
    my $reverse_reads = undef;
    if ($bwa_input =~ /\:|\;|\,|\s+/) {
        my @pair_listing = split(/\:|\;|\,|\s+/, $bwa_input);
        $pair_listing[0] = File::Spec->rel2abs($pair_listing[0]);
        $pair_listing[1] = File::Spec->rel2abs($pair_listing[1]);
        $test_file = $pair_listing[0];
        $forward_reads = qq" <(less $pair_listing[0])";
        $reverse_reads = qq" <(less $pair_listing[1])";
        $bwa_input = qq" ${forward_reads} ${reverse_reads} ";
    } else {
        $test_file = File::Spec->rel2abs($bwa_input);
        $bwa_input = qq" <(less ${test_file}) ";
        $forward_reads = $bwa_input;
    }

    my $jname = qq"bwa_$options->{species}";
    $jname = $options->{jname} if ($options->{jname});
    my $libtype = $options->{libtype};

    my $bwa_dir = qq"outputs/bwa_$options->{species}";
    $bwa_dir = $options->{bwa_dir} if ($options->{bwa_dir});

    my $uncompress_jobid = undef;
    my $index_jobid = undef;
    ## Check that the indexes exist
    my $bwa_reflib = "$options->{libdir}/${libtype}/indexes/$options->{species}.fa";
    my $bwa_reftest = qq"${bwa_reflib}.sa";
    my $index_job;
    if (!-r $bwa_reftest) {
        $index_job = $class->Bio::Adventure::Map::BWA_Index(
            jdepends => $options->{jdepends},
            jprefix => $options->{jprefix} - 1,
            libtype => $libtype,);
        $options->{jdepends} = $index_job->{job_id};
    }

    my $aln_sam = qq"${bwa_dir}/$options->{jbasename}_aln.sam";
    my $mem_sam = qq"${bwa_dir}/$options->{jbasename}_mem.sam";
    my $aln_args = qq"";
    my $mem_args = qq"-M ";
    my $sam_args = qq"";
    my $mem_string = qq!mkdir -p ${bwa_dir}
bwa mem ${mem_args} \\
  -a ${bwa_reflib} \\
  ${bwa_input} \\
  2>${bwa_dir}/bwa.err \\
  1>${mem_sam}
!;
    my $reporter_string = qq"bwa samse ${sam_args} \\
  ${bwa_reflib} \\
  ${bwa_dir}/$options->{jbasename}_aln-forward.sai \\
  ${bwa_input} \\
  2>${bwa_dir}/$options->{jbasename}.samerr \\
  1>${aln_sam}";
    my $aln_string = qq"bwa aln ${aln_args} \\
  ${bwa_reflib} \\
  ${forward_reads} \\
  2>${bwa_dir}/$options->{jbasename}_aln-forward.err \\
  1>${bwa_dir}/$options->{jbasename}_aln-forward.sai";
    if (defined($reverse_reads)) {
        $aln_string = qq"${aln_string}
bwa aln ${aln_args} \\
  ${bwa_reflib} \\
  <(less ${reverse_reads}) \\
  2>${bwa_dir}/$options->{jbasename}_aln-reverse.err \\
  1>${bwa_dir}/$options->{jbasename}_aln-reverse.sai";
        $reporter_string = qq"bwa ${sam_args} \\
  sampe ${bwa_reflib} \\
  ${bwa_dir}/$options->{jbasename}_aln-forward.sai \\
  ${bwa_dir}/$options->{jbasename}_aln-reverse.sai \\
  <(less ${forward_reads}) <(less ${reverse_reads}) \\
  2>${bwa_dir}/$options->{jbasename}.samerr \\
  1>${aln_sam}";
    }

    my $mem_comment = qq!## This is a BWA mem alignment of ${bwa_input} against
## ${bwa_reflib}.
!;
    my $aln_comment = qq!## This is a BWA aln alignment of ${bwa_input} against
## ${bwa_reflib}.
!;
    my $report_comment = qq!## This is a BWA sam report of ${bwa_input} against
## ${bwa_reflib}.
!;

    ## MEM Runs
    my $bwa_job = $class->Submit(
        comment => $mem_comment,
        input => $bwa_input,
        jdepends => $options->{jdepends},
        jname => "bwamem_$options->{species}",
        output => $mem_sam,
        jprefix => $options->{jprefix},
        jstring => $mem_string,
        jmem => '36',
        modules => $options->{modules},
        postscript => $options->{postscript},
        prescript => $options->{prescript},
        jqueue => 'workstation',
    );
    my $mem_sam_job = $class->Bio::Adventure::Convert::Samtools(
        input => $mem_sam,
        jdepends => $bwa_job->{job_id},
        jname => "s2b_mem",
        jmem => '20',
        jprefix => $options->{jprefix} + 1,);
    $bwa_job->{samtools_mem} = $mem_sam_job;

    ## ALN Runs
    my $aln_job = $class->Submit(
        comment => $aln_comment,
        input => $bwa_input,
        jdepends => $mem_sam_job,
        jname => "bwaaln_$options->{species}",
        output => qq"${bwa_dir}/$options->{jbasename}_aln-forward.sai",
        jprefix => $options->{jprefix} + 2,
        jstring => $aln_string,
        jmem => '36',
        modules => $options->{modules},
        postscript => $options->{postscript},
        prescript => $options->{prescript},
        jqueue => 'workstation',
    );
    $bwa_job->{aln} = $aln_job;
    my $rep_job = $class->Submit(
        comment => $report_comment,
        input => $aln_job->{output},
        jdepends => $aln_job->{job_id},
        jname => "bwarep_$options->{species}",
        modules => $options->{modules},
        output => $aln_sam,
        jprefix => $options->{jprefix} + 3,
        jstring => $reporter_string,
        jmem => '36',
        postscript => $options->{postscript},
        prescript => $options->{prescript},
        jqueue => 'workstation',);
    $bwa_job->{reporter} = $rep_job;
    my $aln_sam_job = $class->Bio::Adventure::Convert::Samtools(
        input => $aln_sam,
        jdepends => $rep_job->{job_id},
        jname => "s2b_aln",
        jmem => '20',
        jprefix => $options->{jprefix} + 4,);
    $bwa_job->{samtools_aln} = $aln_sam_job;

    my $mem_htmulti = $class->Bio::Adventure::Count::HT_Multi(
        htseq_id => $options->{htseq_id},
        htseq_input => $mem_sam_job->{output},
        htseq_type => $options->{htseq_type},
        jdepends => $mem_sam_job->{job_id},
        jname => "htmem_${jname}",
        jprefix => $options->{jprefix} + 5,
        mapper => 'bwa',);
    $bwa_job->{htseq_mem} = $mem_htmulti;

    my $aln_htmulti = $class->Bio::Adventure::Count::HT_Multi(
        htseq_id => $options->{htseq_id},
        htseq_input => $aln_sam_job->{output},
        htseq_type => $options->{htseq_type},
        jdepends => $aln_sam_job->{job_id},
        jname => "htaln_${jname}",
        jprefix => $options->{jprefix} + 6,
        mapper => 'bwa',);
    $bwa_job->{htseq_aln} = $aln_htmulti;

    my $bwa_stats = $class->Bio::Adventure::Map::BWA_Stats(
        jdepends => $mem_sam_job->{job_id},
        jname => 'bwastats',
        jprefix => $options->{jprefix} + 7,
        aln_output => $aln_sam_job->{output},
        mem_output => $mem_sam_job->{output},);
    $bwa_job->{stats} = $bwa_stats;

    return($bwa_job);
}

=head2 C<BWA_Index>

Create bwa indexes.

=cut
sub BWA_Index {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        modules => ['bwa'],);
    my $jstring = qq!
if test \! -e "$options->{libdir}/genome/$options->{species}.fa"; then
  ln -sf $options->{libdir}/genome/$options->{species}.fasta $options->{libdir}/genome/$options->{species}.fa
fi
start=\$(pwd)
cd $options->{libdir}/$options->{libtype}/indexes &&
  bwa index $options->{species}.fa \\
  2>bwa_index.err \\
  1>bwa_index.out
cd \$start
!;
    my $comment = qq!## Generating bwa indexes for species: $options->{species} in $options->{libdir}/$options->{libtype}/indexes!;
    my $bwa_index = $class->Submit(
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => "bwaidx",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        modules => $options->{modules},
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($bwa_index);
}

=head2 C<BWA_Stats>

Collect some alignment statistics from bwa.

=cut
sub BWA_Stats {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $aln_input = $options->{aln_output};
    $aln_input = qq"${aln_input}.stats";
    my $mem_input = $options->{mem_output};
    $mem_input = qq"${mem_input}.stats";
    my $stat_output = qq"outputs/bwa_stats.csv";

    my $jname = "bwa_stats";
    $jname = $options->{jname} if ($options->{jname});
    my $jobid = qq"$options->{jbasename}_stats";
    my $count_table = "";
    $count_table = $options->{count_table} if ($options->{count_table});
    my $comment = qq!## This is a stupidly simple job to collect alignment statistics.!;
    my $jstring = qq!
if [ \! -r "${stat_output}" ]; then
    echo "# original reads, reads used, aln-aligned reads, mem-aligned reads, rpm" > ${stat_output}
fi
original_reads=0
if [ -r "outputs/trimomatic_stats.csv" ]; then
    original_reads_tmp=\$(tail -n 1 outputs/trimomatic_stats.csv | awk -F, '{print \$2}')
    original_reads=\${original_reads_tmp:-0}
fi
reads_tmp=\$(grep "^Total reads: " ${aln_input} | awk '{print \$3}' | sed 's/ //g')
reads=\${reads_tmp:-0}
aln_aligned_tmp=\$(grep "^Mapped reads" ${aln_input} | awk '{print \$3}' | sed 's/ .*//g')
aln_aligned=\${aln_aligned_tmp:-0}
mem_aligned_tmp=\$(grep "^Mapped reads" ${mem_input} | awk '{print \$3}' | sed 's/ .*//g')
mem_aligned=\${mem_aligned_tmp:-0}
rpm_tmp=\$(perl -e "printf(1000000 / \${aligned})" 2>/dev/null)
rpm=\${rpm_tmp:-0}
stat_string=\$(printf "$options->{jbasename},%s,%s,%s,%s,%s,${count_table}" "\${original_reads}" "\${reads}" "\${aln_aligned}" "\${mem_aligned}" "\$rpm")
echo "\${stat_string}" >> ${stat_output}!;
    my $stats = $class->Submit(
        comment => $comment,
        input => $aln_input,
        depends => $options->{jdepends},
        jname => $jname,
        jprefix => $options->{jprefix},
        jstring => $jstring,
        cpus => 1,
        jmem => 1,
        jqueue => "throughput",);
    return($stats);
}

=head2 C<Hisat2>

Invoke hisat2

=cut
sub Hisat2 {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input',],
        htseq_type => 'gene',
        htseq_id => 'ID',
        do_htseq => 1,
        jprefix => '40',
        libtype => 'genome',
        modules => ['hisat2', 'samtools', 'htseq'],
        );
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('hisat2-build');
    die("Could not find hisat2 in your PATH.") unless($check);

    if ($options->{species} =~ /:/) {
        my $start_species = $options->{species};
        my @species_lst = split(/:/, $options->{species});
        my @result_lst = ();
        foreach my $sp (@species_lst) {
            print "Invoking hisat2 on ${sp}\n";
            $options->{species} = $sp;
            my $result = $class->Bio::Adventure::Map::Hisat2(%{$options});
            push (@result_lst, $result);
        }
        $options->{species} = $start_species;
        return(@result_lst);
    }

    my $ready;
    if (!$options->{jdepends}) {
        $ready = $class->Check_Input(
            files => $options->{input},);
    }
    my $sleep_time = 3;
    my $hisat_args = '';
    $hisat_args = $options->{hisat_args} if ($options->{hisat_args});

    my $prefix_name = qq"hisat2";
    my $hisat_name = qq"${prefix_name}_$options->{species}_$options->{libtype}";
    my $suffix_name = $prefix_name;
    if ($options->{jname}) {
        $hisat_name .= qq"_$options->{jname}";
        $suffix_name .= qq"_$options->{jname}";
    }

    my $hisat_dir = qq"outputs/$options->{jprefix}hisat2_$options->{species}";
    my $hisat_input = $options->{input};
    my $test_file = "";
    my $paired = 0;
    if ($hisat_input =~ /\:|\;|\,|\s+/) {
        my @pair_listing = split(/\:|\;|\,|\s+/, $hisat_input);
        $paired = 1;
        $pair_listing[0] = File::Spec->rel2abs($pair_listing[0]);
        $pair_listing[1] = File::Spec->rel2abs($pair_listing[1]);
        ## After years of working without problem, suddenly my lesspipe
        ## process substitution pre-filter for arbitrarily compressed files
        ## stopped working and might well give me an aneurysm trying to figure out.
        $hisat_input = qq" -1 <(less $pair_listing[0]) -2 <(less $pair_listing[1]) ";
        ##$hisat_input = qq" -1 $pair_listing[0] -2 $pair_listing[1] ";
        $test_file = $pair_listing[0];
    } else {
        $test_file = File::Spec->rel2abs($hisat_input);
        ## $hisat_input = qq" ${hisat_input} ";
        $hisat_input = qq" <(less ${test_file}) ";
    }

    ## Check that the indexes exist
    my $hisat_reflib = "$options->{libdir}/$options->{libtype}/indexes/$options->{species}";
    my $hisat_reftest = qq"${hisat_reflib}.1.ht2";
    my $hisat_reftestl = qq"${hisat_reflib}.1.ht2l";

    if (!-r $hisat_reftest && !-r $hisat_reftestl) {
        print "Hey! The Indexes do not appear to exist, check this out: ${hisat_reftest}\n";
        sleep(20);
        my $index_job = Bio::Adventure::Map::HT2_Index(
            $class,
            jprefix => $options->{jprefix} - 1,
            jdepends => $options->{jdepends},
            libtype => $options->{libtype},);
        $options->{jdepends} = $index_job->{job_id};
    }
    my $hisat_input_flag = "-q "; ## fastq by default
    $hisat_input_flag = "-f " if (${hisat_input} =~ /\.fasta$/);

    my $cpus = $options->{cpus};
    my $error_file = qq"${hisat_dir}/hisat2_$options->{species}_$options->{libtype}_$options->{jbasename}.err";
    my $comment = qq!## This is a hisat2 alignment of ${hisat_input} against ${hisat_reflib}
!;
    $comment .= qq"## This alignment is using arguments: ${hisat_args}.\n" unless ($hisat_args eq '');
    my $aligned_discordant_filename = qq"${hisat_dir}/$options->{jbasename}_aldis_$options->{species}_$options->{libtype}.fastq.gz";
    my $unaligned_discordant_filename = qq"${hisat_dir}/$options->{jbasename}_unaldis_$options->{species}_$options->{libtype}.fastq.gz";
    my $aligned_concordant_filename = qq"${hisat_dir}/$options->{jbasename}_alcon_$options->{species}_$options->{libtype}.fastq.gz";
    my $unaligned_concordant_filename = qq"${hisat_dir}/$options->{jbasename}_unalcon_$options->{species}_$options->{libtype}.fastq.gz";
    my $sam_filename = qq"${hisat_dir}/$options->{jbasename}_$options->{species}_$options->{libtype}.sam";
    my $jstring = qq!mkdir -p ${hisat_dir}
sleep ${sleep_time}
hisat2 -x ${hisat_reflib} ${hisat_args} \\
  -p ${cpus} \\
  ${hisat_input_flag} ${hisat_input} \\
  --phred$options->{phred} \\
  --un-gz ${unaligned_discordant_filename} \\
  --al-gz ${aligned_discordant_filename} \\
  --un-conc-gz ${unaligned_concordant_filename} \\
  --al-conc-gz ${aligned_concordant_filename} \\
  -S ${sam_filename} \\
  2>${error_file} \\
  1>${hisat_dir}/hisat2_$options->{species}_$options->{libtype}_$options->{jbasename}.out
!;
    ## Example: r1_trimmed_unaligned_concordant_lpanamensis_v36.fastq.1.gz

    my $unaligned_filenames = $unaligned_concordant_filename;
    if ($paired) {
        my $tmp = basename($unaligned_filenames, ('.gz'));
        my $dir = dirname($unaligned_filenames);
        $unaligned_filenames = qq"${dir}/${tmp}.1.gz:${dir}/${tmp}.2.gz";
    }
    my $hisat_job = $class->Submit(
        aligned => $aligned_concordant_filename,
        comment => $comment,
        jdepends => $options->{jdepends},
        input => $hisat_input,
        jname => $hisat_name,
        jstring => $jstring,
        jprefix => $options->{jprefix},
        jmem => 48,
        modules => $options->{modules},
        paired => $paired,
        output => $sam_filename,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        unaligned => $unaligned_filenames,);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');

    ## HT1_Stats also reads the trimomatic output, which perhaps it should not.
    ## my $trim_output_file = qq"outputs/$options->{jbasename}-trimomatic.out";
    my $new_jprefix = qq"$options->{jprefix}_1";
    my $stats = $class->Bio::Adventure::Map::HT2_Stats(
        ht_input => $error_file,
        count_table => qq"$options->{jbasename}.count.xz",
        jdepends => $hisat_job->{job_id},
        jname => qq"hisat2st_${suffix_name}",
        jprefix => $new_jprefix,
        output_dir => $hisat_dir,
        ## trim_input => ${trim_output_file},
        );

    my $sam_jprefix = qq"$options->{jprefix}_2";
    my $sam_jname = qq"s2b_${suffix_name}";
    my $sam_job = $class->Bio::Adventure::Convert::Samtools(
        input => $sam_filename,
        jdepends => $hisat_job->{job_id},
        jprefix => $sam_jprefix,
        jname => $sam_jname,
        paired => $paired,
        species => $options->{species},);

    $new_jprefix = qq"$options->{jprefix}_3";
    my $htseq_input;
    if ($paired == 1) {
        $htseq_input = $sam_job->{paired_output};
    } else {
        $htseq_input = $sam_job->{output};
    }
    my $htmulti;
    if ($options->{do_htseq}) {
        if ($options->{libtype} eq 'rRNA') {
            $htmulti = $class->Bio::Adventure::Count::HTSeq(
                htseq_id => $options->{htseq_id},
                htseq_type => $options->{htseq_type},
                input => $htseq_input,
                jdepends => $sam_job->{job_id},
                jname => $suffix_name,
                jprefix => $new_jprefix,
                libtype => $options->{libtype},
                mapper => 'hisat2',
                paired => $paired,);
        } else {
            $htmulti = $class->Bio::Adventure::Count::HT_Multi(
                htseq_id => $options->{htseq_id},
                htseq_type => $options->{htseq_type},
                input => $htseq_input,
                jdepends => $sam_job->{job_id},
                jname => $suffix_name,
                jprefix => $new_jprefix,
                libtype => $options->{libtype},
                mapper => 'hisat2',
                paired => $paired,);
            $hisat_job->{htseq_job} = $htmulti;
        }
    }  ## End checking if we should do htseq
    return($hisat_job);
}

=head2 C<HT2_Index>

Create a hisat2 index using ${species}.fasta and leaves it in the indexes/
directory.

=cut
sub HT2_Index {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ["species"],
        jprefix => '21',);
    my $libtype = $options->{libtype};
    my $libdir = File::Spec->rel2abs($options->{libdir});
    my $jstring = qq!
if test \! -e "${libdir}/genome/$options->{species}.fa"; then
  ln -sf ${libdir}/genome/$options->{species}.fasta ${libdir}/genome/$options->{species}.fa
fi
if test \! -e "${libdir}/genome/indexes/$options->{species}.fa"; then
  ln -sf ${libdir}/genome/$options->{species}.fasta ${libdir}/genome/indexes/$options->{species}.fa
fi

hisat2-build $options->{libdir}/genome/$options->{species}.fasta \\
  $options->{libdir}/${libtype}/indexes/$options->{species}
!;
    my $comment = qq!## Generating hisat2 indexes for species: $options->{species} in $options->{libdir}/${libtype}/indexes!;
    my $indexer = $class->Submit(
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => qq"ht2idx_$options->{species}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($indexer);
}

=head2 C<HT2_Stats>

Collect alignment statistics from hisat 2.

=cut
sub HT2_Stats {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        output_dir => 'outputs',);
    my $ht_input = $options->{ht_input};
    my $jname = "ht2_stats";
    $jname = $options->{jname} if ($options->{jname});
    my $jobid = qq"$options->{jbasename}_stats";
    my $count_table = "";
    $count_table = $options->{count_table} if ($options->{count_table});
    my $comment = qq!## This is a stupidly simple job to collect alignment statistics.!;
    my $output = "$options->{output_dir}/hisat2_stats.csv";
    my $jstring = qq!
if [ \! -e "${output}" ]; then
    echo "id, original reads, single hits, failed reads, multi-hits, rpm" > ${output}
fi
original_reads_tmp=\$(grep " reads; of these" "${ht_input}" 2>/dev/null | awk '{print \$1}' | sed 's/ //g')
original_reads=\${original_reads_tmp:-0}
one_align_tmp=\$(grep " aligned exactly 1 time" "${ht_input}" | awk '{print \$1}' | sed 's/ .*//g')
one_align=\${one_align_tmp:-0}
failed_tmp=\$(grep " aligned 0 times" "${ht_input}" | tail -n 1 | awk '{print \$1}' | sed 's/ .*//g')
failed=\${failed_tmp:-0}
sampled_tmp=\$(grep " aligned >1 times" "${ht_input}" | awk '{print \$1}' | sed 's/ .*//g')
sampled=\${sampled_tmp:-0}
rpm_tmp=\$(perl -e "printf(1000000 / \$(( \${one_align} + \${sampled} )) ) " 2>/dev/null)
rpm=\${rpm_tmp:-0}
stat_string=\$(printf "$options->{jbasename},%s,%s,%s,%s,%s" "\${original_reads}" "\${one_align}" "\${failed}" "\${sampled}" "\${rpm}")
echo "\$stat_string" >> ${output}!;
    my $stats = $class->Submit(
        comment => $comment,
        input => $ht_input,
        jname => $jname,
        jdepends => $options->{jdepends},
        jprefix => $options->{jprefix},
        jstring => $jstring,
        cpus => 1,
        jmem => 1,
        jqueue => 'throughput',);
    return($stats);
}

=head2 C<Kallisto>

Perform a kallisto transcript quantification.

=cut
sub Kallisto {
    my ($class, %args) = @_;
    my $check = which('kallisto');
    die("Could not find kallisto in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ["species","input"],
    );

    if ($options->{species} =~ /\:/) {
        my @species_lst = split(/:/, $options->{species});
        my @result_lst = ();
        my $start_species = $options->{species};
        foreach my $sp (@species_lst) {
            print "Invoking kallisto on ${sp}\n";
            $options->{species} = $sp;
            my $result = Bio::Adventure::Map::Kallisto($class, %{$options});
            push (@result_lst, $result);
        }
        $options->{species} = $start_species;
        return(@result_lst);
    }

    my $sleep_time = 3;
    my %ka_jobs = ();
    my $ka_depends_on = '';
    my $libtype = 'genome';
    $libtype = $options->{libtype} if ($options->{libtype});
    my $species = $options->{species};

    my $jname = qq"kall_${species}";
    ## $jname = $options->{jname} if ($options->{jname});
    my $ka_args = qq"";
    my $ka_input = $options->{input};
    my $input_name = $ka_input;
    if ($ka_input =~ /\:|\;|\,|\s+/) {
        my @pair_listing = split(/\:|\;|\,|\s+/, $ka_input);
        $ka_args .= " --bias ";
        if ($options->{stranded} != 0) {
            $ka_args .= " --$options->{stranded} ";
        }
        $ka_input = qq" <(less $pair_listing[0]) <(less $pair_listing[1]) ";
        $input_name = $pair_listing[0];
    } else {
        $ka_input = qq" <(less $ka_input) ";
        $ka_args .= " --bias --single -l 40 -s 10 ";
    }

    ## Check that the indexes exist
    my $ka_reflib = "$options->{libdir}/${libtype}/indexes/$options->{species}.idx";
    if (!-r $ka_reflib) {
        my $index_job = $class->Bio::Adventure::Map::Kallisto_Index(
            jdepends => $options->{jdepends},
            libtype => $libtype,);
        $ka_jobs{index} = $index_job;
        $options->{jdepends} = $index_job->{job_id};
    }

    my $outdir = qq"outputs/kallisto_${species}";
    my $error_file = qq"${outdir}/kallisto_${species}.stderr";
    my $output_sam = qq"${outdir}/kallisto_${species}.sam";
    my $output_bam = qq"${outdir}/kallisto_${species}.bam";
    my $output_stats = qq"${outdir}/kallisto_${species}.stats";
    my $sorted_bam = qq"${outdir}/kallisto_${species}-sorted";
    my $comment = qq!## This is a kallisto pseudoalignment of ${ka_input} against
## ${ka_reflib}.
## This jobs depended on: $options->{jdepends}
## Other candidates for making a pretty count table include:
##  perl -F'\\t' -a -n -e 'print "\$F[0] \$F[3]\\n"' ${outdir}/abundance.tsv > ${outdir}/abundance.count
##   awk '{printf("%s %s\\n", \$1, \$4)}' ${outdir}/abundance.tsv > ${outdir}/abundance.count
## The sam->bam conversion is copy/pasted from Sam2Bam() but I figured why start another job
## because kallisto is so fast
!;
    my $dropped_args = qq" --pseudobam ";
    my $jstring = qq!mkdir -p ${outdir} && sleep ${sleep_time} && \\
kallisto quant ${ka_args} \\
  --plaintext -t 4 -b 100 \\
  -o ${outdir} \\
  -i ${ka_reflib} \\
  ${ka_input} \\
  2>${error_file} \\
  1>${output_sam} && \\
  cut -d "	" -f 1,4 ${outdir}/abundance.tsv > ${outdir}/${input_name}_abundance.count && \\
  gzip -9 -f ${outdir}/${input_name}_abundance.count
!;

    ## I am going to stop doing these pseudobam indexes because that is pretty dumb for kallisto to do
    ## It was interesting for the exosome samples, but only because I wanted to compare them against
    ## the full miRNA database.
    my $unused_material = qq!
  samtools view -u -t ${ka_reflib} -S ${output_sam} \\
    2>${output_bam}.err \\
    1>${output_bam} && \\
  samtools sort -l 9 ${output_bam} ${sorted_bam} \\
    2>${sorted_bam}.err \\
    1>${sorted_bam}.out && \\
  rm ${output_bam} && mv ${sorted_bam}.bam ${output_bam} && samtools index ${output_bam} && \\
  bamtools stats -in ${output_bam} 2>${output_stats} 1>&2
!;
    my $kallisto = $class->Submit(
        comment => $comment,
        input => $ka_input,
        depends => $options->{jdepends},
        jname => qq"${jname}",
        jprefix => '30',
        jstring => $jstring,
        jmem => 30,
        prescript => $args{prescript},
        postscript => $args{postscript},);
    return($kallisto);
}

=head2 C<Kallisto_Index

Use kallisto and an annotated_CDS fasta sequence library to create an index.

=cut
sub Kallisto_Index {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        modules => ['kallisto'],
        required => ["species", "genome"],
    );

    my $libtype = $options->{libtype};
    my $genome = File::Spec->rel2abs($options->{genome});
    unless (-r $genome) {
        die("The indexing operation for kallisto will fail because the $options->{species} genome does not exist.")
    }

    my $jstring = qq!
kallisto index -i $options->{libdir}/${libtype}/indexes/$options->{species}.idx ${genome}!;
    my $comment = qq!## Generating kallisto indexes for species: $options->{species} in $options->{libdir}/${libtype}/indexes!;
    my $ka_index = $class->Submit(
        comment => $comment,
        jdepends => $options->{jdepends},
        jstring => $jstring,
        jname => "kalidx",
        jprefix => $options->{jprefix},
        modules => $options->{modules},
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($ka_index);
}

=head2 C<RSEM>

Invoke RSEM.

=cut
sub RSEM {
    my ($class, %args) = @_;
    my $check = which('rsem-prepare-reference');
    die("Could not find RSEM in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input'],
        modules => ['rsem'],);

    if ($options->{species} =~ /\:/) {
        my @species_lst = split(/:/, $options->{species});
        my @result_lst = ();
        my $start_species = $options->{species};
        foreach my $sp (@species_lst) {
            print "Invoking RSEM on ${sp}\n";
            $options->{species} = $sp;
            my $result = Bio::Adventure::Map::RSEM($class, %{$options});
            push (@result_lst, $result);
        }
        $options->{species} = $start_species;
        return(@result_lst);
    }

    my $jbasename = $class->Get_Job_Name();
    my $cds = qq"$options->{libdir}/$options->{libtype}/$options->{species}_cds_nt.fasta";
    my $idx = qq"$options->{libdir}/$options->{libtype}/indexes/rsem/$options->{species}";
    my $test_idx = qq"${idx}.transcripts.fa";
    my $index_job;
    unless (-r $test_idx) {
        print "Need to create the RSEM indexes, invoking RSEM_Index().\n";
        unless (-r $cds) {
            die("RSEM_Index requires a cds fasta file at: ${cds}, create it with a cyoa2 conversion.");
        }
        $index_job = $class->Bio::Adventure::Map::RSEM_Index(
            cds_fasta => $cds,
            index => $idx,
            jdepends => $options->{jdepends},
            jname => qq"rsidx_${jbasename}",
            libtype => $options->{libtype},
            modules => $options->{modules},);
        $options->{jdepends} = $index_job->{job_id};
    }

    my $rsem_input = $options->{input};
    my $test_file = "";
    if ($rsem_input =~ /\:|\;|\,|\s+/) {
        my @pair_listing = split(/\:|\;|\,|\s+/, $rsem_input);
        $rsem_input = qq"--paired-end <(less $pair_listing[0]) <(less $pair_listing[1])";
        $test_file = $pair_listing[0];
    } else {
        $rsem_input = qq" <(less $rsem_input) ";
    }

    my $rsem_dir = qq"outputs/rsem_$options->{species}";
    my $rsem_comment = qq"## This is a rsem invocation script.
";
    my $output_file = "$options->{species}_something.txt";
    my $jstring = qq"mkdir -p ${rsem_dir} && rsem-calculate-expression --bowtie2 \\
  --calc-ci ${rsem_input} \\
  ${idx} \\
  ${rsem_dir}/$options->{species} \\
  2>${rsem_dir}/$options->{species}.err \\
  1>${rsem_dir}/$options->{species}.out
";

    my $rsem = $class->Submit(
        comment => $rsem_comment,
        input => $rsem_input,
        jname => qq"rsem_${jbasename}",
        jdepends => $options->{jdepends},
        jstring => $jstring,
        jprefix => "28",
        mem => 24,
        modules => $options->{modules},
        output => $output_file,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    $rsem->{index_job} = $index_job;
    return($rsem);
}

=item C<RSEM_Index

Use RSEM and an annotated_CDS fasta sequence library to create a transcript index.

=cut
sub RSEM_Index {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ["species",],
        modules => ['rsem']);
    my $comment = qq"## RSEM Index creation.";
    my $jstring = qq!
rsem-prepare-reference --bowtie2 $options->{cds_fasta} $options->{index}
!;
    my $jobid = $class->Submit(
        comment => $comment,
        jdepends => $options->{jdepends},
        jstring => $jstring,
        jname => "rsemidx",
        jprefix => $options->{jprefix},
        modules => $options->{modules},
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($jobid);
}

=head2 C<Salmon>

Perform a salmon quantification of transcript abundances.

=cut
sub Salmon {
    my ($class, %args) = @_;
    my $check = which('salmon');
    die("Could not find salmon in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ["species","input"],
        modules => ['salmon'],);

    if ($options->{species} =~ /\:/) {
        my @species_lst = split(/:/, $options->{species});
        my @result_lst = ();
        my $start_species = $options->{species};
        foreach my $sp (@species_lst) {
            print "Invoking salmon on ${sp}\n";
            $options->{species} = $sp;
            my $result = Bio::Adventure::Map::Salmon($class, %{$options});
            push (@result_lst, $result);
        }
        $options->{species} = $start_species;
        return(@result_lst);
    }

    my $ready = $class->Check_Input(files => $options->{input});
    my $sleep_time = 3;
    my %sa_jobs = ();
    my $libtype = 'genome';
    $libtype = $options->{libtype} if ($options->{libtype});
    my $species = $options->{species};

    my $jname = qq"sal_${species}";
    ## $jname = $options->{jname} if ($options->{jname});
    my $sa_args = qq"";
    my $sa_input = $options->{input};
    my $input_name = $sa_input;
    if ($sa_input =~ /\:|\;|\,|\s+/) {
        my @pair_listing = split(/\:|\;|\,|\s+/, $sa_input);
        $sa_args .= " -1 <(less $pair_listing[0]) -2 <(less $pair_listing[1]) ";
        $input_name = $pair_listing[0];
    } else {
        $sa_args .= " -r <(less $sa_input) ";
    }

    ## Check that the indexes exist
    my $sa_reflib = "$options->{libdir}/${libtype}/indexes/$options->{species}_salmon_index";
    my $index_job;
    if (!-r $sa_reflib) {
        $index_job = $class->Bio::Adventure::Map::Salmon_Index(
            depends => $options->{jdepends},
            libtype => $libtype,);
        $options->{jdepends} = $index_job->{job_id};
    }

    my $outdir = qq"outputs/salmon_${species}";
    my $error_file = qq"${outdir}/salmon_${species}.stderr";
    my $comment = qq!## This is a salmon pseudoalignment of ${sa_input} against
## ${sa_reflib}.
OB## This jobs depended on: $options->{jdepends}
!;
    my $jstring = qq!mkdir -p ${outdir} && sleep ${sleep_time} && \\
salmon quant -i ${sa_reflib} \\
  -l A --gcBias --validateMappings  \\
  ${sa_args} \\
  -o ${outdir} \\
  2>${outdir}/salmon.err 1>${outdir}/salmon.out
!;

    my $salmon = $class->Submit(
        comment => $comment,
        input => $sa_input,
        jdepends => $options->{jdepends},
        jname => qq"${jname}",
        jprefix => '30',
        jstring => $jstring,
        jmem => 48,
        modules => $options->{modules},
        prescript => $args{prescript},
        postscript => $args{postscript},);
    $salmon->{index_job} = $index_job;

    my $stats = $class->Bio::Adventure::Map::Salmon_Stats(
        input => qq"${outdir}/lib_format_counts.json",
        jdepends => $salmon->{job_id},
        jname => "sastats_$options->{species}",
        jprefix => "33",);
    $salmon->{stats_job} = $stats;
    return($salmon);
}

=head2 C<Salmon_Index>

Invoke salmon with an annotated_CDS fasta sequence library to create a
transcript index.

=cut
sub Salmon_Index {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ["species", "genome"],
        modules => ['salmon'],);
    my $libtype = $options->{libtype};
    my $genome = File::Spec->rel2abs($options->{genome});
    unless (-r $genome) {
        die("The indexing operation for salmon will fail because the $options->{species} genome does not exist.")
    }

    my $jstring = qq!
salmon index -t ${genome} -i $options->{libdir}/${libtype}/indexes/$options->{species}_salmon_index!;
    my $comment = qq!## Generating salmon indexes for species: $options->{species} in $options->{libdir}/${libtype}/indexes!;
    my $jobid = $class->Submit(
        comment => $comment,
        jdepends => $options->{jdepends},
        jstring => $jstring,
        jname => "salidx_$options->{species}",
        jmem => 24,
        jprefix => "15",
        modules => $options->{modules},
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($jobid);
}

=head2 C<Salmon_Stats>

Collect some summary statistics from a salmon run.

=cut
sub Salmon_Stats {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $jname = "stats";
    $jname = $options->{jname} if ($options->{jname});
    my $jobid = qq"$options->{jbasename}_stats";
    my $output = "outputs/salmon_stats.csv";
    my $comment = qq!## This is a stupidly simple job to collect salmon alignment statistics.!;
    my $jstring = qq!
if [ \! -r "${output}" ]; then
  echo "basename,species,fragments,assigned,consistent,inconsistent,bias" > ${output}
fi
reads_tmp=\$(grep "^num_compatible" $options->{input} | awk '{print \$3}' | sed 's/^ *//g')
reads=\${reads_tmp:-0}
aligned_tmp=\$(grep "^num_assigned" $options->{input} | awk '{print \$3}' | sed 's/^ *//g')
aligned=\${aligned_tmp:-0}
consistent_tmp=\$(grep "^concordant" $options->{input} | awk '{print \$3}' | sed 's/^ *//g')
consistent=\${consistent_tmp:-0}
inconsistent_tmp=\$(grep "^inconsistent" $options->{input} | awk '{print \$3}' | sed 's/^ *//g')
inconsistent=\${inconsistent_tmp:-0}
bias_tmp=\$(grep "^mapping_bias" $options->{input} | awk '{print \$3}' | sed 's/^ *//g')
bias=\${bias_tmp:-0}
stat_string=\$(printf "$options->{jbasename},$options->{species},%s,%s,%s,%s,%s" "\${reads}" "\${aligned}" "\${consistent}" "\${inconsistent}" "\${bias}")
echo "\$stat_string" >> "${output}"!;
    my $stats = $class->Submit(
        comment => $comment,
        cpus => 1,
        input => $options->{input},
        jname => $jname,
        jdepends => $options->{jdepends},
        jprefix => $args{jprefix},
        jstring => $jstring,
        jmem => 1,
        jqueue => 'throughput',);
    return($stats);
}

=head2 C<STAR>

Invoke STAR for transcript abundances.

=cut
sub STAR {
    my ($class, %args) = @_;
    my $check = which('STAR');
    die("Could not find STAR in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species','input'],
        modules => ['star'],);
    if ($options->{species} =~ /\:/) {
        my @species_lst = split(/:/, $options->{species});
        my @result_lst = ();
        my $start_species = $options->{species};
        foreach my $sp (@species_lst) {
            print "Invoking STAR on ${sp}\n";
            $options->{species} = $sp;
            my $result = $class->Bio::Adventure::Map::STAR(%{$options});
            push (@result_lst, $result);
        }
        return(@result_lst);
    }

    my $ready = $class->Check_Input(
        files => $options->{input},);

    my $sleep_time = 3;
    my $libtype = 'genome';
    $libtype = $options->{libtype} if ($options->{libtype});
    my $species = $options->{species};

    my $jname = qq"star_${species}";
    ## $jname = $options->{jname} if ($options->{jname});
    my $star_inputstring = qq"";
    my $star_input = $options->{input};
    my $input_name = $star_input;
    if ($star_input =~ /\:|\;|\,|\s+/) {
        my @pair_listing = split(/\:|\;|\,|\s+/, $star_input);
        $star_inputstring = qq"$pair_listing[0],$pair_listing[1]";
        $input_name = $pair_listing[0];
    } else {
        $star_inputstring = qq"$star_input";
    }

    ## Check that the indexes exist
    my $star_refdir = "$options->{libdir}/${libtype}/indexes/$options->{species}_star_index";
    my $star_reflib = qq"${star_refdir}/SAindex";
    my $index_job;
    if (!-r $star_reflib) {
        $index_job = $class->Bio::Adventure::Map::STAR_Index(
            jdepends => $options->{jdepends},
            libtype => $libtype,);
        $options->{jdepends} = $index_job->{job_id};
    }

    my $outdir = qq"outputs/star_${species}";
    my $error_file = qq"${outdir}/star_${species}.stderr";
    my $comment = qq!## This is a star pseudoalignment of ${star_input} against
## ${star_reflib}.
## This jobs depended on: $options->{jdepends}
## Currently, this only works with the module star/git_201803
!;
    my $jstring = qq!mkdir -p ${outdir} && sleep ${sleep_time} && \\
STAR \\
  --genomeDir ${star_refdir} \\
  --outFileNamePrefix outputs/star_$options->{species}/${input_name} \\
  --outSAMtype BAM SortedByCoordinate \\
  --outBAMcompression 10 \\
  --chimOutType WithinBAM \\
  --quantMode GeneCounts \\
  --readFilesIn ${star_inputstring} \\
  --readFilesCommand /usr/bin/lesspipe.sh \\
  --runThreadN 6 \\
  2>outputs/star_gene.out 1>&2

STAR \\
  --genomeDir ${star_refdir} \\
  --outFileNamePrefix outputs/star_$options->{species}/${input_name}_tx \\
  --outSAMtype BAM SortedByCoordinate \\
  --outBAMcompression 10 \\
  --chimOutType WithinBAM \\
  --quantMode TranscriptomeSAM \\
  --readFilesIn ${star_inputstring} \\
  --readFilesCommand /usr/bin/lesspipe.sh \\
  --runThreadN 6 \\
  2>outputs/star_tx.out 1>&2

!;

    my $star_job = $class->Submit(
        comment => $comment,
        input => $star_input,
        jdepends => $options->{jdepends},
        jname => qq"${jname}",
        jprefix => '33',
        jstring => $jstring,
        jmem => 96,
        modules => $options->{modules},
        prescript => $args{prescript},
        postscript => $args{postscript},
        jqueue => 'large',);
    $star_job->{index_job} = $index_job;
    return($star_job);
}

=head2 C<STAR_Index>

Create indexes appropriate for STAR.

=cut
sub STAR_Index {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ["species",],
        modules => ['star'],);
    my $comment = qq"## STAR Index creation.";
    my $libtype = 'genome';
    $libtype = $options->{libtype} if ($options->{libtype});
    my $star_refdir = "$options->{libdir}/${libtype}/indexes/$options->{species}_star_index";
    my $jstring = qq!
STAR \\
  --runMode genomeGenerate \\
  --runThreadN 12 \\
  --genomeDir ${star_refdir} \\
  --genomeFastaFiles $options->{libdir}/${libtype}/$options->{species}.fasta \\
  --sjdbGTFfile $options->{libdir}/${libtype}/$options->{species}.gtf \\
  --limitGenomeGenerateRAM 160000000000
!;
    my $jobid = $class->Submit(
        comment => $comment,
        jdepends => $options->{depends},
        jstring => $jstring,
        jname => "staridx",
        jprefix => $options->{jprefix},
        jmem => 180,
        modules => $options->{modules},
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => 'xlarge',);
    return($jobid);
}

=head2 C<Tophat>

Invokes tophat!  It also sorts/indexes the accepted_hits file, collects some
statistics, and passes the hits to htseq-count.

=cut
sub Tophat {
    my ($class, %args) = @_;
    my $check = which('tophat');
    die("Could not find tophat in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ["species", "input", "htseq_type"],
        modules => ['tophat'],);
    if ($options->{species} =~ /\:/) {
        my @species_lst = split(/:/, $options->{species});
        my @result_lst = ();
        my $start_species = $options->{species};
        foreach my $sp (@species_lst) {
            print "Invoking tophat on ${sp}\n";
            $options->{species} = $sp;
            my $result = Bio::Adventure::Map::Tophat($class, %{$options});
            push (@result_lst, $result);
        }
        $options->{species} = $start_species;
        return(@result_lst);
    }

    my $tophat_cpus = 4;
    my $inputs = $options->{input};
    my @in = split(/:/, $inputs);
    $inputs =~ s/:/ /g;
    my $paired = 0;
    if (scalar(@in) > 1) {
        $inputs = qq" <(less $in[0]) <(less $in[1]) ";
        $paired = 1;
    } else {
        $inputs = qq" <(less $in[0]) ";
    }
    my $tophat_args = ' -g 1 --microexon-search --b2-very-sensitive ';
    if ($options->{tophat_args}) {
        $tophat_args = $options->{tophat_args};
    }
    ## $tophat_args .= ' --no-mixed --no-discordant ' if (scalar(@in) > 1);
    ## $tophat_args .= ' ' if (scalar(@in) > 1);

    my $tophat_queue = $options->{jqueue};
    my $tophat_walltime = '18:00:00';
    my $tophat_mem = 8;
    if ($options->{species} eq 'hsapiens' or $options->{species} eq 'mmusculus') {
        $tophat_queue = 'workstation';
        $tophat_walltime =  '144:00:00';
        $tophat_mem = 20;
    }

    my $tophat_dir = qq"outputs/tophat_$options->{species}";
    if ($options->{tophat_dir}) {
        $tophat_dir = $options->{tophat_dir};
    }
    my $libtype = $options->{libtype};
    my $bt_reflib = "$options->{libdir}/${libtype}/indexes/$options->{species}";
    my $bt_reftest = qq"${bt_reflib}.1.bt2";
    my $index_job = undef;
    if (!-r $bt_reftest) {
        print "Did not find the index for $options->{species} at: ${bt_reflib}, indexing now.\n";
        $index_job = Bio::Adventure::Map::BT2_Index(
            $class,
            jdepends => $options->{jdepends},);
        $options->{jdepends} = $index_job->{job_id};
    }
    my $gtf_file = qq"$options->{libdir}/genome/$options->{species}.gtf";
    if (!-r $gtf_file) {
        print "Missing the gtf file for $options->{species}\n";
        print "Using the gff file.\n";
        $gtf_file =~ s/\.gtf/\.gff/;
        ##my $written = $class->Gff2Gtf(gff => "$class->{libdir}/genome/$class->{species}.gff");
    }

    my $spliced = 0;
    if ($options->{spliced}) {
        $spliced = 1;
    }
    my $jname = qq"th_$options->{species}";
    my $jstring = qq!
mkdir -p ${tophat_dir} && tophat ${tophat_args} \\
  -G ${gtf_file} \\
  -p ${tophat_cpus} -o ${tophat_dir} \\
!;
    if ($spliced) {
        $jstring .= qq!  --no-novel-juncs \\
!;
    }
    $jstring .= qq!  $options->{libdir}/genome/indexes/$options->{species} \\
  ${inputs} \\
  2>outputs/tophat.err \\
  1>outputs/tophat.out && \\
 samtools sort -l 9 -n ${tophat_dir}/accepted_hits.bam -o ${tophat_dir}/accepted_sorted.bam && \\
 samtools index ${tophat_dir}/accepted_hits.bam && \\
 samtools sort -l 9 -n ${tophat_dir}/unmapped.bam -o ${tophat_dir}/unmapped_sorted.bam && \\
 samtools index ${tophat_dir}/unmapped.bam
!;
    if ($paired) {
        $jstring .= qq!
if [ -r "${tophat_dir}/accepted_hits.bam" ]; then
  samtools view -b -f 2 ${tophat_dir}/accepted_hits.bam > ${tophat_dir}/accepted_paired.bam && \\
  samtools index ${tophat_dir}/accepted_paired.bam
fi
!;
    }
    my $comment = qq!## I still have no clue what I am doing when I use tophat...
## However, I know that -g 1 will allow only 1 hit in the case of multihits, but randomly place it
## From the manual:  "If there are more alignments with the same score than this
## number, TopHat will randomly report only this many alignments"
## -N 1 will discard anything with >1 mismatch (default is 2)
## -r adjusts the allowable mean distance between the paired reads
## --mate-std-dev sets the deviation of -r
## --microexon-search will tell it to search short exons for reads >=50
!;
    my $tophat = $class->Submit(
        comment => $comment,
        cpus => $tophat_cpus,
        jdepends => $options->{jdepends},
        jname => $options->{jname},
        jprefix => "31",
        jstring => $jstring,
        jmem => $tophat_mem,
        modules => $options->{mdules},
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => $tophat_queue,
        jwalltime => $tophat_walltime,);

    ## Set the input for htseq
    my $accepted = "${tophat_dir}/accepted_hits.bam";
    $accepted = $options->{accepted_hits} if ($options->{accepted_hits});
    my $count_table = "accepted_hits.count";
    $count_table = $options->{count_table} if ($options->{count_table});
    my $htmulti = $class->Bio::Adventure::Count::HT_Multi(
        htseq_input => $accepted,
        htseq_id => $options->{htseq_id},
        htseq_type => $options->{htseq_type},
        jdepends => $tophat->{job_id},
        jname => qq"hts_$options->{species}",
        jprefix => '32',
        mapper => 'tophat',);
    $tophat->{htseq} = $htmulti;
    ## Perform a separate htseq run using only the successfully paired hits
    if ($paired) {
        my $ht_paired = $class->Bio::Adventure::Count::HT_Multi(
            htseq_input => qq"${tophat_dir}/accepted_paired.bam",
            htseq_id => $options->{htseq_id},
            htseq_type => $options->{htseq_type},
            jdepends => $tophat->{job_id},
            jname => qq"htsp_$options->{species}",
            jprefix => '32',
            mapper => 'tophat',);
    }
    ## Tophat_Stats also reads the trimomatic output, which perhaps it should not.
    my $unaccepted = $accepted;
    $unaccepted =~ s/accepted_hits/unmapped/g;
    my $input_read_info = $accepted;
    $input_read_info =~ s/accepted_hits\.bam/prep_reads\.info/g;
    my $stats = $class->Bio::Adventure::Map::Tophat_Stats(
        accepted_input => $accepted,
        count_table => qq"${count_table}.xz",
        jdepends => $tophat->{job_id},
        jname => "tpstats_$options->{species}",
        jprefix => "33",
        prep_input => $input_read_info,
        unaccepted_input => $unaccepted,);
    $tophat->{stats} = $stats;
    return($tophat);
}

=head2 C<Tophat_Stats>

Collect alignment statistics from the accepted_hits.bam/unaligned.bam files
generated by a tophat run.

=cut
sub Tophat_Stats {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $accepted_input = $options->{accepted_input};
    my $accepted_output = qq"${accepted_input}.stats";
    my $unaccepted_input = $options->{unaccepted_input};
    my $unaccepted_output = qq"${unaccepted_input}.stats";
    my $read_info = $options->{prep_input};
    my $jname = "stats";
    $jname = $options->{jname} if ($options->{jname});
    my $jobid = qq"$options->{jbasename}_stats";
    my $count_table = "";
    $count_table = $options->{count_table} if ($options->{count_table});
    my $output = "outputs/tophat_stats.csv";
    my $comment = qq!## This is a stupidly simple job to collect tophat alignment statistics.!;
    my $jstring = qq!
if [ \! -r "${output}" ]; then
  echo "basename,species,original_reads,aligned_reads,failed_reads,rpm,count_table" > ${output}
fi
bamtools stats < "${accepted_input}" \\
    2>${accepted_output} 1>&2 && \\
  bamtools stats < "${unaccepted_input}" \\
    2>${unaccepted_output} 1>&2

original_reads=0
if [ -r "outputs/trimomatic_stats.csv" ]; then
  original_reads_tmp=\$(tail -n 1 outputs/trimomatic_stats.csv | awk -F, '{print \$2}')
  original_reads=\${original_reads_tmp:-0}
fi
reads_tmp=\$(grep "^reads_in " ${read_info} | awk -F= '{print \$2}' | sed 's/ //g')
reads=\${reads_tmp:-0}
aligned_tmp=\$(grep "^Total reads" ${accepted_output} | awk '{print \$3}' | sed 's/ .*//g')
aligned=\${aligned_tmp:-0}
failed_tmp=\$(grep "^Total reads" ${unaccepted_output} | awk '{print \$3}' | sed 's/ .*//g')
failed=\${failed_tmp:-0}
rpm_tmp=\$(perl -e "printf(1000000 / \${aligned})" 2>/dev/null)
rpm=\${rpm_tmp:-0}
stat_string=\$(printf "$options->{jbasename},$options->{species},%s,%s,%s,%s,%s,${count_table}" "\${original_reads}" "\${reads}" "\${aligned}" "\${failed}" "\$rpm")
echo "\$stat_string" >> "${output}"!;
    my $stats = $class->Submit(
        comment => $comment,
        cpus => 1,
        input => $accepted_input,
        jname => $jname,
        jdepends => $options->{jdepends},
        jprefix => $args{jprefix},
        jstring => $jstring,
        jmem => 1,
        jqueue => 'throughput',);
    return($stats);
}

=back

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

    L<bowtie> L<bowtie2> L<tophat> L<bwa> L<kallisto> L<samtools> L<htseq>

=cut

1;
