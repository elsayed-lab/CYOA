package Bio::Adventure::RNASeq_Map;
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

Bio::Adventure::RNASeq_Map - Perform highthroughput sequence alignments with tools like bowtie/tophat/etc

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

    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        bt_type => 'v0M1',
        count => 1,
        libtype => 'genome',
    );
    my $species = $options->{species};
    my $bt_type = $options->{bt_type};
    my $bt_args = $options->{bt_args}->{$bt_type};
    $bt_args = ' --best -v 0 -M 1 ' if (!defined($bt_args));

    my $sleep_time = 3;
    my %bt_jobs = ();
    my $bt_input = $options->{input};
    my $bt_depends_on;
    $bt_depends_on = $options->{depends} if ($options->{depends});
    my $job_basename = $options->{job_basename};

    my $jname = qq"bt${bt_type}_${species}";
    $jname = $options->{jname} if ($options->{jname});
    my $libtype = $options->{libtype};
    my $count = $options->{count};

    my $bt_dir = qq"outputs/bowtie_${species}";
    $bt_dir = $options->{bt_dir} if ($options->{bt_dir});

    my $uncompress_jobid = undef;
    my $index_jobid = undef;
    ## if ($bt_input =~ /\.gz$|\.bz2$|\.xz$/ ) {
    ##    print "The input needs to be uncompressed, doing that now.\n" if ($options->{debug});
    ##    my $uncomp = Bio::Adventure::Compress::Uncompress(
    ##        $class,
    ##        input => $bt_input,
    ##        depends => $bt_depends_on,
    ##        jname => 'uncomp',
    ##        jprefix => '09',
    ##    );
    ##    $bt_input = basename($bt_input, ('.gz', '.bz2', '.xz'));
    ##    $bt_jobs{uncompress} = $uncomp;
    ##    $options = $class->Set_Vars(input => $bt_input);
    ##    $uncompress_jobid = $uncomp->{job_id};
    ##}

    ## Check that the indexes exist
    my $bt_reflib = "$options->{libdir}/${libtype}/indexes/${species}";
    my $bt_reftest = qq"${bt_reflib}.1.ebwt";
    if (!-r $bt_reftest && !$options->{bt1_indexjobs}) {
        $options = $class->Set_Vars(bt1_indexjobs => 1);
        my $index_job = Bio::Adventure::RNASeq_Map::BT1_Index(
            $class,
            depends => $bt_depends_on,
            libtype => $libtype,
        );
        $bt_jobs{index} = $index_job;
        $index_jobid = $index_job->{job_id};
    }

    ## Make a depends string containing the uncompress, indexer, or both
    if (defined($index_jobid) && defined($uncompress_jobid)) {
        $bt_depends_on = qq"${index_jobid}:${uncompress_jobid}";
    } elsif (defined($index_jobid)) {
        $bt_depends_on = $index_jobid;
    } elsif (defined($uncompress_jobid)) {
        $bt_depends_on = $uncompress_jobid;
    } else {
        $bt_depends_on = "";
    }

    my $bowtie_input_flag = "-q"; ## fastq by default
    $bowtie_input_flag = "-f" if ($options->{input} =~ /\.fasta/);

    my $cpus = $options->{cpus};
    my $error_file = qq"${bt_dir}/${job_basename}-${bt_type}.err";
    my $comment = qq!## This is a bowtie1 alignment of ${bt_input} against
## ${bt_reflib} using arguments: ${bt_args}.
!;
    if ($bt_depends_on) {
        $comment .= qq!## This jobs depended on: ${bt_depends_on}.
!;
    }
    my $aligned_filename = qq"${bt_dir}/${job_basename}-${bt_type}_aligned_${species}.fastq";
    my $unaligned_filename = qq"${bt_dir}/${job_basename}-${bt_type}_unaligned_${species}.fastq";
    my $sam_filename = qq"${bt_dir}/${job_basename}-${bt_type}.sam";
    my $jstring = qq!mkdir -p ${bt_dir} && sleep ${sleep_time} && bowtie \\
  ${bt_reflib} \\
  ${bt_args} \\
  -p ${cpus} \\
  ${bowtie_input_flag} <(less ${bt_input}) \\
  --un ${unaligned_filename} \\
  --al ${aligned_filename} \\
  -S ${sam_filename} \\
  2>${error_file} \\
  1>${bt_dir}/${job_basename}-${bt_type}.out
!;

    my $bt_job = $class->Submit(
        aligned => $aligned_filename,
        comment => $comment,
        input => $bt_input,
        depends => $bt_depends_on,
        jname => $jname,
        job_output => $sam_filename,
        jprefix => "10",
        jstring => $jstring,
        postscript => $options->{postscript},
        prescript => $options->{prescript},
        queue => 'workstation',
        unaligned => $unaligned_filename,
    );
    $bt_jobs{bowtie} = $bt_job;

    my $un_comp = Bio::Adventure::Compress::Recompress(
        $class,
        comment => qq"## Compressing the sequences which failed to align against ${bt_reflib} using options ${bt_args}.\n",
        depends => $bt_job->{job_id},
        jname => "xzun",
        jprefix => "11",
        xz_input => "${bt_dir}/${job_basename}-${bt_type}_unaligned_${species}.fastq",
    );
    $bt_jobs{unaligned_compression} = $un_comp;

    my $al_comp = Bio::Adventure::Compress::Recompress(
        $class,
        comment => qq"## Compressing the sequences which successfully aligned against ${bt_reflib} using options ${bt_args}.",
        depends => $bt_job->{job_id},
        jname => "xzal",
        jprefix => "11",
        xz_input => "${bt_dir}/${job_basename}-${bt_type}_aligned_${species}.fastq",
    );
    $bt_jobs{aligned_compression} = $al_comp;

    ## BT1_Stats also reads the trimomatic output, which perhaps it should not.
    my $trim_output_file = qq"outputs/trimomatic_stats.csv";

    my $sam_job = Bio::Adventure::Convert::Samtools(
        $class,
        input => $sam_filename,
        depends => $bt_job->{job_id},
        jname => "s2b_${jname}",
        jprefix => "13",
    );
    $bt_jobs{samtools} = $sam_job;
    $options = $class->Set_Vars(job_output => $sam_job->{job_output});
    my $htmulti;
    if ($count) {
        if ($libtype eq 'rRNA') {
            $htmulti = Bio::Adventure::RNASeq_Count::HTSeq(
                $class,
                htseq_id => $options->{htseq_id},
                htseq_input => $sam_job->{job_output},
                htseq_type => $options->{htseq_type},
                depends => $sam_job->{job_id},
                jname => "ht_${jname}",
                jprefix => '14',
                libtype => $libtype,
                queue => 'workstation',
                suffix => $bt_type,
            );
        } else {
            $htmulti = Bio::Adventure::RNASeq_Count::HT_Multi(
                $class,
                htseq_id => $options->{htseq_id},
                htseq_input => $sam_job->{job_output},
                htseq_type => $options->{htseq_type},
                depends => $sam_job->{job_id},
                jname => "ht_${jname}",
                jprefix => '14',
                libtype => $libtype,
                queue => 'workstation',
                suffix => $bt_type,
            );
            $bt_jobs{htseq} = $htmulti;
        }
    }  ## End if ($count)

    my $stats = Bio::Adventure::RNASeq_Map::BT1_Stats(
        $class, %args,
        bt_input => $error_file,
        bt_type => $bt_type,
        count_table => qq"${job_basename}-${bt_type}.count.xz",
        depends => $bt_job->{job_id},
        jname => "${jname}_stats",
        jprefix => "12",
        trim_input => ${trim_output_file},
    );
    $bt_jobs{stats} = $stats;

    return(\%bt_jobs);
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
        required => ['species', 'input', 'htseq_type'],
        do_htseq => 1,
    );
    my $ready = $class->Check_Input(
        files => $options->{input},
    );
    my $sleep_time = 3;
    my %bt_jobs = ();
    my $libtype = 'genome';
    my $bt_depends_on = "";
    $bt_depends_on = $options->{depends} if ($options->{depends});
    my $bt2_args = $options->{bt2_args};

    my $prefix_name = qq"bt2";
    my $bt2_name = qq"${prefix_name}_$options->{species}";
    my $suffix_name = $prefix_name;
    if ($options->{jname}) {
        $bt2_name .= qq"_$options->{jname}";
        $suffix_name .= qq"_$options->{jname}";
    }

    my $job_basename = $options->{job_basename};
    my $bt_dir = qq"outputs/bowtie2_$options->{species}";
    if ($args{bt_dir}) {
        $bt_dir = $args{bt_dir};
    }
    my $bt_input = $options->{input};
    ##if ($bt_input =~ /\.gz$|\.bz2$|\.xz$/ ) {
    ##    print "The input needs to be uncompressed, doing that now.\n" if ($options->{debug});
    ##    my $uncomp = Bio::Adventure::Compress::Uncompress(
    ##        $class,
    ##        input => $bt_input,
    ##        depends => $bt_depends_on,
    ##    );
    ##    $bt_input =~ s/\:|\;|\,|\s+/ /g;
    ##    $bt_input =~ s/\.gz|\.bz|\.xz//g;
    ##    ##$bt_input = $bt_inputbasename($bt_input, ('.gz', '.bz2', '.xz'));
    ##    $options = $class->Set_Vars(input => $bt_input);
    ##    $bt_depends_on = $uncomp->{job_id};
    ##}

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
    if (!-r $bt_reftest) {
        print "Hey! The Indexes do not appear to exist, check this out: ${bt_reftest}\n";
        sleep(20);
        my $index_job = Bio::Adventure::RNASeq_Map::BT2_Index(
            $class,
            depends => $bt_depends_on,
            libtype => $libtype,
        );
        $bt_jobs{index} = $index_job;
        $bt_depends_on = $index_job->{job_id};
    }
    my $bowtie_input_flag = "-q "; ## fastq by default
    $bowtie_input_flag = "-f " if (${bt_input} =~ /\.fasta$/);

    my $cpus = $options->{cpus};
    my $error_file = qq"${bt_dir}/${job_basename}.err";
    my $comment = qq!## This is a bowtie2 alignment of ${bt_input} against
## ${bt_reflib} using arguments: ${bt2_args}.
## This jobs depended on: ${bt_depends_on}
!;
    my $aligned_filename = qq"${bt_dir}/${job_basename}_aligned_$options->{species}.fastq";
    my $unaligned_filename = qq"${bt_dir}/${job_basename}_unaligned_$options->{species}.fastq";
    my $sam_filename = qq"${bt_dir}/${job_basename}.sam";
    my $jstring = qq!if [[ -e "/scratch0" ]]; then
  scratchdir=\$(mktemp -d "/scratch0/\${USER}.XXXX")
  echo "Working in: \${scratchdir} on \$(hostname)."
  cd "\${scratchdir}" || exit
fi

mkdir -p ${bt_dir} && \\
  sleep ${sleep_time} && \\
  bowtie2 -x ${bt_reflib} ${bt2_args} \\
    -p ${cpus} \\
    ${bowtie_input_flag} ${bt_input} \\
    --un ${unaligned_filename} \\
    --al ${aligned_filename} \\
    -S ${sam_filename} \\
    2>${error_file} \\
    1>${bt_dir}/${job_basename}.out

if [[ -e "/scratch0" ]]; then
  cd $options->{basedir} && \\
    rsync -a "\${scratchdir}/${bt_dir}/" ${bt_dir} && \\
    rm -r "\${scratchdir}"
fi
!;

    my $bt2_job = $class->Submit(
        aligned => $aligned_filename,
        comment => $comment,
        input => $bt_input,
        jname => $bt2_name,
        depends => $bt_depends_on,
        jstring => $jstring,
        jprefix => "15",
        mem => 24,
        output => $sam_filename,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        unaligned => $unaligned_filename,
    );
    $bt_jobs{bowtie} = $bt2_job;

    my $un_comp = Bio::Adventure::Compress::Recompress(
        $class,
        comment => qq"## Compressing the sequences which failed to align against ${bt_reflib} using options ${bt2_args}\n",
        xz_input => "${bt_dir}/${job_basename}_unaligned_$options->{species}.fastq",
        depends => $bt2_job->{job_id},
        jname => "xzun_${suffix_name}",
        jprefix => "16",
        );
    $bt_jobs{unaligned_compression} = $un_comp;

    my $al_comp = Bio::Adventure::Compress::Recompress(
        $class,
        comment => qq"## Compressing the sequences which successfully aligned against ${bt_reflib} using options ${bt2_args}",
        xz_input => "${bt_dir}/${job_basename}_aligned_$options->{species}.fastq",
        jname => "xzal_${suffix_name}",
        jprefix => "17",
        depends => $bt2_job->{job_id},
    );
    $bt_jobs{aligned_compression} = $al_comp;

    ## BT1_Stats also reads the trimomatic output, which perhaps it should not.
    ## my $trim_output_file = qq"outputs/${job_basename}-trimomatic.out";
    my $stats = Bio::Adventure::RNASeq_Map::BT2_Stats(
        $class,
        bt_input => $error_file,
        count_table => qq"${job_basename}.count.xz",
        depends => $bt2_job->{job_id},
        jname => "bt2st_${suffix_name}",
        jprefix => "18",
        ## trim_input => ${trim_output_file},
    );
    my $sam_job = Bio::Adventure::Convert::Samtools(
        $class,
        input => $sam_filename,
        depends => $bt2_job->{job_id},
        jname => "s2b_${suffix_name}",
        jprefix => "19",
    );
    $bt_jobs{samtools} = $sam_job;
    my $htseq_input = $sam_job->{job_output};
    my $htmulti;
    if ($options->{do_htseq}) {
        if ($libtype eq 'rRNA') {
            $htmulti = Bio::Adventure::RNASeq_Count::HTSeq(
                $class,
                htseq_input => $sam_job->{job_output},
                depends => $sam_job->{job_id},
                jname => $suffix_name,
                jprefix => "20",
                libtype => $libtype,
            );
        } else {
            $htmulti = Bio::Adventure::RNASeq_Count::HT_Multi(
                $class,
                htseq_input => $sam_job->{job_output},
                depends => $sam_job->{job_id},
                jname => $suffix_name,
                jprefix => "21",
                libtype => $libtype,
            );
            $bt_jobs{htseq} = $htmulti;
        }
    }
    return(\%bt_jobs);
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
        required => ["species", "input", "htseq_type"],
    );
    my $bt_input = $options->{input};
    my $job_basename = $options->{job_basename};
    my $species = $options->{species};
    my $depends_on = $options->{depends};
    my %bt_types = %{$options->{bt_args}};
    my @jobs = ();
    foreach my $type (keys %bt_types) {
        my $jname = qq"bt${type}_${species}";
        my $job = Bio::Adventure::RNASeq_Map::Bowtie(
            $class,
            input => $bt_input,
            bt_type => $type,
            depends => $depends_on,
            jname => $jname,
            prescript => $args{prescript},
            postscript => $args{postscript},
        );
        push(@jobs, $job);
    }
    return(\@jobs);
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
        required => ["species"],
    );
    my $job_basename = $options->{job_basename};
    $job_basename = qq"rRNA_${job_basename}";
    my $exclude = 0;
    $exclude = $options->{exclude} if ($options->{exclude});
    my $species = $options->{species};
    my $depends_on = $options->{depends};
    my ${bt_dir} = qq"outputs/bowtie_$options->{species}";
    if ($exclude) {
        my $in = $options->{input};
        $options = $class->Set_Vars(
            postscript => qq"mv $in includerrna_${in} && mv ${bt_dir}/${job_basename}-rRNA_unaligned_${species}.fastq $in",
        );
    }
    my $job = Bio::Adventure::RNASeq_Map::Bowtie(
        $class,
        depends => $depends_on,
        jname => qq"btrrna",
        libtype => 'rRNA',
        prescript => $args{prescript},
        postscript => $args{postscript},
    );
    ## Return the basename back to normal so that future tasks don't
    ## get confuseled.
    return($job);
}

=head2 C<BT1_Index>

Create a bowtie1 index using $hpgl->{species}.fasta and leaves it in the
indexes/ directory.

=cut
sub BT1_Index {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args, required => ["species"],
                                   depends => "");
    my $job_basename = $options->{job_basename};
    my $jstring = qq!bowtie-build $options->{libdir}/$options->{libtype}/$options->{species}.fasta \\
  $options->{libdir}/$options->{libtype}/indexes/$options->{species}
!;
    my $comment = qq!## Generating bowtie1 indexes for species: $options->{species} in $options->{libdir}/$options->{libtype}/indexes!;
    my $bt1_index = $class->Submit(
        comment => $comment,
        jname => "bt1idx",
        depends => $options->{depends},
        jstring => $jstring,
        jprefix => "10",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
    );
    return($bt1_index);
}

=head2 C<BT2_Index>

Create a bowtie2 index using ${species}.fasta and leaves it in the indexes/
directory.

=cut
sub BT2_Index {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args, required => ["species"]);
    my $job_basename = $options->{job_basename};
    my $dep = "";
    $dep = $options->{depends};
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
        depends => $dep,
        jname => "bt2idx",
        jprefix => "15",
        jstring => $jstring,
        prescript => $args{prescript},
        postscript => $args{postscript},
    );
    return($indexer);
}

=head2 C<BWA>

Perform a bwa alignment using both the sam(s|p)e and aln algorithms.  It then
converts the output (when appropriate) to sorted/indexed bam and passes them to
htseq.

=cut
sub BWA {
    my ($class, %args) = @_;
    my $check = which('bwa');
    die("Could not find bwa in your PATH.") unless($check);

    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        species => 'lmajor',
        libtype => 'genome',
    );

    my $sleep_time = 3;
    my %bwa_jobs = ();
    my $bwa_input = $options->{input};
    my $bwa_depends_on;
    $bwa_depends_on = $options->{depends} if ($options->{depends});
    my $job_basename = $options->{job_basename};

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
    if (!-r $bwa_reftest) {
        my $index_job = Bio::Adventure::RNASeq_Map::BWA_Index(
            $class,
            depends => $bwa_depends_on,
            libtype => $libtype,
        );
        $bwa_jobs{index} = $index_job;
        $bwa_depends_on = $index_job->{job_id};
    }

    ## Make a depends string containing the uncompress, indexer, or both
    if (defined($index_jobid)) {
        $bwa_depends_on = $index_jobid;
    }

    ## Include logic to deal with paired end/single end reads.
    my $test_file = "";
    my @pair_listing = ();
    my ($forward_reads, $reverse_reads);
    if ($bwa_input =~ /\s+/) {
        ($forward_reads, $reverse_reads) = split(/\s+/, $bwa_input);
    } else {
        $forward_reads = $bwa_input;
    }

    my $aln_sam = qq"${bwa_dir}/${job_basename}_aln.sam";
    my $mem_sam = qq"${bwa_dir}/${job_basename}_mem.sam";
    my $aln_args = qq"";
    my $mem_args = qq"";
    my $sam_args = qq"";
    my $jstring = qq!mkdir -p ${bwa_dir}
bwa mem ${mem_args} -a ${bwa_reflib} <(less ${bwa_input}) \\
  2>${bwa_dir}/bwa.err 1>${mem_sam}
!;
    my $reporter_string = qq"bwa samse ${sam_args} ${bwa_reflib} \\
  ${bwa_dir}/${job_basename}_aln-forward.sai <(less ${bwa_input}) \\
  2>${bwa_dir}/${job_basename}.samerr \\
  1>${aln_sam}";
    my $aln_string = qq"bwa aln ${aln_args} ${bwa_reflib} \\
  ${forward_reads} \\
  2>${bwa_dir}/${job_basename}_aln-forward.err \\
  1>${bwa_dir}/${job_basename}_aln-forward.sai";
    if (defined($reverse_reads)) {
        $aln_string = qq"${aln_string}
bwa aln ${aln_args} ${bwa_reflib} \\
  <(less ${reverse_reads}) \\
  2>${bwa_dir}/${job_basename}_aln-reverse.err \\
  1>${bwa_dir}/${job_basename}_aln-reverse.sai";
        $reporter_string = qq"bwa ${sam_args} sampe ${bwa_reflib} \\
  ${bwa_dir}/${job_basename}_aln-forward.sai ${bwa_dir}/${job_basename}_aln-reverse.sai \\
  <(less ${forward_reads}) <(less ${reverse_reads}) \\
  2>${bwa_dir}/${job_basename}.samerr \\
  1>${aln_sam}";
    }

    my $comment = qq!## This is a BWA alignment of ${bwa_input} against
## ${bwa_reflib}.
## It will perform a separate bwa mem run and bwa aln run.!;

    $jstring = qq!${jstring}
${aln_string}
${reporter_string}
!;
    my $bwa_job = $class->Submit(
        comment => $comment,
        input => $bwa_input,
        depends => $bwa_depends_on,
        jname => "bwa_$options->{species}",
        job_output => qq"${bwa_dir}/${job_basename}_mem.sam",
        jprefix => "20",
        jstring => $jstring,
        postscript => $options->{postscript},
        prescript => $options->{prescript},
        queue => 'workstation',
    );
    $bwa_jobs{bwa} = $bwa_job;

    my $mem_sam_job = Bio::Adventure::Convert::Samtools(
        $class,
        input => $bwa_job->{job_output},
        depends => $bwa_job->{job_id},
        jname => "s2b_mem",
        jprefix => "21",
    );
    $bwa_jobs{samtools_mem} = $mem_sam_job;

    my $aln_sam_job = Bio::Adventure::Convert::Samtools(
        $class,
        input => ${aln_sam},
        depends => $mem_sam_job->{job_id},
        jname => "s2b_aln",
        jprefix => "21",
    );
    $bwa_jobs{samtools_mem} = $aln_sam_job;

    my $mem_htmulti = Bio::Adventure::RNASeq_Count::HT_Multi(
        $class,
        htseq_id => $options->{htseq_id},
        htseq_input => $mem_sam_job->{job_output},
        htseq_type => $options->{htseq_type},
        depends => $mem_sam_job->{job_id},
        jname => "htmem_${jname}",
        jprefix => '22',
    );
    $bwa_jobs{htseq_mem} = $mem_htmulti;

    my $aln_htmulti = Bio::Adventure::RNASeq_Count::HT_Multi(
        $class,
        htseq_id => $options->{htseq_id},
        htseq_input => $aln_sam_job->{job_output},
        htseq_type => $options->{htseq_type},
        depends => $aln_sam_job->{job_id},
        jname => "htaln_${jname}",
        jprefix => '22',
    );
    $bwa_jobs{htseq_aln} = $aln_htmulti;

    my $bwa_stats = Bio::Adventure::RNASeq_Map::BWA_Stats(
        $class,
        depends => $mem_sam_job->{job_id},
        jname => 'bwastats',
        jprefix => "25",
        aln_output => $aln_sam_job->{job_output},
        mem_output => $mem_sam_job->{job_output},
    );
    $bwa_jobs{stats} = $bwa_stats;

    return(\%bwa_jobs);
}

=head2 C<BWA_Stats>

Collect some alignment statistics from bwa.

=cut
sub BWA_Stats {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $job_basename = $options->{job_basename};
    my $aln_input = $options->{aln_output};
    $aln_input = qq"${aln_input}.stats";
    my $mem_input = $options->{mem_output};
    $mem_input = qq"${mem_input}.stats";
    my $stat_output = qq"outputs/bwa_stats.csv";

    my $depends = $options->{depends};
    my $jname = "bwa_stats";
    $jname = $options->{jname} if ($options->{jname});
    my $jobid = qq"${job_basename}_stats";
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
stat_string=\$(printf "${job_basename},%s,%s,%s,%s,%s,${count_table}" "\${original_reads}" "\${reads}" "\${aln_aligned}" "\${mem_aligned}" "\$rpm")
echo "\${stat_string}" >> ${stat_output}!;
    my $stats = $class->Submit(
        comment => $comment,
        input => $aln_input,
        depends => $depends,
        jname => $jname,
        jprefix => $options->{jprefix},
        jstring => $jstring,
        cpus => 1,
        mem => 1,
        queue => "throughput",
        walltime => "00:10:00",
    );
    return($stats);
}

=head2 C<BWA_Index>

Create bwa indexes.

=cut
sub BWA_Index {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $job_basename = $options->{job_basename};
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
        depends => $options->{depends},
        jname => "bwaidx",
        jprefix => "26",
        jstring => $jstring,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
    );
    return($bwa_index);
}

=head2 C<BT2_Stats>

Collects alignment statistics from bowtie 2.

=cut
sub BT2_Stats {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $bt_input = $options->{bt_input};
    my $job_basename = $options->{job_basename};
    my $bt_type = $options->{bt_type};
    my $jname = "bt2_stats";
    $jname = $options->{jname} if ($options->{jname});
    my $jobid = qq"${job_basename}_stats";
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
stat_string=\$(printf "${job_basename},${bt_type},%s,%s,%s,%s,%s" "\${original_reads}" "\${one_align}" "\${failed}" "\${sampled}" "\${rpm}")
echo "\$stat_string" >> ${output}!;
    my $stats = $class->Submit(
        comment => $comment,
        input => $bt_input,
        jname => $jname,
        depends => $options->{depends},
        jprefix => $args{jprefix},
        jstring => $jstring,
        cpus => 1,
        mem => 1,
        queue => "throughput",
        walltime => "00:10:00",
    );
    return($stats);
}

=head2 C<Hisat2>

Invoke hisat2

=cut
sub Hisat2 {
    my ($class, %args) = @_;
    my $check = which('hisat2-build');
    die("Could not find hisat2 in your PATH.") unless($check);

    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input', 'htseq_type'],
        do_htseq => 1,
    );
    my $ready = $class->Check_Input(
        files => $options->{input},
    );
    my $sleep_time = 3;
    my %ht_jobs = ();
    my $libtype = 'genome';
    my $ht_depends_on = "";
    $ht_depends_on = $options->{depends} if ($options->{depends});
    my $ht2_args = '';
    $ht2_args = $options->{ht2_args} if ($options->{ht2_args});

    my $prefix_name = qq"hisat2";
    my $ht2_name = qq"${prefix_name}_$options->{species}";
    my $suffix_name = $prefix_name;
    if ($options->{jname}) {
        $ht2_name .= qq"_$options->{jname}";
        $suffix_name .= qq"_$options->{jname}";
    }

    my $job_basename = $options->{job_basename};
    my $ht_dir = qq"outputs/hisat2_$options->{species}";
    if ($args{ht_dir}) {
        $ht_dir = $args{ht_dir};
    }
    my $ht_input = $options->{input};

    my $test_file = "";
    if ($ht_input =~ /\:|\;|\,|\s+/) {
        my @pair_listing = split(/\:|\;|\,|\s+/, $ht_input);
        $pair_listing[0] = File::Spec->rel2abs($pair_listing[0]);
        $pair_listing[1] = File::Spec->rel2abs($pair_listing[1]);
        $ht_input = qq" -1 <(less $pair_listing[0]) -2 <(less $pair_listing[1]) ";
        $test_file = $pair_listing[0];
    } else {
        $test_file = File::Spec->rel2abs($ht_input);
        $ht_input = qq" <(less ${ht_input}) ";
    }

    ## Check that the indexes exist
    my $ht_reflib = "$options->{libdir}/$options->{libtype}/indexes/$options->{species}";
    my $ht_reftest = qq"${ht_reflib}.1.ht2";
    my $ht_reftestl = qq"${ht_reflib}.1.ht2l";
    if (!-r $ht_reftest && !-r $ht_reftestl) {
        print "Hey! The Indexes do not appear to exist, check this out: ${ht_reftest}\n";
        sleep(20);
        my $index_job = Bio::Adventure::RNASeq_Map::HT2_Index(
            $class,
            depends => $ht_depends_on,
            libtype => $libtype,
        );
        $ht_jobs{index} = $index_job;
        $ht_depends_on = $index_job->{job_id};
    }
    my $hisat_input_flag = "-q "; ## fastq by default
    $hisat_input_flag = "-f " if (${ht_input} =~ /\.fasta$/);

    my $cpus = $options->{cpus};
    my $error_file = qq"${ht_dir}/${job_basename}.err";
    my $comment = qq!## This is a hisat2 alignment of ${ht_input} against
## ${ht_reflib} using arguments: ${ht2_args}.
## This jobs depended on: ${ht_depends_on}
!;
    my $aligned_filename = qq"${ht_dir}/${job_basename}_aligned_$options->{species}.fastq";
    my $unaligned_filename = qq"${ht_dir}/${job_basename}_unaligned_$options->{species}.fastq";
    my $sam_filename = qq"${ht_dir}/${job_basename}.sam";
    my $jstring = qq!if [[ -e "/scratch0" ]]; then
  scratchdir=\$(mktemp -d "/scratch0/\${USER}.XXXX")
  echo "Working in: \${scratchdir} on \$(hostname)."
  cd "\${scratchdir}" || exit
fi

mkdir -p ${ht_dir} && \\
  sleep ${sleep_time} && \\
  hisat2 -x ${ht_reflib} ${ht2_args} \\
    -p ${cpus} \\
    ${hisat_input_flag} ${ht_input} \\
    --phred$options->{phred} \\
    --un ${unaligned_filename} \\
    --al ${aligned_filename} \\
    -S ${sam_filename} \\
    2>${error_file} \\
    1>${ht_dir}/${job_basename}.out

if [[ -e "/scratch0" ]]; then
  cd $options->{basedir} && \\
    rsync -a "\${scratchdir}/${ht_dir}/" ${ht_dir} && \\
    rm -r "\${scratchdir}"
fi
!;

    my $ht2_job = $class->Submit(
        aligned => $aligned_filename,
        comment => $comment,
        input => $ht_input,
        jname => $ht2_name,
        depends => $ht_depends_on,
        jstring => $jstring,
        jprefix => "15",
        mem => 32,
        walltime => '24:00:00',
        output => $sam_filename,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        unaligned => $unaligned_filename,
    );
    $ht_jobs{hisat} = $ht2_job;

    ## HT1_Stats also reads the trimomatic output, which perhaps it should not.
    ## my $trim_output_file = qq"outputs/${job_basename}-trimomatic.out";
    my $stats = Bio::Adventure::RNASeq_Map::HT2_Stats(
        $class,
        ht_input => $error_file,
        count_table => qq"${job_basename}.count.xz",
        depends => $ht2_job->{job_id},
        jname => "ht2st_${suffix_name}",
        jprefix => "18",
        ## trim_input => ${trim_output_file},
    );
    my $sam_job = Bio::Adventure::Convert::Samtools(
        $class,
        input => $sam_filename,
        depends => $ht2_job->{job_id},
        jname => "s2b_${suffix_name}",
        jprefix => "19",
    );
    $ht_jobs{samtools} = $sam_job;
    my $htseq_input = $sam_job->{job_output};
    my $htmulti;
    if ($options->{do_htseq}) {
        if ($libtype eq 'rRNA') {
            $htmulti = Bio::Adventure::RNASeq_Count::HTSeq(
                $class,
                htseq_input => $sam_job->{job_output},
                depends => $sam_job->{job_id},
                jname => $suffix_name,
                jprefix => "20",
                libtype => $libtype,
            );
        } else {
            $htmulti = Bio::Adventure::RNASeq_Count::HT_Multi(
                $class,
                htseq_input => $sam_job->{job_output},
                depends => $sam_job->{job_id},
                jname => $suffix_name,
                jprefix => "21",
                libtype => $libtype,
            );
            $ht_jobs{htseq} = $htmulti;
        }
    }
return(\%ht_jobs);
}

=head2 C<HT2_Index>

Create a hisat2 index using ${species}.fasta and leaves it in the indexes/
directory.

=cut
sub HT2_Index {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args, required => ["species"]);
    my $job_basename = $options->{job_basename};
    my $dep = "";
    $dep = $options->{depends};
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
        depends => $dep,
        jname => "ht2idx",
        jprefix => "15",
        jstring => $jstring,
        prescript => $args{prescript},
        postscript => $args{postscript},
    );
    return($indexer);
}

=head2 C<HT2_Stats>

Collect alignment statistics from hisat 2.

=cut
sub HT2_Stats {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $ht_input = $options->{ht_input};
    my $job_basename = $options->{job_basename};
    my $jname = "ht2_stats";
    $jname = $options->{jname} if ($options->{jname});
    my $jobid = qq"${job_basename}_stats";
    my $count_table = "";
    $count_table = $options->{count_table} if ($options->{count_table});
    my $comment = qq!## This is a stupidly simple job to collect alignment statistics.!;
    my $output = "outputs/hisat2_stats.csv";
    my $jstring = qq!
if [ \! -e "${output}" ]; then
    echo "original reads, single hits, failed reads, multi-hits, rpm" > ${output}
fi
original_reads_tmp=\$(grep " reads; of these" "${ht_input}" 2>/dev/null | awk '{print \$1}' | sed 's/ //g')
original_reads=\${original_reads_tmp:-0}
one_align_tmp=\$(grep " aligned exactly 1 time" "${ht_input}" | awk '{print \$1}' | sed 's/ .*//g')
one_align=\${one_align_tmp:-0}
failed_tmp=\$(grep " aligned 0 times" "${ht_input}" | awk '{print \$1}' | sed 's/ .*//g')
failed=\${failed_tmp:-0}
sampled_tmp=\$(grep " aligned >1 times" "${ht_input}" | awk '{print \$1}' | sed 's/ .*//g')
sampled=\${sampled_tmp:-0}
rpm_tmp=\$(perl -e "printf(1000000 / \$(( \${one_align} + \${sampled} )) ) " 2>/dev/null)
rpm=\${rpm_tmp:-0}
stat_string=\$(printf "${job_basename},%s,%s,%s,%s,%s" "\${original_reads}" "\${one_align}" "\${failed}" "\${sampled}" "\${rpm}")
echo "\$stat_string" >> ${output}!;
    my $stats = $class->Submit(
        comment => $comment,
        input => $ht_input,
        jname => $jname,
        depends => $options->{depends},
        jprefix => $args{jprefix},
        jstring => $jstring,
        cpus => 1,
        mem => 1,
        queue => "throughput",
        walltime => "00:10:00",
    );
    return($stats);
}

=head2 C<Kallisto_Index

Use kallisto and an annotated_CDS fasta sequence library to create an index.

=cut
sub Kallisto_Index {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ["species", "genome"],
    );
    my $job_basename = $options->{job_basename};
    my $dep = "";
    $dep = $options->{depends};
    my $libtype = $options->{libtype};
    my $genome = File::Spec->rel2abs($options->{genome});
    unless (-r $genome) {
        die("The indexing operation for kallisto will fail because the $options->{species} genome does not exist.")
    }

    my $jstring = qq!
kallisto index -i $options->{libdir}/${libtype}/indexes/$options->{species}.idx ${genome}!;
    my $comment = qq!## Generating kallisto indexes for species: $options->{species} in $options->{libdir}/${libtype}/indexes!;
    my $jobid = $class->Submit(
        comment => $comment,
        depends => $dep,
        jstring => $jstring,
        jname => "kalidx",
        jprefix => $options->{jprefix},
        prescript => $options->{prescript},
        postscript => $options->{postscript},
    );
    return($jobid);
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

    ##my $ready = $class->Check_Input(
    ##    files => $options->{input},
    ##);

    my $sleep_time = 3;
    my %ka_jobs = ();
    my $ka_depends_on = '';
    $ka_depends_on = $options->{depends} if ($options->{depends});
    my $libtype = 'genome';
    $libtype = $options->{libtype} if ($options->{libtype});
    my $job_basename = $options->{job_basename};
    my $species = $options->{species};


    my $jname = qq"kall_${species}";
    ## $jname = $options->{jname} if ($options->{jname});
    my $ka_args = qq"";

    my $ka_input = $options->{input};
    ##if ($ka_input =~ /\.gz|\.bz2$|\.xz$/ ) {
    ##    print "The input needs to be uncompressed, doing that now.\n" if ($options->{debug});
    ##    my $ka_input_cleaned = $ka_input;
    ##    $ka_input_cleaned =~ s/\:|\;|\,|\s+/ /g;
    ##    $ka_input_cleaned =~ s/\.gz|\.bz|\.xz//g;
    ##    my @ka_array = split(/ /, $ka_input_cleaned);
    ##    my $ka_short = $ka_array[0];
    ##    my $uncomp = Bio::Adventure::Compress::Uncompress(
    ##        $class,
    ##        input => $ka_input,
    ##        jname => "uncomp_${ka_short}",
    ##        depends => $ka_depends_on,
    ##    );
    ##    $ka_input = $ka_input_cleaned;
    ##    $options = $class->Set_Vars(input => $ka_input);
    ##    $ka_depends_on = $uncomp->{job_id};
    ##}

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
        my $index_job = Bio::Adventure::RNASeq_Map::Kallisto_Index(
            $class,
            depends => $ka_depends_on,
            libtype => $libtype,
        );
        $ka_jobs{index} = $index_job;
        $ka_depends_on = $index_job->{job_id};
    }

    my $outdir = qq"outputs/kallisto_${species}";
    my $error_file = qq"${outdir}/kallisto_${species}.stderr";
    my $output_sam = qq"${outdir}/kallisto_${species}.sam";
    my $output_bam = qq"${outdir}/kallisto_${species}.bam";
    my $output_stats = qq"${outdir}/kallisto_${species}.stats";
    my $sorted_bam = qq"${outdir}/kallisto_${species}-sorted";
    my $comment = qq!## This is a kallisto pseudoalignment of ${ka_input} against
## ${ka_reflib}.
## This jobs depended on: ${ka_depends_on}
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
    my $ka_job = $class->Submit(
        comment => $comment,
        input => $ka_input,
        depends => $ka_depends_on,
        jname => qq"${jname}",
        jprefix => "30",
        jstring => $jstring,
        mem => 30,
        prescript => $args{prescript},
        postscript => $args{postscript},
        queue => "workstation",
    );
    $ka_jobs{kallisto} = $ka_job;

    return(\%ka_jobs);
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
    );
    my $job_basename = $options->{job_basename};
    my $dep = "";
    $dep = $options->{depends};
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
        depends => $dep,
        jstring => $jstring,
        jname => "salidx",
        jprefix => $options->{jprefix},
        prescript => $options->{prescript},
        postscript => $options->{postscript},
    );
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
    );

    my $ready = $class->Check_Input(
        files => $options->{input},
    );

    my $sleep_time = 3;
    my %sa_jobs = ();
    my $sa_depends_on = "";
    $sa_depends_on = $options->{depends} if ($options->{depends});
    my $libtype = 'genome';
    $libtype = $options->{libtype} if ($options->{libtype});
    my $job_basename = $options->{job_basename};
    my $species = $options->{species};


    my $jname = qq"sal_${species}";
    ## $jname = $options->{jname} if ($options->{jname});
    my $sa_args = qq"";

    my $sa_input = $options->{input};
    ##if ($sa_input =~ /\.gz|\.bz2$|\.xz$/ ) {
    ##    print "The input needs to be uncompressed, doing that now.\n" if ($options->{debug});
    ##    my $sa_input_cleaned = $sa_input;
    ##    $sa_input_cleaned =~ s/\:|\;|\,|\s+/ /g;
    ##    $sa_input_cleaned =~ s/\.gz|\.bz|\.xz//g;
    ##    my @sa_array = split(/ /, $sa_input_cleaned);
    ##    my $sa_short = $sa_array[0];
    ##    my $uncomp = Bio::Adventure::Compress::Uncompress(
    ##        $class,
    ##        input => $sa_input,
    ##        jname => "uncomp_${sa_short}",
    ##        depends => $sa_depends_on,
    ##    );
    ##    $sa_input = $sa_input_cleaned;
    ##    $options = $class->Set_Vars(input => $sa_input);
    ##    $sa_depends_on = $uncomp->{job_id};
    ##}

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
    if (!-r $sa_reflib) {
        my $index_job = Bio::Adventure::RNASeq_Map::Salmon_Index(
            $class,
            depends => $sa_depends_on,
            libtype => $libtype,
        );
        $sa_jobs{index} = $index_job;
        $sa_depends_on = $index_job->{job_id};
    }

    my $outdir = qq"outputs/salmon_${species}";
    my $error_file = qq"${outdir}/salmon_${species}.stderr";
    my $comment = qq!## This is a salmon pseudoalignment of ${sa_input} against
## ${sa_reflib}.
## This jobs depended on: ${sa_depends_on}
!;
    my $jstring = qq!mkdir -p ${outdir} && sleep ${sleep_time} && \\
salmon quant -i ${sa_reflib} \\
  -l A --gcBias --validateMappings  \\
  ${sa_args} \\
  -o ${outdir} \\
  2>${outdir}/salmon.err 1>${outdir}/salmon.out
!;

    my $sa_job = $class->Submit(
        comment => $comment,
        input => $sa_input,
        depends => $sa_depends_on,
        jname => qq"${jname}",
        jprefix => "30",
        jstring => $jstring,
        mem => 48,
        prescript => $args{prescript},
        postscript => $args{postscript},
        queue => "workstation",
    );

    my $stats = Bio::Adventure::RNASeq_Map::Salmon_Stats(
        $class,
        input => qq"${outdir}/lib_format_counts.json",
        depends => $sa_job->{job_id},
        jname => "sastats_$options->{species}",
        jprefix => "33",
    );

    $sa_jobs{salmon} = $sa_job;

    return(\%sa_jobs);
}

sub Salmon_Stats {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $depends = $options->{depends};
    my $jname = "stats";
    $jname = $options->{jname} if ($options->{jname});
    my $job_basename = $options->{job_basename};
    my $jobid = qq"${job_basename}_stats";
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
stat_string=\$(printf "${job_basename},$options->{species},%s,%s,%s,%s,%s" "\${reads}" "\${aligned}" "\${consistent}" "\${inconsistent}" "\${bias}")
echo "\$stat_string" >> "${output}"!;
    my $stats = $class->Submit(
        comment => $comment,
        cpus => 1,
        input => $options->{input},
        jname => $jname,
        depends => $depends,
        jprefix => $args{jprefix},
        jstring => $jstring,
        mem => 1,
        queue => "throughput",
        walltime => "00:10:00",
    );
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
        required => ["species","input"],
    );

    my $ready = $class->Check_Input(
        files => $options->{input},
    );

    my $sleep_time = 3;
    my %star_jobs = ();
    my $star_depends_on = "";
    $star_depends_on = $options->{depends} if ($options->{depends});
    my $libtype = 'genome';
    $libtype = $options->{libtype} if ($options->{libtype});
    my $job_basename = $options->{job_basename};
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
    if (!-r $star_reflib) {
        my $index_job = Bio::Adventure::RNASeq_Map::STAR_Index(
            $class,
            depends => $star_depends_on,
            libtype => $libtype,
        );
        $star_jobs{index} = $index_job;
        $star_depends_on = $index_job->{job_id};
    }

    my $outdir = qq"outputs/star_${species}";
    my $error_file = qq"${outdir}/star_${species}.stderr";
    my $comment = qq!## This is a star pseudoalignment of ${star_input} against
## ${star_reflib}.
## This jobs depended on: ${star_depends_on}
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
        depends => $star_depends_on,
        jname => qq"${jname}",
        jprefix => "33",
        jstring => $jstring,
        mem => 96,
        prescript => $args{prescript},
        postscript => $args{postscript},
        queue => "large",
    );
    $star_jobs{salmon} = $star_job;

    return(\%star_jobs);
}

=head2 C<STAR_Index>

Create indexes appropriate for STAR.

=cut
sub STAR_Index {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ["species",],
    );
    my $dep = $options->{depends};
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
        depends => $dep,
        jstring => $jstring,
        jname => "staridx",
        jprefix => $options->{jprefix},
        mem => 180,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => 'xlarge',
        walltime => '20-00:00:00',
    );
    return($jobid);
}

=item C<RSEM_Index

Use RSEM and an annotated_CDS fasta sequence library to create a transcript index.

=cut
sub RSEM_Index {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ["species",],
    );
    my $dep = $options->{depends};
    my $comment = qq"## RSEM Index creation.";
    my $jstring = qq!
rsem-prepare-reference --bowtie2 $options->{cds_fasta} $options->{index}
!;
    my $jobid = $class->Submit(
        comment => $comment,
        depends => $dep,
        jstring => $jstring,
        jname => "rsemidx",
        jprefix => $options->{jprefix},
        prescript => $options->{prescript},
        postscript => $options->{postscript},
    );
    return($jobid);
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
    );

    my $rsem_depends_on = "";
    $rsem_depends_on = $options->{depends} if ($options->{depends});
    my %rsem_jobs = ();
    my $job_basename = $class->Get_Job_Name();
    my $cds = qq"$options->{libdir}/$options->{libtype}/$options->{species}_cds_nt.fasta";
    my $idx = qq"$options->{libdir}/$options->{libtype}/indexes/rsem/$options->{species}";
    my $test_idx = qq"${idx}.transcripts.fa";
    unless (-r $test_idx) {
        print "Need to create the RSEM indexes, invoking RSEM_Index().\n";
        unless (-r $cds) {
            die("RSEM_Index requires a cds fasta file at: ${cds}, create it with a cyoa2 conversion.");
        }
        my $index_job = Bio::Adventure::RNASeq_Map::RSEM_Index(
            $class,
            cds_fasta => $cds,
            index => $idx,
            depends => $rsem_depends_on,
            jname => qq"rsidx_${job_basename}",
            libtype => $options->{libtype},
        );
        $rsem_jobs{index} = $index_job;
        $rsem_depends_on = $index_job->{job_id};
    }

    my $rsem_input = $options->{input};
    ##if ($rsem_input =~ /\.gz$|\.bz2$|\.xz$/ ) {
    ##    print "The input needs to be uncompressed, doing that now.\n" if ($options->{debug});
    ##    my $uncomp = Bio::Adventure::Compress::Uncompress(
    ##        $class,
    ##        input => $rsem_input,
    ##        depends => $rsem_depends_on,
    ##    );
    ##    $rsem_input =~ s/\:|\;|\,|\s+/ /g;
    ##    $rsem_input =~ s/\.gz|\.bz|\.xz//g;
    ##    $options = $class->Set_Vars(input => $rsem_input);
    ##    $rsem_depends_on = $uncomp->{job_id};
    ##}

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

    my $rsem_job = $class->Submit(
        comment => $rsem_comment,
        input => $rsem_input,
        jname => qq"rsem_${job_basename}",
        depends => $rsem_depends_on,
        jstring => $jstring,
        jprefix => "28",
        mem => 24,
        output => $output_file,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        walltime => '36:00:00',
    );
    $rsem_jobs{rsem} = $rsem_job;
    return(\%rsem_jobs);
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
    );
    my $depends = $options->{depends};
    my $tophat_cpus = 4;
    my $inputs = $options->{input};
    my @in = split(/:/, $inputs);
    $inputs =~ s/:/ /g;
    my $number_inputs = scalar(@in);
    if ($number_inputs > 1) {
        $inputs = qq" <(less $in[0]) <(less $in[1]) ";
    } else {
        $inputs = qq" <(less $in[0]) ";
    }
    my $tophat_args = ' -g 1 --microexon-search --b2-very-sensitive ';
    if ($options->{tophat_args}) {
        $tophat_args = $options->{tophat_args};
    }
    ## $tophat_args .= ' --no-mixed --no-discordant ' if (scalar(@in) > 1);
    ## $tophat_args .= ' ' if (scalar(@in) > 1);

    my $tophat_queue = $options->{queue};
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
    my $job_basename = $options->{job_basename};
    my $libtype = $options->{libtype};
    my $bt_reflib = "$options->{libdir}/${libtype}/indexes/$options->{species}";
    my $bt_reftest = qq"${bt_reflib}.1.bt2";
    my $index_job = undef;
    if (!-r $bt_reftest) {
        print "Did not find the index for $options->{species} at: ${bt_reflib}, indexing now.\n";
        $index_job = Bio::Adventure::RNASeq_Map::BT2_Index(
            $class,
            depends => $depends,
        );
        ## Use a colon separated append to make tophat depend on multiple jobs
        ## $depends .= qq":$index_job->{jobid}";
        ## Or just replace the dependency string with this job's and insert it into the stack
        $depends = $index_job->{job_id};
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
    if ($number_inputs > 1) {
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
        comment => ${comment},
        cpus => ${tophat_cpus},
        depends => ${depends},
        jname => ${jname},
        jprefix => "31",
        jstring => ${jstring},
        mem => ${tophat_mem},
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => ${tophat_queue},
        walltime => ${tophat_walltime},
    );

    ## Set the input for htseq
    my $accepted = "${tophat_dir}/accepted_hits.bam";
    $accepted = $options->{accepted_hits} if ($options->{accepted_hits});
    my $count_table = "accepted_hits.count";
    $count_table = $options->{count_table} if ($options->{count_table});
    my $htmulti = Bio::Adventure::RNASeq_Count::HT_Multi(
        $class,
        htseq_input => $accepted,
        htseq_id => $options->{htseq_id},
        htseq_type => $options->{htseq_type},
        depends => $tophat->{job_id},
        jname => qq"hts_$options->{species}",
        jprefix => "32",
    );
    $tophat->{htseq} = $htmulti;
    ## Perform a separate htseq run using only the successfully paired hits
    if ($number_inputs > 1) {
        my $ht_paired = Bio::Adventure::RNASeq_Count::HT_Multi(
            $class,
            htseq_input => qq"${tophat_dir}/accepted_paired.bam",
            htseq_id => $options->{htseq_id},
            htseq_type => $options->{htseq_type},
            depends => $tophat->{job_id},
            jname => qq"htsp_$options->{species}",
            jprefix => "32",
        );
    }
    ## Tophat_Stats also reads the trimomatic output, which perhaps it should not.
    my $unaccepted = $accepted;
    $unaccepted =~ s/accepted_hits/unmapped/g;
    my $input_read_info = $accepted;
    $input_read_info =~ s/accepted_hits\.bam/prep_reads\.info/g;
    my $stats = Bio::Adventure::RNASeq_Map::Tophat_Stats(
        $class,
        accepted_input => $accepted,
        count_table => qq"${count_table}.xz",
        depends => $tophat->{job_id},
        jname => "tpstats_$options->{species}",
        jprefix => "33",
        prep_input => $input_read_info,
        unaccepted_input => $unaccepted,
    );
    $tophat->{stats} = $stats;

    return($tophat);
}

=head2 C<BT1_Stats>

Collect some alignment statistics from bowtie1.

=cut
sub BT1_Stats {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $bt_input = $options->{bt_input};
    my $job_basename = $options->{job_basename};
    my $bt_type = "";
    $bt_type = $options->{bt_type} if ($options->{bt_type});
    my $depends = $options->{depends};
    my $jname = "stats";
    $jname = $options->{jname} if ($options->{jname});
    my $jobid = qq"${job_basename}_stats";
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
stat_string=\$(printf "${job_basename},${bt_type},%s,%s,%s,%s,%s,%s,${count_table}" "\${original_reads}" "\${reads}" "\${one_align}" "\${failed}" "\${sampled}" "\$rpm")
echo "\$stat_string" >> ${stat_output}!;
    my $stats = $class->Submit(
        comment => $comment,
        input => $bt_input,
        depends => $depends,
        jname => $jname,
        jprefix => $args{jprefix},
        jstring => $jstring,
        cpus => 1,
        mem => 1,
        queue => "throughput",
        walltime => "00:10:00",
    );
    return($stats);
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
    my $job_basename = $options->{job_basename};
    my $depends = $options->{depends};

    my $jname = "stats";
    $jname = $options->{jname} if ($options->{jname});
    my $jobid = qq"${job_basename}_stats";
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
stat_string=\$(printf "${job_basename},$options->{species},%s,%s,%s,%s,%s,${count_table}" "\${original_reads}" "\${reads}" "\${aligned}" "\${failed}" "\$rpm")
echo "\$stat_string" >> "${output}"!;
    my $stats = $class->Submit(
        comment => $comment,
        cpus => 1,
        input => $accepted_input,
        jname => $jname,
        depends => $depends,
        jprefix => $args{jprefix},
        jstring => $jstring,
        mem => 1,
        queue => "throughput",
        walltime => "00:10:00",
    );
    return($stats);
}

=back

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

    L<bowtie> L<bowtie2> L<tophat> L<bwa> L<kallisto> L<samtools> L<htseq>

=cut

1;
