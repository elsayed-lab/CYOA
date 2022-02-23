package Bio::Adventure::Map;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use File::Basename;
use File::Copy qw"cp";
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

 Run the OG (except bwa) short read aligner.
 10.1186/gb-2009-10-3-r25

 Perform a bowtie alignment.  Unless instructed otherwise, it will do so with 0
 mismatches and with multi-matches randomly placed 1 time among the
 possibilities (options -v 0 -M 1).  It will then convert the resulting sam alignment
 to a sorted-compressed-indexed bam, count it with htseq-count, compress the various
 output fastq files, and collect a few alignment statistics.

 This requires the arguments: 'input' and 'species'.  The input is likely a
 colon-separated pair of (compressed)fastq files.  The species will be used
 to look for bowtie indexes in ${libdir}/${libtype}/indexes/${species}.

 The argument bt_type(v0M1: no mismatches, randomly place multi-matches in 1 location)
 defines the mismatch and multimatch parameters; I
 pre-defined a few likely option sets for these rather important options.

 The count(1: e.g. yes) argument defines whether htseq-count will be performed.

 libtype(genome: as opposed to rRNA or contaminants etc) defines the type of
 index to search against.

 htseq_type(gene) defines the type of feature to count with htseq-count.  This is effectively
 the third column of a gff file.

 htseq_id(ID: ID is common for eukaryotic organisms, locus_tag is common for bacteria, most
 other species follow their own arbitrary rules) defines the ID type for htseq-count. These
 are the tags in the last column of a gff file.

 jprefix(10): Used all over the place to define the prefix job number.

 modules(bowtie1)

=cut
sub Bowtie {
    my ($class, %args) = @_;
    say('Recall that you can change the bowtie arguments via "bt_type".');
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        bt_type => 'v0M1',
        count => 1,
        libtype => 'genome',
        htseq_type => 'gene',
        htseq_id => 'ID',
        jmem => 12,
        jprefix => '10',
        modules => ['bowtie1'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('bowtie-build');
    die('Could not find bowtie in your PATH.') unless($check);
    my $start_species = $options->{species};
    my $species = $start_species;
    if ($species =~ /\:/) {
        my @species_lst = split(/:/, $species);
        my @result_lst = ();
        foreach my $sp (@species_lst) {
            print "Invoking bowtie on ${sp}\n";
            $class->{species} = $sp;
            my $result = $class->Bio::Adventure::Map::Bowtie(%{$options});
            push (@result_lst, $result);
        }
        $options->{species} = $start_species;
        return(@result_lst);
    }

    my $bt_type = $options->{bt_type};
    my $bt_args = $options->{bt_args}->{$bt_type};
    $bt_args = ' --best -v 0 -M 1 ' if (!defined($bt_args));
    my $bt_input = $options->{input};

    my $paired = 0;
    my $test_file = "";
    if ($bt_input =~ /\:|\;|\,|\s+/) {
        $paired = 1;
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
    ## $jname = $options->{jname} if ($options->{jname});
    my $libtype = $options->{libtype};

    my $count = $options->{count};
    my $bt_dir = qq"outputs/bowtie_${species}";
    $bt_dir = $options->{bt_dir} if ($options->{bt_dir});

    ## Check that the indexes exist
    my $bt_reflib = qq"$options->{libdir}/$options->{libtype}/indexes/$options->{species}";
    my $bt_reftest = qq"${bt_reflib}.1.ebwt";
    my $current_prefix = 10;
    if (defined($options->{jprefix})) {
        $current_prefix = $options->{jprefix};
    }
    my $index_job;
    if (!-r $bt_reftest && !$options->{bt1_indexjobs}) {
        $class->{bt1_indexjobs} = 1;
        my $genome_input = qq"$options->{libdir}/$options->{libtype}/${species}.fasta";
        $index_job = $class->Bio::Adventure::Index::BT1_Index(
            input => $genome_input,
            jdepends => $options->{jdepends},
            jprefix => $current_prefix,
            libtype => $libtype,);
        $options->{jdepends} = $index_job->{jobid};
    }

    my $bowtie_input_flag = "-q"; ## fastq by default
    $bowtie_input_flag = "-f" if ($options->{input} =~ /\.fasta/);
    my $error_file = qq"${bt_dir}/$options->{jbasename}-${bt_type}.stderr";
    my $comment = qq!## This is a bowtie1 alignment of ${bt_input} against
## ${bt_reflib} using arguments: ${bt_args}.!;
    my $aligned_filename = qq"${bt_dir}/$options->{jbasename}-${bt_type}_aligned_${species}.fastq";
    my $unaligned_filename = qq"${bt_dir}/$options->{jbasename}-${bt_type}_unaligned_${species}.fastq";
    my $sam_filename = qq"${bt_dir}/$options->{jbasename}-${bt_type}.sam";
    my $jstring = qq!mkdir -p ${bt_dir}
bowtie \\
  ${bt_reflib} \\
  ${bt_args} \\
  -p $options->{cpus} \\
  ${bowtie_input_flag} ${bt_input} \\
  --un ${unaligned_filename} \\
  --al ${aligned_filename} \\
  -S ${sam_filename} \\
  2>${error_file} \\
  1>${bt_dir}/$options->{jbasename}-${bt_type}.stdout
!;

    my $bt_job = $class->Submit(
        aligned => $aligned_filename,
        comment => $comment,
        input => $bt_input,
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
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

    my $compress_files = qq"${bt_dir}/$options->{jbasename}-${bt_type}_unaligned_${species}.fastq:${bt_dir}/$options->{jbasename}-${bt_type}_aligned_${species}.fastq";
    my $comp = $class->Bio::Adventure::Compress::Recompress(
        comment => '## Compressing the sequences which failed to align against ${bt_reflib} using options ${bt_args}.',
        jdepends => $bt_job->{job_id},
        jname => 'xzun',
        jprefix => $options->{jprefix} + 1,
        jwalltime => '24:00:00',
        input => $compress_files,);
    $bt_job->{compression} = $comp;

    ## BT1_Stats also reads the trimomatic output, which perhaps it should not.
    my $trim_output_file = qq"outputs/trimomatic_stats.csv";
    my $sam_job = $class->Bio::Adventure::Convert::Samtools(
        cpus => 1,
        input => $sam_filename,
        jdepends => $bt_job->{job_id},
        jname => qq"s2b_${jname}",
        jprefix => $options->{jprefix} + 3,
        paired => $paired,
        species => $options->{species},);
    $bt_job->{samtools} = $sam_job;

    my $htmulti;
    if ($count) {
        if ($libtype eq 'rRNA') {
            $htmulti = $class->Bio::Adventure::Count::HTSeq(
                htseq_id => $options->{htseq_id},
                input => $sam_job->{output},
                htseq_type => $options->{htseq_type},
                jdepends => $sam_job->{job_id},
                jname => qq"ht_${jname}",
                jprefix => $options->{jprefix} + 4,
                libtype => $libtype,
                jqueue => 'workstation',
                suffix => $bt_type,
                mapper => 'bowtie1',);
            $bt_job->{rRNA_count} = $htmulti;
        } else {
            $htmulti = $class->Bio::Adventure::Count::HT_Multi(
                htseq_id => $options->{htseq_id},
                input => $sam_job->{output},
                htseq_type => $options->{htseq_type},
                jdepends => $sam_job->{job_id},
                jname => qq"ht_${jname}",
                jprefix => $options->{jprefix} + 4,
                libtype => $libtype,
                jqueue => 'workstation',
                suffix => $bt_type,
                mapper => 'bowtie1',);
            $bt_job->{htseq} = $htmulti;
        }
    }  ## End if ($count)

    my $stats = $class->Bio::Adventure::Metadata::BT1_Stats(
        input => $error_file,
        bt_type => $bt_type,
        count_table => $bt_job->{htseq}->[0]->{output},
        jdepends => $bt_job->{job_id},
        jname => qq"${jname}_stats",
        jprefix => $options->{jprefix} + 5,
        trim_input => ${trim_output_file},);
    $bt_job->{stats} = $stats;

    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($bt_job);
}

=head2 C<Bowtie2>

 Perform a bowtie2 alignment.

 This is pretty much a twin to the Bowtie function above with the
 obvious caveat that it uses bowtie2. It converts the resulting sam
 alignment to a sorted-compressed-indexed bam, count it with
 htseq-count, compress the various output fastq files, and collect a
 few alignment statistics.

=item C<Arguments>

 input(required): colon-separated pair of (compressed)fastq files.
 species(required): used to look for bowtie indexes in
  ${libdir}/${libtype}/indexes/${species}.
 bt_type(v0M1: no mismatches, randomly place multi-matches in 1 location):
  defines the mismatch and multimatch parameters; I pre-defined a few
  likely option sets for these rather important options.
 count(1: e.g. yes): argument defines whether htseq-count will be performed.
 libtype(genome: as opposed to rRNA or contaminants etc): defines the type of
  index to search against.
 htseq_type(gene): defines the type of feature to count with
  htseq-count.  This is effectively the third column of a gff file.
 htseq_id(ID: ID is common for eukaryotic organisms, locus_tag is common for bacteria, most
  other species follow their own arbitrary rules): defines the ID type for htseq-count. These
  are the tags in the last column of a gff file.
 jprefix(10): Prefix for the job/output directory
 modules(bowtie2): Load this environment module.

=cut
sub Bowtie2 {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input',],
        count => 1,
        htseq_type => 'gene',
        htseq_id => 'ID',
        jmem => 28,
        jprefix => 20,
        modules => ['bowtie2']);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('bowtie2-build');
    die('Could not find bowtie2 in your PATH.') unless($check);
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
    my %bt_jobs = ();
    my $libtype = 'genome';
    my $bt2_args = $options->{bt2_args};
    my $prefix_name = qw"bt2";
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
    my $test_file = '';
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
        my $genome_input = qq"$options->{libdir}/$options->{libtype}/$options->{species}.fasta";
        $index_job = $class->Bio::Adventure::Index::BT2_Index(
            input => $genome_input,
            jdepends => $options->{jdepends},
            jprefix => $options->{jprefix} - 1,
            libtype => $libtype,);
        $options->{jdepends} = $index_job->{jobid};
    }
    my $bowtie_input_flag = '-q '; ## fastq by default
    $bowtie_input_flag = '-f ' if (${bt_input} =~ /\.fasta$/);

    my $cpus = $options->{cpus};
    my $error_file = qq"${bt_dir}/$options->{jbasename}.stderr";
    my $comment = qq!## This is a bowtie2 alignment of ${bt_input} against
## ${bt_reflib} using arguments: ${bt2_args}.
!;
    my $aligned_filename = qq"${bt_dir}/$options->{jbasename}_aligned_$options->{species}.fastq";
    my $unaligned_filename = qq"${bt_dir}/$options->{jbasename}_unaligned_$options->{species}.fastq";
    my $sam_filename = qq"${bt_dir}/$options->{jbasename}.sam";
    my $jstring = qq!mkdir -p ${bt_dir}
bowtie2 -x ${bt_reflib} ${bt2_args} \\
    -p ${cpus} \\
    ${bowtie_input_flag} ${bt_input} \\
    --un ${unaligned_filename} \\
    --al ${aligned_filename} \\
    -S ${sam_filename} \\
    2>${error_file} \\
    1>${bt_dir}/$options->{jbasename}.stdout
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
    my $compression_files = qq"${bt_dir}/$options->{jbasename}_unaligned_$options->{species}.fastq:${bt_dir}/$options->{jbasename}_aligned_$options->{species}.fastq";
    my $comp = $class->Bio::Adventure::Compress::Recompress(
        comment => '## Compressing the sequences which failed to align against ${bt_reflib} using options ${bt2_args}.',
        input => $compression_files,
        jdepends => $bt2_job->{job_id},
        jname => qq"xzun_${suffix_name}",
        jprefix => $options->{jprefix} + 1,);
    $bt2_job->{compression} = $comp;
    ## BT1_Stats also reads the trimomatic output, which perhaps it should not.
    ## my $trim_output_file = qq"outputs/$options->{jbasename}-trimomatic.out";
    my $stats = $class->Bio::Adventure::Metadata::BT2_Stats(
        input => $error_file,
        count_table => qq"$options->{jbasename}.count.xz",
        jdepends => $bt2_job->{job_id},
        jname => qq"bt2st_${suffix_name}",
        jprefix => $options->{jprefix} + 3,);
    $bt2_job->{stats} = $stats;
    my $sam_job = $class->Bio::Adventure::Convert::Samtools(
        input => $sam_filename,
        jdepends => $bt2_job->{job_id},
        jname => qq"s2b_${suffix_name}",
        jprefix => $options->{jprefix} + 4,);
    $bt2_job->{samtools} = $sam_job;
    my $htseq_input = $sam_job->{output};
    my $htmulti;
    if ($options->{count}) {
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
    $loaded = $class->Module_Loader(modules => $options->{modules},
        action => 'unload',);
    return($bt2_job);
}

=head2 C<Bowtie_RRNA>

 Search for rRNA using bowtie.

 This function is probably extraneous at this point.  It simply calls Bowtie()
 with the libtype set to 'rRNA' in order to get it to look for ribosomal reads
 instead of the default, genomic reads.

=cut
sub Bowtie_RRNA {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        exclude => 0,
        required => ["species"],);
    my $job_name = qq"rRNA_$options->{jbasename}";
    my $exclude = $options->{exclude};
    my $species = $options->{species};
    my ${bt_dir} = qq"outputs/bowtie_$options->{species}";
    my $job = $class->Bio::Adventure::Map::Bowtie(
        jdepends => $options->{jdepends},
        jname => $job_name,
        libtype => 'rRNA',
        prescript => $args{prescript},
        postscript => $args{postscript},);
    return($job);
}

=head2 C<BT_Multi>

 Run multiple bowtie parameters at once.

 Attempts to run multiple bowtie1 runs for a given species.  One run is performed
 for each of a few parameter sets which are kept in the global 'bt_args' variable.
 and generally include: 0 mismatch, 1 mismatch, 2 mismatches, 1 randomly placed
 hit, 0 randomly placed hits, or default options.

 This should either be removed or modified to work more generally with other tools.

=cut
sub BT_Multi {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        htseq_type => 'gene',
        htseq_id => 'ID',
        required => ['species', 'input'],);
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

=head2 C<BWA>

 The other OG short read aligner!
 10.1093/bioinformatics/btp324

 Perform a bwa alignment using both the sam(s|p)e and aln algorithms.  It then
 converts the output (when appropriate) to sorted/indexed bam and passes them to
 htseq.

=item C<Arguments>

 input(required): likely a colon-separated pair of (compressed)fastq files.
 species(required): will be used to look for bwa indexes in
  ${libdir}/${libtype}/indexes/${species}.
 count(1: e.g. yes): argument defines whether htseq-count will be performed.
 libtype(genome: as opposed to rRNA or contaminants etc): defines the type of
  index to search against.
 htseq_type(gene): defines the type of feature to count with
  htseq-count.  This is effectively the third column of a gff file.
 htseq_id(ID: ID is common for eukaryotic organisms, locus_tag is common for bacteria, most
  other species follow their own arbitrary rules): defines the ID type for htseq-count. These
  are the tags in the last column of a gff file.
 jprefix(30): Prefix for jobname and output directory.
 modules(bwa): load this environment module

=cut
sub BWA {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        count => 1,
        species => 'lmajor',
        libtype => 'genome',
        htseq_type => 'gene',
        htseq_id => 'ID',
        jmem => 24,
        jprefix => 30,
        modules => ['bwa'],);
    my $check = which('bwa');
    die('Could not find bwa in your PATH.') unless($check);

    if ($options->{species} =~ /\:/) {
        my @species_lst = split(/:/, $options->{species});
        my @result_lst = ();
        my $start_species = $options->{species};
        foreach my $sp (@species_lst) {
            print "Invoking bwa on ${sp}\n";
            $options->{species} = $sp;
            my $result = $class->Bio::Adventure::Map::BWA(%{$options});
            push (@result_lst, $result);
        }
        $options->{species} = $start_species;
        return(@result_lst);
    }

    my $bwa_input = $options->{input};
    my $test_file = '';
    my $forward_reads = '';
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

    my $bwa_dir = qq"outputs/$options->{jprefix}bwa_$options->{species}";
    $bwa_dir = $options->{bwa_dir} if ($options->{bwa_dir});

    my $uncompress_jobid = undef;
    my $index_jobid = undef;
    ## Check that the indexes exist
    ## NOTE: BWA Might require the file to be named ...fa instead of ...fasta
    my $bwa_reflib = qq"$options->{libdir}/${libtype}/indexes/$options->{species}.fa";
    my $bwa_reftest = qq"${bwa_reflib}.sa";
    my $index_job;
    if (!-r $bwa_reftest) {
        $index_job = $class->Bio::Adventure::Index::BWA_Index(
            input => $bwa_reflib,
            jdepends => $options->{jdepends},
            jprefix => $options->{jprefix} - 1,
            libtype => $libtype,);
        $options->{jdepends} = $index_job->{job_id};
    }

    my $aln_sam = qq"${bwa_dir}/$options->{jbasename}_aln.sam";
    my $mem_sam = qq"${bwa_dir}/$options->{jbasename}_mem.sam";
    my $aln_args = '';
    my $mem_args = '-M ';
    my $sam_args = '';
    my $mem_string = qq!mkdir -p ${bwa_dir}
bwa mem ${mem_args} \\
  -a ${bwa_reflib} \\
  ${bwa_input} \\
  2>${bwa_dir}/bwa.stderr \\
  1>${mem_sam}
!;
    my $reporter_string = qq"bwa samse ${sam_args} \\
  ${bwa_reflib} \\
  ${bwa_dir}/$options->{jbasename}_aln-forward.sai \\
  ${bwa_input} \\
  2>${bwa_dir}/$options->{jbasename}_samse.stderr \\
  1>${aln_sam}";
    my $aln_string = qq"bwa aln ${aln_args} \\
  ${bwa_reflib} \\
  ${forward_reads} \\
  2>${bwa_dir}/$options->{jbasename}_aln-forward.stderr \\
  1>${bwa_dir}/$options->{jbasename}_aln-forward.sai";
    if (defined($reverse_reads)) {
        $aln_string = qq"${aln_string}
bwa aln ${aln_args} \\
  ${bwa_reflib} \\
  <(less ${reverse_reads}) \\
  2>${bwa_dir}/$options->{jbasename}_aln-reverse.stderr \\
  1>${bwa_dir}/$options->{jbasename}_aln-reverse.sai";
        $reporter_string = qq"bwa ${sam_args} \\
  sampe ${bwa_reflib} \\
  ${bwa_dir}/$options->{jbasename}_aln-forward.sai \\
  ${bwa_dir}/$options->{jbasename}_aln-reverse.sai \\
  <(less ${forward_reads}) <(less ${reverse_reads}) \\
  2>${bwa_dir}/$options->{jbasename}_sam.stderr \\
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
        jname => qq"bwamem_$options->{species}",
        output => $mem_sam,
        jprefix => $options->{jprefix},
        jstring => $mem_string,
        jmem => $options->{jmem},
        modules => $options->{modules},
        postscript => $options->{postscript},
        prescript => $options->{prescript},
        jqueue => 'workstation',
    );
    my $mem_sam_job = $class->Bio::Adventure::Convert::Samtools(
        input => $mem_sam,
        jdepends => $bwa_job->{job_id},
        jname => 's2b_mem',
        jmem => '20',
        jprefix => $options->{jprefix} + 1,);
    $bwa_job->{samtools_mem} = $mem_sam_job;

    ## ALN Runs
    my $aln_job = $class->Submit(
        comment => $aln_comment,
        input => $bwa_input,
        jdepends => $mem_sam_job->{job_id},
        jname => qq"bwaaln_$options->{species}",
        output => qq"${bwa_dir}/$options->{jbasename}_aln-forward.sai",
        jprefix => $options->{jprefix} + 2,
        jstring => $aln_string,
        jmem => $options->{jmem},
        modules => $options->{modules},
        postscript => $options->{postscript},
        prescript => $options->{prescript},
        jqueue => 'workstation',);
    $bwa_job->{aln} = $aln_job;

    my $rep_job = $class->Submit(
        comment => $report_comment,
        input => $aln_job->{output},
        jdepends => $aln_job->{job_id},
        jname => qq"bwarep_$options->{species}",
        modules => $options->{modules},
        output => $aln_sam,
        jprefix => $options->{jprefix} + 3,
        jstring => $reporter_string,
        jmem => $options->{jmem},
        postscript => $options->{postscript},
        prescript => $options->{prescript},
        jqueue => 'workstation',);
    $bwa_job->{reporter} = $rep_job;

    my $aln_sam_job = $class->Bio::Adventure::Convert::Samtools(
        input => $aln_sam,
        jdepends => $rep_job->{job_id},
        jname => 's2b_aln',
        jmem => '20',
        jprefix => $options->{jprefix} + 4,);
    $bwa_job->{samtools_aln} = $aln_sam_job;

    if ($options->{count}) {
        my $mem_htmulti = $class->Bio::Adventure::Count::HT_Multi(
            htseq_id => $options->{htseq_id},
            input => $mem_sam_job->{output},
            htseq_type => $options->{htseq_type},
            jdepends => $mem_sam_job->{job_id},
            jname => qq"htmem_${jname}",
            jprefix => $options->{jprefix} + 5,
            mapper => 'bwa',);
        $bwa_job->{htseq_mem} = $mem_htmulti;
        my $aln_htmulti = $class->Bio::Adventure::Count::HT_Multi(
            htseq_id => $options->{htseq_id},
            input => $aln_sam_job->{output},
            htseq_type => $options->{htseq_type},
            jdepends => $aln_sam_job->{job_id},
            jname => qq"htaln_${jname}",
            jprefix => $options->{jprefix} + 6,
            mapper => 'bwa',);
        $bwa_job->{htseq_aln} = $aln_htmulti;
    }

    my $bwa_stats = $class->Bio::Adventure::Metadata::BWA_Stats(
        jdepends => $mem_sam_job->{job_id},
        jname => 'bwastats',
        jprefix => $options->{jprefix} + 7,
        aln_output => $aln_sam_job->{output},
        mem_output => $mem_sam_job->{output},);
    $bwa_job->{stats} = $bwa_stats;

    return($bwa_job);
}

=head2 C<Hisat2>

 Invoke hisat2!
 10.1038/s41587-019-0201-4

 Hisat2 is currently my favorite aligner.

=item C<Arguments>

 input(required): Colon separated fastq input files.
 species(required): Set the indexes.
 htseq_type('gene'): Define the htseq-count type parameter.
 htseq_id('ID'): Define the htseq-count id attribute parameter.
 count(1): Count after aligning?
 libtype('genome'): Change this for different index types (contaminants/rRNA/genome).
 jmem(24): Expected memory required.
 jprefix('40'): Set the job prefix and output directory.
 modules('hisat2','samtools','htseq','bamtools'): Load these
  environment modules.

=cut
sub Hisat2 {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input',],
        htseq_type => 'gene',
        htseq_id => 'ID',
        count => 1,
        libtype => 'genome',
        jmem => 24,
        jprefix => '40',
        maximum => undef,
        modules => ['hisat2', 'samtools', 'htseq', 'bamtools'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('hisat2-build');
    die('Could not find hisat2 in your PATH.') unless($check);

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
        $ready = $class->Check_Input(files => $options->{input},);
    }
    my $hisat_args = '';
    $hisat_args = $options->{hisat_args} if ($options->{hisat_args});
    $hisat_args .= qq" -k $options->{maximum} " if (defined($options->{maximum}));
    my $prefix_name = 'hisat2';
    $prefix_name .= qq"_k$options->{maximum}" if (defined($options->{maximum}));
    my $hisat_name = qq"${prefix_name}_$options->{species}_$options->{libtype}";
    my $suffix_name = $prefix_name;
    if ($options->{jname}) {
        $hisat_name .= qq"_$options->{jname}";
        $suffix_name .= qq"_$options->{jname}";
    }

    my $hisat_dir = qq"outputs/$options->{jprefix}hisat2_$options->{species}";
    my $hisat_input = $options->{input};
    my $test_file = '';
    my $paired = 0;
    my @pair_listing;
    if ($hisat_input =~ /\:|\;|\,|\s+/) {
        @pair_listing = split(/\:|\;|\,|\s+/, $hisat_input);
        $paired = 1;
        $pair_listing[0] = File::Spec->rel2abs($pair_listing[0]);
        $pair_listing[1] = File::Spec->rel2abs($pair_listing[1]);
        ## After years of working without problem, suddenly my lesspipe
        ## process substitution pre-filter for arbitrarily compressed files
        ## stopped working and might well give me an aneurysm trying to figure out.
        ## It turns out that there is a race condition somewhere which is triggered when
        ## less buffers its output -- so make sure the environment variable 'LESS'
        ## contains --unbuffered
        $hisat_input = qq" -1 <(less $pair_listing[0]) -2 <(less $pair_listing[1]) ";
        $test_file = $pair_listing[0];
    } else {
        $test_file = File::Spec->rel2abs($hisat_input);
        $hisat_input = qq" <(less ${test_file}) ";
    }

    ## Check that the indexes exist
    my $hisat_reflib = qq"$options->{libdir}/$options->{libtype}/indexes/$options->{species}";
    my $hisat_reftest = qq"${hisat_reflib}.1.ht2";
    my $hisat_reftestl = qq"${hisat_reflib}.1.ht2l";
    if (!-r $hisat_reftest && !-r $hisat_reftestl) {
        print "Hey! The Indexes do not appear to exist, check this out: ${hisat_reftest}\n";
        sleep(10);
        my $genome_fasta = qq"$options->{libdir}/$options->{libtype}/$options->{species}.fasta";
        my $index_job = $class->Bio::Adventure::Index::Hisat2_Index(
            input => $genome_fasta,
            jprefix => $options->{jprefix} - 1,
            jdepends => $options->{jdepends},
            libtype => $options->{libtype},);
        ## The following line inserts the indexer into the dependency chain
        ## Do not forget this, it is rather important.
        $options->{jdepends} = $index_job->{job_id};
    }

    my $hisat_input_flag = '-q '; ## fastq by default
    $hisat_input_flag = '-f ' if (${hisat_input} =~ /\.fasta$/);
    my $cpus = $options->{cpus};
    my $error_file = qq"${hisat_dir}/hisat2_$options->{species}_$options->{libtype}_$options->{jbasename}.stderr";
    my $comment = qq!## This is a hisat2 alignment of ${hisat_input} against ${hisat_reflib}
!;
    $comment .= qq"## This alignment is using arguments: ${hisat_args}.\n" unless ($hisat_args eq '');
    my $aligned_discordant_filename = qq"${hisat_dir}/$options->{jbasename}_aldis_$options->{species}_$options->{libtype}.fastq";
    my $unaligned_discordant_filename = qq"${hisat_dir}/$options->{jbasename}_unaldis_$options->{species}_$options->{libtype}.fastq";
    my $aligned_concordant_filename = qq"${hisat_dir}/$options->{jbasename}_alcon_$options->{species}_$options->{libtype}.fastq";
    my $unaligned_concordant_filename = qq"${hisat_dir}/$options->{jbasename}_unalcon_$options->{species}_$options->{libtype}.fastq";
    my $sam_filename = qq"${hisat_dir}/$options->{jbasename}_$options->{species}_$options->{libtype}.sam";
    my $jstring = qq!mkdir -p ${hisat_dir}
hisat2 -x ${hisat_reflib} ${hisat_args} \\
  -p ${cpus} \\
  ${hisat_input_flag} ${hisat_input} \\
  --phred$options->{phred} \\
  --un ${unaligned_discordant_filename} \\
  --al ${aligned_discordant_filename} \\
!;
    ## Record the aligned/unaligned filenames both before and after compression.
    my $aligned_filenames = $aligned_discordant_filename;
    my $aligned_xz = qq"${aligned_filenames}.xz";
    my $unaligned_filenames = $unaligned_discordant_filename;
    my $unaligned_xz = qq"${unaligned_filenames}.xz";
    if ($paired) {
        ## The concordant flag tries to send reads to xxx.1 and xxx.2 even if
        ## the data is not paired.
        $jstring .= qq!  --un-conc ${unaligned_concordant_filename} \\
  --al-conc ${aligned_concordant_filename} \\
!;
        my $aligned_base = basename($aligned_concordant_filename, ('.fastq', '.fasta'));
        my $unaligned_base = basename($unaligned_concordant_filename, ('.fastq', '.fastq'));
        my $the_dirname = dirname($aligned_concordant_filename);
        ## For the set of all aligned filenames, there should now be 3, the original discordant and two more.
        ## For the xz filenames, there should be 2, just the concordant; because these will be used downstream.
        $aligned_filenames .= qq":${the_dirname}/${aligned_base}.1.fastq:${the_dirname}/${aligned_base}.2.fastq";
        $aligned_xz = qq"${the_dirname}/${aligned_base}.1.fastq.xz:${the_dirname}/${aligned_base}.2.fastq.xz";
        $unaligned_filenames .= qq":${the_dirname}/${unaligned_base}.1.fastq:${the_dirname}/${unaligned_base}.2.fastq";
        $unaligned_xz = qq"${the_dirname}/${unaligned_base}.1.fastq.xz:${the_dirname}/${unaligned_base}.2.fastq.xz";
    }
    my $all_filenames = qq"${aligned_filenames}:${unaligned_filenames}";
    $jstring .= qq!  -S ${sam_filename} \\
  2>${error_file} \\
  1>${hisat_dir}/hisat2_$options->{species}_$options->{libtype}_$options->{jbasename}.stdout
!;
    my $hisat_job = $class->Submit(
        comment => $comment,
        jdepends => $options->{jdepends},
        input => $hisat_input,
        jname => $hisat_name,
        jstring => $jstring,
        jprefix => $options->{jprefix},
        jmem => $options->{jmem},
        modules => $options->{modules},
        paired => $paired,
        output => $sam_filename,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        aligned => $aligned_filenames,
        aligned_comp => $aligned_xz,
        unaligned => $unaligned_filenames,
        unaligned_comp => $unaligned_xz,);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    my $xz_jname = qq"xz_${suffix_name}";
    my $comp = $class->Bio::Adventure::Compress::Recompress(
        jname => $xz_jname,
        input => qq"${aligned_filenames}:${unaligned_filenames}",
        jdepends => $hisat_job->{job_id});
    $hisat_job->{compression} = $comp;

    ## Sneak the compression job's ID in place as hisat's.
    $hisat_job->{job_id} = $comp->{job_id};
    ## HT1_Stats also reads the trimomatic output, which perhaps it should not.
    ## my $trim_output_file = qq"outputs/$options->{jbasename}-trimomatic.out";
    my $new_jprefix = qq"$options->{jprefix}_1";
    my $sam_jprefix = qq"$options->{jprefix}_2";
    my $sam_jname = qq"s2b_${suffix_name}";
    my $sam_job = $class->Bio::Adventure::Convert::Samtools(
        cpus => 1,
        input => $sam_filename,
        jdepends => $hisat_job->{job_id},
        jname => $sam_jname,
        jprefix => $sam_jprefix,
        paired => $paired,
        species => $options->{species},);
    $hisat_job->{samtools} = $sam_job;

    $new_jprefix = qq"$options->{jprefix}_3";
    my $htseq_input;
    if ($paired == 1) {
        $htseq_input = $sam_job->{paired_output};
    } else {
        $htseq_input = $sam_job->{output};
    }
    my $htmulti;
    if ($options->{count}) {
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
        }
        $hisat_job->{htseq} = $htmulti;
    }  ## End checking if we should do htseq

    my $htseq_out = '';
    if (ref($hisat_job->{htseq}) eq 'ARRAY') {
        $htseq_out = $hisat_job->{htseq}->[0]->{output};
    } else {
        $htseq_out = $hisat_job->{htseq}->{output};
    }
    my $stats = $class->Bio::Adventure::Metadata::HT2_Stats(
        ht_input => $error_file,
        count_table => $htseq_out,
        jdepends => $hisat_job->{job_id},
        jname => qq"hisat2st_${suffix_name}",
        jprefix => $new_jprefix,
        output_dir => $hisat_dir,);
    $hisat_job->{stats} = $stats;
    ## If this is not the final job in a chain, then make sure
    ## that any jobs queued after it do not start until
    ## samtools/htseq/etc are finished.
    $hisat_job->{job_id} = $stats->{job_id};
    return($hisat_job);
}

=head2 C<Kallisto>

 Perform a kallisto transcript quantification.
 10.1038/nbt.3519

 Kallisto and salmon are my two favorite 'voting' based aligners.

=item C<Arguments>

 input(required): Colon separated fastq input file(s).
 species(required): Define the location of the indexes.
 jmem(24): Expected memory requirements.
 jprefix('46'): Set the output directory/job prefix.
 modules('kallisto'): Load this environment module.

=cut
sub Kallisto {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input'],
        jmem => 24,
        jprefix => '46',
        modules => ['kallisto'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('kallisto');
    die('Could not find kallisto in your PATH.') unless($check);

    if ($options->{species} =~ /\:/) {
        my @species_lst = split(/:/, $options->{species});
        my @result_lst = ();
        my $start_species = $options->{species};
        foreach my $sp (@species_lst) {
            print "Invoking kallisto on ${sp}\n";
            $options->{species} = $sp;
            my $result = $class->Bio::Adventure::Map::Kallisto(%{$options});
            push (@result_lst, $result);
        }
        $options->{species} = $start_species;
        return(@result_lst);
    }

    my %ka_jobs = ();
    my $ka_depends_on = '';
    my $libtype = 'genome';
    $libtype = $options->{libtype} if ($options->{libtype});
    my $species = $options->{species};
    my $jname = qq"kall_${species}";
    ## $jname = $options->{jname} if ($options->{jname});
    my $ka_args = '';
    my $ka_input = $options->{input};
    my $input_name = $ka_input;
    if ($ka_input =~ /\:|\;|\,|\s+/) {
        my @pair_listing = split(/\:|\;|\,|\s+/, $ka_input);
        $ka_args .= ' --bias ';
        if ($options->{stranded} != 0) {
            $ka_args .= " --$options->{stranded} ";
        }
        $ka_input = qq" <(less $pair_listing[0]) <(less $pair_listing[1]) ";
        $input_name = $pair_listing[0];
    } else {
        $ka_input = qq" <(less $ka_input) ";
        $ka_args .= qq' --bias --single -l 40 -s 10 ';
    }

    ## Check that the indexes exist
    my $ka_reflib = qq"$options->{libdir}/${libtype}/indexes/$options->{species}.idx";
    my $index_job;
    if (!-r $ka_reflib) {
        my $transcriptome_fasta = qq"$options->{libdir}/${libtype}/$options->{species}_cds.fasta";
        $index_job = $class->Bio::Adventure::Index::Kallisto_Index(
            input => $transcriptome_fasta,
            jdepends => $options->{jdepends},
            libtype => $libtype,);
        $options->{jdepends} = $index_job->{job_id};
    }

    my $outdir = qq"outputs/$options->{jprefix}kallisto_${species}";
    my $error_file = qq"${outdir}/kallisto_${species}.stderr";
    my $output_sam = qq"${outdir}/kallisto_${species}.sam";
    my $output_bam = qq"${outdir}/kallisto_${species}.bam";
    my $output_stats = qq"${outdir}/kallisto_${species}.stats";
    my $sorted_bam = qq"${outdir}/kallisto_${species}-sorted";
    my $comment = qq!## This is a kallisto pseudoalignment of ${ka_input} against
## ${ka_reflib}.
## Other candidates for making a pretty count table include:
##  perl -F'\\t' -a -n -e 'print "\$F[0] \$F[3]\\n"' ${outdir}/abundance.tsv > ${outdir}/abundance.count
##   awk '{printf("%s %s\\n", \$1, \$4)}' ${outdir}/abundance.tsv > ${outdir}/abundance.count
## The sam->bam conversion is copy/pasted from Sam2Bam() but I figured why start another job
## because kallisto is so fast
!;
    my $dropped_args = ' --pseudobam ';
    my $jstring = qq!mkdir -p ${outdir}
kallisto quant ${ka_args} \\
  --plaintext -t 4 -b 100 \\
  -o ${outdir} \\
  -i ${ka_reflib} \\
  ${ka_input} \\
  2>${error_file} \\
  1>${output_sam} && \\
  cut -d "	" -f 1,4 ${outdir}/abundance.tsv > ${outdir}/${input_name}_abundance.count && \\
  xz -9e -f ${outdir}/${input_name}_abundance.count
!;
    ## Newer kallisto does not seem to do this well anymore...
##  samtools view -u -t ${ka_reflib} -S ${output_sam} \\
##    2>${output_bam}.err \\
##    1>${output_bam} && \\
##  samtools sort -l 9 ${output_bam} ${sorted_bam} \\
##    2>${sorted_bam}.err \\
##    1>${sorted_bam}.out
##mv ${sorted_bam}.bam ${output_bam} && samtools index ${output_bam} && \\
##  bamtools stats -in ${output_bam} 2>${output_stats} 1>&2
##!;
    my $kallisto = $class->Submit(
        comment => $comment,
        input => $ka_input,
        depends => $options->{jdepends},
        jname => qq"${jname}",
        jprefix => '30',
        jstring => $jstring,
        jmem => $options->{jmem},
        output => qq"${outdir}/abundance.tsv",
        count => qq"${outdir}/${input_name}_abundance.count",
        prescript => $args{prescript},
        postscript => $args{postscript},);

    $kallisto->{index} = $index_job;
    return($kallisto);
}

=head2 C<RSEM>

 Invoke RSEM.
 10.1186/1471-2105-12-323

 The most accurate and slow transcript quantification method.

=cut
sub RSEM {
    my ($class, %args) = @_;
    my $check = which('rsem-prepare-reference');
    die('Could not find RSEM in your PATH.') unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input'],
        jmem => 24,
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
        $index_job = $class->Bio::Adventure::Index::RSEM_Index(
            input => $cds,
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
    my $rsem_comment = '## This is a rsem invocation script.';
    my $output_file = qq"$options->{species}_something.txt";
    my $jstring = qq"mkdir -p ${rsem_dir} && rsem-calculate-expression --bowtie2 \\
  --calc-ci ${rsem_input} \\
  ${idx} \\
  ${rsem_dir}/$options->{species} \\
  2>${rsem_dir}/$options->{species}.stderr \\
  1>${rsem_dir}/$options->{species}.stdout
";
    my $rsem = $class->Submit(
        comment => $rsem_comment,
        input => $rsem_input,
        jname => qq"rsem_${jbasename}",
        jdepends => $options->{jdepends},
        jstring => $jstring,
        jprefix => '28',
        mem => $options->{jmem},
        modules => $options->{modules},
        output => $output_file,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    $rsem->{index_job} = $index_job;
    return($rsem);
}

=head2 C<Salmon>

 Perform a salmon quantification of transcript abundances.
 10.1038/nmeth.4197

 My favorite transcript aware quantification method.

=cut
sub Salmon {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input'],
        jmem => 24,
        jprefix => '45',
        modules => ['salmon'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('salmon');
    die('Could not find salmon in your PATH.') unless($check);

    my $depends = $options->{jdepends};
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
    my %sa_jobs = ();
    my $libtype = 'genome';
    $libtype = $options->{libtype} if ($options->{libtype});
    my $species = $options->{species};

    my $jname = qq"sal_${species}";
    ## $jname = $options->{jname} if ($options->{jname});
    my $sa_args = '';
    my $sa_input = $options->{input};
    my $input_name = $sa_input;
    if ($sa_input =~ /\:|\;|\,|\s+/) {
        my @pair_listing = split(/\:|\;|\,|\s+/, $sa_input);
        $sa_args .= qq" -1 <(less $pair_listing[0]) -2 <(less $pair_listing[1]) ";
        $input_name = $pair_listing[0];
    } else {
        $sa_args .= qq" -r <(less $sa_input) ";
    }

    ## Check that the indexes exist
    my $sa_reflib = qq"$options->{libdir}/${libtype}/indexes/$options->{species}_salmon_index";
    my $index_job;
    if (!-r $sa_reflib) {
        my $transcript_file = qq"$options->{libdir}/${libtype}/$options->{species}_cds.fasta";
        $index_job = $class->Bio::Adventure::Index::Salmon_Index(
            input => $transcript_file,
            depends => $options->{jdepends},
            libtype => $libtype,);
        $options->{jdepends} = $index_job->{job_id};
    }

    my $outdir = qq"outputs/$options->{jprefix}salmon_${species}";
    my $error_file = qq"${outdir}/salmon_${species}.stderr";
    my $comment = qq!## This is a salmon pseudoalignment of ${sa_input} against
## ${sa_reflib}.
!;
    my $jstring = qq!mkdir -p ${outdir}
salmon quant -i ${sa_reflib} \\
  -l A --gcBias --validateMappings  \\
  ${sa_args} \\
  -o ${outdir} \\
  2>${error_file} 1>${outdir}/salmon.stdout
!;

    my $salmon = $class->Submit(
        comment => $comment,
        input => $sa_input,
        jdepends => $options->{jdepends},
        jname => qq"${jname}",
        jprefix => '30',
        jstring => $jstring,
        jmem => $options->{jmem},
        modules => $options->{modules},
        output => qq"${outdir}/quant.sf",
        prescript => $args{prescript},
        postscript => $args{postscript},);
    $salmon->{index_job} = $index_job;
    my $stats = $class->Bio::Adventure::Metadata::Salmon_Stats(
        input => qq"${outdir}/lib_format_counts.json",
        jdepends => $salmon->{job_id},
        jname => qq"sastats_$options->{species}",
        jprefix => '33',);
    $salmon->{stats} = $stats;
    return($salmon);
}

=head2 C<STAR>

 Invoke STAR for transcript abundances.
 10.1093/bioinformatics/bts635

=cut
sub STAR {
    my ($class, %args) = @_;
    my $check = which('STAR');
    die('Could not find STAR in your PATH.') unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input'],
        jmem => 48,
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

    my $ready = $class->Check_Input(files => $options->{input});
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
        $star_inputstring = qq"${star_input}";
    }

    ## Check that the indexes exist
    my $star_refdir = qq"$options->{libdir}/${libtype}/indexes/$options->{species}_star_index";
    my $star_reflib = qq"${star_refdir}/SAindex";
    my $index_job;
    if (!-r $star_reflib) {
        my $genome_file = qq"$options->{libdir}/${libtype}/$options->{species}.fasta";
        $index_job = $class->Bio::Adventure::Index::STAR_Index(
            input => $genome_file,
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
    my $jstring = qq!mkdir -p ${outdir}
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
  2>${error_file} 1>${outdir}/star_${species}.stdout

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
  2>>${error_file} 1>>${outdir}/star_${species}.stdout
!;

    my $star_job = $class->Submit(
        comment => $comment,
        input => $star_input,
        jdepends => $options->{jdepends},
        jname => qq"${jname}",
        jprefix => '33',
        jstring => $jstring,
        jmem => $options->{jmem},
        modules => $options->{modules},
        prescript => $args{prescript},
        postscript => $args{postscript},
        jqueue => 'large',);
    $star_job->{index_job} = $index_job;
    return($star_job);
}

=head2 C<Tophat>

 Invokes tophat!
 10.1186/gb-2013-14-4-r36

 It also sorts/indexes the accepted_hits file, collects some
 statistics, and passes the hits to htseq-count.

=cut
sub Tophat {
    my ($class, %args) = @_;
    my $check = which('tophat');
    die('Could not find tophat in your PATH.') unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input', 'htseq_type'],
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
    my $bt_reflib = qq"$options->{libdir}/${libtype}/indexes/$options->{species}";
    my $bt_reftest = qq"${bt_reflib}.1.bt2";
    my $index_job = undef;
    if (!-r $bt_reftest) {
        print "Did not find the index for $options->{species} at: ${bt_reflib}, indexing now.\n";
        my $genome_file = qq"$options->{libdir}/${libtype}/$options->{species}.fasta";
        $index_job = $class->Bio::Adventure::Index::BT2_Index(
            input => $genome_file,
            jdepends => $options->{jdepends},);
        $options->{jdepends} = $index_job->{job_id};
    }
    my $gtf_file = qq"$options->{libdir}/genome/$options->{species}.gtf";
    if (!-r $gtf_file) {
        print "Missing the gtf file for $options->{species}\n";
        print "Using the gff file.\n";
        $gtf_file =~ s/\.gtf/\.gff/;
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
  2>${tophat_dir}/tophat.stderr \\
  1>${tophat_dir}/tophat.stdout && \\
 samtools sort -l 9 -n ${tophat_dir}/accepted_hits.bam -o ${tophat_dir}/accepted_sorted.bam \\
  2>${tophat_dir}/samtools_sort.stderr 1>${tophat_dir}/samtools_sort.stdout && \\
 samtools index ${tophat_dir}/accepted_hits.bam \\
  2>${tophat_dir}/samtools_index.stderr 1>${tophat_dir}/samtools_index.stdout && \\
 samtools sort -l 9 -n ${tophat_dir}/unmapped.bam -o ${tophat_dir}/unmapped_sorted.bam \\
  2>>${tophat_dir}/samtools_sorted.stderr 1>>${tophat_dir}/samtools_sort.stdout && \\
 samtools index ${tophat_dir}/unmapped.bam \\
  2>>${tophat_dir}/samtools_index.stderr 1>>${tophat_dir}/samtools_sort.stdout
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
        jprefix => '31',
        jstring => $jstring,
        jmem => $tophat_mem,
        modules => $options->{mdules},
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => $tophat_queue,
        jwalltime => $tophat_walltime,);

    ## Set the input for htseq
    my $accepted = qq"${tophat_dir}/accepted_hits.bam";
    $accepted = $options->{accepted_hits} if ($options->{accepted_hits});
    my $count_table = 'accepted_hits.count';
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
    my $stats = $class->Bio::Adventure::Metadata::Tophat_Stats(
        accepted_input => $accepted,
        count_table => qq"${count_table}.xz",
        jdepends => $tophat->{job_id},
        jname => qq"tpstats_$options->{species}",
        jprefix => '33',
        prep_input => $input_read_info,
        unaccepted_input => $unaccepted,);
    $tophat->{stats} = $stats;
    return($tophat);
}

=back

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

    L<bowtie> L<bowtie2> L<tophat> L<bwa> L<kallisto> L<samtools> L<htseq>

=cut

1;
