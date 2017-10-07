package Bio::Adventure::RNASeq_Trim;
## LICENSE: gplv2
## ABSTRACT:  Kitty!
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use File::Basename;
use File::Which qw"which";

=head1 NAME
    Bio::Adventure::RNASeq_Trim - Use trimomatic/cutadapt/etc to trim libraries

=head1 SYNOPSIS

    use Bio::Adventure;
    my $hpgl = new Bio::Adventure;
    $hpgl->Cutadapt();

=head2 Methods

=item C<Cutadapt>

    $hpgl->Cutadapt(); will use biopieces/cutadapt to attempt to
    remove sequence adapters from a library.  This is most common used
    by me for ribosome profiling libraries.

=cut
sub Cutadapt {
    my ($class, %args) = @_;
    my $check = which('cutadapt');
    die("Could not find cutadapt in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'type'],
    );
    if ($options->{task}) {
        $options->{type} = $options->{task};
    }
    $options->{type} = lc($options->{type});
    my $type;
    if ($options->{type}) {
        $type = $options->{type};
    } elsif ($args{type}) {
        $type = $args{type};
    } else {
        $type = 'tnseq';
    }

    my $input = $options->{input};
    my $basename = basename($input, @{$options->{suffixes}});
    $basename = basename($basename, @{$options->{suffixes}});
    my $adapter_flags = "";
    my $minlength = 7;
    my $maxlength = 42;
    if ($type eq 'tnseq') {
        $adapter_flags = qq" -a ACAGGTTGGATGATAAGTCCCCGGTCTGACACATC -a ACAGTCCCCGGTCTGACACATCTCCCTAT -a ACAGTCCNCGGTCTGACACATCTCCCTAT ";
        $maxlength = 20;
    } elsif ($type eq 'riboseq') {
        $adapter_flags = qq" -a ACAGGTTGGATGATAAGTCCCCGGTCTGACACATCTCCCTAT -a AGATCGGAAGAGCACACGTCTGAAC -b AGATCGGAAGAGCACACGTCTGAAC ";
        $minlength = 16;
        if ($input =~ /\.csfasta/) {
            $adapter_flags = qq" -a CGCCTTGGCCGTACAGCAGCATATTGGATAAGAGAATGAGGAACCCGGGGCAG -a GCGGAACCGGCATGTCGTCGGGCATAACCCTCTCTTACTCCTTGGGCCCCGTC ";
        }
    } else {
        $adapter_flags = qq" -a ACAGGTTGGATGATAAGTCCCCGGTCTGACACATCTCCCTAT ";
    }

    my $comment = qq!## This script makes use of biopieces and cutadapt to trim away adapters
## and separate the sequence file into a few pieces depending on size
## and adapter status.  It also performs some simple graphs of the data.!;
    my $out_dir = qq"$options->{basedir}/outputs/cutadapt";
    my $type_flag = '';
    my $out_suffix = 'fastq';
    my $input_flags = qq"less ${input} | cutadapt - ";
    if ($input =~ /\.csfasta/) {
        $type_flag = '-c -t --strip-f3';
        $options = $class->Get_Vars(
            args => \%args,
            required => ["qual"],
        );
        ## $input_flags = qq"less ${input} | sed 's/^T//g' | cutadapt - $options->{qual} ";
        $input_flags = qq"less ${input} | cutadapt - $options->{qual} "; ## If we are keeping quality files
        ## $input_flats = qq"less ${input} | cutadapt - ";
        $out_suffix = 'fastq';
    }
    my $output = qq"${basename}-trimmed_ca.${out_suffix}";
    my $job_string = qq!
mkdir -p ${out_dir} && \\
 ${input_flags} ${type_flag} ${adapter_flags} -e 0.1 -n 3 -m ${minlength} -M ${maxlength} \\
  --too-short-output=${out_dir}/${basename}_tooshort.${out_suffix} \\
  --too-long-output=${out_dir}/${basename}_toolong.${out_suffix} \\
  --untrimmed-output=${out_dir}/${basename}_untrimmed.${out_suffix} \\
  2>outputs/cutadapt.err 1>${output}
!;
    my $cutadapt = $class->Submit(
        comment => $comment,
        input => $input,
        job_name => "cutadapt",
        job_prefix => "07",
        job_string => $job_string,
        prescript => $args{prescript},
        postscript => $args{postscript},
        queue => "workstation",
        output => $output,
        wall => "8:00:00",
    );
    if ($class->{type} eq 'tnseq') {
        my $ta_check = Bio::Adventure::TNSeq::TA_Check(
            $class,
            comment => qq"## Check that TAs exist.",
            input => qq"${output}",
            job_name => "tacheck",
            job_depends => $cutadapt->{pbs_id},
            job_prefix => "08",
        );
    }
    my $comp_short = Bio::Adventure::Compress::Recompress(
        $class,
        comment => qq"## Compressing the tooshort sequences.",
        input => qq"${out_dir}/${basename}_tooshort.fastq",
        job_depends => $cutadapt->{pbs_id},
        job_name => "xzcutshort",
        job_prefix => "08",
        queue => "workstation",
        wall => "04:00:00",
    );
    my $comp_long = Bio::Adventure::Compress::Recompress(
        $class,
        comment => qq"## Compressing the toolong sequences.",
        input => qq"${out_dir}/${basename}_toolong.fastq",
        job_depends => $cutadapt->{pbs_id},
        job_name => "xzcutlong",
        job_prefix => "08",
        queue => "workstation",
        wall => "04:00:00",
    );
    my $comp_un = Bio::Adventure::Compress::Recompress(
        $class,
        comment => qq"## Compressing the toolong sequences.",
        input => qq"${out_dir}/${basename}_untrimmed.fastq",
        job_depends => $cutadapt->{pbs_id},
        job_name => "xzuncut",
        job_prefix => "08",
        queue => "workstation",
        wall => "04:00:00",
    );
    my $comp_original = Bio::Adventure::Compress::Recompress(
        $class,
        comment => qq"## Compressing the original sequence.",
        input => qq"$input",
        job_depends => $cutadapt->{pbs_id},
        job_name => "xzorig",
        job_prefix => "08",
        output => qq"sequences/${input}.xz",
        queue => "workstation",
        wall => "04:00:00",
    );
    return($cutadapt);
}

=item C<Trimomatic>

    $hpgl->Trimomatic(); calls the java tool trimomatic to remove
    adapters/low quality sequences.  If $args{input} has a ':' or ','
    then this will assume the input is comprised of two pairwise files
    and will call 'Trimomatic_Pairwise()', otherwise
    'Trimomatic_Single()'.

=cut
sub Trimomatic {
    my ($class, %args) = @_;
    my $check = which('trimomatic');
    die("Could not find trimomatic in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input',],
    );
    my $trim;
    if ($options->{input} =~ /:|\,/) {
        $trim = Bio::Adventure::RNASeq_Trim::Trimomatic_Pairwise($class, %args);
    } else {
        $trim = Bio::Adventure::RNASeq_Trim::Trimomatic_Single($class, %args);
    }
    return($trim);
}

=item C<Trimomatic_Pairwise>

    $hpgl->Trimomatic_Pairwise(); invokes trimomatic with parameters
    suitable for pairwise sequence libraries.

=cut
sub Trimomatic_Pairwise {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input',],
    );
    my $exe = undef;
    my $found_exe = 0;
    my @exe_list = ('trimomatic', 'TrimmomaticPE', 'trimmomatic');
    for my $test_exe (@exe_list) {
        if (which($test_exe)) {
            $exe = $test_exe;
        }
    }
    if (!defined($exe)) {
        die('Unable to find the trimomatic executable.');
    }
    my $input = $options->{input};
    my @input_list = split(/:|\,/, $input);
    if (scalar(@input_list) <= 1) {
        my $ret = Bio::Adventure::RNASeq_Trim::Trimomatic_Single(
            $class,
            input => $input,
        );
        return($ret);
    }
    my $r1 = $input_list[0];
    my $r2 = $input_list[1];
    my @suff = (".fastq",".gz",".xz");

    my $basename = basename($r1, @suff);
    $basename = basename($basename, @suff);
    $basename =~ s/_R1$//g;
    $basename =~ s/_forward$//g;

    my $r1b = basename($r1, @suff);
    $r1b = basename($r1b, @suff);
    my $r2b = basename($r2, @suff);
    $r2b = basename($r2b, @suff);

    my $r1o = qq"${r1b}-trimmed.fastq.gz";
    my $r1op = qq"${r1b}-trimmed_paired.fastq.gz";
    my $r1ou = qq"${r1b}-trimmed_unpaired.fastq.gz";

    my $r2o = qq"${r2b}-trimmed.fastq.gz";
    my $r2op = qq"${r2b}-trimmed_paired.fastq.gz";
    my $r2ou = qq"${r2b}-trimmed_unpaired.fastq.gz";

    $options = $class->Set_Vars(basename => $basename);
    my $output = qq"${r1o}:${r2o}";
    my $comment = qq!## This call to trimomatic removes illumina and epicentre adapters from ${input}.
## It also performs a sliding window removal of anything with quality <25;
## cutadapt provides an alternative to this tool.
## The original sequence data is recompressed and saved in the sequences/ directory.!;
    my $job_string = qq!
## Trimomatic_Pairwise: In case a trimming needs to be redone...
if [[ \! -r "${r1}" ]]; then
  if [[ -r "sequences/${r1b}.fastq.xz" ]]; then
    mv sequences/${r1b}.fastq.xz . && pxz -d ${r1b}.fastq.xz && pigz ${r1b}.fastq &&\\
       mv sequences/${r2b}.fastq.xz . && pxz -d ${r2b}.fastq.xz && pigz ${r2b}.fastq
  else
    echo "Missing files. Did not find ${r1} nor sequences/${r1b}.fastq.xz"
    exit 1
  fi
fi
## Note that trimomatic prints all output and errors to STDERR, so send both to output
${exe} PE -threads 1 -phred33 ${r1} ${r2} ${r1op} ${r1ou} ${r2op} ${r2ou} \\
    ILLUMINACLIP:$options->{libdir}/adapters.fa:2:20:4 SLIDINGWINDOW:4:25 \\
    1>outputs/${basename}-trimomatic.out 2>&1
excepted=\$(grep "Exception" outputs/${basename}-trimomatic.out)
## The following is in case the illumina clipping fails, which it does if this has already been run I think.
if [[ "\${excepted}" \!= "" ]]; then
  trimomatic PE -threads 1 -phred33 ${r1} ${r2} ${r1op} ${r1ou} ${r2op} ${r2ou} SLIDINGWINDOW:4:25 \\
    1>outputs/${basename}-trimomatic.out 2>&1
fi
sleep 10
mv ${r1op} ${r1o} && mv ${r2op} ${r2o}
!;
    ## Example output from trimomatic:
    ## Input Read Pairs: 10000 Both Surviving: 9061 (90.61%) Forward Only Surviving: 457 (4.57%) Reverse Only Surviving: 194 (1.94%) Dropped: 288 (2.88%)
    ## Perhaps I can pass this along to Get_Stats()
    my $trim = $class->Submit(
        comment => $comment,
        input => $input,
        job_name => "trim",
        job_prefix => "05",
        job_string => $job_string,
        output => $output,
        queue => "workstation",
        prescript => $args{prescript},
        postscript => $args{postscript},
        wall => "12:00:00",
    );
    ## Set the input   for future tools to the output from trimming.
    $options = $class->Set_Vars(input => $output);
    my $trim_stats = Bio::Adventure::RNASeq_Trim::Trimomatic_Stats(
        $class,
        basename => $basename,
        job_prefix => "06",
        job_depends => $trim->{pbs_id},
        pairwise => 1,
    );
    $trim->{stats} = $trim_stats;
    return($trim);
}

=item C<Trimomatic_Single>

    $hpgl->Trimomatic_Single(); invokes trimomatic with parameters
    suitable for single-read sequence libraries.

=cut
sub Trimomatic_Single {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input',],
    );
    my $exe = undef;
    my $found_exe = 0;
    my @exe_list = ('trimomatic', 'TrimmomaticSE', 'trimmomatic');
    for my $test_exe (@exe_list) {
        if (which($test_exe)) {
            $exe = $test_exe;
        }
    }
    if (!defined($exe)) {
        die('Unable to find the trimomatic executable.');
    }

    if ($args{interactive}) {
        print "Run with: cyoa --task rnaseq --method trim --input $options->{input}\n";
    }
    my $input = $options->{input};
    my $basename = $input;
    $basename = basename($basename, (".gz"));
    $basename = basename($basename, (".fastq"));
    my $output = qq"${basename}-trimmed.fastq";
    my $comment = qq!## This call to trimomatic removes illumina and epicentre adapters from ${input}.
## It also performs a sliding window removal of anything with quality <25;
## cutadapt provides an alternative to this tool.
## The original sequence data is recompressed and saved in the sequences/ directory.!;
    my $job_string = qq!
## Trimomatic_Single: In case a trimming needs to be redone...
if [[ \! -r "${input}" ]]; then
  if [[ -r "sequences/${input}.xz" ]]; then
    mv sequences/${input}.xz . && pxz -d ${input}.xz
  else
    echo "Missing files. Did not find ${input} nor sequences/${input}.xz"
    exit 1
  fi
fi
## Note that trimomatic prints all output and errors to STDERR, so send both to output
${exe} SE -phred33 ${input} ${output} ILLUMINACLIP:$options->{libdir}/adapters.fa:2:20:4 \\
    SLIDINGWINDOW:4:25 1>outputs/${basename}-trimomatic.out 2>&1
!;
    my $trim = $class->Submit(
        comment => $comment,
        input => $input,
        job_name => "trim",
        job_prefix => "05",
        job_string => $job_string,
        output => $output,
        prescript => $args{prescript},
        postscript => $args{postscript},
        wall => "4:00:00",
    );
    my $trim_stats = Bio::Adventure::RNASeq_Trim::Trimomatic_Stats(
        $class,
        basename => $basename,
        job_depends => $trim->{pbs_id},
        job_prefix => "06",
    );
    ## Set the input for future tools to the output from this trimming operation.
    $trim->{stats} = $trim_stats;
    $options = $class->Set_Vars(input => $output);
    return($trim);
}

=item C<Trimomatic_Stats>

    Collect the trimming statistics from the output file
    'trimomatic.out' and report them in a .csv file by library.

=cut
sub Trimomatic_Stats {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
    );
    my $basename = $options->{basename};
    my $input_file = "outputs/${basename}-trimomatic.out";
    my $job_depends = $options->{job_depends};
    my $job_name = 'trimst';
    $job_name = $args{job_name} if ($args{job_name});
    my $comment = qq!## This is a stupidly simple job to collect trimomatic statistics!;
    my $stat_output = qq"outputs/trimomatic_stats.csv";
    my $job_string = qq!
if [ \! -r ${stat_output} ]; then
  echo "total_reads,surviving_reads,dropped_reads" > ${stat_output}
fi
total_reads_tmp=\$(grep "^Input Reads" ${input_file} | awk '{print \$3}')
total_reads=\${total_reads_tmp:-0}
surviving_reads_tmp=\$(grep "^Input Reads" ${input_file} | awk '{print \$5}')
surviving_reads=\${surviving_reads_tmp:-0}
dropped_reads_tmp=\$(grep "^Input Reads" ${input_file} | awk '{print \$8}')
dropped_reads=\${dropped_reads_tmp:-0}

stat_string=\$(printf "${basename},%s,%s,%s" "\${total_reads}" "\${surviving_reads}" "\${dropped_reads}")
echo "\$stat_string" >> ${stat_output}
!;
    if ($args{pairwise}) {
        ## The output looks a bit different for pairwise input:
        ## Input Read Pairs: 10000 Both Surviving: 9061 (90.61%) Forward Only Surviving: 457 (4.57%) Reverse Only Surviving: 194 (1.94%) Dropped: 288 (2.88%)
        $job_string = qq!
if [ \! -r ${stat_output} ]; then
  echo "total_reads,surviving_both,surviving_forward,surviving_reverse,dropped_reads" > ${stat_output}
fi
total_reads_tmp=\$(grep "^Input Read Pairs" ${input_file} | awk '{print \$4}')
total_reads=\${total_reads_tmp:-0}
surviving_both_tmp=\$(grep "^Input Read Pairs" ${input_file} | awk '{print \$7}')
surviving_both=\${surviving_both_tmp:-0}
surviving_forward_tmp=\$(grep "^Input Read Pairs" ${input_file} | awk '{print \$12}')
surviving_forward=\${surviving_forward_tmp:-0}
surviving_reverse_tmp=\$(grep "^Input Read Pairs" ${input_file} | awk '{print \$17}')
surviving_reverse=\${surviving_reverse_tmp:-0}
dropped_reads_tmp=\$(grep "^Input Read Pairs" ${input_file} | awk '{print \$20}')
dropped_reads=\${dropped_reads_tmp:-0}

stat_string=\$(printf "${basename},%s,%s,%s,%s,%s" "\${total_reads}" "\${surviving_both}" "\${surviving_forward}" "\${surviving_reverse}" "\${dropped_reads}")
echo "\$stat_string" >> ${stat_output}
!;
    }
    my $stats = $class->Submit(
        comment => $comment,
        cpus => 1,
        input => $input_file,
        job_depends => $job_depends,
        job_name => $job_name,
        job_prefix => $args{job_prefix},
        job_string => $job_string,
        mem => 1,
        queue => "throughput",
        wall => "00:10:00",
    );
    return($stats);
}

=back

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

    L<Bio::Tools::Run::StandAloneBlast>

=cut

1;
