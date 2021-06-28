package Bio::Adventure::Trim;
## LICENSE: gplv2
## ABSTRACT:  Kitty!
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use File::Basename;
use File::ShareDir ':ALL';
use File::Temp qw":POSIX";
use File::Which qw"which";

=head1 NAME

Bio::Adventure::Trim - Use trimomatic/cutadapt/etc to trim libraries

=head1 SYNOPSIS

This file is responsible for invoking the various sequence trimmers.

=head1 METHODS

=head2 C<Cutadapt>

Use biopieces/cutadapt to attempt to remove sequence adapters from a library.
This is most common used by me for ribosome profiling libraries.

=cut
sub Cutadapt {
    my ($class, %args) = @_;
    my $check = which('cutadapt');
    die('Could not find cutadapt in your PATH.') unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'type'],
        minlength => 8,
        maxlength => 42,
        maxerr => 0.1,
        maxremoved => 3,
	arbitrary => undef,
	left => undef,
	right => undef,
	either => undef,
    );
    my %start_options = %{$options};
    my $job_basename = $class->Get_Job_Name();

    if ($options->{task}) {
        $options->{type} = $options->{task};
    }
    $options->{type} = lc($options->{type});
    my $type = 'rnaseq';
    if ($options->{type}) {
        $type = $options->{type};
    } else {
        $type = 'tnseq';
    }

    my $arbitrary = '';
    if (defined($options->{arbitrary})) {
        $arbitrary = $options->{arbitrary};
    }

    my $input = $options->{input};
    my $basename = basename($input, @{$options->{suffixes}});
    $basename = basename($basename, @{$options->{suffixes}});
    my $adapter_flags = "";
    if ($type eq 'tnseq') {
        $adapter_flags = qq" -a ACAGGTTGGATGATAAGTCCCCGGTCTGACACATC -a ACAGTCCCCGGTCTGACACATCTCCCTAT -a ACAGTCCNCGGTCTGACACATCTCCCTAT ";
        $options->{maxlength} = 20;
    } elsif ($type eq 'riboseq') {
        $adapter_flags = qq" -a ACAGGTTGGATGATAAGTCCCCGGTCTGACACATCTCCCTAT -a AGATCGGAAGAGCACACGTCTGAAC -b AGATCGGAAGAGCACACGTCTGAAC ";
        $options->{minlength} = 16;
        if ($input =~ /\.csfasta/) {
            $adapter_flags = qq" -a CGCCTTGGCCGTACAGCAGCATATTGGATAAGAGAATGAGGAACCCGGGGCAG -a GCGGAACCGGCATGTCGTCGGGCATAACCCTCTCTTACTCCTTGGGCCCCGTC ";
        }
    } else {
        $adapter_flags = qq" -a ACAGGTTGGATGATAAGTCCCCGGTCTGACACATCTCCCTAT -a AGATCGGAAGAGCACACGTCTGAAC -b AGATCGGAAGAGCACACGTCTGAAC ";
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
        $input_flags = qq"less ${input} | cutadapt - $options->{qual} "; ## If we are keeping quality files
        $out_suffix = 'fastq';
    }
    my $output = qq"${basename}-trimmed_ca.${out_suffix}";
    my $jstring = qq!
mkdir -p ${out_dir} && \\
 ${input_flags} \\
  ${type_flag} ${adapter_flags} ${arbitrary} \\
  -e $options->{maxerr} -n $options->{maxremoved} \\
  -m $options->{minlength} -M $options->{maxlength} \\
  --too-short-output=${out_dir}/${basename}_tooshort.${out_suffix} \\
  --too-long-output=${out_dir}/${basename}_toolong.${out_suffix} \\
  --untrimmed-output=${out_dir}/${basename}_untrimmed.${out_suffix} \\
  2>outputs/cutadapt.err 1>${output}
!;
    my $cutadapt = $class->Submit(
        comment => $comment,
        input => $input,
        jname => qq"ca_${job_basename}",
        jprefix => "07",
        jstring => $jstring,
        prescript => $args{prescript},
        postscript => $args{postscript},
        queue => "workstation",
        output => $output,
        walltime => "8:00:00",
    );
    if ($type eq 'tnseq') {
        my $ta_check = Bio::Adventure::TNSeq::TA_Check(
            $class,
            comment => qq"## Check that TAs exist.",
            input => qq"${output}",
            jname => qq"tach_${job_basename}",
            depends => $cutadapt->{job_id},
            jprefix => "08",
        );
    }
    my $comp_short = Bio::Adventure::Compress::Recompress(
        $class,
        comment => qq"## Compressing the tooshort sequences.",
        xz_input => qq"${out_dir}/${basename}_tooshort.fastq",
        depends => $cutadapt->{job_id},
        jname => "xzcutshort_${job_basename}",
        jprefix => "08",
        queue => "workstation",
        walltime => "04:00:00",
    );
    my $comp_long = Bio::Adventure::Compress::Recompress(
        $class,
        comment => qq"## Compressing the toolong sequences.",
        xz_input => qq"${out_dir}/${basename}_toolong.fastq",
        depends => $cutadapt->{job_id},
        jname => "xzcutlong_${job_basename}",
        jprefix => "08",
        queue => "workstation",
        walltime => "04:00:00",
    );
    my $comp_un = Bio::Adventure::Compress::Recompress(
        $class,
        comment => qq"## Compressing the toolong sequences.",
        xz_input => qq"${out_dir}/${basename}_untrimmed.fastq",
        depends => $cutadapt->{job_id},
        jname => "xzuncut_${job_basename}",
        jprefix => "08",
        queue => "workstation",
        walltime => "04:00:00",
    );
    my $comp_original = Bio::Adventure::Compress::Recompress(
        $class,
        comment => qq"## Compressing the original sequence.",
        xz_input => $input,
        depends => $cutadapt->{job_id},
        jname => "xzorig_${job_basename}",
        jprefix => "08",
        queue => "workstation",
        walltime => "04:00:00",
    );
    my $comp_output = Bio::Adventure::Compress::Recompress(
        $class,
        comment => qq"## Compressing the output sequence.",
        xz_input => $output,
        depends => $cutadapt->{job_id},
        jname => "xzout_${job_basename}",
        jprefix => "08",
        queue => "workstation",
        walltime => "04:00:00",
    );

    $class->{options} = \%start_options;
    return($cutadapt);
}

sub Racer {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        length => 1000000,
        required => ['input',],
    );
    my %start_options = %{$options};
    my $job_basename = $class->Get_Job_Name();
    my $input = $options->{input};
    my @input_list = split(/:|\,/, $input);
    my @suffixes = (".fastq",".gz",".xz");
    my @base_list = ();
    for my $in (@input_list) {
        my $shorted = basename($in, @suffixes);
        push(@base_list, basename($shorted, @suffixes));
    }
    my $comment = qq!## This calls RACER to try to remove
## arbitrary errors in sequencing data.
!;
    my $jstring = qq"";

    foreach my $c (0 .. $#input_list) {
        my $name = File::Temp::tempnam('.', 'racer');
        my $output = qq"$base_list[$c]-corrected.fastq";
        $jstring .= "zcat $input_list[$c] > ${name}.fastq &&
  RACER \\
  ${name}.fastq \\
  ${output} \\
  $options->{length} \\
  2>outputs/racer.out 1>&2 &&
  rm ${name}.fastq
gzip -9 ${output}
";
    }

    my $racer = $class->Submit(
        comment => $comment,
        input => $input,
        jname => "racer_${job_basename}",
        jprefix => "07",
        jstring => $jstring,
        queue => "workstation",
        cpus => 4,
        mem => 30,
        prescript => $args{prescript},
        postscript => $args{postscript},
        walltime => "12:00:00",
    );
    ## Set the input   for future tools to the output from trimming.
    return($racer);
}

=head2 C<Trimomatic>

Call trimomatic to remove adapters/low quality sequences.  If $args{input} has a
':' or ',' then this will assume the input is comprised of two pairwise files
and will call 'Trimomatic_Pairwise()', otherwise 'Trimomatic_Single()'.

=cut
sub Trimomatic {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input',],
    );
    my %start_options = %{$options};
    my $trim;
    if ($options->{input} =~ /:|\,/) {
        $trim = Bio::Adventure::Trim::Trimomatic_Pairwise($class, %args);
    } else {
        $trim = Bio::Adventure::Trim::Trimomatic_Single($class, %args);
    }
    $class->{options} = \%start_options;
    return($trim);
}

=head2 C<Trimomatic_Pairwise>

Invoke trimomatic with parameters suitable for pairwise sequence libraries.

=cut
sub Trimomatic_Pairwise {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input',],
        quality => '20',
    );
    my %start_options = %{$options};
    my $job_basename = $class->Get_Job_Name();
    my $exe = undef;
    my $found_exe = 0;
    my @exe_list = ('trimomatic PE', 'TrimmomaticPE', 'trimmomatic PE');
    for my $test_exe (@exe_list) {
        my @executable_list = split(/\s+/, $test_exe);
        my $executable = $executable_list[0];
        if (which($executable)) {
            $exe = $test_exe;
        }
    }
    if (!defined($exe)) {
        die('Unable to find the trimomatic executable.');
    }

    my $adapter_file = module_file('Bio::Adventure', 'genome/adapters.fa');
    my $input = $options->{input};
    my @input_list = split(/:|\,/, $input);
    if (scalar(@input_list) <= 1) {
        my $ret = Bio::Adventure::Trim::Trimomatic_Single(
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
    my $reader = qq"${r1} ${r2}";
    if ($r1 =~ /\.xz$/) {
        $reader = qq"<(less ${r1}) <(less ${r2})";
    }
    my $r1o = qq"${r1b}-trimmed.fastq.gz";
    my $r1op = qq"${r1b}-trimmed_paired.fastq.gz";
    my $r1ou = qq"${r1b}-trimmed_unpaired.fastq.gz";

    my $r2o = qq"${r2b}-trimmed.fastq.gz";
    my $r2op = qq"${r2b}-trimmed_paired.fastq.gz";
    my $r2ou = qq"${r2b}-trimmed_unpaired.fastq.gz";

    my $leader_trim = "";
    if ($options->{task} eq 'dnaseq') {
        $leader_trim = 'HEADCROP:20 LEADING:3 TRAILING:3';
    }

    $options = $class->Set_Vars(basename => $basename);
    my $output = qq"${r1o}:${r2o}";
    my $comment = qq!## This call to trimomatic removes illumina and epicentre adapters from ${input}.
## It also performs a sliding window removal of anything with quality <25;
## cutadapt provides an alternative to this tool.
## The original sequence data is recompressed and saved in the sequences/ directory.!;
    my $jstring = qq!
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
${exe} \\
  -threads 1 \\
  -phred33 \\
  ${reader} \\
  ${r1op} ${r1ou} \\
  ${r2op} ${r2ou} \\
  ${leader_trim} ILLUMINACLIP:${adapter_file}:2:$options->{quality}:10:2:keepBothReads \\
  SLIDINGWINDOW:4:$options->{quality} MINLEN:40 \\
  1>outputs/${basename}-trimomatic.out 2>&1
excepted=\$(grep "Exception" outputs/${basename}-trimomatic.out)
## The following is in case the illumina clipping fails, which it does if this has already been run I think.
if [[ "\${excepted}" \!= "" ]]; then
  ${exe} \\
    -threads 1 \\
    -phred33 \\
    ${reader} \\
    ${r1op} ${r1ou} \\
    ${r2op} ${r2ou} \\
    ${leader_trim} SLIDINGWINDOW:4:25 MINLEN:50\\
    1>outputs/${basename}-trimomatic.out 2>&1
fi
sleep 10
mv ${r1op} ${r1o} && mv ${r2op} ${r2o}
ln -s ${r1o} r1_trimmed.fastq.gz
ln -s ${r2o} r2_trimmed.fastq.gz
!;
    ## Example output from trimomatic:
    ## Input Read Pairs: 10000 Both Surviving: 9061 (90.61%) Forward Only Surviving: 457 (4.57%) Reverse Only Surviving: 194 (1.94%) Dropped: 288 (2.88%)
    ## Perhaps I can pass this along to Get_Stats()
    my $trim = $class->Submit(
        comment => $comment,
        input => $input,
        jname => "trim_${job_basename}",
        jprefix => "05",
        jstring => $jstring,
        output => $output,
        queue => "workstation",
        cpus => 3,
        mem => 40,
        prescript => $args{prescript},
        postscript => $args{postscript},
        walltime => "12:00:00",
    );
    ## Set the input   for future tools to the output from trimming.
    $options = $class->Set_Vars(input => $output);
    my $trim_stats = Bio::Adventure::Trim::Trimomatic_Stats(
        $class,
        basename => $basename,
        jprefix => "06",
        jname => "trst_${job_basename}",
        depends => $trim->{job_id},
        pairwise => 1,
    );
    $trim->{stats} = $trim_stats;
    $class->{options} = \%start_options;
    return($trim);
}

=head2 C<Trimomatic_Single>

Invoke trimomatic with parameters suitable for single-read sequence libraries.

=cut
sub Trimomatic_Single {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input',],
    );
    my %start_options = %{$options};
    my $exe = undef;
    my $found_exe = 0;
    my @exe_list = ('trimomatic SE', 'TrimmomaticSE', 'trimmomatic SE');
    for my $test_exe (@exe_list) {
        my @executable_list = split(/\s+/, $test_exe);
        my $executable = $executable_list[0];
        if (which($executable)) {
            $exe = $test_exe;
        }
    }
    if (!defined($exe)) {
        die('Unable to find the trimomatic executable.');
    }
    my $adapter_file = module_file('Bio::Adventure', 'genome/adapters.fa');

    if ($args{interactive}) {
        print "Run with: cyoa --task rnaseq --method trim --input $options->{input}\n";
    }

    my $leader_trim = "";
    if ($options->{task} eq 'dnaseq') {
        $leader_trim = 'HEADCROP:20 LEADING:3 TRAILING:3';
    }

    my $input = $options->{input};
    my $basename = $input;
    $basename = basename($basename, (".gz"));
    $basename = basename($basename, (".fastq"));
    my $job_basename = $class->Get_Job_Name();
    my $output = qq"${basename}-trimmed.fastq.gz";
    my $comment = qq!## This call to trimomatic removes illumina and epicentre adapters from ${input}.
## It also performs a sliding window removal of anything with quality <25;
## cutadapt provides an alternative to this tool.
## The original sequence data is recompressed and saved in the sequences/ directory.!;
    my $jstring = qq!
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
${exe} \\
  -phred33 \\
  ${input} \\
  ${output} \\
  ${leader_trim} ILLUMINACLIP:${adapter_file}:2:30:10 \\
  SLIDINGWINDOW:4:25 MINLEN:50 \\
  1>outputs/${basename}-trimomatic.out 2>&1
!;
    my $trim = $class->Submit(
        comment => $comment,
        input => $input,
        jname => "trim_${job_basename}",
        jprefix => "05",
        jstring => $jstring,
        output => $output,
        prescript => $args{prescript},
        postscript => $args{postscript},
        walltime => "4:00:00",
    );
    my $trim_stats = Bio::Adventure::Trim::Trimomatic_Stats(
        $class,
        basename => $basename,
        depends => $trim->{job_id},
        jname => "trst_${job_basename}",
        jprefix => "06",
    );
    $trim->{stats} = $trim_stats;
    $class->{options} = \%start_options;
    return($trim);
}

=head2 C<Trimomatic_Stats>

Collect the trimming statistics from the output file 'trimomatic.out' and report
them in a .csv file by library.

=cut
sub Trimomatic_Stats {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
    );
    ## Dereferencing the options to keep a safe copy.
    my %start_options = %{$options};
    my $basename = $options->{basename};
    my $input_file = "outputs/${basename}-trimomatic.out";
    my $depends = $options->{depends};
    my $jname = 'trimst';
    $jname = $args{jname} if ($args{jname});
    my $comment = qq!## This is a stupidly simple job to collect trimomatic statistics!;
    my $stat_output = qq"outputs/trimomatic_stats.csv";
    my $jstring = qq!
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
        $jstring = qq!
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
        depends => $depends,
        jname => $jname,
        jprefix => $args{jprefix},
        jstring => $jstring,
        mem => 1,
        queue => "throughput",
        walltime => "00:10:00",
    );
    $class->{options} = \%start_options;
    return($stats);
}

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

L<trimomatic> L<cutadapt>

=cut

1;
