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
use File::Path qw"make_path";
use File::ShareDir qw":ALL";
use File::Temp qw":POSIX";
use File::Which qw"which";

=head1 NAME

 Bio::Adventure::Trim - Use trimomatic/cutadapt/etc to trim libraries

=head1 SYNOPSIS

 This file is responsible for invoking the various sequence trimmers.

=head1 METHODS

=head2 C<Cutadapt>

 Invoke cutadapt on a pile of sequence.
 10.14806/ej.17.1.200

 Use biopieces/cutadapt to attempt to remove sequence adapters from a library.
 This is most common used by me for ribosome profiling libraries.

=cut
sub Cutadapt {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        arbitrary => undef,
        maxlength => 42,
        maxerr => 0.1,
        maxremoved => 3,
        minlength => 8,
        modules => ['cutadapt'],
        left => undef,
        right => undef,
        either => undef,
        type => 'rnaseq',
        jmem => 12,
        jwalltime => '48:00:00',
        jprefix => '12',);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $job_name = $class->Get_Job_Name();
    my $inputs = $class->Get_Paths($options->{input});
    my $check = which('cutadapt');
    die('Could not find cutadapt in your PATH.') unless($check);
    my $minlength = $options->{minlength};
    my $maxlength = $options->{maxlength};
    my $arbitrary = '';
    if (defined($options->{arbitrary})) {
        $arbitrary = $options->{arbitrary};
    }

    my $input = $options->{input};
    my @suffixes = split(/\,/, $options->{suffixes});
    my $basename = basename($input, @suffixes);
    $basename = basename($basename, @suffixes);
    my $adapter_flags = "";

    if ($options->{type} eq 'old_tnseq') {
        $adapter_flags = ' -a ACAGGTTGGATGATAAGTCCCCGGTCTGACACATC -a ACAGTCCCCGGTCTGACACATCTCCCTAT -a ACAGTCCNCGGTCTGACACATCTCCCTAT ';
        $maxlength = 20;
    } elsif ($options->{type} eq 'riboseq') {
        $adapter_flags = ' -a ACAGGTTGGATGATAAGTCCCCGGTCTGACACATCTCCCTAT -a AGATCGGAAGAGCACACGTCTGAAC -b AGATCGGAAGAGCACACGTCTGAAC ';
        $minlength = 16;
        $maxlength = 42;
    } else {
        $minlength = 30;
        $maxlength = undef;
        $adapter_flags = ' -a ACAGGTTGGATGATAAGTCCCCGGTCTGACACATCTCCCTAT -a AGATCGGAAGAGCACACGTCTGAAC -b AGATCGGAAGAGCACACGTCTGAAC ';
    }

    my $comment = qq!## This script makes use of biopieces and cutadapt to trim away adapters
## and separate the sequence file into a few pieces depending on size
## and adapter status.  It also performs some simple graphs of the data.!;
    my $out_dir = qq"outputs/$options->{jprefix}cutadapt";
    my $stderr = qq"${out_dir}/cutadapt.stderr";
    my $stdout = qq"${out_dir}/cutadapt.stdout";
    my $type_flag = '';
    my $out_suffix = 'fastq';
    my $input_flags = qq"less ${input} | cutadapt - ";
    if ($input =~ /\.csfasta/) {
        $type_flag = '-c -t --strip-f3';
        $options = $class->Get_Vars(
            args => \%args,
            required => ['qual'],
        );
        $input_flags = qq"less ${input} | cutadapt - $options->{qual} "; ## If we are keeping quality files
        $out_suffix = 'fastq';
    }
    my $minarg = '';
    my $too_short = undef;
    if ($minlength) {
        $too_short = qq"${out_dir}/${basename}_too_short.${out_suffix}";
        $minarg = qq" -m ${minlength} --too-short-output=${too_short}";
    }
    my $maxarg = '';
    my $too_long = undef;
    if ($maxlength) {
        $too_long = qq"${out_dir}/${basename}_too_long.${out_suffix}";
        $maxarg = qq" -M ${maxlength} --too-long-output=${too_long}";
    }
    my $maxerr = '';
    if ($options->{maxerr}) {
        $maxerr = qq" -e $options->{maxerr}";
    }
    my $output = qq"${out_dir}/${basename}-trimmed_ca.${out_suffix}";
    my $compressed_out = qq"${output}.xz";
    my $jstring = qq!## I am putting the args in an odd order in case some are not defined.
mkdir -p ${out_dir}
${input_flags} \\
  ${type_flag} ${adapter_flags} ${arbitrary} \\
  -o ${output} ${minarg} \\
  -n $options->{maxremoved} ${maxerr} ${maxarg} \\
  --untrimmed-output=${out_dir}/${basename}_untrimmed.${out_suffix} \\
  2>${stderr} \\
  1>${stdout}
xz -9e -f ${output}
xz -9e -f ${out_dir}/${basename}_untrimmed.${out_suffix}
!;
    if (defined($too_short)) {
        $jstring .= qq!
xz -9e -f ${too_short}
!;
    }
    if (defined($too_long)) {
        $jstring .= qq!
xz -9e -f ${too_long}
!;
    }

    my $cutadapt = $class->Submit(
        comment => $comment,
        input => $input,
        jmem => $options->{jmem},
        jname => qq"cutadapt_${job_name}",
        jprefix => $options->{jprefix},
        jqueue => 'workstation',
        jstring => $jstring,
        jwalltime => '8:00:00',
        modules => $options->{modules},
        stderr => $stderr,
        stdout => $stdout,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        output => $compressed_out,);
    if ($options->{type} eq 'old_tnseq') {
        my $ta_check = $class->Bio::Adventure::TNSeq::TA_Check(
            comment => '## Check that TAs exist.',
            input => $compressed_out,
            jdepends => $cutadapt->{job_id},
            jname => qq"tacheck_${job_name}",
            jprefix => $options->{jprefix} + 1,);
        $cutadapt->{tacheck} = $ta_check;
    }
    return($cutadapt);
}

=head2 C<Racer>

 Use the RACER command from hitec to correct sequencer-based errors.
 10.1093/bioinformatics/btq653

=cut
sub Racer {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        compress => 1,
        length => 1000000,
        jmem => 24,
        jprefix => '10',
        jwalltime => '40:00:00',
        modules => ['hitec'],
        required => ['input', ],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $job_name = $class->Get_Job_Name();
    my $input = $options->{input};
    my @input_list = split(/:|\,/, $input);
    my @suffixes = ('.fastq', '.gz', '.xz');
    my $decompress_input = 0;
    my @base_list = ();
    for my $in (@input_list) {
        $decompress_input = 1 unless ($in =~ /\.fastq$/);
        my $shorted = basename($in, @suffixes);
        push(@base_list, basename($shorted, @suffixes));
    }
    my $comment = qq!## This calls RACER to try to remove
## arbitrary errors in sequencing data.
!;
    my $jstring = qq'';
    my $output_dir = qq"outputs/$options->{jprefix}racer";
    my @created = make_path($output_dir);
    my @output_files;
    my $stdout = qq"${output_dir}/racer.stdout";
    my $stderr = qq"${output_dir}/racer.stderr";
    foreach my $c (0 .. $#input_list) {
        my $name = File::Temp::tempnam($output_dir, 'racer');
        my $output = qq"${output_dir}/$base_list[$c]-corrected.fastq";
        if ($decompress_input) {
            $jstring .= qq"less $input_list[$c] > ${name}.fastq\n";
        } else {
            $jstring .= qq"ln -sf \$(pwd)/$input_list[$c] ${name}.fastq\n";
        }
        $jstring .= qq"RACER \\
  ${name}.fastq \\
  ${output} \\
  $options->{length} \\
  1>>${stdout} \\
  2>>${stderr}
rm ${name}.fastq
echo \"Finished correction of $input_list[$c].\" >> ${stdout}
";
        if ($options->{compress}) {
            $jstring .= qq"xz -9e -f ${output}
";
            push(@output_files, "${output}.xz");
        } else {
            push(@output_files, "${output}");
        }
    }

    my $output = '';
    for my $o (@output_files) {
        $output .= qq"${o}:";
    }
    $output =~ s/\:$//g;

    my $racer = $class->Submit(
        comment => $comment,
        input => $input,
        jcpu => 4,
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
        jname => "racer_${job_name}",
        jprefix => $options->{jprefix},
        jqueue => 'workstation',
        jstring => $jstring,
        jwalltime => '12:00:00',
        modules => $options->{modules},
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        output => $output,
        stdout => $stdout,
        stderr => $stderr);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($racer);
}

=head2 C<Trimomatic>

 Invoke trimomatic to get sequencing data ready to play with.
 10.1093/bioinformatics/btu170

 Call trimomatic to remove adapters/low quality sequences.  If $args{input} has a
 ':' or ',' then this will assume the input is comprised of two pairwise files
 and will call 'Trimomatic_Pairwise()', otherwise 'Trimomatic_Single()'.

=cut
sub Trimomatic {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        compress => 1,
        jcpu => 4,
        jmem => 24,
        jprefix => '01',
        jwalltime => '36:00:00',
        length => 50,
        modules => ['trimomatic'],
        quality => '20',
        required => ['input',],);
    my $trim;
    if ($options->{input} =~ /:|\,/) {
        $trim = $class->Bio::Adventure::Trim::Trimomatic_Pairwise(%{$options});
    } else {
        $trim = $class->Bio::Adventure::Trim::Trimomatic_Single(%{$options});
    }
    return($trim);
}

=head2 C<Trimomatic_Pairwise>

 Invoke trimomatic with parameters suitable for pairwise sequence libraries.

=cut
sub Trimomatic_Pairwise {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        compress => 1,
        jcpu => 4,
        jmem => 24,
        jprefix => '01',
        jwalltime => '48:00:00',
        length => 50,
        modules => ['trimomatic'],
        quality => '20',
        required => ['input',],);
    my $loaded = $class->Module_Loader(modules => $options->{modules}, verbose => 1);
    my $output_dir = qq"outputs/$options->{jprefix}trimomatic";
    my $job_name = $class->Get_Job_Name();
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

    my $adapter_file = dist_file('Bio-Adventure', 'genome/adapters.fa');
    my $input = $options->{input};
    my @input_list = split(/:|\,/, $input);
    if (scalar(@input_list) <= 1) {
        my $ret = $class->Bio::Adventure::Trim::Trimomatic_Single(input => $input,);
        return($ret);
    }
    my $r1 = $input_list[0];
    my $r2 = $input_list[1];
    my @suff = ('.fastq', '.gz', '.xz');

    my $basename = basename($r1, @suff);
    $basename = basename($basename, @suff);
    $basename =~ s/_R1$//g;
    $basename =~ s/_forward$//g;

    my $r1b = basename($r1, @suff);
    $r1b = basename($r1b, @suff);
    my $r2b = basename($r2, @suff);
    $r2b = basename($r2b, @suff);
    my $reader = qq"${r1} ${r2}";
    if ($r1 =~ /\.fastq\.xz$/) {
        $reader = qq"<(less ${r1}) <(less ${r2})";
    }
    my $r1o = qq"${r1b}-trimmed.fastq";
    my $r1op = qq"${r1b}-trimmed_paired.fastq";
    my $r1ou = qq"${r1b}-trimmed_unpaired.fastq";

    my $r2o = qq"${r2b}-trimmed.fastq";
    my $r2op = qq"${r2b}-trimmed_paired.fastq";
    my $r2ou = qq"${r2b}-trimmed_unpaired.fastq";

    my $leader_trim = '';
    if ($options->{task} eq 'dnaseq') {
        $leader_trim = 'HEADCROP:20 LEADING:3 TRAILING:3';
    }
    my $suffix_trim = '';
    if (defined($options->{arbitrary}) &&
        $options->{arbitrary} =~ /^CROP/) {
        $suffix_trim = $options->{arbitrary};
    }

    my $output = qq"${r1o}:${r2o}";
    my $output_unpaired = qq"${r1ou}:${r2ou}";
    my $compress_string = '';
    if ($options->{compress}) {
        $output = qq"${r1o}.xz:${r2o}.xz";
        $output_unpaired = qq"${r1ou}.xz:${r2ou}.xz";
        $compress_string = qq"
## Recompress the unpaired reads, this should not take long.
xz -9e -f ${r1ou}
xz -9e -f ${r2ou}
## Recompress the paired reads.
xz -9e -f ${r1o}
xz -9e -f ${r2o}
ln -sf ${r1o}.xz r1_trimmed.fastq.xz
ln -sf ${r2o}.xz r2_trimmed.fastq.xz
";
    }

    my $comment = qq!## This call to trimomatic removes illumina and epicentre adapters from ${input}.
## It also performs a sliding window removal of anything with quality <25;
## cutadapt provides an alternative to this tool.
## The original sequence data is recompressed and saved in the sequences/ directory.!;
    my $stdout = qq"${output_dir}/${basename}-trimomatic.stdout";
    my $stderr = qq"${output_dir}/${basename}-trimomatic.stderr";
    my $jstring = qq!mkdir -p ${output_dir}
## Note that trimomatic prints all output and errors to STDERR, so send both to output
${exe} \\
  -threads 1 \\
  -phred33 \\
  ${reader} \\
  ${r1op} ${r1ou} \\
  ${r2op} ${r2ou} \\
  ${leader_trim} ILLUMINACLIP:${adapter_file}:2:$options->{quality}:10:2:keepBothReads ${suffix_trim} \\
  SLIDINGWINDOW:4:$options->{quality} MINLEN:$options->{length} \\
  1>${stdout} \\
  2>${stderr}
excepted=\$( { grep "Exception" "${output_dir}/${basename}-trimomatic.stdout" || test \$? = 1; } )
## The following is in case the illumina clipping fails, which it does if this has already been run I think.
if [[ "\${excepted}" \!= "" ]]; then
  ${exe} \\
    -threads 1 \\
    -phred33 \\
    ${reader} \\
    ${r1op} ${r1ou} \\
    ${r2op} ${r2ou} \\
    ${leader_trim} SLIDINGWINDOW:4:25 MINLEN:$options->{length} \\
    1>${stdout} \\
    2>${stderr}
fi
sleep 10
mv ${r1op} ${r1o}
mv ${r2op} ${r2o}
!;

    if ($options->{compress}) {
        $jstring .= qq"${compress_string}\n";
    }

    ## Example output from trimomatic:
    ## Input Read Pairs: 10000 Both Surviving: 9061 (90.61%) Forward Only Surviving: 457 (4.57%) Reverse Only Surviving: 194 (1.94%) Dropped: 288 (2.88%)
    ## Perhaps I can pass this along to Get_Stats()
    $loaded = $class->Module_Loader(modules => $options->{modules});
    print "-------------------
TESTME: About to submit first trim job.
-----------------\n";
    my $trim = $class->Submit(
        args => \%args,
        comment => $comment,
        input => $input,
        jcpu => $options->{jcpu},
        jmem => $options->{jmem},
        jname => qq"trim_${job_name}",
        jprefix => $options->{jprefix},
        jqueue => 'workstation',
        jstring => $jstring,
        jwalltime => $options->{jwalltime},
        length => $options->{length},
        modules => $options->{modules},
        output => $output,
        output_unpaired => $output_unpaired,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        stdout => $stdout,
        stderr => $stderr);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    my $new_prefix = qq"$options->{jprefix}_1";
    my $trim_stats = $class->Bio::Adventure::Metadata::Trimomatic_Stats(
        basename => $basename,
        jdepends => $trim->{job_id},
        jcpu => 1,
        jprefix => $new_prefix,
        jname => "trst_${job_name}",
        jwalltime => '00:03:00',
        pairwise => 1,
        input => $stderr,
        output_dir => $output_dir,);
    $trim->{stats} = $trim_stats;
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
        compress => 1,
        jmem => 24,
        jcpu => 4,
        jprefix => '01',
        jwalltime => '48:00:00',
        length => 50,
        modules => ['trimomatic'],
        quality => '20',
        required => ['input',],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
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
    my $output_dir = qq"outputs/$options->{jprefix}trimomatic";
    if (!defined($exe)) {
        die('Unable to find the trimomatic executable.');
    }
    my $adapter_file = dist_file('Bio-Adventure', 'genome/adapters.fa');
    my $leader_trim = "";
    if (defined($options->{task}) && $options->{task} eq 'dnaseq') {
        $leader_trim = 'HEADCROP:20 LEADING:3 TRAILING:3';
    }

    my $input = $options->{input};
    my $basename = $input;
    $basename = basename($basename, ('.gz'));
    $basename = basename($basename, ('.fastq'));
    my $job_name = $class->Get_Job_Name();
    my $output = qq"${basename}-trimmed.fastq";
    my $comment = qq!## This call to trimomatic removes illumina and epicentre adapters from ${input}.
## It also performs a sliding window removal of anything with quality <25;
## cutadapt provides an alternative to this tool.
## The original sequence data is recompressed and saved in the sequences/ directory.!;
    my $stdout = qq"${output_dir}/${basename}-trimomatic.stdout";
    my $stderr = qq"${output_dir}/${basename}-trimomatic.stderr";
    my $jstring = qq!mkdir -p ${output_dir}
## Note that trimomatic prints all output and errors to STDERR, so send both to output
${exe} \\
  -phred33 \\
  <(less ${input}) \\
  ${output} \\
  ${leader_trim} ILLUMINACLIP:${adapter_file}:2:30:10 \\
  SLIDINGWINDOW:4:25 MINLEN:$options->{length} \\
  1>${stdout} 2>${stderr}
xz -9e -f ${output}
ln -sf ${output}.xz r1_trimmed.fastq.xz
!;
    my $trim = $class->Submit(
        comment => $comment,
        input => $input,
        jcpu => $options->{jcpu},
        jmem => $options->{jmem},
        jname => qq"trim_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jwalltime => $options->{jwalltime},
        length => $options->{length},
        modules => $options->{modules},
        output => $output,
        stderr => $stderr,
        stdout => $stdout,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    my $trim_stats = $class->Bio::Adventure::Metadata::Trimomatic_Stats(
        basename => $basename,
        input => $stderr,
        jcpu => 1,
        jdepends => $trim->{job_id},
        jname => qq"trst_${job_name}",
        jprefix => '06',
        jwalltime => '00:03:00',
        output_dir => $output_dir,
    );
    $trim->{stats} = $trim_stats;
    return($trim);
}

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

L<trimomatic> L<cutadapt>

=cut

1;
