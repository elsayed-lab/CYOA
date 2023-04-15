package Bio::Adventure::QA;
## LICENSE: gplv2
## ABSTRACT:  Kitty!
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use feature 'try';
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use File::Basename;
use File::Which qw"which";

=head1 NAME

 Bio::Adventure::QA - Use trimomatic/cutadapt/etc to trim libraries

=head1 SYNOPSIS

 The functions here invoke various quality assurance programs.

=head1 METHODS

=head2 C<Biopieces_Graph>

 Biopieces is an older but fun toolkit.
 https://github.com/maasha/biopieces

 Reads in a fastq file and uses biopieces to make some graphs describing the
 sequences therein.

=cut
sub Biopieces_Graph {
    my ($class, %args) = @_;
    my $check = which('read_fastq');
    die("Could not find biopieces in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jmem => 8,
        jprefix => '02',
        modules => ['biopieces'],);
    my $input = $options->{input};
    my $jname = $class->Get_Job_Name();
    $jname = qq"biop_${jname}";
    my @inputs = split(/\,|\:|\;/, $input);
    my $comment = '## This script uses biopieces to draw some simple graphs of the sequence.';
    my $output_dir = 'output/biopieces';
    my $stdout = qq"${output_dir}/biopieces.stdout";
    my $stderr = qq"${output_dir}/biopieces.stderr";
    my $bp;
    if (scalar(@inputs) > 1) { ## multiple comma/colon/semicolon inputs were provided.
        foreach my $in (@inputs) {
            my $short_in = basename($in, ('.gz','.fastq'));
            $short_in = basename($short_in, ('.gz','.fastq'));
            my $jstring = qq!
## Do not forget that _only_ the last command in a biopieces string is allowed to have the -x.
mkdir -p outputs/biopieces
less ${in} | read_fastq -i - -e base_$options->{phred} |\\
 plot_scores -T 'Quality Scores' -t svg -o outputs/biopieces/${short_in}_quality_scores.svg |\\
 plot_nucleotide_distribution -T 'NT. Distribution' -t svg -o outputs/biopieces/${short_in}_ntdist.svg |\\
 plot_lendist -T 'Length Distribution' -k SEQ_LEN -t svg -o outputs/biopieces/${short_in}_lendist.svg |\\
 analyze_gc |\\
     bin_vals -b 20 -k 'GC%' |\\
     plot_distribution -k 'GC%_BIN' -t svg -o outputs/biopieces/${short_in}_gcdist.svg |\\
 analyze_gc |\\
     mean_vals -k 'GC%' -o outputs/biopieces/${short_in}_gc.txt |\\
 count_records -o outputs/biopieces/${short_in}_count.txt -x
!;
            $bp = $class->Submit(
                input => $in,
                comment => $comment,
                jdepends => $options->{jdepends},
                jmem => $options->{jmem},
                jname => qq"${jname}_${in}",
                jprefix => $options->{jprefix},
                jqueue => 'long',
                jstring => $jstring,
                modules => $options->{modules},
                stderr => $stderr,
                stdout => $stdout,
                prescript => $args{prescript},
                postscript => $args{postscript},);
        }
    } else { ## A single input was provided
        my $jstring = qq!
## Do not forget that _only_ the last command in a biopieces string is allowed to have the -x.
mkdir -p outputs/biopieces
less ${input} | read_fastq -i - -e base_33 |\\
 plot_scores -T 'Quality Scores' -t svg -o outputs/biopieces/${jname}_quality_scores.svg |\\
 plot_nucleotide_distribution -T 'NT. Distribution' -t svg -o outputs/biopieces/${jname}_ntdist.svg |\\
 plot_lendist -T 'Length Distribution' -k SEQ_LEN -t svg -o outputs/biopieces/${jname}_lendist.svg |\\
 analyze_gc |\\
     bin_vals -b 20 -k 'GC%' |\\
     plot_distribution -k 'GC%_BIN' -t svg -o outputs/biopieces/${jname}_gcdist.svg |\\
 analyze_gc |\\
     mean_vals -k 'GC%' -o outputs/biopieces/${jname}_gc.txt |\\
 count_records -o outputs/biopieces/${jname}_count.txt -x
!;
        $bp = $class->Submit(
            input => $input,
            comment => $comment,
            jdepends => $options->{jdepends},
            jmem => $options->{jmem},
            jname => 'biop',
            jprefix => $options->{jprefix},
            jstring => $jstring,
            modules => $options->{modules},
            stderr => $stderr,
            stdout => $stdout,
            prescript => $options->{prescript},
            postscript => $options->{postscript},);
    }
    return($bp);
}

=head2 C<Fastqc>

 Invoke fastqc on some sequence files.
 https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

=cut
sub Fastqc {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        filtered => 'unfiltered',
        jprefix => '01',
        modules => ['fastqc'],
        required => ['input',],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});

    my $job_name = $class->Get_Job_Name();
    my $input_paths = $class->Get_Paths($options->{input});
    my $dirname = $input_paths->[0]->{dirname};
    my $jname = qq"fqc_${job_name}";
    if (defined($dirname)) {
        $jname .= qq"_${dirname}";
    }
    my $outdir = qq"outputs/$options->{jprefix}fastqc";

    my $fastqc_job;
    my $input_file_string = '';
    my $subshell = 0;
    my $modified_input = undef;
    if (scalar(@{$input_paths}) > 1) {
        my @inputs;
        for my $element (@{$input_paths}) {
            push(@inputs, $element->{fullpath});
        }
        for my $in (@inputs) {
            $modified_input = basename($in, ('.gz', '.bz2', '.xz')) unless(defined($modified_input));
            $input_file_string = qq"$input_file_string ${in} ";
            if ($in =~ /\.xz$|\.bz2$/) {
                $subshell = 1;
                $input_file_string = qq"$input_file_string <(less $in) ";
            }
        }
    } elsif ($options->{input} =~ /\.xz$|\.bz2$/) {
        $modified_input = basename($options->{input}, ('.gz', '.bz2', '.xz')) unless ($modified_input);
        $subshell = 1;
        $input_file_string = qq" <(less $options->{input}) ";
    } else {
        $modified_input = basename($options->{input}, ('.gz')) unless ($modified_input);
        $input_file_string = qq" $options->{input} ";
    }
    $modified_input = basename($modified_input, ('.fastq'));
    $modified_input = qq"${modified_input}_fastqc";
    my $final_output = qq"${outdir}/${modified_input}";
    my $txtfile = qq"${final_output}/fastqc_data.txt";
    my $summary = qq"${final_output}/summary.txt";

    ## This is where fastqc should put its various output files
    ## with one important exception: It tries to be smart and take its own basename of the input
    ## but since I am using a bash subshell <(), it will get /dev/fd/xxx and so put the outputs into
    ## outputs/${jprefix}fastqc/xxx_fastqc...

    my $stdout = qq"${outdir}/${jname}-$options->{filtered}_fastqc.stdout";
    my $stderr = qq"${outdir}/${jname}-$options->{filtered}_fastqc.stderr";
    my $jstring = qq!mkdir -p ${outdir}
fastqc --extract \\
  -o ${outdir} \\
  ${input_file_string} \\
  2>${stderr} \\
  1>${stdout}
if [ "\$?" -ne "0" ]; then
  echo "fastqc failed, this is not considered fatal for a pipeline."
  exit 0
fi
!;

    if ($subshell) {
        $jstring .= qq!
## Note that if this is using a subshell, fastqc will assume that the inputs
## are /dev/fd/xx (usually 63 or 64).
## We can likely cheat and get the subshell fd with this:
badname=\$(basename <(env))
echo \${badname}
## with the caveat that this is subject to race conditions if a bunch of other things are
## creating subshells on this host.
mv ${outdir}/\${badname}_fastqc.html ${outdir}/${modified_input}.html 2>/dev/null
echo "move finished with: $?"
mv ${outdir}/\${badname}_fastqc.zip ${outdir}/${modified_input}.zip 2>/dev/null
echo "move finished with: $?"
mv \$(/bin/ls -d ${outdir}/\${badname}_fastqc) ${outdir}/${modified_input} 2>/dev/null
echo "move finished with: $?"
!;
    }

    my $comment = qq!## This FastQC run is against $options->{filtered} data and is used for
## an initial estimation of the overall sequencing quality.!;
    my $fqc = $class->Submit(
        comment => $comment,
        jcpu => 8,
        jname => $jname,
        jprefix => $options->{jprefix},
        jqueue => 'throughput',
        jstring => $jstring,
        modules => $options->{modules},
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        output => qq"$options->{jprefix}fastqc.html",
        summary => $summary,
        txtfile => $txtfile,
        stdout => $stdout,
        stderr => $stderr);

    my $stats = $class->Bio::Adventure::Metadata::Fastqc_Stats(
        input => $fqc->{txtfile},
        jcpu => 1,
        jmem => 1,
        jwalltime => '00:03:00',
        jdepends => $fqc->{job_id},);

    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($fqc);
}

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

L<fastqc> L<biopieces>

=cut

1;
