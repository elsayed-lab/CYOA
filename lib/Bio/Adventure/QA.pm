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
        jprefix => '02',
        modules => ['biopieces'],);
    my $input = $options->{input};
    my $jname = $class->Get_Job_Name();
    $jname = qq"biop_${jname}";
    my @inputs = split(/\,|\:|\;/, $input);
    my $comment = qq!## This script uses biopieces to draw some simple graphs of the sequence.!;
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
                comment => $comment,
                input => $in,
                jdepends => $options->{jdepends},
                jname => qq"${jname}_${in}",
                jprefix => $options->{jprefix},
                jqueue => 'long',
                jstring => $jstring,
                modules => $options->{modules},
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
            comment => $comment,
            input => $input,
            jdepends => $options->{jdepends},
            jname => 'biop',
            jprefix => $options->{jprefix},
            jstring => $jstring,
            modules => $options->{modules},
            prescript => $options->{prescript},
            postscript => $options->{postscript},);
    }
    return($bp);
}

=head2 C<Fastqc>

Invoke fastqc on some sequence files.

=cut
sub Fastqc {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input',],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $fastqc_job;
    if ($options->{input} =~ /:|\,/) {
        $fastqc_job = Bio::Adventure::QA::Fastqc_Pairwise($class, %args);
    } else {
        $fastqc_job = Bio::Adventure::QA::Fastqc_Single($class, %args);
    }
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($fastqc_job);
}

=head2 C<Fastqc_Pairwise>

Invoke fastq with options suitable for pairwise sequence data, separating the
forward and reverse read outputs.

=cut
sub Fastqc_Pairwise {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        filtered => 'unfiltered',
        jprefix => '01',
        modules => ['fastqc'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});

    my @input_list = split(/:|\,|\s+/, $options->{input});
    my ($r1, $r2);
    if (scalar(@input_list) <= 1) {
        my $ret = $class->Bio::Adventure::QA::Fastqc_Single(%args);
        return($ret);
    } else {
        $r1 = $input_list[0];
        $r2 = $input_list[1];
    }

    my $job_name = $class->Get_Job_Name();
    my $input_paths = $class->Get_Paths($options->{input});
    my $dirname = $input_paths->[0]->{dirname};
    my $jname = qq"fqc_${job_name}";
    if (defined($dirname)) {
        $jname .= qq"_${dirname}";
    }
    my $outdir = qq"outputs/$options->{jprefix}fastqc";

    ## This is where fastqc should put its various output files
    ## with one important exception: It tries to be smart and take its own basename of the input
    ## but since I am using a bash subshell <(), it will get /dev/fd/xxx and so put the outputs into
    ## outputs/${jprefix}fastqc/xxx_fastqc...
    my $modified_inputname = basename($r1, (".fastq.gz",".fastq.xz", ".fastq")) . "_fastqc";
    my $final_output = qq"${outdir}/${modified_inputname}";

    my $jstring = qq!mkdir -p ${outdir} && \\
  fastqc --extract -o ${outdir} <(less ${r1}) <(less ${r2}) \\
  2>outputs/${jname}-$options->{filtered}_fastqc.out 1>&2

## Note that because I am using a subshell, fastqc will assume that the inputs
## are /dev/fd/xx (usually 63 or 64).
## We can likely cheat and get the subshell fd with this:
badname=\$(basename <(env))
echo \${badname}
## with the caveat that this is subject to race conditions if a bunch of other things are
## creating subshells on this host.
mv ${outdir}/\${badname}_fastqc.html ${outdir}/${modified_inputname}.html
mv ${outdir}/\${badname}_fastqc.zip ${outdir}/${modified_inputname}.zip
mv \$(/bin/ls -d ${outdir}/\${badname}_fastqc) ${outdir}/${modified_inputname}

!;
    my $comment = qq!## This FastQC run is against $options->{filtered} data and is used for
## an initial estimation of the overall sequencing quality.!;
    my $fqc = $class->Submit(
        comment => $comment,
        cpus => 8,
        jname => $jname,
        jprefix => $options->{jprefix},
        jqueue => 'throughput',
        jstring => $jstring,
        modules => $options->{modules},
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        output => qq"$options->{jprefix}fastqc.html",);
    return($fqc);
}

=head2 C<Fastqc_Single>

Invoke fastqc with options suitable for a single-ended reads.

=cut
sub Fastqc_Single {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        filtered => 'unfiltered',
        jprefix => '01',
        modules => ['fastqc'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $job_name = $class->Get_Job_Name();
    my $input_paths = $class->Get_Paths($options->{input});
    my $jname = qq"fqc_${job_name}_$input_paths->{dirname}";
    my $outdir = qq"outputs/$options->{jprefix}fastqc";

    ## This is where fastqc should put its various output files
    ## with one important exception: It tries to be smart and take its own basename of the input
    ## but since I am using a bash subshell <(), it will get /dev/fd/xxx and so put the outputs into
    ## outputs/${jprefix}fastqc/xxx_fastqc...
    my $modified_inputname = basename($options->{input}, (".fastq.gz",".fastq.xz", ".fastq")) . "_fastqc";
    my $final_output = qq"${outdir}/${modified_inputname}";

    my $jstring = qq!mkdir -p ${outdir} && \\
  fastqc -q --extract -o ${outdir} <(less $options->{input}) \\
  2>outputs/${jname}-$options->{filtered}_fastqc.out 1>&2
## Note that because I am using a subshell, fastqc will assume that the inputs
## are /dev/fd/xx (usually 63 or 64).
## We can likely cheat and get the subshell fd with this:
badname=\$(basename <(env))
echo \${badname}
## with the caveat that this is subject to race conditions if a bunch of other things are
## creating subshells on this host.
mv ${outdir}/\${badname}_fastqc.html ${outdir}/${modified_inputname}.html
mv ${outdir}/\${badname}_fastqc.zip ${outdir}/${modified_inputname}.zip
mv \$(/bin/ls -d ${outdir}/\${badname}_fastqc) ${outdir}/${modified_inputname}
!;
    my $comment = qq!## This FastQC run is against $options->{filtered} data and is used for
## an initial estimation of the overall sequencing quality.!;
    my $fqc = $class->Submit(
        comment => $comment,
        cpus => 8,
        jname => $jname,
        jprefix => $options->{jprefix},
        jqueue => 'throughput',
        jstring => $jstring,
        modules => $options->{modules},
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        output => qq"$options->{jprefix}fastqc.html",);
    $outdir .= "/" . basename($options->{input}, (".fastq.gz",".fastq.xz", ".fastq")) . "_fastqc";
    my $newname = qq"fqcstats_${job_name}_$input_paths->{dirname}";
    my $fqc_stats = $class->Bio::Adventure::Metadata::Fastqc_Stats(
        input => $final_output,
        jdepends => $fqc->{job_id},
        jname => $jname,
        jprefix => $options->{jprefix} + 1,);
    $fqc->{stats} = $fqc_stats;
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
