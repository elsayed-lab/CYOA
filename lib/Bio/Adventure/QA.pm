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
        my $r1 = $input_list[0];
        my $r2 = $input_list[1];
    }

    my $job_name = $class->Get_Job_Name();
    my $input_paths = $class->Get_Paths($options->{input});
    my $jname = qq"fqc_${job_name}_$input_paths->{dirname}";
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
    my $fqc_stats = $class->Bio::Adventure::QA::Fastqc_Stats(
        input => $final_output,
        jdepends => $fqc->{job_id},
        jname => $jname,
        jprefix => $options->{jprefix} + 1,);
    $fqc->{stats} = $fqc_stats;
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($fqc);
}

=head2 C<Fastqc_Stats>

Collect some information from the fastqc output files and present them in a
simple-to-read csv file.

=cut
sub Fastqc_Stats {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jname => 'fqcst',
        paired => 1,);
    ## Dereferencing the options to keep a safe copy.
    my $jname = $options->{jname};
    ## Dereferencing the options to keep a safe copy.
    my $input_file = qq"$options->{input}/fastqc_data.txt";
    if ($options->{paired}) {
        $input_file = qq"$options->{input}/fastqc_data.txt";
    }
    my $stat_output = qq"outputs/fastqc_stats.csv";
    if ($options->{direction}) {
        $stat_output = qq"outputs/fastqc_$options->{direction}_stats.csv";
        $jname = qq"$options->{jname}_$options->{direction}";
    }
    my $comment = qq!## This is a stupidly simple job to collect alignment statistics.!;
    my $jstring = qq!
if [ \! -r "${stat_output}" ]; then
  echo "name,total_reads,poor_quality,per_quality,per_base_content,per_sequence_gc,per_base_n,per_seq_length,over_rep,adapter_content,kmer_content" > $stat_output
fi
total_reads_tmp=\$(grep "^Total Sequences" ${input_file} | awk -F '\\\\t' '{print \$2}')
total_reads=\${total_reads_tmp:-0}
poor_quality_tmp=\$(grep "^Sequences flagged as poor quality" ${input_file} | awk -F '\\\\t' '{print \$2}')
poor_quality=\${poor_quality_tmp:-0}
per_quality_tmp=\$(grep "Per base sequence quality" ${input_file} | awk -F '\\\\t' '{print \$2}')
per_quality=\${per_quality_tmp:-0}
per_base_content_tmp=\$(grep "Per base sequence content" ${input_file} | awk -F '\\\\t' '{print \$2}')
per_base_content=\${per_base_content_tmp:-0}
per_sequence_gc_tmp=\$(grep "Per sequence GC content" ${input_file} | awk -F '\\\\t' '{print \$2}')
per_sequence_gc=\${per_sequence_gc_tmp:-0}
per_base_n_tmp=\$(grep "Per base N content" ${input_file} | awk -F '\\\\t' '{print \$2}')
per_base_n=\${per_base_n_tmp:-0}
per_seq_length_tmp=\$(grep "Sequence Length Distribution" ${input_file} | awk -F '\\\\t' '{print \$2}')
per_seq_length=\${per_seq_length_tmp:-0}
over_rep_tmp=\$(grep "Overrepresented sequences" ${input_file} | awk -F '\\\\t' '{print \$2}')
over_rep=\${over_rep_tmp:-0}
adapter_content_tmp=\$(grep "Adapter Content" ${input_file} | awk -F '\\\\t' '{print \$2}')
adapter_content=\${adapter_content_tmp:-0}
kmer_content_tmp=\$(grep "Kmer Content" ${input_file} | awk -F '\\\\t' '{print \$2}')
kmer_content=\${kmer_content_tmp:-0}

stat_string=\$(printf "$options->{jname},%s,%s,%s,%s,%s,%s,%s,%s,%s,%s" "\${total_reads}" "\${poor_quality}" "\${per_quality}" "\${per_base_content}" "\${per_sequence_gc}" "\${per_base_n}" "\${per_seq_length}" "\${over_rep}" "\${adapter_content}" "\${kmer_content}")
echo "\$stat_string" >> ${stat_output}
!;
    my $stats = $class->Submit(
        comment => $comment,
        cpus => 1,
        input => $input_file,
        jdepends => $options->{jdepends},
        jmem => 1,
        jname => $jname,
        jprefix => $options->{jprefix},
        jqueue => 'throughput',
        jstring => $jstring,
        jwalltime => '00:1:00',
        output => qq"${stat_output}",);
    ## Added to return the state of the system to what it was
    ## before we messed with the options.
    return($stats);
}

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

L<fastqc> L<biopieces>

=cut

1;
