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
    );
    my $input = $options->{input};
    my $bp_depends_on;
    my $job_basename = $class->Get_Job_Name();
    my $jname = qq"biop_${job_basename}";
    $bp_depends_on = $options->{depends};
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
                depends => $bp_depends_on,
                jname => qq"${jname}_${in}",
                jprefix => "02",
                jstring => $jstring,
                prescript => $args{prescript},
                postscript => $args{postscript},
                queue => "long",
            );
        }
    } else {                    ## A single input was provided
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
            depends => $bp_depends_on,
            jname => "biop",
            jprefix => "02",
            jstring => $jstring,
            prescript => $options->{prescript},
            postscript => $options->{postscript},
            queue => "workstation",
        );
    }
    return($bp);
}

=head2 C<Fastqc>

Invoke fastqc on some sequence files.

=cut
sub Fastqc {
    my ($class, %args) = @_;
    try {
        use Env::Modulecmd { load  => 'fastqc' };
    }
    my $check = which('fastqc');
    die("Could not find fastqc in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input',],
    );
    my $fastqc_job;
    if ($options->{input} =~ /:|\,/) {
        $fastqc_job = Bio::Adventure::QA::Fastqc_Pairwise($class, %args);
    } else {
        $fastqc_job = Bio::Adventure::QA::Fastqc_Single($class, %args);
    }
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
        type => 'unfiltered',
    );
    my $type = $options->{type};
    my $input = $options->{input};
    my @input_list = split(/:|\,/, $input);
    if (scalar(@input_list) <= 1) {
        my $ret = Bio::Adventure::QA::Fastqc_Single($class, %args);
        return($ret);
    }
    my $r1 = $input_list[0];
    my $r2 = $input_list[1];
    my $basename = basename($r1, (".fastq", ".gz", ".xz"));
    $basename = basename($basename, (".fastq", ".gz", ".xz"));
    $basename = basename($basename, (".fastq", ".gz", ".xz"));
    my $output_suffix = "forward";
    if ($basename =~ /_R1/ or $basename =~ /forward/) {
        $output_suffix = "forward";
    } elsif ($basename =~ /_R2/ or $basename =~ /reverse/) {
        $output_suffix = "reverse";
    } else {
        $output_suffix = "all";
    }
    $basename =~ s/\_R1$//g;
    $basename =~ s/_forward//g;
    my $outdir = qq"outputs/fastqc";
    my $jstring = qq!mkdir -p ${outdir} && \\
  fastqc --extract -o ${outdir} ${r1} ${r2} \\
  2>${outdir}.out 1>&2
!;
    my $comment = qq!## This FastQC run is against ${type} data and is used for
## an initial estimation of the overall sequencing quality.!;
    my $job_basename = $class->Get_Job_Name();
    my $jname = qq"fqc_${job_basename}";
    my $fqc = $class->Submit(
        comment => $comment,
        cpus => 8,
        jname => $jname,
        jprefix => "00",
        jstring => $jstring,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => "throughput",
    );
    my $forward_indir = qq"${outdir}/${r1}_fastqc";
    $forward_indir =~ s/\.fastq//g;
    my $reverse_indir = qq"${outdir}/${r2}_fastqc";
    $reverse_indir =~ s/\.fastq//g;
    my $fsf = Bio::Adventure::QA::Fastqc_Stats(
        $class,
        direction => 'forward',
        indir => $forward_indir,
        depends => $fqc->{pbs_id},
        jname => $basename,
        jprefix => "01",
        paired => 1,
    );
    my $fsr = Bio::Adventure::QA::Fastqc_Stats(
        $class,
        direction => 'reverse',
        indir => $reverse_indir,
        depends => $fqc->{pbs_id},
        jname => $basename,
        jprefix => "01",
        paired => 1,
    );
    $fqc->{stats_forward} = $fsf;
    $fqc->{stats_reverse} = $fsr;
    return($fqc);
}

=head2 C<Fastqc_Single>

Invoke fastqc with options suitable for a single-ended reads.

=cut
sub Fastqc_Single {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        filtered => 'unfiltered',
    );
    my $outdir = qq"outputs/fastqc";
    my $jstring = qq!mkdir -p ${outdir} &&\\
  fastqc -q --extract -o ${outdir} $options->{input} \\
  2>outputs/$options->{jname}-$options->{filtered}_fastqc.out 1>&2
!;
    my $comment = qq!## This FastQC run is against $options->{filtered} data and is used for
## an initial estimation of the overall sequencing quality.!;
    my $job_basename = $class->Get_Job_Name();
    my $jname = qq"fqc_${job_basename}";
    my $fqc = $class->Submit(
        comment => $comment,
        cpus => 8,
        jname => $jname,
        jprefix => "00",
        jstring => $jstring,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => "throughput",
    );
    $outdir .= "/" . basename($options->{input}, (".fastq.gz",".fastq.xz", ".fastq")) . "_fastqc";
    my $fqc_stats = Bio::Adventure::QA::Fastqc_Stats(
        $class,
        indir => $outdir,
        depends => $fqc->{pbs_id},
        jname => $options->{jname},
        jprefix => "01",
    );
    $fqc->{stats} = $fqc_stats;
    return($fqc);
}

=head2 C<Fastqc_Stats>

Collect some information from the fastqc output files and present them in a
simple-to-read csv file.

=cut
sub Fastqc_Stats {
    my ($class, %args) = @_;
    my $opt = $class->Get_Vars(
        args => \%args,
        jname => 'fqcst',
        depends => '',
        paired => 1,
    );
    ## Dereferencing the options to keep a safe copy.
    my %start_options = %{$opt};
    my $jname = $opt->{jname};
    my $input_file = qq"$opt->{indir}/fastqc_data.txt";
    if ($opt->{paired}) {
        $input_file = qq"$opt->{indir}/fastqc_data.txt";
    }
    my $depends = $opt->{depends};
    my $stat_output = qq"outputs/fastqc_stats.csv";
    if ($opt->{direction}) {
        $stat_output = qq"outputs/fastqc_$opt->{direction}_stats.csv";
        $jname = qq"${jname}_$opt->{direction}";
    }
    my $comment = qq!## This is a stupidly simple job to collect alignment statistics.!;
    my $jstring = qq!
if [ \! -r $stat_output ]; then
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

stat_string=\$(printf "${jname},%s,%s,%s,%s,%s,%s,%s,%s,%s,%s" "\${total_reads}" "\${poor_quality}" "\${per_quality}" "\${per_base_content}" "\${per_sequence_gc}" "\${per_base_n}" "\${per_seq_length}" "\${over_rep}" "\${adapter_content}" "\${kmer_content}")
echo "\$stat_string" >> ${stat_output}
!;
    my $stats = $class->Submit(
        comment => $comment,
        cpus => 1,
        input => $input_file,
        depends => $depends,
        jname => $jname,
        jprefix => $opt->{jprefix},
        jstring => $jstring,
        mem => 1,
        queue => "throughput",
        walltime => "00:10:00",
    );
    ## Added to return the state of the system to what it was
    ## before we messed with the options.
    $class->{options} = \%start_options;
    return($stats);
}

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

L<fastqc> L<biopieces>

=cut

1;
