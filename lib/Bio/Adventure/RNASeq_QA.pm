package Bio::Adventure::RNASeq_QA;
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

    Bio::Adventure::RNASeq_QA - Use trimomatic/cutadapt/etc to trim libraries

=head1 SYNOPSIS

    use Bio::Adventure;
    my $hpgl = new Bio::Adventure;
    $hpgl->Fastqc();

=head2 Methods

=over 4

=item C<Biopieces_Graph>

    Reads in a fastq file and uses biopieces to make some graphs
    describing the sequences therein.

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
less ${in} | read_fastq -i - -e base_33 |\\
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

sub Fastqc {
    my ($class, %args) = @_;
    my $check = which('fastqc');
    die("Could not find fastqc in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input',],
    );
    my $fastqc_job;
    if ($options->{input} =~ /:|\,/) {
        $fastqc_job = Bio::Adventure::RNASeq_QA::Fastqc_Pairwise($class, %args);
    } else {
        $fastqc_job = Bio::Adventure::RNASeq_QA::Fastqc_Single($class, %args);
    }
    return($fastqc_job);
}

=item C<Fastqc_Pairwise>

    Fastqc_Pairwise()

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
        my $ret = Bio::Adventure::RNASeq_QA::Fastqc_Single($class, %args);
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
    my $jname = qq"fcq_${job_basename}";
    my $fqc = $class->Submit(
        comment => $comment,
        cpus => 8,
        jname => $jname,
        jprefix => "00",
        jstring => $jstring,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => "workstation",
    );
    my $forward_indir = qq"${outdir}/${r1}_fastqc";
    $forward_indir =~ s/\.fastq//g;
    my $reverse_indir = qq"${outdir}/${r2}_fastqc";
    $reverse_indir =~ s/\.fastq//g;
    my $fsf = Bio::Adventure::RNASeq_QA::Fastqc_Stats(
        $class,
        direction => 'forward',
        indir => $forward_indir,
        depends => $fqc->{pbs_id},
        jname => $basename,
        jprefix => "01",
        paired => 1,
    );
    my $fsr = Bio::Adventure::RNASeq_QA::Fastqc_Stats(
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

=item C<Fastqc_Single>

    $hpgl->Fastqc_Single()

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
    my $jname = qq"fcq_${job_basename}";
    my $fqc = $class->Submit(
        comment => $comment,
        cpus => 8,
        jname => $jname,
        jprefix => "00",
        jstring => $jstring,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => "workstation",
    );
    $outdir .= "/" . basename($options->{input}, (".fastq.gz",".fastq.xz", ".fastq")) . "_fastqc";
    my $fqc_stats = Bio::Adventure::RNASeq_QA::Fastqc_Stats(
        $class,
        indir => $outdir,
        depends => $fqc->{pbs_id},
        jname => $options->{jname},
        jprefix => "01",
    );
    $fqc->{stats} = $fqc_stats;
    return($fqc);
}

=item C<Fastqc_Stats>

    $hpgl->Fastqc_Stats()

=cut
sub Fastqc_Stats {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        jname => 'fqcst',
        depends => '',
        paired => 1,
    );
    my $jname = $options->{jname};
    my $input_file = qq"$options->{indir}/fastqc_data.txt";
    if ($options->{paired}) {
        $input_file = qq"$options->{indir}/fastqc_data.txt";
    }
    my $depends = $options->{depends};
    my $stat_output = qq"outputs/fastqc_stats.csv";
    if ($options->{direction}) {
        $stat_output = qq"outputs/fastqc_$options->{direction}_stats.csv";
        $jname = qq"${jname}_$options->{direction}";
    }
    my $comment = qq!## This is a stupidly simple job to collect alignment statistics.!;
    my $jstring = qq!
if [ \! -r $stat_output ]; then
  echo "name,total_reads,poor_quality,per_quality,per_base_content,per_sequence_gc,per_base_n,per_seq_length,over_rep,adapter_content,kmer_content" > $stat_output
fi
total_reads_tmp=\$(grep "^Total Sequences" $input_file | awk -F '\\\\t' '{print \$2}')
total_reads=\${total_reads_tmp:-0}
poor_quality_tmp=\$(grep "^Sequences flagged as poor quality" $input_file | awk -F '\\\\t' '{print \$2}')
poor_quality=\${poor_quality_tmp:-0}
per_quality_tmp=\$(grep "Per base sequence quality" $input_file | awk -F '\\\\t' '{print \$2}')
per_quality=\${per_quality_tmp:-0}
per_base_content_tmp=\$(grep "Per base sequence content" $input_file | awk -F '\\\\t' '{print \$2}')
per_base_content=\${per_base_content_tmp:-0}
per_sequence_gc_tmp=\$(grep "Per sequence GC content" $input_file | awk -F '\\\\t' '{print \$2}')
per_sequence_gc=\${per_sequence_gc_tmp:-0}
per_base_n_tmp=\$(grep "Per base N content" $input_file | awk -F '\\\\t' '{print \$2}')
per_base_n=\${per_base_n_tmp:-0}
per_seq_length_tmp=\$(grep "Sequence Length Distribution" $input_file | awk -F '\\\\t' '{print \$2}')
per_seq_length=\${per_seq_length_tmp:-0}
over_rep_tmp=\$(grep "Overrepresented sequences" $input_file | awk -F '\\\\t' '{print \$2}')
over_rep=\${over_rep_tmp:-0}
adapter_content_tmp=\$(grep "Adapter Content" $input_file | awk -F '\\\\t' '{print \$2}')
adapter_content=\${adapter_content_tmp:-0}
kmer_content_tmp=\$(grep "Kmer Content" $input_file | awk -F '\\\\t' '{print \$2}')
kmer_content=\${kmer_content_tmp:-0}

stat_string=\$(printf "${jname},%s,%s,%s,%s,%s,%s,%s,%s,%s,%s" "\${total_reads}" "\${poor_quality}" "\${per_quality}" "\${per_base_content}" "\${per_sequence_gc}" "\${per_base_n}" "\${per_seq_length}" "\${over_rep}" "\${adapter_content}" "\${kmer_content}")
echo "\$stat_string" >> $stat_output
!;
    my $stats = $class->Submit(
        comment => $comment,
        cpus => 1,
        input => $input_file,
        depends => $depends,
        jname => $jname,
        jprefix => $options->{jprefix},
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

    L<Bio::Adventure>

=cut

1;
