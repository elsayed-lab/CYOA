package CYOA;

=head1 NAME

    CYOA::RNASeq_QA - Use trimomatic/cutadapt/etc to trim libraries

=head1 SYNOPSIS

    use CYOA;
    my $hpgl = new CYOA;
    $hpgl->Fastqc();

=head2 Methods

=over 4

=item C<Biopieces_Graph>

    Reads in a fastq file and uses biopieces to make some graphs
    describing the sequences therein.

=cut
sub Biopieces_Graph {
    my $me = shift;
    my %args = @_;
    $me->Check_Options(['input']);
    my $input = $me->{input};
    my $bp_depends_on;
    my $basename = $me->{basename};
    $bp_depends_on = $args{depends} if ($args{depends});
    my @inputs = split(/\,|\:|\;/, $input);
    my $comment = qq!## This script uses biopieces to draw some simple graphs of the sequence.!;
    my $bp;
    if (scalar(@inputs) > 1) { ## multiple comma/colon/semicolon inputs were provided.
        foreach my $in (@inputs) {
            my $short_in = basename($in, ('.gz','.fastq'));
            $short_in = basename($short_in, ('.gz','.fastq'));
            my $job_string = qq!
mkdir -p outputs/biopieces
less ${in} | read_fastq -i - -e base_33 |\\
 plot_scores -T \'Quality Scores\' -t svg -o outputs/biopieces/${short_in}_quality_scores.svg |\\
 plot_nucleotide_distribution -T 'NT. Distribution' -t svg -o outputs/biopieces/${short_in}_ntdist.svg |\\
 plot_lendist -x -T 'Length Distribution' -k SEQ_LEN -t svg -o outputs/biopieces/${short_in}_lendist.svg |\\
 analyze_gc | bin_vals -b 5 -k GC\% | plot_distribution -k GC\%_BIN -t svg -o outputs/biopieces/${short_in}_gc_dist.svg |\\
 analyze_gc | mean_vals -k GC\% -o outputs/biopieces/${short_in}_gc.txt |\\
 count_records -o outputs/biopieces/${short_in}_count.txt -x
!;
            $bp = $me->Qsub(job_name => "biop",
                               depends => $bp_depends_on,
                               job_string => $job_string,
                               comment => $comment,
                               input => $in,
                               prescript => $args{prescript},
                               postscript => $args{postscript},
                );
        }
    } else {  ## A single input was provided
        my $job_string = qq!
less ${input} | read_fastq -i - -e base_33 |\\
 plot_scores -T \'Quality Scores\' -t svg -o outputs/biopieces/${basename}_quality_scores.svg |\\
 plot_nucleotide_distribution -T 'NT. Distribution' -t svg -o outputs/biopieces/${basename}_ntdist.svg |\\
 plot_lendist -x -T 'Length Distribution' -k SEQ_LEN -t svg -o outputs/biopieces/${basename}_lendist.svg |\\
 analyze_gc | bin_vals -b 5 -k GC\% | plot_distribution -k GC\%_BIN -t svg -o outputs/biopieces/${basename}_gc_dist.svg |\\
 analyze_gc | mean_vals -k GC\% -o outputs/biopieces/${basename}_gc.txt |\\
 count_records -o outputs/biopieces/${basename}_count.txt -x
!;
        $bp = $me->Qsub(job_name => "biop",
                        depends => $bp_depends_on,
                        job_string => $job_string,
                        comment => $comment,
                        input => $input,
                        prescript => $args{prescript},
                        postscript => $args{postscript},
            );
    }
    return($bp);
}

sub Fastqc {
    my $me = shift;
    my %args = @_;
    $me->Check_Options(['input',]);
    my $fastqc_job;
    if ($me->{input} =~ /:|\,/) {
        $fastqc_job = $me->Fastqc_Pairwise(%args);
    } else {
        $fastqc_job = $me->Fastqc_Single(%args);
    }
    return($fastqc_job);
}

=item C<Fastqc_Pairwise>

    Fastqc_Pairwise()

=cut
sub Fastqc_Pairwise {
    my $me = shift;
    my %args = @_;
    my $type = $args{type};
    $type = "unfiltered" unless ($type);
    my $input = $me->{input};
    my @input_list = split(/:|\,/, $input);
    if (scalar(@input_list) <= 1) {
        my $ret = $me->Fastqc_Single(@_);
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
    $me->{basename} = $basename;
    my $outdir = qq"outputs/fastqc";
    my $job_string = qq!mkdir -p ${outdir} && \\
  fastqc --extract -o ${outdir} ${r1} ${r2} \\
  2>${outdir}.out 1>&2
!;
    my $comment = qq!## This FastQC run is against ${type} data and is used for
## an initial estimation of the overall sequencing quality.!;
    my $fqc_jobid = qq"${basename}_fqc";
    my $fqc = $me->Qsub(job_name => "fqc",
                          qsub_cpus => 8,
                          qsub_queue => "workstation",
                          job_string => $job_string,
                          comment => $comment,
			  prescript => $args{prescript},
			  postscript => $args{postscript},
        );
    my $forward_indir = qq"${outdir}/${r1}_fastqc";
    $forward_indir =~ s/\.fastq//g;
    my $reverse_indir = qq"${outdir}/${r2}_fastqc";
    $reverse_indir =~ s/\.fastq//g;
    my $fqc_stats_foward = $me->Fastqc_Stats(basename => $basename,
                                             indir => $forward_indir,
                                             depends => $fqc->{pbs_id},
                                             paired => 1,
                                             direction => 'forward');
    my $fqc_stats_reverse = $me->Fastqc_Stats(basename => $basename,
                                              indir => $reverse_indir,
                                              depends => $fqc->{pbs_id},
                                              paired => 1,
                                              direction => 'reverse');
    $fqc->{stats_forward} = $fcq_stats_forward;
    $fqc->{stats_reverse} = $fcq_stats_reverse;
    return($fqc);
}

=item C<Fastqc_Single>

    $hpgl->Fastqc_Single()

=cut
sub Fastqc_Single {
    my $me = shift;
    my %args = @_;
    my $filtered = $args{filtered};
    $filtered = "unfiltered" unless ($filtered);
    my $input = $me->{input};
    my $basename = $me->{basename};
    my $outdir = qq"$me->{basename}/outputs/${basename}-${filtered}_fastqc";
    my $job_string = qq!mkdir -p ${outdir} &&\\
  fastqc --extract -o ${outdir} ${input} \\
  2>outputs/${basename}-${filtered}_fastqc.out 1>&2
!;
    my $comment = qq!## This FastQC run is against ${filtered} data and is used for
## an initial estimation of the overall sequencing quality.!;
    my $fqc_jobid = qq"${basename}_fqc";
    my $fqc = $me->Qsub(job_name => "fqc",
                          qsub_cpus => 8,
                          qsub_queue => "workstation",
                          job_string => $job_string,
                          comment => $comment,
			  prescript => $args{prescript},
			  postscript => $args{postscript},
        );
    my $fqc_stats = $me->Fastqc_Stats(basename => $basename, indir => $outdir, depends => $fqc->{pbs_id});
    $fqc->{stats} = $fqc_stats;
    return($fqc);
}

=item C<Fastqc_Stats>

    $hpgl->Fastqc_Stats()

=cut
sub Fastqc_Stats {
    my $me = shift;
    my %args = @_;
    my $basename = $me->{basename};
    my $input_file = qq"$args{indir}/fastqc_data.txt";
    if ($args{paired}) {
        $input_file = qq"$args{indir}/fastqc_data.txt";
    }
    my $depends = "";
    $depends = $args{depends} if ($args{depends});
    my $job_name = "fqcst";
    $job_name = $args{job_name} if ($args{job_name});
    my $stat_output = qq"outputs/fastqc_stats.csv";
    if ($args{direction}) {
        $stat_output = qq"outputs/fastqc_$args{direction}_stats.csv";
        $job_name = qq"${job_name}_$args{direction}";
    }
    my $comment = qq!## This is a stupidly simple job to collect alignment statistics.!;
    my $job_string = qq!
if [ \! -r $stat_output ]; then
  echo "total_reads,poor_quality,per_quality,per_base_content,per_sequence_gc,per_base_n,per_seq_length,over_rep,adapter_content,kmer_content" > $stat_output
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

stat_string=\$(printf "${basename},%s,%s,%s,%s,%s,%s,%s,%s,%s,%s" "\${total_reads}" "\${poor_quality}" "\${per_quality}" "\${per_base_content}" "\${per_sequence_gc}" "\${per_base_n}" "\${per_seq_length}" "\${over_rep}" "\${adapter_content}" "\${kmer_content}")
echo "\$stat_string" >> $stat_output
!;
    my $stats = $me->Qsub(job_name => $job_name,
                          depends => $depends,
                          qsub_queue => "throughput",
                          qsub_cpus => 1,
                          qsub_mem => 1,
                          qsub_wall => "00:10:00",
                          job_string => $job_string,
                          input => $input_file,
                          comment => $comment,
        );
    return($stats);
}

=back

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

    L<CYOA>

=cut

1;
