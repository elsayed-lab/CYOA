package HPGL;

=head2
    Biopieces_Graph()
    Reads in a fastq file and uses biopieces to make some graphs
    describing the sequences therein.
=cut
sub Biopieces_Graph {
    my $me = shift;
    my %args = @_;
    my $input = $me->{input};
    my $bp_depends_on;
    my $basename = $me->{basename};
    $bp_depends_on = $args{depends} if ($args{depends});
    my @inputs = split(/\,/, $input);
    my $comment = qq!## This script uses biopieces to draw some simple graphs of the sequence.!;
    my $job_string = qq!
mkdir -p biopieces
xzcat -f ${input} | read_fastq -i - -e base_33 |\\
 plot_scores -T \'Quality Scores\' -t svg -o biopieces/${basename}_quality_scores.svg |\\
 plot_nucleotide_distribution -T 'NT. Distribution' -t svg -o biopieces/${basename}_ntdist.svg |\\
 plot_lendist -T 'Length Distribution' -k SEQ_LEN -t svg -o biopieces/${basename}_lendist.svg |\\
!;
    my $bp = $me->Qsub(job_name => "biop",
                       depends => $bp_depends_on,
                       job_string => $job_string,
                       comment => $comment,
                       input => $input,
		       prescript => $args{prescript},
		       postscript => $args{postscript},
	);
    return($bp);
}


=head2
    Fastqc_Pairwise()
=cut
sub Fastqc_Pairwise {
    my $me = shift;
    my %args => @_;
    my $type = $args{type};
    $type = "unfiltered" unless ($type);
    my $input = $me->{input};
    my @input_list = split(/\:|\,/, $input);
    if (scalar(@input_list) <= 1) {
        my $ret = $me->Fastqc_Single(@_);
        return($ret);
    }
    my $r1 = $input_list[0];
    my $r2 = $input_list[1];
    my $basename = basename($r1, (".fastq"));
    $basename =~ s/\_R1$//g;
    $me->{basename} = $basename;
    my $job_string = qq!mkdir -p ${basename}-${type}_fastqc &&\\
  fastqc --extract -o ${basename}-${type}_fastqc ${r1} ${r2} \\
  2>outputs/${basename}-${type}_fastqc.out 1>&2
!;
    my $comment = qq!## This FastQC run is against ${type} data and is used for
## an initial estimation of the overall sequencing quality.!;
    my $fqc_jobid = qq"${basename}_fqc";
    my $jobid = $me->Qsub(job_name => "fqc",
                          qsub_cpus => 8,
                          qsub_queue => "workstation",
                          job_string => $job_string,
                          comment => $comment,
			  prescript => $args{prescript},
			  postscript => $args{postscript},
        );
}

=head2
    Fastqc_Single()
=cut
sub Fastqc_Single {
    my $me = shift;
    my %args = @_;
    my $filtered = $args{filtered};
    $filtered = "unfiltered" unless ($filtered);
    my $input = $me->{input};
    my $basename = $me->{basename};
    my $outdir = qq"${basename}-${filtered}_fastqc";
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
    return($fqc);
}

sub Fastqc_Stats {
    my $me = shift;
    my %args = @_;
    my $basename = $me->{basename};
    my $input_file = qq"$args{indir}/${basename}_fastqc/fastqc_data.txt";
    my $depends = "";
    $depends = $args{depends} if ($args{depends});
    my $job_name = "fqcst";
    $job_name = $args{job_name} if ($args{job_name});
    my $comment = qq!
## This is a stupidly simple job to collect alignment statistics.
!;
    my $job_string = qq!
total_reads_tmp=\$(grep "^Total Sequences" $input_file | awk -F '\\t' '{print \$2}')
total_reads=\${total_reads_tmp:-0}
poor_quality_tmp=\$(grep "^Sequences flagged as poor quality" $input_file | awk -F '\\t' '{print \$2}')
poor_quality=\${poor_quality_tmp:-0}
per_quality_tmp=\$(grep "Per base sequence quality" $input_file | awk -F '\\t' '{print \$2}')
per_quality=\${per_quality_tmp:-0}
per_base_content_tmp=\$(grep "Per base sequence content" $input_file | awk -F '\\t' '{print \$2}')
per_base_content=\${per_base_content_tmp:-0}
per_sequence_gc_tmp=\$(grep "Per sequence GC content" $input_file | awk -F '\\t' '{print \$2}')
per_sequence_gc=\${per_sequence_gc_tmp:-0}
per_base_n_tmp=\$(grep "Per base N content" $input_file | awk -F '\\t' '{print \$2}')
per_base_n=\${per_base_n_tmp:-0}
per_seq_length_tmp=\$(grep "Sequence Length Distribution" $input_file | awk -F '\\t' '{print \$2}')
per_seq_length=\${per_seq_length_tmp:-0}
over_rep_tmp=\$(grep "Overrepresented sequences" $input_file | awk -F '\\t' '{print \$2}')
over_rep=\${over_rep_tmp:-0}
adapter_content_tmp=\$(grep "Adapter Content" $input_file | awk -F '\\t' '{print \$2}')
adapter_content=\${adapter_content_tmp:-0}
kmer_content_tmp=\$(grep "Kmer Content" $input_file | awk -F '\\t' '{print \$2}')
kmer_content=\${kmer_content_tmp:-0}

stat_string=\$(printf "${basename},%s,%s,%s,%s,%s,%s,%s,%s,%s,%s" "\${total_reads}" "\${poor_quality}" "\${per_quality}" "\${per_base_content}" "\${per_sequence_gc}" "\${per_base_n}" "\${per_seq_length}" "\${over_rep}" "\${adapter_content}" "\${kmer_content}")
echo "\$stat_string" >> stats/fastqc_stats.csv
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

1;
