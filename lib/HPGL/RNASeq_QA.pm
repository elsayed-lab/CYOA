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
  fastqc -o ${basename}-${type}_fastqc ${r1} ${r2} \\
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
    my %args => @_;
    my $type = $args{type};
    $type = "unfiltered" unless ($type);
    my $input = $me->{input};
    my $basename = $me->{basename};
    my $job_string = qq!mkdir -p ${basename}-${type}_fastqc &&\\
  fastqc -o ${basename}-${type}_fastqc ${input} \\
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


1;
