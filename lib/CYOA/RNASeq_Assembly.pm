package CYOA;
use common::sense;
use autodie;

=head1 NAME

    CYOA::RNASeq_Assembly - Perform some denovo assembly tasks

=head1 SYNOPSIS

    use CYOA;
    my $hpgl = new CYOA;
    $hpgl->Trinity();

=head2 Methods

=over 4

=item C<Trinity>

    $hpgl->Trinity() submits a trinity denovo sequence assembly and
    runs its default post-processing tools.

=cut
sub Trinity {
    my $me = shift;
    my %args = @_;
    my $min_length = 600;
    $min_length = $args{min_length} if (defined($args{min_length}});
    my %trin_jobs = ();
    my $trin_depends_on;
    my $basename = $me->{basename};
    my $inputs = $me->{input};
    my @in = split(/\:/, $inputs);
    $inputs =~ s/\:/ /g;
    my $comment = qq!## This is a trinity submission script
!;
    my $job_string = qq!Trinity --seqType fq --min_contig_length ${min_length} --normalize_reads --trimmomatic --max_memory 90G --left $in[0] --right $in[1] --CPU 6 2>trinity.output 1>&2
!;
    my $trin_job = $me->Qsub(job_name => "trin",
                             depends => $trin_depends_on,
                             job_string => $job_string,
                             input => $inputs,
                             job_prefix => "45",
                             comment => $comment,
                             output => qq"trinity.output",
                             prescript => $args{prescript},
                             postscript => $args{postscript},
                             qsub_queue => "large",
                             qsub_cpus => 6,
                             qsub_mem => 90,
                             qsub_wall => "144:00:00",);

    my $rsem_job = $me->Trinity_Post(depends => $trin_job->{pbs_id},
                                     job_name => "trin_rsem",
                                     rsem_input => qq"trinity_out_dir/Trinity.fasta",);

    return($trin_job);
}

sub Trinity_Post {
    my $me = shift;
    my %args = @_;
    my $basename = $me->{basename};
    my $inputs = $me->{input};
    my @in = split(/\:/, $inputs);
    $inputs =~ s/\:/ /g;
    my $rsem_input = $args{rsem_input};
    my $job_name = qq"trin_rsem";
    $job_name = $args{job_name} if ($args{job_name});

    my $trinity_dir = qx"dirname Trinity";
    my $comment = qq!## This is a trinity submission script for rsem.
!;
    my $job_string = qq!
${trinity_dir}/utils/align_and_estimate_abundance.pl --output_dir trinity_out_dir/align_estimate.out --transcripts ${rsem_input} --seqType fq --left $in[0] --right $in[1] --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference 2>trinity_align_estimate.output 1>&2
${trinity_dir}/utils/filter_fasta_by_rsem_values.pl --rsem_output trinity_out_dir/align_estimate.out/RSEM.isoforms.results --fasta trinity_out_dir/Trinity.fasta --output filter_fasta_by_rsem.out 2>trinity_filter_rsem.output 1>&2
${trinity_dir}/utils/SAM_nameSorted_to_uniq_count_stats.pl trinity_out_dir/bowtie_out/bowtie_out.nameSorted.bam 2>>./count_stats.output 1>&2
!;
    my $trinpost_job = $me->Qsub(job_name => "trinpost",
                                 depends => $args{depends},
                                 job_string => $job_string,
                                 input => $inputs,
                                 comment => $comment,
                                 job_prefix => "46",
                                 output => qq"trinitypost.output",
                                 prescript => $args{prescript},
                                 postscript => $args{postscript},
                                 qsub_queue => "long",
                                 qsub_wall => "144:00:00",);
}


=back

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

    L<bowtie> L<bowtie2> L<tophat> L<bwa> L<kallisto> L<samtools> L<htseq>

=cut

1;
