package Bio::Adventure::RNASeq_Assembly;
## LICENSE: gplv2
## ABSTRACT:  Kitty!
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use File::Which qw"which";

=head1 NAME

    Bio::Adventure::RNASeq_Assembly - Perform some denovo assembly tasks

=head1 SYNOPSIS

    use Bio::Adventure;
    my $hpgl = new Bio::Adventure;
    $hpgl->Trinity();

=head2 Methods

=over 4

=item C<Trinity>

    $hpgl->Trinity() submits a trinity denovo sequence assembly and
    runs its default post-processing tools.

=cut
sub Trinity {
    my ($class, %args) = @_;
    my $check = which('Trinity');
    die("Could not find trinity in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        min_length => 600,
        required => ['input'],
    );
    my $min_length = $options->{min_length};
    my %trin_jobs = ();
    my $trin_depends_on;
    my $basename = $options->{basename};
    my $inputs = $options->{input};
    my @in = split(/\:/, $inputs);
    $inputs =~ s/\:/ /g;
    my $comment = qq!## This is a trinity submission script
!;
    my $jstring = qq!Trinity --seqType fq --min_contig_length ${min_length} --normalize_reads \\
   --trimmomatic --max_memory 90G --left $in[0] --right $in[1] --CPU 6 2>trinity.output 1>&2
!;
    my $trin_job = $class->Submit(
        cpus => 6,
        comment => $comment,
        input => $inputs,
        depends => $trin_depends_on,
        jname => "trin",
        jprefix => "45",
        jstring => $jstring,
        mem => 90,
        output => qq"trinity.output",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => "large",
        walltime => "144:00:00",
    );
    my $rsem_job = Bio::Adventure::RNASeq_Assembly::Trinity_Post(
        $class,
        %args,
        depends => $trin_job->{pbs_id},
        jname => "trin_rsem",
        rsem_input => qq"trinity_out_dir/Trinity.fasta",
    );

    return($trin_job);
}

sub Trinity_Post {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        jname => "trin_rsem",
    );
    my $basename = $options->{basename};
    my $inputs = $options->{input};
    my @in = split(/\:/, $inputs);
    $inputs =~ s/\:/ /g;
    my $rsem_input = $options->{rsem_input};
    my $jname = $options->{jname};

    my $trinity_dir = qx"dirname Trinity";
    my $comment = qq!## This is a trinity submission script for rsem.
!;
    my $jstring = qq!
${trinity_dir}/utils/align_and_estimate_abundance.pl --output_dir trinity_out_dir/align_estimate.out \\
  --transcripts ${rsem_input} --seqType fq --left $in[0] --right $in[1] --est_method RSEM \\
  --aln_method bowtie --trinity_mode --prep_reference 2>trinity_align_estimate.output 1>&2
${trinity_dir}/utils/filter_fasta_by_rsem_values.pl \\
  --rsem_output trinity_out_dir/align_estimate.out/RSEM.isoforms.results \\
  --fasta trinity_out_dir/Trinity.fasta \\
  --output filter_fasta_by_rsem.out 2>trinity_filter_rsem.output 1>&2
${trinity_dir}/utils/SAM_nameSorted_to_uniq_count_stats.pl \\
  trinity_out_dir/bowtie_out/bowtie_out.nameSorted.bam 2>>./count_stats.output 1>&2
!;
    my $trinpost_job = $class->Submit(
        comment => $comment,
        input => $inputs,
        depends => $options->{depends},
        jname => "trinpost",
        jprefix => "46",
        jstring => $jstring,
        output => qq"trinitypost.output",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => "long",
        walltime => "144:00:00",
    );
}


=back

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

    L<bowtie> L<bowtie2> L<tophat> L<bwa> L<kallisto> L<samtools> L<htseq>

=cut

1;
