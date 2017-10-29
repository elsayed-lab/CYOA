package Bio::Adventure::RNASeq_Assembly;
## LICENSE: gplv2
## ABSTRACT:  Kitty!
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use File::Basename;
use File::Spec;
use File::Which qw"which";

=head1 NAME

    Bio::Adventure::RNASeq_Assembly - Perform some denovo assembly tasks

=head1 SYNOPSIS

    use Bio::Adventure;
    my $hpgl = new Bio::Adventure;
    $hpgl->Trinity();

=head2 Methods

=over 4

=item C<Transdecoder>

    $hpgl->Transdecoder() submits a trinity denovo sequence assembly and
    runs its default post-processing tools.

=cut
sub Transdecoder {
    my ($class, %args) = @_;
    my $check = which('TransDecoder.LongOrfs');
    die("Could not find transdecoder in your PATH.") unless($check);
    my $transdecoder_exe_dir = dirname($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
    );
    my $transdecoder_input = File::Spec->rel2abs($options->{input});
    my $job_basename = $class->Get_Job_Name();
    my $output_dir = qq"outputs/trinity_${job_basename}";
    my $comment = qq!## This is a transdecoder submission script
!;
    my $jstring = qq!mkdir -p ${output_dir} && cd ${output_dir} && \\
  TransDecoder.LongOrfs -t ${transdecoder_input} \\
    2>${output_dir}/transdecoder_longorfs_${job_basename}.err \\
    1>${output_dir}/transdecoder_longorfs_${job_basename}.out
TransDecoder.Predict -t ${transdecoder_input} \\
  2>${output_dir}/transdecoder_predict_${job_basename}.err \\
  1>${output_dir}/transdecoder_predict_${job_basename}.out
${transdecoder_exe_dir}/util/cdna_alignment_orf_to_genome_orf.pl \\
  ${output_dir}/transcripts.fasta.transdecoder.gff3 \\
  ${output_dir}/transcripts.gff3 \\
  ${transdecoder_input} \\
  2>${output_dir}/cdna_alignment_${job_basename}.err \\
  1>${output_dir}/transcripts.fasta.transdecoder.genome.gff3
!;
    my $transdecoder_job = $class->Submit(
        comment => $comment,
        jname => "transdecoder_${job_basename}",
        jprefix => "47",
        jstring => $jstring,
        output => qq"outputs/transdecoder.sbatchout",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => "workstation",
    );
    return($transdecoder_job);
}

=item C<Trinotate>

    $hpgl->Trinotate() submits a trinity denovo sequence assembly and
    runs its default post-processing tools.

=cut
sub Trinotate {
    my ($class, %args) = @_;
    my $check = which('Trinotate');
    die("Could not find trinotate in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
    );
    my $job_basename = $class->Get_Job_Name();
    my $output_dir = qq"outputs/trinity_${job_basename}";
    my $trinotate_exe_dir = dirname($check);
    my $input_dir = dirname($options->{input});
    my $input_name = basename($options->{input});
    my $comment = qq!## This is a trinotate submission script
!;
    my $jstring = qq!cd ${input_dir}
${trinotate_exe_dir}/auto/autoTrinotate.pl \\
  --Trinotate_sqlite ${trinotate_exe_dir}/sample_data/Trinotate.boilerplate.sqlite \\
  --transcripts $input_name \\
  --gene_to_trans_map Trinity.fasta.gene_trans_map \\
  --conf ${trinotate_exe_dir}/auto/conf.txt \\
  --CPU 6 \\
  2>trinotate_${job_basename}.err \\
  1>trinotate_${job_basename}.out
!;
    my $trinotate_job = $class->Submit(
        cpus => 6,
        comment => $comment,
        jname => "trinotate_${job_basename}",
        jprefix => "48",
        jstring => $jstring,
        output => qq"trinotate.sbatchout",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => "workstation",
    );
    return($trinotate_job);
}

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
        contig_length => 600,
        required => ['input'],
    );
    my $job_basename = $class->Get_Job_Name();
    my %trin_jobs = ();
    my $trin_depends_on;
    my $output_dir = qq"outputs/trinity_${job_basename}";
    my $input_string = "";
    if ($options->{input} =~ /\:|\;|\,|\s+/) {
        my @in = split(/\:|\;|\,|\s+/, $options->{input});
        $input_string = qq"--left $in[0] --right $in[1] ";
    } else {
        $input_string = qq"--single $options->{input} ";
    }
    my $comment = qq!## This is a trinity submission script
!;
    my $jstring = qq!mkdir -p ${output_dir} && \\
  Trinity --seqType fq --min_contig_length $options->{contig_length} --normalize_reads \\
    --trimmomatic --max_memory 90G --CPU 6 \\
    --output ${output_dir} \\
    ${input_string} \\
    2>${output_dir}/trinity_${job_basename}.err \\
    1>${output_dir}/trinity_${job_basename}.out
!;
    my $trin_job = $class->Submit(
        cpus => 6,
        comment => $comment,
        depends => $trin_depends_on,
        jname => "trin_${job_basename}",
        jprefix => "45",
        jstring => $jstring,
        mem => 96,
        output => qq"outputs/trinity.sbatchout",
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
        input => $options->{input},
    );
    my $trinotate_job = Bio::Adventure::RNASeq_Assemly::Trinotate(
        $class,
        %args,
        depends => $trin_job->{pbs_id},
        jname => "trinotate",
        input => qq"${output_dir}/Trinity.fasta",
    );
    my $jobs = {
        trinity => $trin_job,
        trinity_post => $rsem_job,
        trinotate => $trinotate_job,
    };
    return($jobs);
}

sub Trinity_Post {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jname => "trin_rsem",
    );
    my $job_basename = $class->Get_Job_Name();
    my $trinity_out_dir = qq"outputs/trinity_${job_basename}";

    my $rsem_input = qq"${trinity_out_dir}/Trinity.fasta";
    my $trinity_path = which('Trinity');
    my $trinity_exe_dir = dirname($trinity_path);
    my $input_string = "";
    if ($options->{input} =~ /\:|\;|\,|\s+/) {
        my @in = split(/\:|\;|\,|\s+/, $options->{input});
        $input_string = qq"--left $in[0] --right $in[1] ";
    } else {
        $input_string = qq"--single $options->{input} ";
    }

    my $comment = qq!## This is a trinity post-processing submission script.
!;
    my $jstring = qq!
cd ${trinity_out_dir}
${trinity_exe_dir}/util/TrinityStats.pl Trinity.fasta \\
  2>${trinity_out_dir}/trinpost_stats.err \\
  1>${trinity_out_dir}/trinpost_stats.out &

${trinity_exe_dir}/util/align_and_estimate_abundance.pl \\
  --output_dir align_estimate.out \\
  --transcripts ${rsem_input} \\
  --seqType fq \\
  ${input_string} \\
  --est_method RSEM \\
  --aln_method bowtie \\
  --trinity_mode --prep_reference \\
  2>trinpost_align_estimate.err \\
  1>trinpost_align_estimate.out

${trinity_exe_dir}/util/abundance_estimates_to_matrix.pl \\
  --est_method RSEM \\
  --gene_trans_map Trinity.fasta.gene_trans_map \\
  align_estimate.out/RSEM.isoform.results \\
  2>trinpost_estimate_to_matrix.err \\
  1>trinpost_estimate_to_matrix.out

${trinity_exe_dir}/util/SAM_nameSorted_to_uniq_count_stats.pl \\
  bowtie_out/bowtie_out.nameSorted.bam \\
  2>trinpost_count_stats.err \\
  1>trinpost_count_stats.out
!;
    my $trinpost_job = $class->Submit(
        comment => $comment,
        input => $options->{input},
        depends => $options->{depends},
        jname => "trinpost_$job_basename",
        jprefix => "46",
        jstring => $jstring,
        mem => 90,
        output => qq"outputs/trinitypost.sbatchout",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => "large",
        walltime => "144:00:00",
##        queue => "long",
##        walltime => "144:00:00",
    );
    return($trinpost_job);
}


=back

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

    L<bowtie> L<bowtie2> L<tophat> L<bwa> L<kallisto> L<samtools> L<htseq>

=cut

1;
