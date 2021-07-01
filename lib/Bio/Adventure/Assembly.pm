package Bio::Adventure::Assembly;
## LICENSE: gplv2
## ABSTRACT:  Kitty!
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use Cwd qw"abs_path getcwd";
use File::Basename;
use File::Spec;
use File::Which qw"which";
use Parse::CSV;

=head1 NAME

Bio::Adventure::Assembly - Perform some de-novo assembly tasks.

=head1 SYNOPSIS

Invocations for a few assemblers: Abyss, Shovill/Unicycler (spades), trinity.

=head1 METHODS

=head2 C<Abyss>

This defaults to a k=41 abyss-pe or abyss-se assembly, depending on
how many input files provided.  Ideally, it should at least do a
little work to try to optimize k.  Abyss is the first program I have
ever seen which uses make as an interpreter.

=cut
sub Abyss {
    my ($class, %args) = @_;
    ## abyss-pe k=41 name=EAb01 in="r1_trimmed-corrected.fastq.gz r2_trimmed-corrected.fastq.gz"
    my $check = which('abyss-pe');
    die("Could not find abyss in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        k => 41,
    );
    my $job_basename = $class->Get_Job_Name();
    my %abyss_jobs = ();
    my $outname = basename(getcwd());
    my $output_dir = qq"outputs/abyss_${outname}";
    my $k_string = qq"k=$options->{k} ";
    my $name_string = qq"name=${outname} ";
    my $input_string = "";
    my $executable = "abyss-pe";
    if ($options->{input} =~ /\:|\;|\,|\s+/) {
        my @in = split(/\:|\;|\,|\s+/, $options->{input});
        my $r1 = abs_path($in[0]);
        my $r2 = abs_path($in[1]);
        $input_string = qq!in="${r1} ${r2}" !;
    } else {
        my $r1 = abs_path($options->{input});
        $input_string = qq!in=${r1} "!;
        $executable = "abyss-se"
    }
    my $comment = qq!## This is a abyss submission script
!;
    my $jstring = qq!start=\$(pwd)
mkdir -p ${output_dir}
cd ${output_dir}
rm -f ./*
${executable} -C \$(pwd) \\
    ${k_string} ${name_string} \\
    ${input_string} \\
    2>abyss_${outname}.err \\
    1>abyss_${outname}.out
cd \${start}
!;
    my $abyss_job = $class->Submit(
        cpus => 6,
        comment => $comment,
        jname => "abyss_${job_basename}",
        jprefix => "46",
        jstring => $jstring,
        mem => 30,
        output => qq"outputs/abyss.sbatchout",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => "workstation",
        walltime => "4:00:00",
    );
    return($abyss_job);
}

=head2 C<Trinity>

$hpgl->Trinity() submits a trinity denovo sequence assembly and runs its default
post-processing tools.

=over

=item I<input> * File(s) containing reads to be assembled by trinity.

=item I<contig_length> (600) Minimum length for a contig to keep.

=back

=head3 C<Invocation>

> cyoa --task assembly --method trinity --input forward.fastq.gz:reverse.fastq.gz

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
    my $rsem_job = Bio::Adventure::Assembly::Trinity_Post(
        $class,
        %args,
        depends => $trin_job->{pbs_id},
        jname => "trin_rsem",
        input => $options->{input},
    );
    my $trinotate_job = Bio::Adventure::Annotation::Trinotate(
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

=head2 C<Trinity_Post>

Perform some of the post-processing tools provided by trinity.

=over

=item I<input> * Input fasta from trinity.

=back

=head3 C<Invocation>

> cyoa --task assembly --method trinitypost --input trinity.fasta

=cut
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
start=\$(pwd)
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

cd \${start}
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


=head2 C<Velvet>

Submit sequences for a generic assembly by velvet and pass it to ragoo.

=over

=item I<input> * Input reads passed to velveth.

=item I<species> * Species name for passing to ragoo.py.

=back

=head3 C<Invocation>

> cyoa --task assembly --method velvet --input forward.fastq.gz:reverse.fastq.gz

=cut
sub Velvet {
    my ($class, %args) = @_;
    my $check = which('velveth');
    die("Could not find velvet in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        kmer => 31,
        required => ['input', 'species'],
    );
    my $job_basename = $class->Get_Job_Name();
    my %velvet_jobs = ();
    my $output_dir = qq"outputs/velvet_${job_basename}";
    my $input_string = "";
    if ($options->{input} =~ /\:|\;|\,|\s+/) {
        my @in = split(/\:|\;|\,|\s+/, $options->{input});
        $input_string = qq" -fastq -short2 -separate <(less $in[0]) <(less $in[1])";
    } else {
        $input_string = qq" -fastq -short <(less $options->{input})";
    }
    my $comment = qq!## This is a velvet submission script
!;
    my $jstring = qq!mkdir -p ${output_dir} && \\
  velveth ${output_dir} $options->{kmer} \\
    $input_string \\
    2>${output_dir}/velveth_${job_basename}.err \\
    1>${output_dir}/velveth_${job_basename}.out
  velvetg ${output_dir} \\
    -exp_cov auto -cov_cutoff auto \\
    2>${output_dir}/velvetg_${job_basename}.err \\
    1>${output_dir}/velvetg_${job_basename}.out
  new_params=\$(velvet-estimate-exp_cov.pl ${output_dir}/stats.txt \|
    grep velvetg parameters \|
    sed 's/velvetg parameters: //g')
  ##velvetg ${output_dir} \${new_params} -read_trkg yes -amos_file yes \\
  ##  2>${output_dir}/second_velvetg.txt 2>&1
  ragoo.py \\
    ${output_dir}/configs.fa \\
    $options->{libdir}/$options->{libtype}/$options->{species}.fasta
!;
    my $velvet_job = $class->Submit(
        cpus => 6,
        comment => $comment,
        jname => "velveth_${job_basename}",
        jprefix => "46",
        jstring => $jstring,
        mem => 30,
        output => qq"outputs/velvet.sbatchout",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => "workstation",
        walltime => "4:00:00",
    );
    return($velvet_job);
}

=head2 C<Shovill>

Use shovill to perform/optimize a spades assembly.  Shovill has a lot
of interesting options, this only includes a few at the moment.

=cut
sub Shovill {
    my ($class, %args) = @_;
    my $check = which('shovill');
    die("Could not find shovill in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        depth => 20,
        arbitrary => '',
    );
    my $job_basename = $class->Get_Job_Name();
    my %shovill_jobs = ();
    my $outname = basename(getcwd());
    my $output_dir = qq"outputs/shovill_${outname}";
    my $input_string = "";
    if ($options->{input} =~ /\:|\;|\,|\s+/) {
        my @in = split(/\:|\;|\,|\s+/, $options->{input});
        $input_string = qq" -R1 $in[0] -R2 $in[1]";
    } else {
        $input_string = qq" -R1 $options->{input}";
    }
    my $comment = qq!## This is a shovill submission script
!;
    my $jstring = qq!mkdir -p ${output_dir} && \\
  shovill $options->{arbitrary} --force --keepfiles --depth $options->{depth} \\
     --outdir ${output_dir} \\
    $input_string \\
    2>${output_dir}/shovill_${outname}.err \\
    1>${output_dir}/shovill_${outname}.out
!;
    my $shovill_job = $class->Submit(
        cpus => 6,
        comment => $comment,
        jname => "shovill_${job_basename}",
        jprefix => "46",
        jstring => $jstring,
        mem => 30,
        output => qq"outputs/shovill.sbatchout",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => "workstation",
        walltime => "4:00:00",
    );
    return($shovill_job);
}


=head2 C<Unicycler>

Use unicycler to assemble bacterial/viral reads.

=cut
sub Unicycler {
    my ($class, %args) = @_;
    my $check = which('unicycler');
    die("Could not find unicycler in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        depth => 20,
        mode => 'normal',
        arbitrary => '',
    );
    my $job_basename = $class->Get_Job_Name();
    my %unicycler_jobs = ();
    my $outname = basename(getcwd());
    my $output_dir = qq"outputs/unicycler_${outname}";
    my $input_string = "";
    if ($options->{input} =~ /\:|\;|\,|\s+/) {
        my @in = split(/\:|\;|\,|\s+/, $options->{input});
        $input_string = qq" -1 $in[0] -2 $in[1]";
    } else {
        $input_string = qq" -1 $options->{input}";
    }
    my $comment = qq!## This is a unicycler submission script
!;
    my $jstring = qq!mkdir -p ${output_dir}
unicycler --mode $options->{mode} $options->{arbitrary} \\
  ${input_string} \\
  -o ${output_dir} \\
  2>${output_dir}/unicycler_${outname}.err \\
  1>${output_dir}/unicycler_${outname}.out
!;
    my $unicycler_job = $class->Submit(
        cpus => 6,
        comment => $comment,
        jname => "unicycler_${job_basename}",
        jprefix => "46",
        jstring => $jstring,
        mem => 30,
        output => qq"outputs/unicycler.sbatchout",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => "workstation",
        walltime => "4:00:00",
    );
    return($unicycler_job);

}

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

L<trinity> L<trinotate> L<transdecoder> L<velvet> L<ragoo>

=cut

1;
