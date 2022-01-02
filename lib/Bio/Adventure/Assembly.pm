package Bio::Adventure::Assembly;
## LICENSE: gplv2
## ABSTRACT:  Kitty!
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use feature 'try';
use warnings qw"all";
no warnings 'experimental::try';
use Moo;
extends 'Bio::Adventure';

use Cwd qw"abs_path getcwd cwd";
use File::Basename;
use File::Spec;
use File::Which qw"which";
use Parse::CSV;

use Bio::Tools::Run::Alignment::StandAloneFasta;

no warnings 'experimental::try';

=head1 NAME

 Bio::Adventure::Assembly - Perform some de-novo assembly tasks.

=head1 SYNOPSIS

 Invocations for a few assemblers: Abyss, Shovill/Unicycler (spades), trinity.

=head1 METHODS

=head2 C<Abyss>

 Abyss: doi:10.1101/gr.214346.116 is a de novo assembler which uses short reads.

 This defaults to a k=41 abyss-pe or abyss-se assembly, depending on
 how many input files provided.  Ideally, it should at least do a
 little work to try to optimize k.  Abyss is the first program I have
 ever seen which uses make as an interpreter.

=item C<Arguments>

 input(required): one or more fastq files to assemble.
 k(41): Kmer size for the de Bruijn graph.
 jmem(12): Memory allocation on the cluster.
 modules('abyss'): Environment module to load.

=cut
sub Abyss {
    my ($class, %args) = @_;
    ## abyss-pe k=41 name=EAb01 in="r1_trimmed-corrected.fastq.gz r2_trimmed-corrected.fastq.gz"
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        k => 41,
        jmem => 12,
        modules => ['abyss'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('abyss-pe');
    die('Could not find abyss in your PATH.') unless($check);

    my $job_name = $class->Get_Job_Name();
    my %abyss_jobs = ();
    my $outname = basename(cwd());
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
    my $comment = '## This is a abyss submission script.';
    my $jstring = qq!start=\$(pwd)
mkdir -p ${output_dir}
cd ${output_dir}
rm -f ./*
${executable} -C \$(pwd) \\
    ${k_string} ${name_string} \\
    ${input_string} \\
    2>abyss_${outname}.stderr \\
    1>abyss_${outname}.stdout
cd \${start}
!;
    my $abyss = $class->Submit(
        cpus => 6,
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => "abyss_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => $options->{jmem},
        modules => $options->{modules},
        output => qq"${output_dir}/${outname}.fasta",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => 'workstation',);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($abyss);
}

=head2 C<Assembly_Coverage>

 Use hisat2 and bbmap's pileup script to calculate coverate on a per-contig basis.

 Given an assembly and a pile of reads, this will return the coverage
 on a per-nucleotide/per-contig basis, hopefully providing some metrics
 which will allow one to discriminate the quality of the contigs.  This
 is particularly notable when working with phage assemblies where one
 might expect prophage sequence to reside in the bacterial host at
 detectable and assemblable(is this a word?) levels, but at coverages
 which are notably lower than the actual viral sequences.  E.g. this
 should provide in effect a metric of the amount of bacterial
 contamination which snuck past previous filters.

=item C<Arguments>

 input(required): Filename(s) of corrected/filtered reads which were
  used to assemble the genome.
 library(required): Filename of the assembly.
 jmem(12): Memory allocated on the cluster.
 jprefix(14): job name/output directory prefix.
 modules('hisat', 'bbmap'): Environment modules loaded to run this.

=cut
sub Assembly_Coverage {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'library'],
        ## input is the corrected/filtered reads, library is the assembly
        jmem => 12,
        jprefix => 14,
        modules => ['hisat', 'bbmap'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});

    my $job_name = $class->Get_Job_Name();

    my $outname = basename(cwd());
    my $output_dir = qq"outputs/$options->{jprefix}assembly_coverage_${outname}";
    my $input_string = '';
    if ($options->{input} =~ /\:|\;|\,|\s+/) {
        my @in = split(/\:|\;|\,|\s+/, $options->{input});
        my $r1 = abs_path($in[0]);
        my $r2 = abs_path($in[1]);
        $input_string = qq" -1 <(less ${r1}) -2 <(less ${r2}) ";
    } else {
        my $r1 = abs_path($options->{input});
        $input_string = qq"-1 <(less ${r1}) ";
    }
    my $comment = qq!## This is a script to remap the reads against an assembly
## and calculate the coverage by contig.
!;
    my $jstring = qq!start=\$(pwd)
mkdir -p ${output_dir}
hisat2-build $options->{library} ${output_dir}/coverage_test \\
  2>${output_dir}/coverage_hisat_index.stderr 1>${output_dir}/coverage_hisat_index.stdout
hisat2 -x ${output_dir}/coverage_test -q \\
  ${input_string} -S ${output_dir}/coverage.sam \\
  2>${output_dir}/coverage_hisat.stderr 1>${output_dir}/coverage_hisat.stdout
pileup.sh in=${output_dir}/coverage.sam \\
  out=${output_dir}/coverage.tsv \\
  basecov=${output_dir}/base_coverage.tsv \\
  covwindow=100 \\
  k=19 \\
  overwrite=true \\
  2>${output_dir}/pileup.stderr 1>${output_dir}/pileup.stdout
samtools view -u -t $options->{library} \\
  -S ${output_dir}/coverage.sam -o ${output_dir}/coverage.bam \\
  2>${output_dir}/coverage_samtools.stderr 1>${output_dir}/coverage_samtools.stdout
rm ${output_dir}/coverage.sam
samtools sort -l 9 ${output_dir}/coverage.bam \\
  -o ${output_dir}/coverage_sorted.bam \\
  2>>${output_dir}/coverage_samtools.stderr 1>>${output_dir}/coverage_samtools.stdout
mv ${output_dir}/coverage_sorted.bam ${output_dir}/coverage.bam
samtools index ${output_dir}/coverage.bam \\
  2>>${output_dir}/coverage_samtools.stderr 1>>${output_dir}/coverage_samtools.stdout
!;
    my $coverage = $class->Submit(
        cpus => 6,
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => qq"coverage_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => $options->{jmem},
        jqueue => 'workstation',
        jwalltime => '4:00:00',
        modules => $options->{modules},
        output => qq"${output_dir}/coverage.txt",
        output_bam => qq"${output_dir}/coverage.bam",
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($coverage);
}

=head2 C<Collect_Assembly>

 Collect the final files created by an assembly.

 My little assembly pipeline generates quite a pile of
 files/directories.  Choosing the appropriate final outputs can be a
 bit daunting.  This function defines a few likely candidates and
 copies them to a single working directory.

=item C<Arguments>

 input_fsa(''): .fsa assembly to copy.
 input(''): I do not think this is used any longer.
 input_genbank(''): Output from the genbank file creator.
 input_stripped(''): Output from the genbank creator with some features
  stripped out.
 input_tsv(''): Tsv output from the genbank creator.
 input_xlsx(''): Xlsx output from the same.
 input_faa(''): Amino acid fasta file from the same.
 input_cds(''): CDS fasta file from same.
 jmem(2): Memory allocated on the cluster.
 jprefix(81): Default prefix for the job/directory.
 jname('collect'): Default jobname.

=cut
sub Collect_Assembly {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        input_fsa => '',
        input => '',
        input_genbank => '',
        input_tsv => '',
        output => '',
        jmem => 2,
        jprefix => 81,
        jname => 'collect',);
    my $output_dir = qq"outputs/$options->{jprefix}$options->{jname}";
    my @output_files = split(/:/, $options->{output});
    my $jstring = qq!
mkdir -p ${output_dir}
cp $options->{input_fsa} ${output_dir}
cp $options->{input} ${output_dir}
cp $options->{input_genbank} ${output_dir}
cp $options->{input_tsv} ${output_dir}
!;
    for my $o (@output_files) {
        $jstring .= qq"cp ${o} ${output_dir}
";
    }

    my $collect = $class->Submit(
        jdepends => $options->{jdepends},
        jname => $options->{jname},
        jstring => $jstring,
        jprefix => $options->{jprefix},
        jmem => $options->{jmem},
        output => $output_dir,);
    return($collect);
}

=head2 C<Unicycler_Filter_Worker>

 This does the work for Unicycler_Filter_Depth().

 Filter (for the moment only) a unicycler assembly by the ratio of the
 highest depth observed vs. the depth of each contig.  If that ratio
 falls below options->{coverage} then that contig should be dropped.
 The remaining contigs should be depth normalized to <=1. If there is
 only 1 contig, just return it with depth set to 1.  This should be
 trivially improved to handle other assembly methods by using the
 coverage calculation script above.

=item C<Arguments>

 input(required): Fasta assembly from unicycler to test.
 output(''): Set of filtered contigs.
 coverage(0.2): Minimal relative coverage ratio allowed.
 output_log(''): Location to write the log of the tasks performed.

=cut
sub Unicycler_Filter_Worker {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input',],
        coverage => 0.2,
        output_log => '',
        output => '',);
    my $paths = $class->Get_Paths($options->{output});
    my $log = FileHandle->new(">$options->{output_log}");
    my $input_contigs = Bio::SeqIO->new(-file => $options->{input}, -format => 'Fasta');
    my $output_contigs = Bio::SeqIO->new(-file => ">$options->{output}", -format => 'Fasta');
    my $log_string = qq"Starting depth coverage filter of $options->{input}.
Any contigs with a coverage ratio vs. the highest coverage of < $options->{coverage} will be dropped.
Writing filtered contigs to $options->{output}
";
    print $log $log_string;
    my $number_contigs = 0;
    ## The highest and lowest coverage values should begin at improbable values.
    my $highest_coverage = 0;
    my $lowest_coverage = 1000000;
    my $initial_coverage_data = {};
    my $final_coverage_data = {};
  INITIAL: while (my $contig_seq = $input_contigs->next_seq()) {
      my $number_contigs++;
      my $current_id = $contig_seq->id();
      my $current_desc = $contig_seq->desc();
      my ($contig_length_string, $contig_coverage_string, $circular) = split(/\s+/, $current_desc);
      my ($contig_length, $contig_coverage) = $current_desc =~ /^length=(.*)? depth=(.*)?x/;
      if ($contig_coverage > $highest_coverage) {
          $highest_coverage = $contig_coverage;
      }
      if ($contig_coverage < $lowest_coverage) {
          $lowest_coverage = $contig_coverage;
      }
      $initial_coverage_data->{$current_id} = {
          desc => $current_desc,
          coverage => $contig_coverage,
          length => $contig_length,
          circular => $circular,
          sequence => $contig_seq->seq(),
      };
  } ## Finished iterating over the input contigs, collected coverage and sequences.
    $log_string = qq"The range of observed coverages is ${lowest_coverage} <= x <= ${highest_coverage}\n";
    print $log $log_string;
    $input_contigs->close();

    ## Now write out the new data.
  WRITELOOP: foreach my $input_id (sort keys %{$initial_coverage_data}) {
      my $input_datum = $initial_coverage_data->{$input_id};
      my $normalized_coverage = $input_datum->{coverage} / $highest_coverage;
      if ($normalized_coverage < $options->{coverage}) {
          print $log "Skipping ${input_id}, its normalized coverage is: ${normalized_coverage}\n";
          next WRITELOOP;
      }
      $final_coverage_data->{$input_id} = $normalized_coverage;
      my $new_desc = qq"length=$input_datum->{length} coverage=${normalized_coverage}x";
      if (defined($input_datum->{circular})) {
          $new_desc .= qq" $input_datum->{circular}";
      }
      $log_string = qq"Writing ${input_id} with normalized coverage: ${normalized_coverage}\n";
      print $log $log_string;
      my $output_seq = Bio::Seq->new(-seq => $input_datum->{sequence},
                                     -desc => $new_desc,
                                     -id => $input_id);
      $output_contigs->write_seq($output_seq);
  }

    $log->close();
    $output_contigs->close();
    return($final_coverage_data);
}

=head2 C<Unicycler_Filter_Depth>

 Parse unicycler contig headers, extract relative coverage, and filter.

 Filter (for the moment only) a unicycler assembly by the ratio of the
 highest depth observed vs. the depth of each contig.  If that ratio
 falls below options->{coverage} then that contig should be dropped.
 The remaining contigs should be depth normalized to <=1. If there is
 only 1 contig, just return it with depth set to 1.  This should be
 trivially improved to handle other assembly methods by using the
 coverage calculation script above.

=item C<Arguments>

 input(required): Fasta assembly from unicycler to test.
 output(''): Set of filtered contigs.
 coverage(0.2): Minimal relative coverage ratio allowed.
 output_log(''): Location to write the log of the tasks performed.

=cut
sub Unicycler_Filter_Depth {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        coverage => 0.2,  ## The ratio of each sequence's coverage / the maximum coverage observed.
        jmem => 4,
        jprefix => '13',
        output => 'final_assembly.fasta',);
    my $job_name = $class->Get_Job_Name();
    my $outname = basename(cwd());
    my $output_dir = qq"outputs/$options->{jprefix}filter_depth";
    my $output_log = qq"${output_dir}/filter_depth.log";
    my $output = qq"${output_dir}/$options->{output}";
    my $comment = '## This is a submission script for a depth filter.';
    my $jstring = qq!
use Bio::Adventure;
use Bio::Adventure::Assembly;
my \$result = Bio::Adventure::Assembly::Unicycler_Filter_Worker(\$h,
  coverage => '$options->{coverage}',
  input => '$options->{input}',
  output => '${output}',
  output_log => '${output_log}',);
!;
    my $depth_filtered = $class->Submit(
        jdepends => $options->{jdepends},
        comment => $comment,
        jmem => $options->{jmem},
        jname => qq"filter_depth_${job_name}",
        jprefix => $options->{jprefix},
        jqueue => 'workstation',
        jstring => $jstring,
        language => 'perl',
        output => $output,
        output_log => $output_log,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        shell => '/usr/bin/env perl',);
    return($depth_filtered);
}

=head2 C<Shovill>

 Perform a shovill/spades assembly. https://github.com/tseemann/shovill

 I am learning that spades is the basis for more than a couple of
 assemblers. Shovill is one of them, and very feature-full.  Shovill
 has a lot of interesting options, this only includes a few at the
 moment.

=item C<Arguments>

 input(required): Filtered/corrected/trimmed reads.
 depth(40): Used to filter over/under represented contigs.
 jmem(12): Expected memory usage.
 jprefix('13'): Added to the beginning of the job name.
 arbitrary(''): Add arbitrary arguments here.
 modules('shovill'): Use this environment module.

=cut
sub Shovill {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        depth => 40,
        jmem => 12,
        jprefix => '13',
        arbitrary => '',
        modules => ['shovill',]);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('shovill');
    die('Could not find shovill in your PATH.') unless($check);
    my $job_name = $class->Get_Job_Name();
    my $outname = basename(cwd());
    my $output_dir = qq"outputs/$options->{jprefix}shovill";
    my $input_string = "";
    if ($options->{input} =~ /\:|\;|\,|\s+/) {
        my @in = split(/\:|\;|\,|\s+/, $options->{input});
        $input_string = qq" -R1 $in[0] -R2 $in[1]";
    } else {
        $input_string = qq" -R1 $options->{input}";
    }
    my $comment = '## This is a shovill submission script.';
    my $jstring = qq!mkdir -p ${output_dir}
shovill $options->{arbitrary} --force --keepfiles --depth $options->{depth} \\
   --outdir ${output_dir} \\
  $input_string \\
  2>${output_dir}/shovill_${outname}.stderr \\
  1>${output_dir}/shovill_${outname}.stdout
if [[ -f ${output_dir}/contigs.fa ]]; then
  mv ${output_dir}/contigs.fa ${output_dir}/final_assembly.fasta
elif [[ -f ${output_dir}/spades.fasta.uncorrected ]]; then
  cp ${output_dir}/spades.fasta.uncorrected ${output_dir}/final_assembly.fasta
else
  mv ${output_dir}/spades.fasta ${output_dir}/final_assembly.fasta
fi
!;
    my $shovill_job = $class->Submit(
        cpus => 6,
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => qq"shovill_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => $options->{jmem},
        jqueue => 'workstation',
        modules => $options->{modules},
        output => qq"${output_dir}/final_assembly.fasta",
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($shovill_job);
}


=head2 C<Trinity>

 Run the trinity de-novo transcriptome assembler.  doi:10.1038/nbt.1883.

 This function invokes trinity along with its default set of
 post-processing tools.

=item C<Arguments>

 input(required): Fastq input files, trimmed and presumably corrected.
 contig_length(600): Minimum contig length for a transcript.
 jmem(80): Trinity is quite memory intensive.
 jprefix(60): Prefix for job name.
 modules('trinity'): Use this environment module.

=cut
sub Trinity {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        contig_length => 600,
        jmem => 80,
        jprefix => '60',
        modules => ['trinity'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('Trinity');
    die('Could not find trinity in your PATH.') unless($check);

    my $job_name = $class->Get_Job_Name();
    my $output_dir = qq"outputs/$options->{jprefix}trinity_${job_name}";
    my $input_string = '';
    if ($options->{input} =~ /\:|\;|\,|\s+/) {
        my @in = split(/\:|\;|\,|\s+/, $options->{input});
        $input_string = qq"--left $in[0] --right $in[1] ";
    } else {
        $input_string = qq"--single $options->{input} ";
    }
    my $comment = '## This is a trinity submission script.';
    my $jstring = qq!mkdir -p ${output_dir} && \\
  Trinity --seqType fq --min_contig_length $options->{contig_length} --normalize_reads \\
    --trimmomatic --max_memory 90G --CPU 6 \\
    --output ${output_dir} \\
    ${input_string} \\
    2>${output_dir}/trinity_${job_name}.stderr \\
    1>${output_dir}/trinity_${job_name}.stdout
!;
    my $trinity = $class->Submit(
        cpus => 6,
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => qq"$options->{jprefix}trin_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => $options->{jmem},
        jqueue => 'large',
        modules => $options->{modules},
        output => qq"${output_dir}/Trinity.tsv",
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    my $rsem = $class->Bio::Adventure::Assembly::Trinity_Post(
        %args,
        jdepends => $trinity->{job_id},
        jname => qq"$options->{jprefix}_1trin_rsem",
        input => $options->{input},);
    my $trinotate = $class->Bio::Adventure::Annotation::Trinotate(
        %args,
        jdepends => $trinity->{job_id},
        jname => qq"$options->{jprefix}_2trinotate",
        input => qq"${output_dir}/Trinity.fasta",);
    $trinity->{rsem_job} = $rsem;
    $trinity->{trinotate_job} = $trinotate;
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($trinity);
}

=head2 C<Trinity_Post>

 Perform some of the post-processing tools provided by trinity.

 Trinity comes with some scripts to grade the quality of the
 transcripts it generated.  This invokes a few of them.

=item C<Arguments>

 input(required): The fasta assembly from trinity.
 jmem(24): Expected memory requirement.
 jname(trin_rsem): Name on the cluster.
 jprefix('61'): Job prefix.
 modules('rsem'): RSEM is used to quanify the transcripts.

=cut
sub Trinity_Post {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jmem => 24,
        jname => 'trin_rsem',
        jprefix => '61',
        modules => ['rsem'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $job_name = $class->Get_Job_Name();
    my $trinity_out_dir = qq"outputs/trinity_${job_name}";

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

    my $comment = '## This is a trinity post-processing submission script.';
    my $jstring = qq!
start=\$(pwd)
cd ${trinity_out_dir}
${trinity_exe_dir}/util/TrinityStats.pl Trinity.fasta \\
  2>${trinity_out_dir}/trinpost_stats.stderr \\
  1>${trinity_out_dir}/trinpost_stats.stdout &

${trinity_exe_dir}/util/align_and_estimate_abundance.pl \\
  --output_dir align_estimate.out \\
  --transcripts ${rsem_input} \\
  --seqType fq \\
  ${input_string} \\
  --est_method RSEM \\
  --aln_method bowtie \\
  --trinity_mode --prep_reference \\
  2>trinpost_align_estimate.stderr \\
  1>trinpost_align_estimate.stdout

${trinity_exe_dir}/util/abundance_estimates_to_matrix.pl \\
  --est_method RSEM \\
  --gene_trans_map Trinity.fasta.gene_trans_map \\
  align_estimate.out/RSEM.isoform.results \\
  2>trinpost_estimate_to_matrix.stderr \\
  1>trinpost_estimate_to_matrix.stdout

${trinity_exe_dir}/util/SAM_nameSorted_to_uniq_count_stats.pl \\
  bowtie_out/bowtie_out.nameSorted.bam \\
  2>trinpost_count_stats.stderr \\
  1>trinpost_count_stats.stdout

cd \${start}
!;
    my $trinpost = $class->Submit(
        comment => $comment,
        input => $options->{input},
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
        jname => qq"$options->{jprefix}trinpost_${job_name}",
        jprefix => $options->{jprefix},
        jqueue => 'large',
        jstring => $jstring,
        modules => $options->{modules},
        output => qq"${trinity_out_dir}/RSEM.isoform.results",
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($trinpost);
}


=head2 C<Unicycler>

 Use unicycler to assemble bacterial/viral reads. doi:10.1371/journal.pcbi.1005595

 Yet another Spades-based assembler, this one has some nice pre and
 post-processing tools written into it, which means that I do not need
 to figure out the arguments for pilon and stuff!

=item C<Arguments>

 input(required): Fastq input, it need not be filtered, trimmed,
  corrected.
 arbitrary(''): Add extra arguments here.
 depth(20): Include an extra depth filter(I think this is not
  currently used)
 min_length(1000): Drop contigs less than this size.
 mode('bold'): Unicycler has three modes, fast, normal, bold.  Bold
  tries hard to reduce the assembly graph to a single contig at the
  potential cost of correctness.
 jmem(24): Expected memory usage.
 jprefix('13'): Prefix for the job/output name.
 modules('trimomatic','bowtie2','spades','unicycer'): The environment
  module dependencies.

=cut
sub Unicycler {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        arbitrary => '',
        depth => 20,
        min_length => 1000,
        mode => 'bold',
        jmem => 24,
        jprefix => '13',
        modules => ['trimomatic', 'bowtie2', 'spades', 'unicycler'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('unicycler');
    die('Could not find unicycler in your PATH.') unless($check);

    my $job_name = $class->Get_Job_Name();
    my $outname = basename(cwd());
    my $output_dir = qq"outputs/$options->{jprefix}unicycler";
    my $input_string = '';
    my $ln_string = '';
    if ($options->{input} =~ /\:|\;|\,|\s+/) {
        my @in = split(/\:|\;|\,|\s+/, $options->{input});
        $input_string = qq" -1 ${output_dir}/r1.fastq.gz -2 ${output_dir}/r2.fastq.gz";
        $ln_string = qq"less $in[0] | gzip > ${output_dir}/r1.fastq.gz
less $in[1] | gzip > ${output_dir}/r2.fastq.gz
";
    } else {
        $input_string = qq" -1 ${output_dir}/r1.fastq.gz";
        $ln_string = qq"less $options->{input} | gzip > ${output_dir}/r1.fastq.gz
";
    }
    my $comment = '## This is a unicycler submission script.';
    my $jstring = qq!mkdir -p ${output_dir}
${ln_string}
unicycler $options->{arbitrary} \\
  --mode $options->{mode} \\
  --min_fasta_length $options->{min_length} \\
  ${input_string} \\
  -o ${output_dir} \\
  2>${output_dir}/unicycler_${outname}.stderr \\
  1>${output_dir}/unicycler_${outname}.stdout
mv ${output_dir}/assembly.fasta ${output_dir}/${outname}_final_assembly.fasta
rm -f r1.fastq.gz r2.fastq.gz
ln -sf ${output_dir}/${outname}_final_assembly.fasta unicycler_assembly.fasta
!;
    my $unicycler = $class->Submit(
        jdepends => $options->{jdepends},
        cpus => 6,
        comment => $comment,
        jmem => $options->{jmem},
        jname => qq"unicycler_${job_name}",
        jprefix => $options->{jprefix},
        jqueue => 'workstation',
        jstring => $jstring,
        modules => $options->{modules},
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        output => qq"${output_dir}/${outname}_final_assembly.fasta",
        output_gfa => qq"${output_dir}/assembly.gfa",
        output_log => qq"${output_dir}/unicycler.log",);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($unicycler);
}

=head2 C<Velvet>

 Perform a velvet assembly and pass it to ragoo.
 velvet: doi:10.1371/journal.pone.0008407
 ragoo: doi:10.1186/s13059-019-1829-6

 Yet another kmer-based assembler.

=item C<Arguments>

 input(required): Fastq reads, presumably filtered/corrected/trimmed.
 species(''): Used by ragoo when provided.
 kmer(31): Kmer length for velvet.
 jmem(24): Expected memory usage.
 modules(velvet): Environment modules used.

=item C<Invocation>

 > cyoa --task assembly --method velvet --input forward.fastq.gz:reverse.fastq.gz

=cut
sub Velvet {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        species => '',
        kmer => 31,
        jmem => 24,
        modules => ['velvet'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('velveth');
    die('Could not find velvet in your PATH.') unless($check);
    my $job_name = $class->Get_Job_Name();
    my $output_dir = qq"outputs/velvet_${job_name}";
    my $input_string = "";
    if ($options->{input} =~ /\:|\;|\,|\s+/) {
        my @in = split(/\:|\;|\,|\s+/, $options->{input});
        $input_string = qq" -fastq -short2 -separate <(less $in[0]) <(less $in[1])";
    } else {
        $input_string = qq" -fastq -short <(less $options->{input})";
    }
    my $comment = '## This is a velvet submission script.';
    my $jstring = qq!mkdir -p ${output_dir} && \\
  velveth ${output_dir} $options->{kmer} \\
    $input_string \\
    2>${output_dir}/velveth_${job_name}.stderr \\
    1>${output_dir}/velveth_${job_name}.stdout
  velvetg ${output_dir} \\
    -exp_cov auto -cov_cutoff auto \\
    2>${output_dir}/velvetg_${job_name}.stderr \\
    1>${output_dir}/velvetg_${job_name}.stdout
  new_params=\$(velvet-estimate-exp_cov.pl ${output_dir}/stats.txt \|
    grep velvetg parameters \|
    sed 's/velvetg parameters: //g')
  ##velvetg ${output_dir} \${new_params} -read_trkg yes -amos_file yes \\
  ##  2>${output_dir}/second_velvetg.txt 2>&1
!;
    if ($options->{species}) {
        $jstring .= qq!
  ragoo.py \\
    ${output_dir}/configs.fa \\
    $options->{libdir}/$options->{libtype}/$options->{species}.fasta
!;
    }
    my $velvet = $class->Submit(
        cpus => 6,
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => qq"velveth_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => $options->{jmem},
        modules => $options->{modules},
        output => qq"$output_dir/Sequences",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => 'workstation',);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($velvet);
}


=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

L<trinity> L<trinotate> L<transdecoder> L<velvet> L<ragoo>

=cut

1;
