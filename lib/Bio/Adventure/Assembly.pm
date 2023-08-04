package Bio::Adventure::Assembly;
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

=over

=item C<Arguments>

 input(required): one or more fastq files to assemble.
 k(41): Kmer size for the de Bruijn graph.
 jmem(12): Memory allocation on the cluster.
 modules('abyss'): Environment module to load.

=item C<Invocation>

> cyoa --task assemble --method abyss --input r1.fastq.xz:r2.fastq.xz

=cut
sub Abyss {
    my ($class, %args) = @_;
    ## abyss-pe k=41 name=EAb01 in="r1_trimmed-corrected.fastq.gz r2_trimmed-corrected.fastq.gz"
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        k => 41,
        jmem => 12,);
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
        comment => $comment,
        jcpu => 6,
        jdepends => $options->{jdepends},
        jname => "abyss_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => $options->{jmem},
        output => qq"${output_dir}/${outname}.fasta",
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($abyss);
}

=back

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

=over

=item C<Arguments>

 input(required): Filename(s) of corrected/filtered reads which were
  used to assemble the genome.
 library(required): Filename of the assembly.
 jmem(12): Memory allocated on the cluster.
 jprefix(14): job name/output directory prefix.
 modules('hisat2', 'bbmap'): Environment modules loaded to run this.

=item C<Invocation>

> cyoa --task assem --method assemblyc --input r1.fastq.xz:r2.fastq.xz --library assembly.fasta

=cut
sub Assembly_Coverage {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'library'],
        ## input is the corrected/filtered reads, library is the assembly
        jmem => 18,
        jprefix => 14,);
    my $job_name = $class->Get_Job_Name();
    my $outname = basename(cwd());
    my $output_dir = qq"outputs/$options->{jprefix}assembly_coverage_${outname}";
    my $input_string = '';
    if ($options->{input} =~ /\:|\;|\,|\s+/) {
        my @in = split(/\:|\;|\,|\s+/, $options->{input});
        my $r1 = abs_path($in[0]);
        my $r2 = abs_path($in[1]);
        $input_string = qq" -1 <(less ${r1}) -2 <(less ${r2}) ";
        if ($r1 =~ /\.fastq$/) {
            $input_string = qq" -1 ${r1} -2 ${r2} ";
        }
    } else {
        my $r1 = abs_path($options->{input});
        $input_string = qq"-1 <(less ${r1}) ";
        if ($r1 =~ /\.fastq$/) {
            $input_string = qq" -1 ${r1} ";
        }
    }
    my $comment = qq!## This is a script to remap the reads against an assembly
## and calculate the coverage by contig.
!;
    my $stdout = qq"${output_dir}/coverage.stdout";
    my $stderr = qq"${output_dir}/coverage.stderr";
    my $jstring = qq!start=\$(pwd)
mkdir -p ${output_dir}
hisat2-build $options->{library} ${output_dir}/coverage_test \\
  2>${stderr} \\
  1>${stdout}
hisat2 -x ${output_dir}/coverage_test -q \\
  ${input_string} -S ${output_dir}/coverage.sam \\
  2>>${stderr} \\
  1>>${stdout}
pileup.sh in=${output_dir}/coverage.sam \\
  out=${output_dir}/coverage.tsv \\
  basecov=${output_dir}/base_coverage.tsv \\
  covwindow=100 \\
  k=19 \\
  overwrite=true \\
  2>>${stderr} \\
  1>>${stdout}
samtools view -u -t $options->{library} \\
  -S ${output_dir}/coverage.sam -o ${output_dir}/coverage.bam \\
  2>>${stderr} \\
  1>>${stdout}
rm ${output_dir}/coverage.sam
samtools sort -l 9 ${output_dir}/coverage.bam \\
  -o ${output_dir}/coverage_sorted.bam \\
  2>>${stderr} \\
  1>>${stdout}
mv ${output_dir}/coverage_sorted.bam ${output_dir}/coverage.bam
samtools index ${output_dir}/coverage.bam \\
  2>>${stderr} \\
  1>>${stdout}
!;
    my $coverage = $class->Submit(
        comment => $comment,
        jcpu => 6,
        jdepends => $options->{jdepends},
        jname => qq"coverage_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => $options->{jmem},
        jwalltime => '4:00:00',
        output => qq"${output_dir}/coverage.txt",
        output_bam => qq"${output_dir}/coverage.bam",
        output_tsv => qq"${output_dir}/coverage.tsv",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        stdout => $stdout,
        stderr => $stderr);
    return($coverage);
}

=back

=head2 C<Collect_Assembly>

 Collect the final files created by an assembly.

 My little assembly pipeline generates quite a pile of
 files/directories.  Choosing the appropriate final outputs can be a
 bit daunting.  This function defines a few likely candidates and
 copies them to a single working directory.

=over

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

=item C<Invocation>

 This is really only invoked at the end of a pipeline.

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
        jcpu => 1,
        jdepends => $options->{jdepends},
        jname => $options->{jname},
        jstring => $jstring,
        jprefix => $options->{jprefix},
        jmem => $options->{jmem},
        output => $output_dir,);
    return($collect);
}

=back

=head2 C<Unicycler_Filter_Depth>

 Parse unicycler contig headers, extract relative coverage, and filter.

 Filter (for the moment only) a unicycler assembly by the ratio of the
 highest depth observed vs. the depth of each contig.  If that ratio
 falls below options->{coverage} then that contig should be dropped.
 The remaining contigs should be depth normalized to <=1. If there is
 only 1 contig, just return it with depth set to 1.  This should be
 trivially improved to handle other assembly methods by using the
 coverage calculation script above.

=over

=item C<Arguments>

 input(required): Fasta assembly from unicycler to test.
 output(''): Set of filtered contigs.
 coverage(0.2): Minimal relative coverage ratio allowed.
 output_log(''): Location to write the log of the tasks performed.

=item C<Invocation>

> cyoa --task ass --method unicyclerfil --input unicycler.fasta

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
        jcpu => 1,
        jdepends => $options->{jdepends},
        comment => $comment,
        jmem => $options->{jmem},
        jname => qq"filter_depth_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        language => 'perl',
        output => $output,
        output_log => $output_log,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($depth_filtered);
}

=back

=head2 C<Unicycler_Filter_Worker>

 This does the work for Unicycler_Filter_Depth().

 Filter (for the moment only) a unicycler assembly by the ratio of the
 highest depth observed vs. the depth of each contig.  If that ratio
 falls below options->{coverage} then that contig should be dropped.
 The remaining contigs should be depth normalized to <=1. If there is
 only 1 contig, just return it with depth set to 1.  This should be
 trivially improved to handle other assembly methods by using the
 coverage calculation script above.

=over

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
        coverage => 0.2,
        jcpu => 1,
        output_log => '',
        output => '',
        required => ['input',],);

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
    my $highest_coverage = 1;
    my $lowest_coverage = 1000000;
    my $initial_coverage_data = {};
    my $final_coverage_data = {};
  INITIAL: while (my $contig_seq = $input_contigs->next_seq()) {
      my $number_contigs++;
      my $current_id = $contig_seq->id();
      my $current_desc = $contig_seq->desc();
      my ($contig_length, $contig_coverage);
      my $circular = 1;
      if ($current_desc =~ /cov=/) {
          ## Then this was done with shovill.
          my ($contig_id, $contig_length_string, $contig_coverage_string, $stuff) = split(/\s+/, $current_desc);
          ($contig_length, $contig_coverage) = $current_desc =~ /len=(\d+)? cov=(.* )/;
          ($contig_coverage, $stuff) = split(/ /, $contig_coverage);
      } elsif ($current_desc =~ /depth=/) {
          ## This assembly was via spades/unicycler
          my ($contig_length_string, $contig_coverage_string, $circular) = split(/\s+/, $current_desc);
          ($contig_length, $contig_coverage) = $current_desc =~ /^length=(.*)? depth=(.*)?x/;
      } else {
          ## If we cannot figure out what to do, just set the length to 100k and coverage to 1.
          print $log "Unable to determine the assembler used, not filtering.\n";
          ($contig_length, $contig_coverage) = (100000, 1);
      }
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

=back

=head2 C<Shovill>

 Perform a shovill/spades assembly. https://github.com/tseemann/shovill

 I am learning that spades is the basis for more than a couple of
 assemblers. Shovill is one of them, and very feature-full.  Shovill
 has a lot of interesting options, this only includes a few at the
 moment.

=over

=item C<Arguments>

 input(required): Filtered/corrected/trimmed reads.
 depth(40): Used to filter over/under represented contigs.
 jmem(12): Expected memory usage.
 jprefix('13'): Added to the beginning of the job name.
 arbitrary(''): Add arbitrary arguments here.
 modules('shovill'): Use this environment module.

=item C<Invocation>

> cyoa --task ass --method shov --input r1.fastq.xz:r2.fastq.xz

=cut
sub Shovill {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        depth => 40,
        jmem => 12,
        jprefix => '13',
        arbitrary => '',);
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
        comment => $comment,
        jcpu => $options->{jcpu},
        jdepends => $options->{jdepends},
        jname => qq"shovill_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => $options->{jmem},
        output => qq"${output_dir}/final_assembly.fasta",
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($shovill_job);
}

=back

=head2 C<Trinity>

 Run the trinity de-novo transcriptome assembler.  doi:10.1038/nbt.1883.

 This function invokes trinity along with its default set of
 post-processing tools.

=over

=item C<Arguments>

 input(required): Fastq input files, trimmed and presumably corrected.
 contig_length(600): Minimum contig length for a transcript.
 jmem(80): Trinity is quite memory intensive.
 jprefix(60): Prefix for job name.
 modules('trinity'): Use this environment module.

=item C<Invocation>

> cyoa --task ass --method trinity --input assembly.fasta

=cut
sub Trinity {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        clean => 1,
        contig_length => 600,
        jcpu => 8,
        jmem => 80,
        jprefix => '60',);
    my $job_name = $class->Get_Job_Name();
    my $output_dir = qq"outputs/$options->{jprefix}trinity_${job_name}";
    my $input_string = '';
    if ($options->{input} =~ /\:|\;|\,|\s+/) {
        my @in = split(/\:|\;|\,|\s+/, $options->{input});
        $input_string = qq"--left <(less $in[0]) --right <(less $in[1]) ";
    } else {
        $input_string = qq"--single <(less $options->{input}) ";
    }
    my $trim_flag = '';
    if ($options->{trim}) {
        $trim_flag = '--trimmomatic';
    }
    my $arbitrary_args = $class->Passthrough_Args(arbitrary => $options->{arbitrary});
    my $comment = '## This is a trinity submission script.';
    my $jstring = qq!mkdir -p ${output_dir}
Trinity --seqType fq --min_contig_length $options->{contig_length} --normalize_reads ${arbitrary_args} \\
  ${trim_flag} --max_memory $options->{jmem}G --CPU $options->{jcpu} \\
  --output ${output_dir} \\
  ${input_string} \\
  2>${output_dir}/trinity_${job_name}.stderr \\
  1>${output_dir}/trinity_${job_name}.stdout
!;
    if ($options->{clean}) {
        $jstring .= qq"
rm -rf chrysalis insilico_read_normalization read_partitions \\
      __salmon_filt.chkpts salmon_outdir Trinity.tmp.fasta.salmon.idx \\
      scaffolding_entries.sam *.cmds *.ok *.fa *.finished *.timing
";
    }
    my $trinity = $class->Submit(
        comment => $comment,
        jcpu => $options->{jcpu},
        jdepends => $options->{jdepends},
        jname => qq"trin_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => $options->{jmem},
        output => qq"${output_dir}/Trinity.tsv",
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    my $rsem = $class->Bio::Adventure::Assembly::Trinity_Post(
        %args,
        jcpu => $options->{jcpu},
        jdepends => $trinity->{job_id},
        jname => qq"1trin_rsem",
        input => $options->{input},);
    my $trinotate = $class->Bio::Adventure::Annotation::Trinotate(
        %args,
        jdepends => $trinity->{job_id},
        jname => qq"2trinotate",
        input => qq"${output_dir}/Trinity.fasta",);
    $trinity->{rsem_job} = $rsem;
    $trinity->{trinotate_job} = $trinotate;
    return($trinity);
}

=back

=head2 C<Trinity_Post>

 Perform some of the post-processing tools provided by trinity.

 Trinity comes with some scripts to grade the quality of the
 transcripts it generated.  This invokes a few of them.

=over

=item C<Arguments>

 input(required): The fasta assembly from trinity.
 jmem(24): Expected memory requirement.
 jname(trin_rsem): Name on the cluster.
 jprefix('61'): Job prefix.
 modules('rsem'): RSEM is used to quanify the transcripts.

=item C<Invocation>

> cyoa --task ass --method trinitypo --input trinity_assembly.fasta

=cut
sub Trinity_Post {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jmem => 24,
        jname => 'trin_rsem',
        jprefix => '61',);
    my $job_name = $class->Get_Job_Name();
    my $trinity_out_dir = dirname($options->{input});
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
  2>trinpost_stats.stderr \\
  1>trinpost_stats.stdout &

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
        jstring => $jstring,
        output => qq"${trinity_out_dir}/RSEM.isoform.results",
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($trinpost);
}

=back

=head2 C<Unicycler>

 Use unicycler to assemble bacterial/viral reads. doi:10.1371/journal.pcbi.1005595

 Yet another Spades-based assembler, this one has some nice pre and
 post-processing tools written into it, which means that I do not need
 to figure out the arguments for pilon and stuff!

=over

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

=item C<Invocation>

> cyoa --task ass --method unicycler --input r1.fastq.xz:r2.fastq.xz

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
        jprefix => '13',);
    my $job_name = $class->Get_Job_Name();
    my $outname = basename(cwd());
    my $output_dir = qq"outputs/$options->{jprefix}unicycler";
    my $input_string = '';
    my $shovill_input = '';
    my $ln_string = '';
    my $backup_string = '';
    if ($options->{input} =~ /\:|\;|\,|\s+/) {
        my @in = split(/\:|\;|\,|\s+/, $options->{input});
        $input_string = qq" -1 ${output_dir}/r1.fastq.gz -2 ${output_dir}/r2.fastq.gz";
        $shovill_input = qq" --R1 ${output_dir}/r1.fastq.gz --R2 ${output_dir}/r2.fastq.gz ";
        $ln_string = qq"less $in[0] > ${output_dir}/r1.fastq
less $in[1] > ${output_dir}/r2.fastq\n";
        if ($in[0] =~ /\.fastq$/) {
            $input_string = qq" -1 ${output_dir}/r1.fastq -2 ${output_dir}/r2.fastq ";
            $shovill_input = qq" --R1 ${output_dir}/r1.fastq --R2 ${output_dir}/r2.fastq ";
            $ln_string = qq"cp $in[0] ${output_dir}/r1.fastq
cp $in[1] ${output_dir}/r2.fastq\n";
        }
    } else {
        $input_string = qq" -1 ${output_dir}/r1.fastq.gz ";
        $shovill_input = qq" --R1 ${output_dir}/r1.fastq.gz ";
        $ln_string = qq" less $options->{input} > ${output_dir}/r1.fastq ";
        if ($options->{input} =~ /\.fastq$/) {
            $input_string = qq" -1 ${output_dir}/r1.fastq ";
            $shovill_input = qq" --R1 ${output_dir}/r1.fastq ";
            $ln_string = qq"cp $options->{input} ${output_dir}/r1.fastq\n";
        }
    }
    my $comment = '## This is a unicycler submission script.';
    my $stdout = qq"${output_dir}/unicycler_${outname}.stdout";
    my $stderr = qq"${output_dir}/unicycler_${outname}.stderr";
    my $jstring = qq!mkdir -p ${output_dir}
${ln_string}
if unicycler $options->{arbitrary} \\
  --mode $options->{mode} \\
  --min_fasta_length $options->{min_length} \\
  ${input_string} \\
  -o ${output_dir} \\
  2>${stderr} \\
  1>${stdout}; then
  echo "unicycler passed."
else
  if shovill --nocorr --force \\
    --outdir ${output_dir} \\
    ${shovill_input} \\
    --minlen $options->{min_length} \\
    2>>${stderr} \\
    1>>${stdout}; then
    echo "Second unicycler attempt passed."
    mv ${output_dir}/contigs.fa ${output_dir}/assembly.fasta
  else
    echo "Both unicycler attempts failed."
    exit 1
  fi
fi

mv ${output_dir}/assembly.fasta ${output_dir}/${outname}_final_assembly.fasta
rm -f r1.fastq.gz r2.fastq.gz r1.fastq r2.fastq
ln -sf ${output_dir}/${outname}_final_assembly.fasta unicycler_assembly.fasta
!;
    my $unicycler = $class->Submit(
        comment => $comment,
        cpus => $options->{jcpu},
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
        jname => qq"unicycler_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        output => qq"${output_dir}/${outname}_final_assembly.fasta",
        output_gfa => qq"${output_dir}/assembly.gfa",
        output_log => qq"${output_dir}/unicycler.log",
        stdout => $stdout,
        stderr => $stderr);
    return($unicycler);
}

=back

=head2 C<Velvet>

 Perform a velvet assembly and pass it to ragoo.
 velvet: doi:10.1371/journal.pone.0008407
 ragoo: doi:10.1186/s13059-019-1829-6

 Yet another kmer-based assembler.

=over

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
        jmem => 24,);
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
    { grep velvetg parameters || test \$? = 1; } \|
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
        comment => $comment,
        jcpu => $options->{jcpu},
        jdepends => $options->{jdepends},
        jname => qq"velveth_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => $options->{jmem},
        output => qq"$output_dir/Sequences",
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($velvet);
}

=back

=head1 AUTHOR - atb

  Email  <abelew@gmail.com>

=head1 SEE ALSO

L<trinity> L<trinotate> L<transdecoder> L<velvet> L<ragoo>

=cut

1;
