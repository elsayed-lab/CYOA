=head1 NAME

    Bio::Adventure::Arguments - a place to document how I handle arguments in this code.

=head1 SYNOPSIS

    There are multiple ways to set parameters used by the various
    tools in Bio::Adventure and cyoa, and (too) many parameters
    available.

    If the CLI is invoked without the options 'task' and 'method', it
    creates a menu-driven interface which attempts to guide the user
    through the process of submitting jobs.  Any double dashed
    keywords will get set via ReadLine, thus:

    > cyoa --task pipe --meth phagea --host_filter 1 \
       --input r1.fastq.xz:r2.fastq.xz

    will pick up the argument 'phagea' and see that it has only one
    possible match in the set of known TODOs (phageassembly) and
    therefore flag that as something TO DO.  Similarly, I suspect one
    could just type --ho and ReadLine will see that 'host_filter' is
    the only argument which starts with 'ho' and therefore set it.

    When performed via the API, arguments sent to functions directly
    get precedence over those set for an object:

    my $cyoa = Bio::Adventure->new(species => 'mtuberculosis_rv', --host_filter 1);
    $cyoa->Hisat2(input => 'r1.fastq.xz');
    $cyoa->Hisat2(species => 'hg38_100');

    In this instance, the constructor set the species, but later it
    was set for a specific invocation for hisat; since host_filter
    is 1, the second invocation of hisat will run against the
    unaligned reads from the first invocation.  Hmm, now that I typed
    that sentence, I realized that 'host_filter' is a particularly
    poor name and will be changed to 'filter'.

=head1 DESCRIPTION

    As the above suggests, I attempted to make this as flexible as
    possible while still following the guidelines of L<Moo>.

    The following attempts to write down the purpose of the various
    arguments currently defined, categorized by type of analysis.

=over 4

=item C<General Options>: The options in this section are used all over the place.

  task - Top-level task to perform.  Most of the time this can just be
    filled in arbitrarily with no effect.  However, some tools will
    take note of it and adjust their parameters accordingly.  Most
    notably, when trimming one may set task to 'tnseq', 'assemble', or
    'riboseq' in order to set some options which I think are more
    appropriate for the task at hand.
  method - Which function to actually call?  The set of available
    methods are defined in L<Bio::Adventure::Get_TODOs()> and
    included below.
  basedir - Base directory from which all other directories inherit. (cwd)
  input - Generic input file, used by most functions.
  libdir(${HOME}/libraries) - The location of the various fasta/gff
    and index files.
  libtype(genome) - The type of library for which to search.  This is
    almost always 'genome', and will most likely be set to rRNA when
    filtering out rRNA sequences from data, or perhaps 'contaminants'.
  output - Generic output file, used by most functions.
  species - The species used for the given analysis, this will be used
    to look in the I<libdir> for a genome in the files
    I<species>.fasta, a gff file in the same directory with the same
    prefix, indices for the various mappers/aligners with the same
    prefix within the I<libdir>/I<libtype>/indices directory.  There
    are a few other places this gets used, notably for codon
    adaptation index calculations or species-specific hmm libraries.
    It should also be used for any peptide libraries when performing
    comet/DIA/etc searches for proteomics data.


  raw_dir - Location of the raw data when it might lie elsewhere.
  type - Used in a few places, genomic vs ribosomal, tnseq vs rnaseq, etc.
  suffixes - Set of suffixes to remove from files when looking for basenames. (.fastq, .gz, .fasta, .sam, .bam, .count, .csfast, .qual)
  debug - Print debugging information? (false)
  help - Print help information.

=item C<Cluster Specific Options>

  cluster - Use a computing cluster or just invoke scripts one at a time?
  jdepends - A string describing the set of dependencies for a job.
  jname - A machine-readable name for each job.
  jprefix - A prefix number for each job, to make reading successions easier.
  jstring - The string comprising the actual job to be run.
  job_id - The job ID returned by Torque/Slurm for a job.
  jmem - Memory in Gb to request from the cluster.
  jwalltime - Amount of walltime to request from the cluster. (10:00:00)
  queue_array_string - When running an array of jobs, this holds the requisite formatted string.
  language - Choose if the cluster options should be appropriate for a shell or perl script. (bash)
  qsub_args - Default arguments when using torque. (-j oe -V -m n)
  sbatch_depends - Dependency string on a slurm cluster. (afterok:)
  sbatch_dependsarray - Dependency string on a slurm for array jobs. (afterok:)
  qsub_depends - Dependency string on a torque cluster. (depend=afterok:)
  qsub_dependsarray - Dependency string on a torque cluster for array jobs. (depend=afterokarray:)
  queue - Default qos/queue to send jobs. (workstation)
  queues - Set of available queues. (throughput, workstation, long, large)
  shell - Default shell to request when using torque/slurm. (/usr/bin/bash)

=item C<Alignment Arguments>

  align_blast_format(5) - Which blast output format should be used?  Some
    have more robust parsers than others.  5 is blastxml.
  align_jobs(40) - When splitting up a large number of alignments into
    smaller pieces, this defines how many concurrent alignments to perform.
  align_parse(1) - When turned on, the collected alignments will be
    parsed into a relatively simple table describing how many hits
    were observed for each query and some of the statistics for each
    alignment as a set of TSV columns.
  best_only(0) - When this is turned on, then parsed alignments will
    have only 1 line per query sequence, describing the best hit.
  blast_params( -e 10 ) - Set my favorite blast parameters.  This may
    warrant moving into the catch-all 'arbitrary' argument.
  blast_tool(blastn) - Pick the appropriate blast executable here.  I
    think I added a little logic which should warn the user if one is
    using a blastp for nucleotide sequence or vice-versa?  At the very
    least, my blast invocations should remind the forgetful (me) of
    which methods to use in which situations.

  fasta_format - Which fastq output format should we use? 0 is tabular. (0)


  peptide - Explicitly tell blast if this is (or not) a peptide search. (F)

  library - Choose the name of the blast/fasta database to search against.
  evalue - E-value cutoff when performing blast/fasta searches.  This is not currently entirely used. (1)
  identity - Percent identical cutoff when doing blast/fasta.  Not used as it should be. (70)
  fasta_args - Default arguments when invoking fasta36 programs. (-b 20 -d 20)
  fasta_tool - Which fasta tool to invoke? (ggsearch36)

=item C<Mapping Arguments>

  mapper - Set the default mapping program. (salmon)
  bt_args - Some default arguments when using bowtie. (--very-sensitive -L 14)
  bt2_args - Similar for bowtie2.
  bt_type - Choose a specific type of bowtie mapping to perform. (v0M1)
  btmulti - Or perform bowtie with all of the chosen options in bt_args? (false)
  stranded - Choose options for stranded libraries (primarily kallisto) (false)

=item C<Counting Arguments>

  htseq_args - Some default arguments for htseq. (--order=name --idattr=gene_id --minaqual=10
     --type exon --stranded=yes --mode=union)
  htseq_stranded - Use htseq in stranded mode? (no)
  htseq_type - Which gff type (3rd column) to use in htseq? (exon)
  htseq_id - Which gff ID (tag in the last column) to use in htseq? (gene_id)
  mi_genome - Database of miRNA sequences for counting small RNA sequences.
  mature_fasta - Database of miRNA mature sequences when counting small RNA sequences.
  mirbase_data - Location of the full miRbase data when counting small RNA.

=item C<Conversion Arguments>

  bamfile - Where to put output bam data?  Used in a few places.
  feature_type - Primarily used when looking at gff files, ergo htseq and converting gff to other formats.
  taxid - Used for making gaf files.
  species - When converting from sam->bam, use this to pull the appropriate genome fasta file.
  paired - Explicitly set whether a sam file is from paired data and therefore filter unpaired reads.

=item C<Preparation Arguments>

  csv_file - A csv file containing the metadata for multiple samples. (all_samples.csv)
  sampleid - Chosen ID for a given sample.

=item C<TNSeq Arguments>

  index_file - File containing a key/value pairing of sample IDs and TNSeq indexes. (indexes.txt)
  runs - How many gumbel runs to apply in essentiality? (1000)

=item C<Riboseq Arguments>

  orientation - When looking at reads vs. starts/ends, count start->end or end->start? (start)
  riboasite,psite,esite - Count A/P/E site reads. (1,1,1)
  ribominsize/maxsize - Range for riboseq sizes. (24,36)
  ribominpos/maxpos - Relative positions when counting vs. start/stop codons. (-30, 30)
  riboanchor - Anchor reads with respect to the start/end. (start)
  ribosizes - Set of sizes to count. (25,26,27,28,29,30,31,32,33,34)

=item C<Trimming Arguments>

  maxlength - Maximum length to keep when trimming sequence data (usually TNSeq). (42)
  minlength - Minimum length to keep when trimming sequence data (usually TNSeq). (8)
  phred - Score range to work in, also used for some QC tasks. (33)
  qual - cutadapt specific quality string.

=item C<QC Arguments>

  paired - Primarily Fastqc, boolean check if the data is paired. (false)

=item C<Phylogeny Arguments>

  outgroup - Name of a chosen phylogenetic outgroup.
  starting_tree - Name of a starting tree for improving phylogenetic trees.

=item C<SNP Arguments>

  trim - Trim raw sequence with options more appropriate for SNP calling.
  varfilter - Filter variants by variance. (true)
  vcf_cutoff - Minimum reads for a variant call. (10)
  vcf_minpct - Minimum percent agreement for a variant call. (0.8)

=item C<Methods>  If anyone is interested, this works because I set
  every key of this hash to a GetOpt::Long parameter definition string
  which ends in a '+'.  That tells GetOpt::Long to increment this
  element of the hash reference provided, which is called 'todo'.
  Every key of the todo hash is a function name; so when GetOpt::Long
  decides that the method chosen matches any of the argument strings
  below, the todo value for the associated function name gets set to 1
  (or 2, etc, if one were to set the same option multiple times).
  Thus, when cyoa looks for a list of tasks to perform, it reads the
  todo list and runs any function it finds defined.

  Finally, note that one need type only as much of the following as is unique.

  "abricate+" - Abricate is a Torsten Seeman method for searching
    resistance genes and other features which may be of concern.
  "abyss+" - Abyss is a fascinating makefile-based assembler.
  "angsdfilter+" - Angsd is a pretty neat toolkit for population
    genetics data.  I do not understand more than 20% of what it can do.
  "annotatephage+" - Pass an assembly to a pipeline of tasks intended
    to help annotate phage genomes.  It essentially runs the back-half
    of the tasks in the 'phage_pipeline_outline.png' picture in 'share/'.
  "aragorn+" - This tool searches for likely tRNA sequences in a fasta file.
  "assemblycoverage+" - Calculate the coverage of every base in an assembly.
  "biopieces+" - Run a few biopieces tools to get a summary of a data set.
  "blastmerge+" - Merge together a set of Blast searches.
  "blastparse+" - Parse a table of Blast searches into a relatively simple table.
  "blastsplitalign+" - Split a query dataset into multiple pieces and
    concurrently search them with blast.
  "bowtie+" - Run bowtie1 against a genome.
  "bowtierrna+" - Run bowtie with some changes for rRNA (mostly for ribosome profiling).
  "bt2+" - Run bowtie2 against a genome.
  "btmulti+" - Run multiple parameter sets of bowtie.
  "bwa+" - Run bwa mappings against a genome.
  "caical+" - Calculate the codon adaptation indices for a set of genes against a species.
  "calibrate+" - (ribosome profiling) Given a set of calibrations of read lengths vs. the
    known starts/ends of every gene in a genome, set the beginning of
    the A/P/E sites of all reads on a per-length basis.
  "cgview+" - Use some default parameters to run cgview and visualize
    a circular genome.  Note that cgview can do a lot more than I have
    currently implemented.
  "classifyphage+" - This uses a set of (t)blastx searches of an
    unknown phage assembly vs. all of the ICTV reference genomes in
    the hopes of finding a likely classification.
  "collectassembly+" - The assembly pipelines produce a ridiculous
    number of outputs, this collects the most likely ones into one directory.
  "countstates+" - Count the putative state (wrt the ribosome) of all reads in a riboseq data set.
  "concat+" - Concatenate a set of blast/fasta searches.
  "cutadapt+" - Trim high-throughput data with cutadapt.
  "essentialitytas+" - Count up the TAs in a mariner-based TNSeq
    experiment into a .wig file appropriate for the TRANSIT methods.
  "extendkraken+" - Add sequences to an existing kraken2 database.
  "extracttrinotate+" - Summarize the peculiar outputs from trinotate.
  "splitalignfasta+" - Split up a fasta36 search into multiple pieces.
  "fastado+" - Perform an arbitrary fasta36 search.
  "fastaparse+" - Parse the output of a fasta36 search into something
    easier to read.
  "fastqct+" - Generate the fun pictures from fastqc.
  "fastqdump+" - Download SRA sequences via fastqdump.
  "featureextract+" - Extract features by keyword from a genbank file.
  "filterdepth+" - Filter an assembly by sequencing depth.
  "filterkraken+" - Filter out putative host-reads from an assembly
    via the most overrepresented species observed by kraken2.
  "freebayes+" - Quantify potential variants with freebayes.
  "gb2gff+" - Dump a genbank file into a set of fasta/gff files.
  "gff2fasta+" - Separate fasta-containing gff files.
  "glimmer+" - Run glimmer using the recommendations from the pdf manual.
  "glimmersingle+" - Run glimmer without training, essentially just
    getting long ORFs.
  "graphreads+" - My first attempt to plot reads!
  "gubbins+" - Use gubbins to perform some phylogenetic analyses (TBH
    I have never used this and just wrote it to help someone else, who
    told me it worked well)
  "gumbel+" - Run the pre-TRANSIT gumbel bayesian essentiality
    classifier on a wig file.
  "hisat+" - My favorite mapper, hisat2.
  "htmulti+" - Run htseq on an alignment using a few different feature types.
  "htseq+" - Run a single htseq invocation on an alignment.
  "indexsomethingsomething+" - Run the indexer for mappers/blast/etc.
  "interproscan+" - Run interproscan on a pile of (translated) sequences.
  "jellyfish+" - Invoke jellyfish on a sequence database with a few
    k-values and gather some summary statistics from the result.
  "kallisto+" - Quantify reads with kallisto.
  "kraken+" - Classify reads with kraken2.
  "mash+" - Extract pairwise distance information for a set of
    sequences.  I do not yet understand how mash works at all.
  "mergeannotations+" - Assembly specific, bring together annotations
    from multiple sources (interpro, trinotate, abricate, etc)
  "mergecds+" - Assembly specific, bring together CDS feature
    predictions from multiple sources (glimmer, prodigal, phanotate,
    prokka (which yes uses prodigal but with different parameters)).
  "mimap+" - Attempt to extract the mature miRNA species which are
    observed in a high-throughput dataset.
  "mpileup+" - Use samtools' mpileup, vcftools, and friends to look
    for variants.
  "phageterm+" - Search for DTR regions with phageterm.
  "phanotate+" - Hunt phage CDS features with phanotate.
  "phastaf+" - Another Torsten Seeman production, classify (pro)phages via diamond.
  "prodigal+" - Look for CDS features and optionally create a training set.
  "parseblast+" - Parse a blast search into a hopefully easier to read table.
  "parsebcf+" - Given a bcf file, create some tables and a modified
    copy of the original genome.
  "phagepromoter+" - Search phage assemblies for promoter sequences.
  "posttrinity+" - Run the various trinity post-assembly analyses
    suggested by their documentation (rsem, etc).
  "prokka+" - Yet another Torsten Seeman production, for automated annotation.
  "racer+" - Correct sequencer-based errors in a high-throughput data
    set.  I wonder if this would work for long read sequencing?  I am
    guessing no.
  "recatalog+" - Use a little bioperl to count up every RE site in an assembly.
  "resfinder+" - Search an assembly for drug resistance mutations/genes.
  "rhopredict+" - Search an assembly for rho termination sequences.
  "rgi+" - Another drug resistance search tool, it also does some
    other searches which I forget.
  "rsem+" - This provides a more thorough gene-quantification strategy
    than just counting hits with htseq.  I really think we should use
    it more, but I do not fully understand it and so tend to eschew it.
  "runessentiality+" - Run various essentiality tools (I think this is
    where I put my invocation of the full TNSeq TRANSIT tools).
  "sam2bam+" - Convert a sam alignment to sorted/indexed bam, along
    with some potential filtered versions.
  "salmon+" - I think of salmon as the big sister to kallisto.
  "terminasereorder+" - Reorder an assembly to put a putative
    terminase gene first.
  "shovill+" - Shovill provides a way to optimize the parameters used
    when performing a spades assembly.
  "slsearch+" - Perform an explicit search for spliced leader
    sequences using pattern matching.
  "snippy+" - This is another variant search tool, but I wrote it for
    someone else and have never actually learned how it works.
  "snpsearch+" - Run mpileup and parse the results.
  "snpratio+" - Parse the results from either freebayes or mpileup.
  "snpgenome+" - Write a new genome from a set of variants.
  "sortindexes+" - A poor-man's demultiplexer for the oddly-designed
    TNSeq libraries from Yoann.  Given the new library protocols, this
    is likely no longer needed.
  "splitalign+" - Split up an input and align it using fasta/blast
  "star+" - Use the STAR mapper instead of hisat/bowtie/bwa/salmon/kallisto/etc.
  "tacheck+" - Filter a mariner-based TNSeq dataset for correctly
    placed TAs.
  "tophat+" - Superceded by hisat, but still here, run tophat2.
  "tpp+" - Invoke the TNSeq TRANSIT preprocessing methods (e.g. bwa).
  "trainprodigal+" - Create a prodigal training set.
  "transdecoder+" - Extract CDS sequences with transdecoder (usually
    for trinity/trinotate)
  "trimomatic+" - My favorite sequence trimmer!
  "trinity+" - De-novo short-read assembler for transcriptome data.
  "trinotate+" - Run trinotate, optionally with alternate databases.
  "trnascan+" - One of a few(couple?) tRNA sequence search strategies.
  "unicycler+" - Similar to shovill, except this is optimized for
    bacterial/circular genomes. (I use it for phage sequences).
  "velvet+" - Yet another assembler, if I recall this one has some
    nice methods for guiding an assembly.
  "vienna+" - Perform a series of MFE/secondary structure calculations.
  "rosalindplus+" - Attempt to ensure that the rosalind(+) strand
    contains the majority of ORFs in an assembly.
  "psomethingsomething+" - Some(most,all?) of the pipelines start with a 'p'.

=back

=cut

1;
