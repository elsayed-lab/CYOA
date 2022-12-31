# NAME

    Bio::Adventure - A perl library to make preprocessing high-throughput data easier.

# SYNOPSIS

      Bio::Adventure.pm and the associated script cyoa seek to make it easier to create
      and submit template scripts to a computing cluster.

    A few likely commandline invocations and their purposes include:

    > cyoa --task pipe --method phage --input r1.fastq.xz:r2.fastq.xz

    This passes the paired end, xz compressed reads to a phage assembly
    and annotation pipeline, diagrammed below.[1]

    > cyoa --task pipe --method rnaseq --input r1.fastq.xz:r2.fastq.xz \
      --gff_type gene --gff_id gene_id --host_filter 1 --species hg38_100:lpanamensis_v50

    This performs a default (hisat2) rnaseq pipeline which trims the reads, runs hisat twice,
    first against the hg38_100 human assembly followed by Leishmania panamensis version 50 from
    the tritrypdb.  This of course assumes a gff/fasta file in the library directory, but will
    check for the indexes and create them if needed.  When host_filter is on, it will pass only
    the unmapped reads from hg38 to the panamensis mapping.  It converts the sam to
    sorted/indexed bam and generates a few alternate bam files (no mismatches, only paired reads),
    passes them to htseq-count and freebayes for variant searching.

    > cyoa --task map --method salmon --species mm38_100 \
      --input r1_trimmed.fastq.xz:r2_trimmed.fastq.xz

    This invokes salmon on trimmed reads using the mouse ensembl release.

    > cyoa --task align --method blast --library nr --input few_thousand_sequences.fasta

    This will split the fasta file into multiple pieces and pass them to blast on a cluster
    against the NR database and return back some summary tables.

    > cyoa --task align --method fasta --library database.fasta --input query.fasta

    This does much the same, but uses a fasta database and fasta36.

    > cyoa --task tnseq --method essentiality --input processed_data.wig --species spyogenes

    Invoke the DeJesus essentiality tool on a processed wig file.


    There are many more methods, tools, and command line options.  One way to explore
    them is to pull up the text menu interface via just `cyoa`.
    This is also an API for invoking these methods:

      use Bio::Adventure;
      my $hpgl = new Bio::Adventure;
      $hpgl->{species} = 'mmusculus';  ## There had better be a mmusculus.[fasta&gff] in $hpgl->{libdir}
      my $hpgl = new Bio::Adventure(species => 'scerevisiae', libdir => '/home/bob/libraries');
      $hpgl->{input} = 'test.fastq';  ## Also via with script.pl --input=test.fastq

      ## Run Fastqc on untrimmed sequence.
      my $bp = $hpgl->Fastqc()
      ## Generates test-trimmed.fastq from test.fastq
      my $trim = $hpgl->Trimomatic();
      ## Graph statistics from test-trimmed.fastq
      my $bp = $hpgl->Biopieces_Graph(depends => $trim->{pbs_id});
      ## Run bowtie1, convert the output to sorted/indexed bam, and run htseq-count against mmusculus
      ## bowtie1 outputs go (by default) into bowtie_out/
      my $bt = $hpgl->Bowtie();
      ## Run htseq-count with a different gff file, use 'mRNA' as the 3rd column of the gff, and 'locus_tag' as the identifier in the last column.
      my $ht = $hpgl->HTSeq(depends => $bt->{pbs_id}, input => $bt->{output}, htseq_gff => 'other.gff', htseq_type => 'mRNA', htseq_identifier => 'locus_tag'
      ## Run tophat
      my $tp = $hpgl->TopHat(depends => $trim->{pbs_id});
      ## or BWA
      my $bw = $hpgl->BWA(depends => $trim->{pbs_id});

# DESCRIPTION

    This library should write out cluster compatible job files for
    torque or slurm clusters and submit them to appropriate queues.
    It should also collect the outputs and clean up the mess.

## The TODO list

    The following list is the set of potential tasks addressable by cyoa categorized
    somewhat arbitrarily.  Each of the names in parentheses is one methodname which
    may be passed to cyoa either in full or abbreviated to its least unique string.

    1.  Alignment: (TODO: add blat, hmmer, tcoffee, etc)
        * (blastsplit): Split the input sequence into subsets and align with blast.
        * (fastasplit): Split the input sequence into subsets and align with fasta36.
        * (concat): Merge split searches into a single set of results.
        * (fastaparse): Parse fasta36 output into a reasonably simple table of hits.
        * (blastparse): Parse blast output into a reasonably simple table of hits.
        * (fastamerge): Merge and Parse fasta36 output into a reasonably simple table of hits.
        * (blastmerge): Merge and Parse blast output into a reasonably simple table of hits.
    2.  Annotation (TODO: there are a bunch of missing tasks here)
        * (abricate): Search for Resistance genes across databases.
        * (classify_virus): Use ICTV data to classify viral sequences/contigs.
        * (extend_kraken): Extend a kraken2 database with some new sequences.
        * (interproscan): Use interproscan to analyze ORFs.
        * (kraken2): Taxonomically classify reads.
        * (phageterm): Invoke phageterm to hunt for likely phage ends.
        * (phastaf): Invoke phastaf to attempt classifying phage sequence.
        * (rhopredict): Search for rho terminators.
        * (terminasereorder): Reorder an assembly based on the results of a blast search.
        * (prokka): Invoke prokka to annotate a genome.
        * (resfinder): Search for antimicrobial resistance genes.
        * (rgi): Search for resistance genes with genecards and peptide fasta input.
        * (mergeannotations): Merge annotations into a genbank file.
    3. Assembly (TODO: Some of these probably should be elsewhere?, also Spades)
        * (abyss): Run abyss to create a new assembly.'
        * (assemblycoverage): Calculate Assembly coverage across contigs.
        * (extract_trinotate): Extract the most likely hits from Trinotate.
        * (filterdepth): Filter contigs based on sequencing depth.
        * (transdecoder):  Run transdecoder on a putative transcriptome.
        * (trinotate): Perform de novo transcriptome annotation with trinotate.
        * (trinity): Perform de novo transcriptome assembly with trinity.
        * (trinitypost): Perform post assembly analyses with trinity.
        * (velvet): Perform de novo genome assembly with velvet.
        * (unicycler): Perform de novo assembly with unicycler/spades.
        * (shovill): Perform the shovill pre/post processing with spades.
        * (terminasereorder): Reorder an existing assembly to the most likely terminase.
        * (collect_assembly): Copy files generated by an assembly into one directory.
        * (rosalind_plus): Make sure the Rosalind strand has the most CDS.
    4. Conversion (I think I have some phylogenetic converters which should be added, notably for cactus)
        * (sam2bam): Convert a sam mapping to compressed/sorted/indexed bam file(s).
        * (gb2gff): Convert a genbank flat file to gff/fasta files.
        * (gff2fasta): Convert a gff file to a fasta file.
    5. Counting (should caical be here?)
        * (htseq): Count mappings with htseq.
        * (htmulti): Use different option sets for counting with htseq across feature types.
        * (jellyfish): Perform a kmer count of some data along with some post-processing.
        * (mash): Use mash to count pairwise distances among sequences.
        * (mimap): Count (im)mature miRNA species in high-throughput data.
        * (countstates): Count ribosome positions in ribosome profiling data.
        * (slsearch): Count frequency of spliced leader (or an arbitrary) sequences.
    6.  Feature Prediction (TODO: the crispr-cas9 finder)
        * (aragorn): Search for tRNAs with aragorn.'
        * (glimmer): Look for ORFs in bacterial/viral sequence.
        * (glimmersingle): Use glimmer to search for ORFs without training.
        * (phanotate): Look for ORFs in bacterial/viral sequence.
        * (prodigal): Look for ORFs in bacterial/viral sequence.
        * (rhopredict): Search for rho terminators.
        * (trnascan): Search for tRNA genes with trnascan.
        * (trainprodigal): Train prodgial using sequences from a species/strain.
    7. Mapping
        * (bowtie): Map trimmed reads with bowtie1 and count with htseq.
        * (bt2): Map trimmed reads with bowtie2 and count with htseq.
        * (hisat): Map trimmed reads with hisat2 and count with htseq.
        * (btmulti): Map trimmed reads and count using multiple bowtie1 option sets.
        * (bwa): Map reads with bwa and count with htseq.
        * (kallisto): Pseudo-align and count reads using kallisto.
        * (salmon): Pseudo-align and count reads using salmon.
        * (star): Pseudo-align and count reads using STAR.
        * (mimap): Attempt to map reads explicitly to mature miRNA species.
        * (rrnabowtie): Map rRNA reads using bowtie1.
        * (rsem): Quantify reads using rsem.
        * (tophat): Map reads using tophat2 and count with htseq.
    7. Indexing
        * Basically, take the above and prefix it with 'index' to make the various indices.
        * (extend_kraken): Extend a kraken2 database with some new sequences.
    8. Phage-specific analyses (TODO, some are missing, some are redundant)
        * (filterkraken): Filter out host sequences using a kraken report.
        * (classifyphage): Use ICTV data to classify a phage assembly.
        * (phageterm): Invoke phageterm to search for terminal repeats.
        * (phagepromoter): Search for phage promoters.
        * (phastaf): Search for phage regions in arbitrary assemblies.
        * (restriction): Search for restriction sites.
        * (terminasereorder): Reorder a phage assembly to the putative terminase.
        * (caical): Calculate codon adaptation index of a phage vs. some host.
        * (recatalog): Create a catalog of restriction enzyme hits.
    9. Phylogeny (TODO: Add phylip etc and the MSA methods)
        * (gubbins): Run Gubbins with an input msa.
    10. Population Genetics
        * (angsd): Get the set of filtered variants given a list of bam files.
    11. Preparation (TODO: I have some entrez downloaders now which probably would be good here)
        * (fastqdump): Download data from the SRA
    12. QA
        * (biopieces): Use biopieces to graph some metrics of the data.
        * (cutadapt): Perform adapter trimming with cutadapt.
        * (fastqc): Use fastqc to check the overall quality of the raw data.
        * (racer): Perform sequence correction with hitec/RACER.
        * (trimomatic): Perform adapter trimming with Trimomatic.
    13. Ribosome profiling
        * (biopieces): Make some plots of the demultiplexed/trimmed reads.
        * (cutadapt): Use cutadapt to remove riboseq adapters.
        * (rrnabowtie): Filter away rRNA reads using bowtie1.
        * (btmulti): Use bowtie1 to find putative ribosomal positions.
        * (calibrate): Empirically calibrate the positions of the a/p/e sites of the ribosomes.
        * (countstates): Count the positions of a/p/e/etc sites across the transcriptome.
        * (graphreads): Plot the coverage of ribosomes across the genome by ORF.
    14. Variant searching (TODO: Maybe add a separate deduplication step?  Add the variant categorizer)
        * (trim): Trim sequence with an additional rule to remove the first 10 nucleotides.
        * (bwa): Map reads with bwa and count with htseq.
        * (bowtie): Map trimmed reads with bowtie1 and count with htseq.
        * (bt2): Map trimmed reads with bowtie2 and count with htseq.
        * (freebayes): Use freebayes to create a vcf file and filter it.
        * (hisat): Map trimmed reads with hisat2 and count with htseq.
        * (snpsearch): Use mpileup to create a vcf file and filter it. (bam input)
        * (snpratio): Count the variant positions by position and create a new genome. (bcf input)
        * (snp): Perform alignments and search for variants with mpileup. (fastq input)
        * (snippy): Invoke snippy. (fastq and genbank inputs)
    15. TNSeq (TODO: Add the transit runner)
        * (sortindex): Demultiplex raw reads based on the peculiar TNSeq indexes.
        * (cutadapt): Use cutadapt to remove the odd tnseq adapters.
        * (tacheck): Make certain that all reads have a leading or terminal TA.
        * (biopieces): Make some plots of the demultiplexed/trimmed reads.
        * (essentialityta): Count the hits/TA in preparation for essentiality.
        * (runessentiality): Run the essentiality suite of tools.
        * (gumbel): Run the essentiality suite of tools on the ta counts.
        * (tpp): Run the transit preprocessing script.
    16.  Visualization (TODO: Clean up my circos invocation and add it)
        * (cgview): Invoke cgview to visualize a genome.

# Installation

* perl Build.PL
* Build installdeps
* Build install

The only likely dependency which causes trouble is Bio::DB::Sam
because it looks for a samtools < version 1.0 libbam.a

# Footnotes

1. [share/phage_annotation_pipeline.png](share/phage_annotation_pipeline.png)

# AUTHOR - atb

Email abelew@gmail.com

# SEE ALSO

    L<Bio::Seq> L<autodie> L<local::lib> L<Bio::Adventure::RNASeq_Aligners>
    L<Bio::Adventure::RNASeq_Assembly> L<Bio::Adventure::RNASeq_Count> L<Bio::Adventure::RNASeq_Aligners> L<Bio::Adventure::RNASeq_QA>
    L<Bio::Adventure::RNASeq_Trim> L<Bio::Adventure::Align_Blast> L<Bio::Adventure::Align_Fasta> L<Bio::Adventure::Align>
    L<Bio::Adventure::Cleanup> L<Bio::Adventure::Compress> L<Bio::Adventure::Convert> L<Bio::Adventure::PBS>
    L<Bio::Adventure::Prepare> L<Bio::Adventure::Riboseq> L<Bio::Adventure::SeqMisc> L<Bio::Adventure::SNP> L<Bio::Adventure::Status>
    L<Bio::Adventure::TNSeq>
