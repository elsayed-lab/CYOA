package Bio::Adventure::Config;
use Exporter;
@EXPORT = qw"Get_Modules";
@ISA = qw"Exporter";

sub Get_Modules {
    my ($package, $filename, $line, $subroutine, $hasargs,
        $wantarray, $evaltext, $is_require, $hints, $bitmask, $hinthash) = caller(1);
    my %module_data = (
        'Abricate' => {
            modules => ['any2fasta', 'abricate', 'blast', 'blastdb'],
            exe => 'abricate' },
        'Abyss' => {
            modules => 'abyss' },
        'Angsd_Filter' => {
            modules => 'angsd' },
        'Aragorn' => {
            modules => 'aragorn', exe => 'aragorn' },
        'Assembly_Coverage' => {
            modules => ['hisat2', 'samtools', 'bbmap', 'cyoa'], },
        'Bacphlip' => {
            modules => 'bacphlip', exe => 'bacphilp' },
        'Biopieces_Graph' => {
            modules => ['biopieces'] },
        'Bowtie' => {
            modules => 'bowtie1' },
        'Bowtie2' => {
            modules => 'bowtie2' },
        'BT1_Index' => {
            modules => ['bowtie1'] },
        'BT2_Index' => {
            modules => ['bowtie2'] },
        'BWA' => {
            modules => 'bwa' },
        'BWA_Index' => {
            modules => ['bwa'] },
        'Caical' => {
            modules => ['caical', 'cyoa'] },
        'Casfinder' => {
            modules => ['casfinder',], exe => 'casfinder.sh' },
        'Casoff' => {
            modules => 'casoff', exe => 'cas-offinder' },
        'CGView' => {
            modules => ['cgview'], },
        'Check_Blastdb' => {
            modules => ['blast'], },
        'Classify_Phage' => {
            modules => ['blast', 'blastdb', 'cyoa' ], exe => 'tblastx' },
        'Cogent' => {
            modules => ['cogent'], exe => 'cogent'},
        'Collect_Assembly' => { modules => 'cyoa' },
        'Consolidate_TAs' => {
            modules => 'cyoa', },
        'Cutadapt' => {
            modules => ['cutadapt', 'cyoa'], exe => 'cutadapt'},
        'Essentiality_TAs' => {
            modules => 'cyoa', },
        'Extend_Kraken_DB' => {
            modules => ['kraken'], exe => ['kraken2'] },
        'Fastqc' => {
            modules => ['fastqc'] },
        'Fastq_Dump' => {
            modules => ['sra', 'cyoa'], },
        'Fastp' => {
            modules => ['fastp'], exe => 'fastp' },
        'Filter_Host_Kraken' => {
            modules => ['kraken', 'cyoa', 'hisat2', 'htseq', 'samtools'], },
        'Freebayes_SNP_Search' => {
            modules => ['gatk', 'freebayes', 'samtools', 'bcftools', 'vcftools'],
            exe => ['gatk', 'freebayes'], },
        'Gb2Gff' => { modules => 'cyoa' },
        'Generate_Samplesheet' => {
            modules => 'R' },
        'Glimmer' => { modules => ['glimmer', 'cyoa'], exe => 'glimmer3' },
        'Glimmer_Single' => {
            modules => ['glimmer', 'cyoa'], exe => 'glimmer3' },
        'Gubbins' => { modules => 'gubbins' },
        'Guess_Strand' => { modules => 'cyoa' },
        'Hisat2' => {
            modules => 'hisat2', exe => 'hisat2'},
        'Hisat2_Index' => {
            modules => ['hisat2'], exe => ['hisat2'], },
        'HT_Multi' => {
            modules => 'htseq' },
        'HTSeq' => {
            modules => 'htseq', exe => 'htseq-count' },
        'Interpro_Long2Wide' => { modules => 'cyoa' },
        'Interproscan' => {
            modules => ['interproscan', 'cyoa'], exe => 'interproscan.sh' },
        'Jellyfish' => {
            modules => ['jellyfish', 'cyoa'], exe => 'jellyfish' },
        'Kallisto' => {
            modules => 'kallisto', exe => 'kallisto' },
        'Kallisto_Index' => {
            modules => ['kallisto'], exe => ['kallisto'], },
        'Kraken' => {
            modules => 'kraken', exe => 'kraken2', },
        'Kraken_Best_Hit' => {
            modules => ['kraken', 'cyoa'] },
        'Make_Blast_Job' => {
            modules => ['blastdb', 'blast', 'cyoa'], exe => 'tblastx' },
        'Make_Codon_Table' => {
            modules => ['cyoa'] },
        'Mash' => {
            modules => 'mash', exe => 'mash' },
        'Merge_Annotations' => {
            modules => ['ncbi_tools/6.1', 'cyoa'], exe => 'tbl2asn' },
        'Merge_CDS_Predictions' => {
            modules => ['ncbi_tools/6.1', 'cyoa'], exe => 'tbl2asn' },
        'Merge_Parse_Blast' => {
            modules => ['cyoa'], },
        'Mpileup' => {
            modules => 'samtools' },
        'Mpileup_SNP_Search' => {
            modules => ['libgsl/2.7.1', 'libhts/1.13', 'samtools/1.13', 'bcftools', 'vcftools'],
            exe => ['samtools', 'bcftools'], },
        'OrthoFinder' => {
            modules => ['orthofinder', 'cyoa'], exe => ['orthofinder'] },
        'OrthoMCL_Pipeline' => {
            modules => ['orthomcl'], exe => ['orthomclPairs'] },
        'Phagepromoter' => {
            modules => 'phagepromoter' },
        'Phageterm' => {
            modules => ['phageterm', 'cyoa'], exe => 'PhageTerm.py' },
        'Phanotate' => {
            modules => ['trnascan', 'phanotate'], exe => 'phanotate.py' },
        'Phastaf' => {
            modules => ['phastaf', 'cyoa'] },
        'Prodigal' => {
            modules => ['prodigal', 'cyoa'] },
        'Prokka' => {
            modules => ['prokka', 'blast'], exe => 'prokka'},
        'Racer' => {
            modules => ['hitec'], exe => ['RACER'], },
        'Resfinder' => { modules => 'resfinder', exe => 'run_resfinder.py' },
        'Restriction_Catalog' => { modules => 'cyoa' },
        'Rgi' => {
            modules => ['kma', 'jellyfish', 'bowtie', 'bwa', 'diamond', 'rgi'], },
        'Rho_Predict' => { modules => 'rhotermpredict' },
        'RNAFold_Windows' => { modules => ['vienna', 'cyoa' ],
                               exe => 'RNAFold' },
        'Rosalind_Plus' => { modules => ['prodigal', 'cyoa'] },
        'RSEM' => {
            modules => 'rsem', exe => 'rsem-prepare-reference' },
        'RSEM_Index' => {
            modules => ['rsem', 'bowtie2'], exe => ['bowtie2'] },
        'Run_Essentiality' => {
            modules => 'essentiality' },
        'Run_Parse_Blast' => {
            modules => ['blast', 'cyoa'], exe => 'blastp' },
        'Salmon' => {
            modules => 'salmon', exe => 'salmon', },
        'Salmon_Index' => {
            modules => ['salmon'], exe => ['salmon']},
        'Samtools' => {
            modules => ['samtools', 'bamtools'], exe => 'samtools' },
        'Shovill' => {
            modules => 'shovill', exe => 'shovill' },
        'SLSearch' => {
            modules => 'cyoa' },
        'Snippy' => {
            modules => ['snippy', 'gubbins', 'fasttree'],
            exe => ['snippy', 'gubbins'], },
        'SNP_Ratio' => {
            modules => ['freebayes', 'libgsl', 'libhts', 'gatk',
                        'samtools', 'bcftools', 'vcftools', 'cyoa'],
            exe => ['samtools', 'cyoa', 'freebayes'] },
        'SNP_Ratio_Intron' => {
            modules => ['freebayes', 'libgsl', 'libhts', 'gatk',
                        'samtools', 'bcftools', 'vcftools', 'cyoa'], },
        'SNP_Ratio_Worker' => {
            modules => ['freebayes', 'libgsl', 'libhts', 'samtools',
                        'bcftools', 'vcftools', 'cyoa'], },
        'SNP_Ratio_Intron_Worker' => {
            modules => ['freebayes', 'libgsl', 'libhts', 'samtools',
                        'bcftools', 'vcftools', 'cyoa'], },
        'Sort_Indexes' => { modules => 'cyoa' },
        'Split_Align_Blast' => {
            modules => ['blast', 'blastdb', 'cyoa'], exe => 'blastn' },
        'STAR' => {
            modules => 'star' },
        'STAR_Index' => {
            modules => ['star'] },
        'TA_Check' => {
            modules => 'cyoa' },
        'Terminase_ORF_Reorder' => {
            modules => ['fasta', 'blast', 'blastdb', 'cyoa' ], },
        'Train_Prodigal' => { modules => 'prodigal' },
        'Tophat' => {
            modules => 'tophat' },
        'Transdecoder' => { modules => 'transdecoder' },
        'Transit_TPP' => {
            modules => ['bwa', 'transit', 'htseq' ], exe => 'tpp' },
        'Transposonpsi' => {
            modules => 'transposonpsi', exe => 'transposonPSI.pl' },
        'Trimomatic_Pairwise' => {  ## I can set the exe for single vs pairwise now
            modules => ['trimomatic'], },
        'Trimomatic_Single' => {  ## I can set the exe for single vs pairwise now
            modules => ['trimomatic'], },
        'Trinity' => {
            modules => 'trinity', exe => 'Trinity', },
        'Trinity_Post' => {
            modules => ['rsem', 'cyoa'], },
        'Trinotate' => {
            modules => ['divsufsort', 'transdecoder', 'blast', 'blastdb',
                        'signalp/4.1', 'hmmer', 'tmhmm', 'rnammer', 'trinotate', ],
            exe => ['autoTrinotate.pl', 'Trinotate'] },
        'tRNAScan' => {
            modules => ['infernal', 'trnascan'], exe => 'trnascan' },
        'Unicycler' => {
            modules => ['trimomatic', 'bowtie2', 'spades', 'unicycler',
                        'flash', 'shovill', 'bwa', 'pilon'],
            exe => 'unicycler',
        },
        'Unicycler_Filter_Depth' => {
            modules => 'cyoa', },
        'Velvet' => {
            modules => 'velvet', exe => 'velveth' },
        'Xref_Crispr' => {
            modules => 'cyoa', },
        );
    my @function_array = split(/::/, $subroutine);
    my $function = $function_array[$#function_array];
    my $datum = $module_data{$function};
    if ($datum) {
        ## Sometimes I don't bother to make the module list
        ## a list, but just a scalar.
        my $module_type = ref($datum->{modules});
        if (!$module_type) {
            my @mod_lst = ();
            push(@mod_lst, $datum->{modules});
            $datum->{modules} = \@mod_lst;
        }
    } else {
        print "The function: $function does not appear to have defined modules.\n";
    }
    return(%{$datum});
}

1;
