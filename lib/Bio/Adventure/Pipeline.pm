package Bio::Adventure::Pipeline;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use Cwd qw"abs_path getcwd cwd";
use File::Basename;
use File::Spec;
use File::Which qw"which";
use Time::HiRes qw"sleep";

=head1 NAME

Bio::Adventure::Pipeline - Example implementations of a few pipelines.

=head1 SYNOPSIS

use Bio::adventure;
my $cyoa = Bio::Adventure->new();
$cyoa->Phage_Assemble()

=head2 Introduction

This function provides a few sets of predefined orders for running a few analyses.
Currently the only ones likely to work are the Phage and Bacterial assemblies.  The
others were written for a much earlier version of this and have not been updated in
a very long time.

=cut
sub Annotate_Assembly {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input',],);
    my $prefix = sprintf("%02d", 10);
    my $final_locustag = basename(cwd());
    my $filtered_dir = dirname($options->{input});
    my $last_job = '';

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning phastaf.\n";
    my $phastaf = $class->Bio::Adventure::Phage::Phastaf(
        input => $options->{input},
        jname => 'phastaf',
        jprefix => $prefix,);
    $last_job = $phastaf->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nPerforming initial prokka annotation.\n";
    my $prokka = $class->Bio::Adventure::Annotation::Prokka(
        input => $options->{input},
        locus_tag => $final_locustag,
        jdepends => $last_job,
        jname => 'prokka',
        jprefix => $prefix,);
    $last_job = $prokka->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning prodigal to get RBSes.\n";
    my $prodigal = $class->Bio::Adventure::Feature_Prediction::Prodigal(
        input => $prokka->{output_assembly},
        jdepends => $last_job,
        jprefix => $prefix,);
    $last_job = $prodigal->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning Jellyfish on the assembly.\n";
    my $jelly = $class->Bio::Adventure::Count::Jellyfish(
        input => $prokka->{output_assembly},
        jdepends => $last_job,
        jname => 'jelly',
        jprefix => $prefix,);
    $last_job = $jelly->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning aragorn on the assembly to search for tmRNAs.\n";
    my $aragorn = $class->Bio::Adventure::Feature_Prediction::Aragorn(
        input => $prokka->{output_assembly},
        jdepends => $last_job,
        jname => 'aragorn',
        jprefix => $prefix,);
    $last_job = $aragorn->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nInvoking trinotate.\n";
    my $trinotate = $class->Bio::Adventure::Annotation::Trinotate(
        input => $prokka->{output_cds},
        jdepends => $last_job,
        jname => 'trinotate',
        jprefix => $prefix,);
    $last_job = $trinotate->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nSearching for resistance genes with abricate.\n";
    my $abricate = $class->Bio::Adventure::Resistance::Abricate(
        input => $prokka->{output_cds},
        jdepends => $last_job,
        jname => 'abricate',
        jprefix => $prefix,);
    $last_job = $abricate->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning interproscan.\n";
    my $interpro = $class->Bio::Adventure::Annotation::Interproscan(
        input => $prokka->{output_peptide},
        jdepends => $last_job,
        jname => 'interproscan',
        jprefix => $prefix,);
    $last_job = $interpro->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nMerging annotation files.\n";
    my $merge = $class->Bio::Adventure::Metadata::Merge_Annotations(
        evalue => undef,
        input_fsa => $prokka->{output_assembly},
        input_genbank => $prokka->{output_genbank},
        input_tsv => $prokka->{output_tsv},
        input_abricate => $abricate->{output},
        input_aragorn => $aragorn->{output},
        input_classifier => '',
        input_interpro => $interpro->{output_tsv},
        input_prodigal => $prodigal->{output},
        input_trinotate => $trinotate->{output},
        jdepends => $last_job,
        jname => 'mergeannot',
        jprefix => $prefix,);
    $last_job = $merge->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning cgview.\n";
    my $cgview = $class->Bio::Adventure::Visualization::CGView(
        input => $merge->{output_gbk},
        jdepends => $last_job,
        jprefix => $prefix);
    $last_job = $cgview->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning Vienna RNAfold on the assembly.\n";
    my $vienna = $class->Bio::Adventure::Structure::RNAFold_Windows(
        input => $prokka->{output_assembly},
        jdepends => $last_job,
        jname => 'vienna',
        jprefix => $prefix,);
    $last_job = $vienna->{job_id};
    sleep(0.2);

    my $ret = {
        phastaf => $phastaf,
        prokka => $prokka,
        prodigal => $prodigal,
        jellyfish => $jelly,
        aragorn => $aragorn,
        trinotate => $trinotate,
        abricate => $abricate,
        interproscan => $interpro,
        merge => $merge,
        cgview => $cgview,
        vienna => $vienna,
    };
    return($ret)
}

sub Annotate_Phage {
    my ($class, %args) = @_;
    print "This annotation pipeline starts with the output of unicycler.\n";
    my $options = $class->Get_Vars(
        args => \%args,);

    my $prefix = sprintf("%02d", 7);
    my $final_locustag = basename(cwd());
    my $assembly_output = qq"outputs/07unicycler/${final_locustag}_final_assembly.fasta";
    my $filtered_reads = qq"outputs/05filter_kraken_host/r1_host_filtered.fastq.xz:outputs/05filter_kraken_host/r2_host_filtered.fastq.xz";
    my $last_job;

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nDepth filtering initial assembly.\n";
    my $depth_filtered = $class->Bio::Adventure::Assembly::Filter_Depth(
        jdepends => $last_job,
        input => $assembly_output,
        jprefix => $prefix,
        jname => 'depth_filter',);
    $last_job = $depth_filtered->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning phastaf.\n";
    my $phastaf = $class->Bio::Adventure::Phage::Phastaf(
        jdepends => $last_job,
        input => $depth_filtered->{output},
        jprefix => $prefix,
        jname => 'phastaf',);
    $last_job = $phastaf->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning virus ICTV classifier.\n";
    my $ictv = $class->Bio::Adventure::Phage::Classify_Phage(
        jdepends => $last_job,
        input => $depth_filtered->{output},
        jprefix => $prefix,
        jname => 'ictv',);
    $last_job = $ictv->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nSetting the Rosalind strand to the one with the most ORFs.\n";
    my $rosalindplus = $class->Bio::Adventure::Annotation::Rosalind_Plus(
        jdepends => $last_job,
        input => $depth_filtered->{output},
        jprefix => $prefix,
        jname => 'rosalindplus',);
    $last_job = $rosalindplus->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nUsing Phageterm to reorient the assembly to the terminii.\n";
    my $phageterm = $class->Bio::Adventure::Phage::Phageterm(
        jdepends => $last_job,
        input => $filtered_reads,
        jprefix => $prefix,
        jname => 'phageterm',
        library => $rosalindplus->{output},);
    $last_job = $phageterm->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nPerforming reordering by putative terminase if phageterm failed.\n";
    my $termreorder = $class->Bio::Adventure::Phage::Terminase_ORF_Reorder(
        fasta_tool => 'fastx36',
        input => $phageterm->{output},
        library => 'terminase',
        species => 'phages',
        test_file => $phageterm->{test_file},
        jprefix => $prefix,
        jname => 'reorder',
        jdepends => $last_job,);
    $last_job = $termreorder->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nPerforming initial prokka annotation.\n";
    my $prokka = $class->Bio::Adventure::Annotation::Prokka(
        jdepends => $last_job,
        input => $termreorder->{output},
        jprefix => $prefix,
        jname => 'prokka',
        locus_tag => $final_locustag,);
    $last_job = $prokka->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning prodigal to get RBSes.\n";
    my $prodigal = $class->Bio::Adventure::Feature_Prediction::Prodigal(
        jdepends => $last_job,
        ## Use output_assembly to avoid having too much clutter on the sequence ID line.
        input => $prokka->{output_assembly},
        jprefix => $prefix,);
    $last_job = $prodigal->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning glimmer for less stringent ORFs.\n";
    my $glimmer = $class->Bio::Adventure::Feature_Prediction::Glimmer_Single(
        jdepends => $last_job,
        input => $prokka->{output_assembly},
        jprefix => $prefix,);
    $last_job = $glimmer->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning phanotate.\n";
    my $phanotate = $class->Bio::Adventure::Feature_Prediction::Phanotate(
        jdepends => $last_job,
        input => $prokka->{output_assembly},
        jprefix => $prefix,);
    $last_job = $phanotate->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nMerging the CDS predictions and writing an initial gbk file.\n";
    my $cds_merge = $class->Bio::Adventure::Annotation_Genbank::Merge_CDS_Predictions(
        jdepends => $last_job,
        input => $prokka->{output_genbank},
        input_glimmer => $glimmer->{output},
        input_phanotate => $phanotate->{output},
        input_prodigal => $prodigal->{output_gff},
        jprefix => $prefix);
    $last_job = $cds_merge->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning Jellyfish on the assembly.\n";
    my $jelly = $class->Bio::Adventure::Count::Jellyfish(
        jdepends => $last_job,
        input => $cds_merge->{output_fsa},
        jprefix => $prefix,
        jname => 'jelly',);
    $last_job = $jelly->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning aragorn on the assembly to search for tmRNAs.\n";
    my $aragorn = $class->Bio::Adventure::Feature_Prediction::Aragorn(
        jdepends => $last_job,
        input => $cds_merge->{output_fsa},
        jprefix => $prefix,
        jname => 'aragorn',);
    $last_job = $aragorn->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nInvoking trinotate.\n";
    my $trinotate = $class->Bio::Adventure::Annotation::Trinotate(
        jdepends => $last_job,
        input => $cds_merge->{output_cds},
        jprefix => $prefix,
        jname => 'trinotate',
        config => 'phage.txt',);
    $last_job = $trinotate->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nSearching for resistance genes with abricate.\n";
    my $abricate = $class->Bio::Adventure::Resistance::Abricate(
        jdepends => $last_job,
        jprefix => $prefix,
        jname => 'abricate',
        input => $cds_merge->{output_cds},);
    $last_job = $abricate->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning interproscan.\n";
    my $interpro = $class->Bio::Adventure::Annotation::Interproscan(
        jdepends => $last_job,
        jprefix => $prefix,
        jname => 'interproscan',
        input => $cds_merge->{output_faa},);
    $last_job = $interpro->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nMerging annotation files.\n";
    my $merge = $class->Bio::Adventure::Metadata::Merge_Annotations(
        jdepends => $last_job,
        input_fsa => $cds_merge->{output_fsa},
        input_genbank => $cds_merge->{output_gbk},
        input_tsv => $cds_merge->{output_tsv},
        input_abricate => $abricate->{output},
        input_aragorn => $aragorn->{output},
        input_classifier => $ictv->{output},
        input_interpro => $interpro->{output_tsv},
        input_phageterm => $phageterm->{output_dtr},
        input_prodigal => $prodigal->{output},
        input_trinotate => $trinotate->{output},
        jprefix => $prefix,
        jname => 'mergeannot',);
    $last_job = $merge->{job_id};
    sleep(0.2);

    ## An extra invocation of merge_annotations which will not modify the final gbk file.
    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nMerging annotation files a second time.\n";
    my $merge2 = $class->Bio::Adventure::Metadata::Merge_Annotations(
        jdepends => $last_job,
        input_fsa => $cds_merge->{output_fsa},
        input_genbank => $cds_merge->{output_gbk},
        input_tsv => $cds_merge->{output_tsv},
        input_abricate => $abricate->{output},
        input_aragorn => $aragorn->{output},
        input_classifier => $ictv->{output},
        input_interpro => $interpro->{output_tsv},
        input_phageterm => $phageterm->{output_dtr},
        input_prodigal => $prodigal->{output},
        input_trinotate => $trinotate->{output},
        jprefix => $prefix,
        jname => 'mergeannot2',
        keep_genes => 0,
        locus_tag => 0,
        evalue => undef,);
    sleep(0.2);
    $last_job = $merge2->{job_id};

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning cgview.\n";
    my $cgview = $class->Bio::Adventure::Visualization::CGView(
        jdepends => $last_job,
        jprefix => $prefix,
        input => $merge->{output_gbk},);
    $last_job = $cgview->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning Vienna RNAfold on the assembly.\n";
    my $vienna = $class->Bio::Adventure::Structure::RNAFold_Windows(
        jdepends => $last_job,
        input => $merge2->{output_gbk},
        jprefix => $prefix,
        jname => 'vienna',);
    $last_job = $cds_merge->{job_id};
    sleep(0.2);

}

sub Riboseq {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species']);
    my $prefix = sprintf("%02d", 0);
    my $last_job;

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning FastQC.\n";
    my $fastqc = $class->Bio::Adventure::QA::Fastqc(
        input => $options->{input},
        jprefix => $prefix,
        jname => 'fastqc',);
    $last_job = $fastqc->{job_id};
    sleep(0.2);

    print "\nRunning cutadapt.\n";
    my $cutadapt = $class->Bio::Adventure::Trim::Cutadap(
        input => $options->{input},
        task => 'riboseq',
        jdepends => $last_job,
        jprefix => $prefix,
        jname => 'cutadapt',);
    $last_job = $cutadapt->{job_id};
    sleep(0.2);

    print "\nRunning Biopieces.\n";
    my $pieces = $class->Bio::Adventure::QA::Biopieces_Graph(
        input => $cutadapt->{output},
        jdepends => $last_job,
        jprefix => $prefix,
        jname => 'biopieces',);
    $last_job = $pieces->{job_id};
    sleep(0.2);

    print "\nFiltering rRNA reads.\n";
    my $rrna_filter = $class->Bio::Adventure::Map::Bowtie_RRNA(
        input => $cutadapt->{output},
        species => $options->{species},
        jdepends => $last_job,
        jprefix => $prefix,
        jname => 'rrna_filter',);
    $last_job = $rrna_filter->{job_id};
    sleep(0.2);

    print "\nMapping remaining reads.\n";
    my $map_reads = $class->Bio::Adventure::Map::Bowtie(
        input => $rrna_filter->{output},
        species => $options->{species},
        jdepends => $last_job,
        jprefix => $prefix,
        jname => 'map_reads',);
    $last_job = $map_reads->{job_id};
    sleep(0.2);

    my $ret = {
        fastqc => $fastqc,
        cutadapt => $cutadapt,
        biopieces => $pieces,
        rrna => $rrna_filter,
        bt => $map_reads};
    return($ret);
}

sub Bowtie {
    my ($class, %args) = @_;
    $args{aligner} = 'bowtie';
    my $rnaseq_jobs = $class->RNAseq(%args);
    return($rnaseq_jobs);
}

sub Bowtie2 {
    my ($class, %args) = @_;
    $args{aligner} = 'bowtie2';
    my $rnaseq_jobs = $class->RNAseq(%args);
    return($rnaseq_jobs);
}

sub BWA {
    my ($class, %args) = @_;
    $args{aligner} = 'bwa';
    my $rnaseq_jobs = $class->RNAseq(%args);
    return($rnaseq_jobs);
}

sub Hisat {
    my ($class, %args) = @_;
    $args{aligner} = 'hisat';
    my $rnaseq_jobs = $class->RNAseq(%args);
    return($rnaseq_jobs);
}

sub Kallisto {
    my ($class, %args) = @_;
    $args{aligner} = 'kallisto';
    my $rnaseq_jobs = $class->RNAseq(%args);
    return($rnaseq_jobs);
}

sub RNAseq {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        htseq_type => 'gene',
        htseq_id => 'ID',
        mapper => 'hisat2',);
    my $prefix = sprintf("%02d", 1);
    my $final_locustag = basename(cwd());

    print "Starting fastqc.\n";
    my $fastqc = $class->Bio::Adventure::QA::Fastqc(
        jprefix => $prefix,);
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nStarting trimmer.\n";
    my $trim = $class->Bio::Adventure::Trim::Trimomatic(
        input => $fastqc->{input},
        jprefix => $prefix,
        jname => 'trimomatic',);
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nStarting mapper.\n";
    my $mapper;
    if ($options->{mapper} eq 'bowtie1') {
        $mapper = $class->Bio::Adventure::Map::Bowtie(
            jdepends => $trim->{job_id},
            input => $trim->{output},
            jprefix => $prefix,
            species => $options->{species},
            htseq_type => $options->{htseq_type},
            htseq_id => $options->{htseq_id},);
    } elsif ($options->{mapper} eq 'bowtie2') {
        $mapper = $class->Bio::Adventure::Map::Bowtie2(
            jdepends => $trim->{job_id},
            input => $trim->{output},
            jprefix => $prefix,
            species => $options->{species},
            htseq_type => $options->{htseq_type},
            htseq_id => $options->{htseq_id},);
    } elsif ($options->{mapper} eq 'bwa') {
        $mapper = $class->Bio::Adventure::Map::BWA(
            jdepends => $trim->{job_id},
            input => $trim->{output},
            jprefix => $prefix,
            species => $options->{species},
            htseq_type => $options->{htseq_type},
            htseq_id => $options->{htseq_id},);
    } elsif ($options->{mapper} eq 'kallisto') {
        $mapper = $class->Bio::Adventure::Map::Kallisto(
            jdepends => $trim->{job_id},
            input => $trim->{output},
            jprefix => $prefix,
            species => $options->{species},);
    } elsif ($options->{mapper} eq 'salmon') {
        $mapper = $class->Bio::Adventure::Map::Salmon(
            jdepends => $trim->{job_id},
            input => $trim->{output},
            jprefix => $prefix,
            species => $options->{species},);
    } else {
        $mapper = $class->Bio::Adventure::Map::Hisat(
            jdepends => $trim->{job_id},
            input => $trim->{output},
            jprefix => $prefix,
            species => $options->{species},
            htseq_type => $options->{htseq_type},
            htseq_id => $options->{htseq_id},);
    }
    $mapper->{trim} = $trim;
    return($mapper);
}

sub TNseq {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species']);
    my $prefix = sprintf("%02d", 0);
    my $last_job;

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning FastQC.\n";
    my $fastqc = $class->Bio::Adventure::QA::Fastqc(
        input => $options->{input},
        jprefix => $prefix,
        jname => 'fastqc',);
    $last_job = $fastqc->{job_id};
    sleep(0.2);

    print "\nRunning Cutadapt with TNSeq adapters.\n";
    my $cutadapt = $class->Bio::Adventure::Trim::Cutadapt(
        input => $options->{input},
        task => 'tnseq',
        jprefix => $prefix,
        jname => 'cutadapt',);
    $last_job = $cutadapt->{job_id};
    sleep(0.2);

    print "\nRunning Biopieces.\n";
    my $pieces = $class->Bio::Adventure::QA::Biopieces_Graph(
        input => $cutadapt->{output},
        jprefix => $prefix,
        jname => 'biopieces',);
    $last_job = $pieces->{job_id};
    sleep(0.2);

    print "\nMapping remaining reads.\n";
    my $map_reads = $class->Bio::Adventure::Map::Bowtie(
        input => $cutadapt->{output},
        species => $options->{species},
        jdepends => $last_job,
        jprefix => $prefix,
        jname => 'map_reads',);
    $last_job = $map_reads->{job_id};
    sleep(0.2);

    my $ret = {
        fastqc => $fastqc,
        cutadapt => $cutadapt,
        biopieces => $pieces,
        bt => $map_reads, };
    return($ret);
}

=head2 C<Assemble>

This should probably be renamed to 'Bacterial_Assemble()' or something.  It is pretty
bacteria-specific.  With that in mind, it should do the following:

1.  Take one or two fastq files as input.
2.  Trim them with trimomatic.
3.  Run fastqc on them.
4.  Use hitec's RACER command for error correction.
5.  Classify the reads with Kraken2.
6.  Create an initial assembly with Unicycler.
7.  Search for (pro)phage regions with phastaf.
8.  Generate an initial prokka annotation.

The final outputs are ready for the Annotate_Assembly() function,
but I chose to keep that separate so that the user may check it out
before committing to a potentially long series of blast searches.

=cut
sub Assemble {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],);
    my $final_locustag = basename(cwd());

    my $prefix = sprintf("%02d", 1);
    print "\nStarting trimmer.\n";
    my $trim = $class->Bio::Adventure::Trim::Trimomatic(
        input => $options->{input},
        jprefix => $prefix,
        jname => 'trimomatic',);
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nStarting fastqc.\n";
    my $fastqc = $class->Bio::Adventure::QA::Fastqc(
        input => $trim->{input},
        jmem => 3,
        jnice => 100,
        jprefix => $prefix,);
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nStarting read correction.\n";
    my $correct = $class->Bio::Adventure::Trim::Racer(
        jdepends => $trim->{job_id},
        input => $trim->{output},
        jprefix => $prefix,
        jname => 'racer',);
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nClassifying sequences with Kraken2.\n";
    my $kraken = $class->Bio::Adventure::Annotation::Kraken(
        jdepends => $correct->{job_id},
        input => $correct->{output},
        jprefix => $prefix,
        jname => 'kraken',
        library => 'viral',);
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nCreating initial assembly with Unicycler.\n";
    my $assemble = $class->Bio::Adventure::Assembly::Unicycler(
        jdepends => $correct->{job_id},
        input => $correct->{output},
        jprefix => $prefix,
        jname => 'unicycler',);
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning phastaf.\n";
    my $phastaf = $class->Bio::Adventure::Phage::Phastaf(
        jdepends => $assemble->{job_id},
        input => $assemble->{output},
        jprefix => $prefix,
        jname => 'phastaf',);
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nPerforming initial prokka annotation.\n";
    my $prokka = $class->Bio::Adventure::Annotation::Prokka(
        jdepends => $assemble->{job_id},
        input => $assemble->{output},
        jprefix => $prefix,
        jname => 'prokka',
        locus_tag => $final_locustag,);
    sleep(0.2);
}

=head2 C<Phage_Assemble>

This function is an extension of the bacterial assembly methods above, but since phage
genomes are nice and small, it throws a few more tools into the mix.  Therefore, it does
the following:

1.  Take one or two fastq files as input.
2.  Trim them with trimomatic.
3.  Run fastqc on them.
4.  Use hitec's RACER command for error correction.
5.  Classify the reads with Kraken2's standard(and/or bacterial) database
6.  Reads the Kraken output to find an appropriate bacterial host species.
7.  Uses hisat2 to filter out the host reads using the species from #6.
8.  Classify the reads with Kraken2 with an extended viral database.
9.  Create an unicycler assembly.
10. Filter it by depth coverage with the assumption that sometimes prophage sequences sneak past #7.
11. Run phastaf as an initial taxonomic classification (also to see if it is a nonsense chimera).
12. Uses tblastx against the curated set of ICTV viral references to generate another taxonomic classification.
13. Counts up the +/- strand ORFs to make the Rosalind strand have the most.
14. Runs phageterm to reorganize the genome if there is a detectable DTR.
15. If #14 is a no, search for a terminase and reorganize the genome accordingly.
16. Create an initial prokka annotation.
17. Run prodigal separately to search for Shine-Dalgarnos.
18. Run glimmer3 separately for a less stringent set of putative CDS sequences.
19. Rewrite the prokka files with information from #17-18.
20. Perform a kmer count of the assembly with jellyfish.
21. Run Trinotate on the assembly (supplemented with a hand-curated viral blast database).
22. Run interproscan on the assembly.
23. Run abricate on the assembly.
24. Merge the results from #16-23 into a new, more verbose genbank file with confidence cues.
25. Run cgview to make a fun picture of the assembly (I want to replace this with circos).
26. Redo #24 but without the notes regarding the annotation confidence.

=cut
sub Phage_Assemble {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        host_species => '');
    my $prefix = sprintf("%02d", 0);
    my $final_locustag = basename(cwd());

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nStarting trimmer.\n";
    my $trim = $class->Bio::Adventure::Trim::Trimomatic(
        jmem => 4,
        jname => 'trimomatic',
        jprefix => $prefix,);
    my $last_job = $trim->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nStarting fastqc.\n";
    my $fastqc = $class->Bio::Adventure::QA::Fastqc(
        input => $trim->{input},
        jmem => 3,
        jnice => 100,
        jprefix => $prefix,);
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nStarting read correction.\n";
    my $correct = $class->Bio::Adventure::Trim::Racer(
        jdepends => $last_job,
        input => $trim->{output},
        jmem => 3,
        jprefix => $prefix,);
    $last_job = $correct->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nClassifying sequences with Kraken using the standard database.\n";
    my $kraken_std = $class->Bio::Adventure::Annotation::Kraken(
        jdepends => $last_job,
        input => $correct->{output},
        jprefix => $prefix,
        jname => 'krakenstd',
        library => 'standard',);
    $last_job = $kraken_std->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nUsing Kraken results to filter likely host sequences.\n";
    my $filter = $class->Bio::Adventure::Phage::Filter_Host_Kraken(
        jdepends => $last_job,
        input => $kraken_std->{output},
        input_fastq => $correct->{output},
        jprefix => $prefix,
        jname => 'krakenfilter',);
    $last_job = $filter->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nClassifying sequences with Kraken using the viral database.\n";
    my $kraken = $class->Bio::Adventure::Annotation::Kraken(
        jdepends => $last_job,
        input => $filter->{output},
        jprefix => $prefix,
        jname => 'kraken',
        library => 'viral',);
    $last_job = $kraken->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nCreating initial assembly with Unicycler.\n";
    my $assemble = $class->Bio::Adventure::Assembly::Unicycler(
        jdepends => $last_job,
        input => $filter->{output},
        jprefix => $prefix,
        jname => 'unicycler',);
    $last_job = $assemble->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nDepth filtering initial assembly.\n";
    my $depth_filtered = $class->Bio::Adventure::Assembly::Filter_Depth(
        jdepends => $last_job,
        input => $assemble->{output},
        jprefix => $prefix,
        jname => 'depth_filter',);
    $last_job = $depth_filtered->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning phastaf.\n";
    my $phastaf = $class->Bio::Adventure::Phage::Phastaf(
        jdepends => $last_job,
        input => $depth_filtered->{output},
        jprefix => $prefix,
        jname => 'phastaf',);
    $last_job = $phastaf->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning virus ICTV classifier.\n";
    my $ictv = $class->Bio::Adventure::Phage::Classify_Phage(
        jdepends => $last_job,
        input => $depth_filtered->{output},
        jprefix => $prefix,
        jname => 'ictv',);
    $last_job = $ictv->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nSetting the Rosalind strand to the one with the most ORFs.\n";
    my $rosalindplus = $class->Bio::Adventure::Annotation::Rosalind_Plus(
        jdepends => $last_job,
        input => $depth_filtered->{output},
        jprefix => $prefix,
        jname => 'rosalindplus',);
    $last_job = $rosalindplus->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nUsing Phageterm to reorient the assembly to the terminii.\n";
    my $phageterm = $class->Bio::Adventure::Phage::Phageterm(
        jdepends => $last_job,
        input => $filter->{output},
        jprefix => $prefix,
        jname => 'phageterm',
        library => $rosalindplus->{output},);
    $last_job = $phageterm->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nChecking coverage.\n";
    my $coverage = $class->Bio::Adventure::Assembly::Assembly_Coverage(
        jdepends => $last_job,
        input => $filter->{output},
        library => $assemble->{output},
        jprefix => $prefix,
        jname => 'coverage',);
    $last_job = $coverage->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nPerforming reordering by putative terminase if phageterm failed.\n";
    my $termreorder = $class->Bio::Adventure::Phage::Terminase_ORF_Reorder(
        jdepends => $last_job,
        fasta_tool => 'fastx36',
        input => $phageterm->{output},
        jprefix => $prefix,
        jname => 'reorder',
        library => 'terminase',
        species => 'phages',
        test_file => $phageterm->{test_file},);
    $last_job = $termreorder->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nPerforming initial prokka annotation.\n";
    my $prokka = $class->Bio::Adventure::Annotation::Prokka(
        jdepends => $last_job,
        input => $termreorder->{output},
        jprefix => $prefix,
        jname => 'prokka',
        locus_tag => $final_locustag,);
    $last_job = $prokka->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning prodigal to get RBSes.\n";
    my $prodigal = $class->Bio::Adventure::Feature_Prediction::Prodigal(
        jdepends => $last_job,
        ## Use output_assembly to avoid having too much clutter on the sequence ID line.
        input => $prokka->{output_assembly},
        jprefix => $prefix,);
    $last_job = $prodigal->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning glimmer for less stringent ORFs.\n";
    my $glimmer = $class->Bio::Adventure::Feature_Prediction::Glimmer_Single(
        jdepends => $last_job,
        input => $prokka->{output_assembly},
        jprefix => $prefix,);
    $last_job = $glimmer->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning phanotate.\n";
    my $phanotate = $class->Bio::Adventure::Feature_Prediction::Phanotate(
        jdepends => $last_job,
        input => $prokka->{output_assembly},
        jprefix => $prefix,);
    $last_job = $phanotate->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nMerging the CDS predictions and writing an initial gbk file.\n";
    my $cds_merge = $class->Bio::Adventure::Annotation_Genbank::Merge_CDS_Predictions(
        jdepends => $last_job,
        input => $prokka->{output_genbank},
        input_glimmer => $glimmer->{output},
        input_phanotate => $phanotate->{output},
        input_prodigal => $prodigal->{output_gff},
        jprefix => $prefix);
    $last_job = $cds_merge->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning Jellyfish on the assembly.\n";
    my $jelly = $class->Bio::Adventure::Count::Jellyfish(
        jdepends => $last_job,
        input => $cds_merge->{output_fsa},
        jprefix => $prefix,
        jname => 'jelly',);
    $last_job = $jelly->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning aragorn on the assembly to search for tmRNAs.\n";
    my $aragorn = $class->Bio::Adventure::Feature_Prediction::Aragorn(
        jdepends => $last_job,
        input => $cds_merge->{output_fsa},
        jprefix => $prefix,
        jname => 'aragorn',);
    $last_job = $aragorn->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nInvoking trinotate.\n";
    my $trinotate = $class->Bio::Adventure::Annotation::Trinotate(
        jdepends => $last_job,
        input => $cds_merge->{output_cds},
        jprefix => $prefix,
        jname => 'trinotate',
        config => 'phage.txt',);
    $last_job = $trinotate->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nSearching for resistance genes with abricate.\n";
    my $abricate = $class->Bio::Adventure::Resistance::Abricate(
        jdepends => $last_job,
        jprefix => $prefix,
        jname => 'abricate',
        input => $cds_merge->{output_cds},);
    $last_job = $abricate->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning interproscan.\n";
    my $interpro = $class->Bio::Adventure::Annotation::Interproscan(
        jdepends => $last_job,
        jprefix => $prefix,
        jname => 'interproscan',
        input => $cds_merge->{output_faa},);
    $last_job = $interpro->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nMerging annotation files.\n";
    my $merge = $class->Bio::Adventure::Metadata::Merge_Annotations(
        jdepends => $last_job,
        input_fsa => $cds_merge->{output_fsa},
        input_genbank => $cds_merge->{output_gbk},
        input_tsv => $cds_merge->{output_tsv},
        input_abricate => $abricate->{output},
        input_aragorn => $aragorn->{output},
        input_classifier => $ictv->{output},
        input_interpro => $interpro->{output_tsv},
        input_phageterm => $phageterm->{output_dtr},
        input_prodigal => $prodigal->{output},
        input_trinotate => $trinotate->{output},
        jprefix => $prefix,
        jname => 'mergeannot',);
    $last_job = $merge->{job_id};
    sleep(0.2);

    ## An extra invocation of merge_annotations which will not modify the final gbk file.
    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nMerging annotation files a second time.\n";
    my $merge2 = $class->Bio::Adventure::Metadata::Merge_Annotations(
        jdepends => $last_job,
        input_fsa => $cds_merge->{output_fsa},
        input_genbank => $cds_merge->{output_gbk},
        input_tsv => $cds_merge->{output_tsv},
        input_abricate => $abricate->{output},
        input_aragorn => $aragorn->{output},
        input_classifier => $ictv->{output},
        input_interpro => $interpro->{output_tsv},
        input_phageterm => $phageterm->{output_dtr},
        input_prodigal => $prodigal->{output},
        input_trinotate => $trinotate->{output},
        jprefix => $prefix,
        jname => 'mergeannot2',
        suffix => 'stripped',
        keep_genes => 0,
        locus_tag => 0,
        evalue => 0,);
    sleep(0.2);
    $last_job = $merge2->{job_id};

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning cgview.\n";
    my $cgview = $class->Bio::Adventure::Visualization::CGView(
        jdepends => $last_job,
        jprefix => $prefix,
        input => $merge->{output_gbk},);
    $last_job = $cgview->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning Vienna RNAfold on the assembly.\n";
    my $vienna = $class->Bio::Adventure::Structure::RNAFold_Windows(
        jdepends => $last_job,
        input => $cds_merge->{output},
        jprefix => $prefix,
        jname => 'vienna',);
    $last_job = $vienna->{job_id};
    sleep(0.2);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nCollecting output files.\n";
    my $collect = $class->Bio::Adventure::Metadata::Collect_Assembly(
        jdepends => $last_job,
        input_fsa => $merge->{output_fsa},
        input => $merge->{output_gbk},
        input_stripped => $merge2->{output_gbk},
        input_tsv => $merge->{output_tsv},
        input_xlsx => $merge->{output_xlsx},
        input_faa => $cds_merge->{output_faa},
        input_cds => $cds_merge->{output_cds},
        jprefix => $prefix,
        jname => 'collect',);
    $last_job = $collect->{job_id};
    sleep(0.2);

    my $ret = {
        trim => $trim,
        fastqc => $fastqc,
        correction => $correct,
        kraken_standard => $kraken_std,
        host_filter => $filter,
        kraken_viral => $kraken,
        assembly => $assemble,
        depth_filter => $depth_filtered,
        phastaf => $phastaf,
        ictv => $ictv,
        rosalindplus => $rosalindplus,
        phageterm => $phageterm,
        terminase_reorder => $termreorder,
        prokka => $prokka,
        prodigal => $prodigal,
        glimmer => $glimmer,
        cds_merge => $cds_merge,
        aragorn => $aragorn,
        jellyfish => $jelly,
        trinotate => $trinotate,
        abricate => $abricate,
        interproscan => $interpro,
        cgview => $cgview,
        merge_qualities => $merge,
        merge_unmodified => $merge2,
    };
    return($ret)
}

1;
