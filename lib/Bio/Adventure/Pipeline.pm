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
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nPerforming initial prokka annotation.\n";
    my $prokka = $class->Bio::Adventure::Annotation::Prokka(
        input => $options->{input},
        locus_tag => $final_locustag,
        jdepends => $last_job,
        jname => 'prokka',
        jprefix => $prefix,);
    $last_job = $prokka->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning prodigal to get RBSes.\n";
    my $prodigal = $class->Bio::Adventure::Feature_Prediction::Prodigal(
        input => $prokka->{output_assembly},
        jdepends => $last_job,
        jprefix => $prefix,);
    $last_job = $prodigal->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning Jellyfish on the assembly.\n";
    my $jelly = $class->Bio::Adventure::Count::Jellyfish(
        input => $prokka->{output_assembly},
        jdepends => $last_job,
        jname => 'jelly',
        jprefix => $prefix,);
    $last_job = $jelly->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning aragorn on the assembly to search for tmRNAs.\n";
    my $aragorn = $class->Bio::Adventure::Feature_Prediction::Aragorn(
        input => $prokka->{output_assembly},
        jdepends => $last_job,
        jname => 'aragorn',
        jprefix => $prefix,);
    $last_job = $aragorn->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning relaxed tRNAscan to search for tmRNAs.\n";
    my $trnascan = $class->Bio::Adventure::Feature_Prediction::tRNAScan(
        input => $prokka->{output_assembly},
        jdepends => $last_job,
        jname => 'trnascan',
        arbitrary => ' -r ',
        suffix => 'relaxed',
        jprefix => $prefix,);
    $last_job = $trnascan->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nInvoking trinotate.\n";
    my $trinotate = $class->Bio::Adventure::Annotation::Trinotate(
        input => $prokka->{output_cds},
        jdepends => $last_job,
        jname => 'trinotate',
        jprefix => $prefix,);
    $last_job = $trinotate->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nSearching for resistance genes with abricate.\n";
    my $abricate = $class->Bio::Adventure::Resistance::Abricate(
        input => $prokka->{output_cds},
        jdepends => $last_job,
        jname => 'abricate',
        jprefix => $prefix,);
    $last_job = $abricate->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning interproscan.\n";
    my $interpro = $class->Bio::Adventure::Annotation::Interproscan(
        input => $prokka->{output_peptide},
        jdepends => $last_job,
        jname => 'interproscan',
        jprefix => $prefix,);
    $last_job = $interpro->{job_id};
    sleep($options->{jsleep});

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
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning cgview.\n";
    my $cgview = $class->Bio::Adventure::Visualization::CGView(
        input => $merge->{output_gbk},
        jdepends => $last_job,
        jprefix => $prefix);
    $last_job = $cgview->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning Vienna RNAfold on the assembly.\n";
    my $vienna = $class->Bio::Adventure::Structure::RNAFold_Windows(
        input => $prokka->{output_assembly},
        jdepends => $last_job,
        jname => 'vienna',
        jprefix => $prefix,);
    $last_job = $vienna->{job_id};
    sleep($options->{jsleep});

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
    my $depth_filtered = $class->Bio::Adventure::Assembly::Unicycler_Filter_Depth(
        jdepends => $last_job,
        input => $assembly_output,
        jprefix => $prefix,
        jname => 'depth_filter',);
    $last_job = $depth_filtered->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning phastaf.\n";
    my $phastaf = $class->Bio::Adventure::Phage::Phastaf(
        jdepends => $last_job,
        input => $depth_filtered->{output},
        jprefix => $prefix,
        jname => 'phastaf',);
    $last_job = $phastaf->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning virus ICTV classifier.\n";
    my $ictv = $class->Bio::Adventure::Phage::Classify_Phage(
        jdepends => $last_job,
        input => $depth_filtered->{output},
        jprefix => $prefix,
        jname => 'ictv',);
    $last_job = $ictv->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nSetting the Rosalind strand to the one with the most ORFs.\n";
    my $rosalindplus = $class->Bio::Adventure::Annotation::Rosalind_Plus(
        jdepends => $last_job,
        input => $depth_filtered->{output},
        jprefix => $prefix,
        jname => 'rosalindplus',);
    $last_job = $rosalindplus->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nUsing Phageterm to reorient the assembly to the terminii.\n";
    my $phageterm = $class->Bio::Adventure::Phage::Phageterm(
        jdepends => $last_job,
        input => $filtered_reads,
        jprefix => $prefix,
        jname => 'phageterm',
        library => $rosalindplus->{output},);
    $last_job = $phageterm->{job_id};
    sleep($options->{jsleep});

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
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nPerforming initial prokka annotation.\n";
    my $prokka = $class->Bio::Adventure::Annotation::Prokka(
        jdepends => $last_job,
        input => $termreorder->{output},
        jprefix => $prefix,
        jname => 'prokka',
        locus_tag => $final_locustag,);
    $last_job = $prokka->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning prodigal to get RBSes.\n";
    my $prodigal = $class->Bio::Adventure::Feature_Prediction::Prodigal(
        jdepends => $last_job,
        ## Use output_assembly to avoid having too much clutter on the sequence ID line.
        input => $prokka->{output_assembly},
        jprefix => $prefix,);
    $last_job = $prodigal->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning glimmer for less stringent ORFs.\n";
    my $glimmer = $class->Bio::Adventure::Feature_Prediction::Glimmer_Single(
        jdepends => $last_job,
        input => $prokka->{output_assembly},
        jprefix => $prefix,);
    $last_job = $glimmer->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning phanotate.\n";
    my $phanotate = $class->Bio::Adventure::Feature_Prediction::Phanotate(
        jdepends => $last_job,
        input => $prokka->{output_assembly},
        jprefix => $prefix,);
    $last_job = $phanotate->{job_id};
    sleep($options->{jsleep});

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
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning Jellyfish on the assembly.\n";
    my $jelly = $class->Bio::Adventure::Count::Jellyfish(
        jdepends => $last_job,
        input => $cds_merge->{output_fsa},
        jprefix => $prefix,
        jname => 'jelly',);
    $last_job = $jelly->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning aragorn on the assembly to search for tmRNAs.\n";
    my $aragorn = $class->Bio::Adventure::Feature_Prediction::Aragorn(
        jdepends => $last_job,
        input => $cds_merge->{output_fsa},
        jprefix => $prefix,
        jname => 'aragorn',);
    $last_job = $aragorn->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nInvoking trinotate.\n";
    my $trinotate = $class->Bio::Adventure::Annotation::Trinotate(
        jdepends => $last_job,
        input => $cds_merge->{output_cds},
        jprefix => $prefix,
        jname => 'trinotate',
        config => 'phage.txt',);
    $last_job = $trinotate->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nSearching for resistance genes with abricate.\n";
    my $abricate = $class->Bio::Adventure::Resistance::Abricate(
        jdepends => $last_job,
        jprefix => $prefix,
        jname => 'abricate',
        input => $cds_merge->{output_cds},);
    $last_job = $abricate->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning interproscan.\n";
    my $interpro = $class->Bio::Adventure::Annotation::Interproscan(
        jdepends => $last_job,
        jprefix => $prefix,
        jname => 'interproscan',
        input => $cds_merge->{output_faa},);
    $last_job = $interpro->{job_id};
    sleep($options->{jsleep});

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
    sleep($options->{jsleep});

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
    sleep($options->{jsleep});
    $last_job = $merge2->{job_id};

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning cgview.\n";
    my $cgview = $class->Bio::Adventure::Visualization::CGView(
        jdepends => $last_job,
        jprefix => $prefix,
        input => $merge->{output_gbk},);
    $last_job = $cgview->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning Vienna RNAfold on the assembly.\n";
    my $vienna = $class->Bio::Adventure::Structure::RNAFold_Windows(
        jdepends => $last_job,
        input => $merge2->{output_gbk},
        jprefix => $prefix,
        jname => 'vienna',);
    $last_job = $cds_merge->{job_id};
    sleep($options->{jsleep});

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
    sleep($options->{jsleep});

    print "\nRunning cutadapt.\n";
    my $cutadapt = $class->Bio::Adventure::Trim::Cutadap(
        input => $options->{input},
        task => 'riboseq',
        jdepends => $last_job,
        jprefix => $prefix,
        jname => 'cutadapt',);
    $last_job = $cutadapt->{job_id};
    sleep($options->{jsleep});

    print "\nRunning Biopieces.\n";
    my $pieces = $class->Bio::Adventure::QA::Biopieces_Graph(
        input => $cutadapt->{output},
        jdepends => $last_job,
        jprefix => $prefix,
        jname => 'biopieces',);
    $last_job = $pieces->{job_id};
    sleep($options->{jsleep});

    print "\nFiltering rRNA reads.\n";
    my $rrna_filter = $class->Bio::Adventure::Map::Bowtie_RRNA(
        input => $cutadapt->{output},
        species => $options->{species},
        jdepends => $last_job,
        jprefix => $prefix,
        jname => 'rrna_filter',);
    $last_job = $rrna_filter->{job_id};
    sleep($options->{jsleep});

    print "\nMapping remaining reads.\n";
    my $map_reads = $class->Bio::Adventure::Map::Bowtie(
        input => $rrna_filter->{output},
        species => $options->{species},
        jdepends => $last_job,
        jprefix => $prefix,
        jname => 'map_reads',);
    $last_job = $map_reads->{job_id};
    sleep($options->{jsleep});

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

sub Process_RNAseq {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        host_filter => 0,
        gff_type => 'gene',
        gff_tag => 'ID',
        intron => 0,
        mapper => 'hisat2',);
    my $prefix = sprintf("%02d", 0);
    my $cwd_name = basename(cwd());
    my @jobs = ();

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Starting trimmer.\n";
    my $trim = $class->Bio::Adventure::Trim::Trimomatic(
        input => $options->{input},
        jprefix => $prefix,
        jname => 'trimomatic',);
    push(@jobs, $trim);
    my $last_job = $trim->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Starting fastqc.\n";
    my $fastqc = $class->Bio::Adventure::QA::Fastqc(
        input => $trim->{input},
        jnice => 100,
        jprefix => $prefix,);
    push(@jobs, $fastqc);
    $last_job = $fastqc->{job_id};
    sleep($options->{jsleep});

    ## Have some logic to handle a first species which will be used to filter
    ## any following species provided.
    my @species_list = split(/:/, $options->{species});
    my @type_list = split(/:/, $options->{gff_type});
    my @id_list = split(/:/, $options->{gff_tag});

    my $first_species = shift @species_list;
    my $first_type = shift @type_list;
    my $first_id = shift @id_list;

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Performing initial mapping against ${first_species}.\n";
    my $first_map = $class->Bio::Adventure::Map::Hisat2(
        jdepends => $last_job,
        input => $trim->{output},
        species => $first_species,
        gff_type => $first_type,
        gff_tag => $first_id,
        jprefix => $prefix,);
    $last_job = $first_map->{job_id};
    push(@jobs, $first_map);
    sleep($options->{jsleep});
    my $last_sam_job = $first_map->{samtools}->{job_id};

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Performing freebayes search against ${first_species}.\n";
    my $first_snp = $class->Bio::Adventure::SNP::Freebayes_SNP_Search(
        jdepends => $last_sam_job,
        input => $first_map->{samtools}->{paired_output},
        species => $first_species,
        gff_type => $first_type,
        gff_tag => $first_id,
        intron => $options->{intron},
        jprefix => $prefix,);
    push(@jobs, $first_snp);
    sleep($options->{jsleep});

    if (scalar(@species_list) > 0) {
        my $c = 0;
        for my $sp (@species_list) {
            my $nth_species = $sp;
            my $nth_type = $type_list[$c];
            $nth_type = $first_type unless (defined($nth_type));
            my $nth_id = $id_list[$c];
            $nth_id = $first_id unless (defined($nth_id));

            ## Handle if we want to host-filter the data
            my $nth_map;
            $prefix = sprintf("%02d", ($prefix + 1));
            if ($options->{host_filter}) {
                print "\n${prefix}: Performing additional mapping against ${nth_species} with filtering.\n";
                $nth_map = $class->Bio::Adventure::Map::Hisat2(
                    jdepends => $last_job,
                    input => $first_map->{unaligned_comp},
                    species => $nth_species,
                    gff_type => $nth_type,
                    gff_tag => $nth_id,
                    jprefix => $prefix,);
                $last_sam_job = $nth_map->{samtools}->{job_id};
            } else {
                print "\n${prefix}: Performing additional mapping against ${nth_species} without filtering.\n";
                $nth_map = $class->Bio::Adventure::Map::Hisat2(
                    jdepends => $last_sam_job,
                    input => $trim->{output},
                    species => $nth_species,
                    gff_type => $nth_type,
                    gff_tag => $nth_id,
                    jprefix => $prefix,);
            }
            push(@jobs, $nth_map);
            $last_sam_job = $nth_map->{samtools}->{job_id};
            sleep($options->{jsleep});

            $prefix = sprintf("%02d", ($prefix + 1));
            print "\n${prefix}: Performing freebayes search against ${nth_species}.\n";
            my $nth_snp = $class->Bio::Adventure::SNP::Freebayes_SNP_Search(
                jdepends => $last_sam_job,
                input => $nth_map->{samtools}->{paired_output},
                species => $nth_species,
                gff_type => $nth_type,
                gff_tag => $nth_id,
                intron => $options->{intron},
                jprefix => $prefix,);
            $last_job = $nth_map->{job_id};
            push(@jobs, $nth_snp);
            sleep($options->{jsleep});
            $c++;
        } ## End iterating over extra species
    } ## End checking for extra species
    return(\@jobs);
}

sub RNAseq {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        gff_type => 'gene',
        gff_tag => 'ID',
        mapper => 'hisat2',);
    my $prefix = sprintf("%02d", 1);
    my $final_locustag = basename(cwd());

    print "Starting fastqc.\n";
    my $fastqc = $class->Bio::Adventure::QA::Fastqc(
        jprefix => $prefix,);
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nStarting trimmer.\n";
    my $trim = $class->Bio::Adventure::Trim::Trimomatic(
        input => $fastqc->{input},
        jprefix => $prefix,
        jname => 'trimomatic',);
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nStarting mapper.\n";
    my $mapper;
    if ($options->{mapper} eq 'bowtie1') {
        $mapper = $class->Bio::Adventure::Map::Bowtie(
            jdepends => $trim->{job_id},
            input => $trim->{output},
            jprefix => $prefix,
            species => $options->{species},
            gff_type => $options->{gff_type},
            gff_tag => $options->{gff_tag},);
    } elsif ($options->{mapper} eq 'bowtie2') {
        $mapper = $class->Bio::Adventure::Map::Bowtie2(
            jdepends => $trim->{job_id},
            input => $trim->{output},
            jprefix => $prefix,
            species => $options->{species},
            gff_type => $options->{gff_type},
            gff_tag => $options->{gff_tag},);
    } elsif ($options->{mapper} eq 'bwa') {
        $mapper = $class->Bio::Adventure::Map::BWA(
            jdepends => $trim->{job_id},
            input => $trim->{output},
            jprefix => $prefix,
            species => $options->{species},
            gff_type => $options->{gff_type},
            gff_tag => $options->{gff_tag},);
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
            gff_type => $options->{gff_type},
            gff_tag => $options->{gff_tag},);
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
    sleep($options->{jsleep});

    print "\nRunning Cutadapt with TNSeq adapters.\n";
    my $cutadapt = $class->Bio::Adventure::Trim::Cutadapt(
        input => $options->{input},
        task => 'tnseq',
        jprefix => $prefix,
        jname => 'cutadapt',);
    $last_job = $cutadapt->{job_id};
    sleep($options->{jsleep});

    print "\nRunning Biopieces.\n";
    my $pieces = $class->Bio::Adventure::QA::Biopieces_Graph(
        input => $cutadapt->{output},
        jprefix => $prefix,
        jname => 'biopieces',);
    $last_job = $pieces->{job_id};
    sleep($options->{jsleep});

    print "\nMapping remaining reads.\n";
    my $map_reads = $class->Bio::Adventure::Map::Bowtie(
        input => $cutadapt->{output},
        species => $options->{species},
        jdepends => $last_job,
        jprefix => $prefix,
        jname => 'map_reads',);
    $last_job = $map_reads->{job_id};
    sleep($options->{jsleep});

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
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nStarting fastqc.\n";
    my $fastqc = $class->Bio::Adventure::QA::Fastqc(
        input => $trim->{input},
        jmem => 3,
        jnice => 100,
        jprefix => $prefix,);
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nStarting read correction.\n";
    my $correct = $class->Bio::Adventure::Trim::Racer(
        jdepends => $trim->{job_id},
        input => $trim->{output},
        jprefix => $prefix,
        jname => 'racer',);
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nClassifying sequences with Kraken2.\n";
    my $kraken = $class->Bio::Adventure::Annotation::Kraken(
        jdepends => $correct->{job_id},
        input => $correct->{output},
        jprefix => $prefix,
        jname => 'kraken',
        library => 'viral',);
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nCreating initial assembly with Unicycler.\n";
    my $assemble = $class->Bio::Adventure::Assembly::Unicycler(
        jdepends => $correct->{job_id},
        input => $correct->{output},
        jprefix => $prefix,
        jname => 'unicycler',);
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning phastaf.\n";
    my $phastaf = $class->Bio::Adventure::Phage::Phastaf(
        jdepends => $assemble->{job_id},
        input => $assemble->{output},
        jprefix => $prefix,
        jname => 'phastaf',);
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nPerforming initial prokka annotation.\n";
    my $prokka = $class->Bio::Adventure::Annotation::Prokka(
        jdepends => $assemble->{job_id},
        input => $assemble->{output},
        jprefix => $prefix,
        jname => 'prokka',
        locus_tag => $final_locustag,);
    sleep($options->{jsleep});
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
        host_filter => 1,
        host_species => '',
        jsleep => 5,
        required => ['input'],);
    my $prefix = sprintf("%02d", 0);
    my $final_locustag = basename(cwd());

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Starting trimmer.\n";
    my $trim = $class->Bio::Adventure::Trim::Trimomatic(
        compress => 0,
        input => $options->{input},
        jmem => 4,
        jname => 'trimomatic',
        jprefix => $prefix,);
    my $last_job = $trim->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Starting fastqc.\n";
    my $fastqc = $class->Bio::Adventure::QA::Fastqc(
        input => $trim->{output},
        jmem => 3,
        jnice => 100,
        jprefix => $prefix,);
    $last_job = $fastqc->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Starting read correction.\n";
    my $correct = $class->Bio::Adventure::Trim::Racer(
        compress => 0,
        input => $trim->{output},
        jdepends => $last_job,
        jmem => 3,
        jprefix => $prefix,);
    $last_job = $correct->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Classifying sequences with Kraken using the standard database.\n";
    my $kraken_std = $class->Bio::Adventure::Annotation::Kraken(
        jdepends => $last_job,
        input => $correct->{output},
        jprefix => $prefix,
        jname => 'krakenstd',
        library => 'standard',);
    $last_job = $kraken_std->{job_id};
    sleep($options->{jsleep});

    my $filter = {};
    if ($options->{host_filter}) {
        $prefix = sprintf("%02d", ($prefix + 1));
        print "\n${prefix}: Using Kraken results to filter likely host sequences.\n";
        $filter = $class->Bio::Adventure::Phage::Filter_Host_Kraken(
            input => $kraken_std->{output},
            input_fastq => $correct->{output},
            jdepends => $last_job,
            jprefix => $prefix,
            jname => 'krakenfilter',);
        $last_job = $filter->{job_id};
        sleep($options->{jsleep});
    } else {
        print "\n: Not performing host filter.\n";
        $filter = $correct;
    }

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Classifying sequences with Kraken using the viral database.\n";
    my $kraken = $class->Bio::Adventure::Annotation::Kraken(
        jdepends => $last_job,
        input => $filter->{output},
        jprefix => $prefix,
        jname => 'kraken',
        library => 'viral',);
    $last_job = $kraken->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Creating initial assembly with Unicycler.\n";
    my $assemble = $class->Bio::Adventure::Assembly::Unicycler(
        jdepends => $last_job,
        input => $filter->{output},
        jprefix => $prefix,
        jname => 'unicycler',);
    $last_job = $assemble->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Depth filtering initial assembly.\n";
    my $depth_filtered = $class->Bio::Adventure::Assembly::Unicycler_Filter_Depth(
        jdepends => $last_job,
        input => $assemble->{output},
        jprefix => $prefix,
        jname => 'depth_filter',);
    $last_job = $depth_filtered->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Running phastaf.\n";
    my $phastaf = $class->Bio::Adventure::Phage::Phastaf(
        jdepends => $last_job,
        input => $depth_filtered->{output},
        jprefix => $prefix,
        jname => 'phastaf',);
    $last_job = $phastaf->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Running virus ICTV classifier.\n";
    my $ictv = $class->Bio::Adventure::Phage::Classify_Phage(
        jdepends => $last_job,
        input => $depth_filtered->{output},
        jprefix => $prefix,
        jname => 'ictv',);
    $last_job = $ictv->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Setting the Rosalind strand to the one with the most ORFs.\n";
    my $rosalindplus = $class->Bio::Adventure::Annotation::Rosalind_Plus(
        jdepends => $last_job,
        input => $depth_filtered->{output},
        jprefix => $prefix,
        jname => 'rosalindplus',);
    $last_job = $rosalindplus->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Using Phageterm to reorient the assembly to the terminii.\n";
    my $phageterm = $class->Bio::Adventure::Phage::Phageterm(
        jdepends => $last_job,
        input => $filter->{output},
        jprefix => $prefix,
        jname => 'phageterm',
        library => $rosalindplus->{output},);
    $last_job = $phageterm->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Checking coverage.\n";
    my $coverage = $class->Bio::Adventure::Assembly::Assembly_Coverage(
        jdepends => $last_job,
        input => $filter->{output},
        library => $assemble->{output},
        jprefix => $prefix,
        jname => 'coverage',);
    $last_job = $coverage->{job_id};
    sleep($options->{jsleep});

    print "\nCompressing raw filtered files.\n";
    my $compress_filtered = $class->Bio::Adventure::Compress::Compress(
        jdepends => $last_job,
        input => $filter->{output},
        jprefix => $prefix,);
    $last_job = $compress_filtered->{job_id};
    sleep($options->{jsleep});

    print "\nCompressing corrected fastq files.\n";
    my $compress_corrected = $class->Bio::Adventure::Compress::Compress(
        jdepends => $last_job,
        input => $correct->{output},
        jprefix => $prefix,);
    $last_job = $compress_corrected->{job_id};
    sleep($options->{jsleep});

    print "\nCompressing trimmed files.\n";
    my $compress_trimmed = $class->Bio::Adventure::Compress::Compress(
        input => $trim->{output},
        jdepends => $last_job,
        jprefix => $prefix,);
    $last_job = $compress_trimmed->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Performing reordering by putative terminase if phageterm failed.\n";
    my $termreorder = $class->Bio::Adventure::Phage::Terminase_ORF_Reorder(
        fasta_tool => 'fastx36',
        input => $phageterm->{output},
        jdepends => $last_job,
        jname => 'reorder',
        jprefix => $prefix,
        library => 'terminase',
        species => 'phages',
        test_file => $phageterm->{test_file},);
    $last_job = $termreorder->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Performing initial prokka annotation.\n";
    my $prokka = $class->Bio::Adventure::Annotation::Prokka(
        input => $termreorder->{output},
        jdepends => $last_job,
        jname => 'prokka',
        jprefix => $prefix,
        locus_tag => $final_locustag,);
    $last_job = $prokka->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Running prodigal to get RBSes.\n";
    my $prodigal = $class->Bio::Adventure::Feature_Prediction::Prodigal(
        input => $prokka->{output_assembly},
        jdepends => $last_job,
        ## Use output_assembly to avoid having too much clutter on the sequence ID line.
        jprefix => $prefix,);
    $last_job = $prodigal->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Running glimmer for less stringent ORFs.\n";
    my $glimmer = $class->Bio::Adventure::Feature_Prediction::Glimmer_Single(
        input => $prokka->{output_assembly},
        jdepends => $last_job,
        jprefix => $prefix,);
    $last_job = $glimmer->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Running phanotate.\n";
    my $phanotate = $class->Bio::Adventure::Feature_Prediction::Phanotate(
        input => $prokka->{output_assembly},
        jdepends => $last_job,
        jprefix => $prefix,);
    $last_job = $phanotate->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Merging the CDS predictions and writing an initial gbk file.\n";
    my $cds_merge = $class->Bio::Adventure::Annotation_Genbank::Merge_CDS_Predictions(
        input => $prokka->{output_genbank},
        input_glimmer => $glimmer->{output},
        input_phanotate => $phanotate->{output},
        input_prodigal => $prodigal->{output_gff},
        jdepends => $last_job,
        jprefix => $prefix);
    $last_job = $cds_merge->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Running Jellyfish on the assembly.\n";
    my $jelly = $class->Bio::Adventure::Count::Jellyfish(
        input => $cds_merge->{output_fsa},
        jdepends => $last_job,
        jprefix => $prefix,
        jname => 'jelly',);
    $last_job = $jelly->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Running aragorn on the assembly to search for tmRNAs.\n";
    my $aragorn = $class->Bio::Adventure::Feature_Prediction::Aragorn(
        input => $cds_merge->{output_fsa},
        jdepends => $last_job,
        jprefix => $prefix,
        jname => 'aragorn',);
    $last_job = $aragorn->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Running relaxed tRNAscan to search for tmRNAs.\n";
    my $trnascan = $class->Bio::Adventure::Feature_Prediction::tRNAScan(
        arbitrary => ' -r ',
        input => $prokka->{output_assembly},
        jdepends => $last_job,
        jname => 'trnascan',
        jprefix => $prefix,
        suffix => 'relaxed',);
    $last_job = $trnascan->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Invoking trinotate.\n";
    my $trinotate = $class->Bio::Adventure::Annotation::Trinotate(
        config => 'phage.txt',
        input => $cds_merge->{output_cds},
        jdepends => $last_job,
        jprefix => $prefix,
        jname => 'trinotate',);
    $last_job = $trinotate->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Searching for resistance genes with abricate.\n";
    my $abricate = $class->Bio::Adventure::Resistance::Abricate(
        input => $cds_merge->{output_cds},
        jdepends => $last_job,
        jprefix => $prefix,
        jname => 'abricate',);
    $last_job = $abricate->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Running interproscan.\n";
    my $interpro = $class->Bio::Adventure::Annotation::Interproscan(
        input => $cds_merge->{output_faa},
        jdepends => $last_job,
        jprefix => $prefix,
        jname => 'interproscan',);
    $last_job = $interpro->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Merging annotation files.\n";
    my $merge = $class->Bio::Adventure::Metadata::Merge_Annotations(
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
        jdepends => $last_job,
        jprefix => $prefix,
        jname => 'mergeannot',);
    $last_job = $merge->{job_id};
    sleep($options->{jsleep});

    ## An extra invocation of merge_annotations which will not modify the final gbk file.
    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Merging annotation files a second time.\n";
    my $merge2 = $class->Bio::Adventure::Metadata::Merge_Annotations(
        evalue => 0,
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
        jdepends => $last_job,
        jprefix => $prefix,
        jname => 'mergeannot2',
        keep_genes => 0,
        locus_tag => 0,
        suffix => 'stripped',);
    sleep($options->{jsleep});
    $last_job = $merge2->{job_id};

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Running cgview.\n";
    my $cgview = $class->Bio::Adventure::Visualization::CGView(
        input => $merge->{output_gbk},
        jdepends => $last_job,
        jprefix => $prefix,);
    ## Not setting last_job, allowing the next jobs to skip past this.
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Running Vienna RNAfold on the assembly.\n";
    my $vienna = $class->Bio::Adventure::Structure::RNAFold_Windows(
        input => $merge->{output_fsa},
        jdepends => $last_job,
        jprefix => $prefix,
        jname => 'vienna',);
    ## Not setting last_job, allowing the next jobs to skip past this.
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Searching for Restriction Sites.\n";
    my $re_search = $class->Bio::Adventure::Phage::Restriction_Catalog(
        input => $merge->{output_fsa},
        jdepends => $last_job,
        jprefix => $prefix,
        jname => 'restrict',);
    ## Not setting last_job, allowing the next jobs to skip past this.
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Running caical against the assumed host.\n";
    my $caical = $class->Bio::Adventure::Phage::Caical(
        input => $cds_merge->{output_cds},
        jdepends => $last_job,
        jprefix => $prefix,
        jname => 'restrict',
        species => 'host_species.txt',);
    ## Not setting last_job, allowing the next jobs to skip past this.
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Searching for Phage promoters.\n";
    my $phagepromoter = $class->Bio::Adventure::Feature_Prediction::Phagepromoter(
        input => $merge->{output_fsa},
        jdepends => $last_job,
        jprefix => $prefix,
        jname => 'phagepromoter',);
    ## Not setting last_job, allowing the next jobs to skip past this.
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Searching for rho terminators.\n";
    my $rhopredict = $class->Bio::Adventure::Feature_Prediction::Rho_Predict(
        input => $merge->{output_fsa},
        jdepends => $last_job,
        jprefix => $prefix,
        jname => 'rhopredict',);
    ## Not setting last_job, allowing the next jobs to skip past this.
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Using bacphlip to classify this phage.\n";
    my $bacphlip = $class->Bio::Adventure::Phage::Bacphlip(
        input => $merge->{output_fsa},
        jdepends => $last_job,
        jprefix => $prefix,
        jname => 'bacphlip',);
    $last_job = $bacphlip->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Collecting output files.\n";
    my $collect = $class->Bio::Adventure::Metadata::Collect_Assembly(
        input_fsa => $merge->{output_fsa},
        input_genbank => $merge->{output_gbk},
        input_stripped => $merge2->{output_gbk},
        input_tbl => $merge->{output_tbl},
        input_tbl_stripped => $merge2->{output_tbl},
        input_tsv => $merge->{output_tsv},
        input_xlsx => $merge->{output_xlsx},
        input_faa => $cds_merge->{output_faa},
        input_cds => $cds_merge->{output_cds},
        jdepends => $last_job,
        jprefix => $prefix,
        jname => 'collect',);
    $last_job = $collect->{job_id};
    sleep($options->{jsleep});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\n${prefix}: Cleaning fastq files.\n";
    my $clean_phage = $class->Bio::Adventure::Cleanup::Cleanup_Phage_Assembly(
        jname => 'cleanphage',);
    $last_job = $clean_phage->{job_id};
    sleep($options->{jsleep});

    my $ret = {
        '01trim' => $trim,
        '02fastqc' => $fastqc,
        '03correction' => $correct,
        '04kraken_standard' => $kraken_std,
        '05host_filter' => $filter,
        '06kraken_viral' => $kraken,
        '07assembly' => $assemble,
        '08depth_filter' => $depth_filtered,
        '09phastaf' => $phastaf,
        '10ictv' => $ictv,
        '11rosalindplus' => $rosalindplus,
        '12phageterm' => $phageterm,
        '13coverage' => $coverage,
        '14terminase_reorder' => $termreorder,
        '15prokka' => $prokka,
        '16prodigal' => $prodigal,
        '17glimmer' => $glimmer,
        '18phanotate' => $phanotate,
        '19cds_merge' => $cds_merge,
        '20jellyfish' => $jelly,
        '21aragorn' => $aragorn,
        '22trnascan' => $trnascan,
        '23trinotate' => $trinotate,
        '24abricate' => $abricate,
        '25interproscan' => $interpro,
        '26merge_qualities' => $merge,
        '27merge_unmodified' => $merge2,
        '28cgview' => $cgview,
        '29rnafold' => $vienna,
        '30research' => $re_search,
        '31caical' => $caical,
        '32phagepromoter' => $phagepromoter,
        '33rhopredict' => $rhopredict,
        '34bacphlip' => $bacphlip,
        '35collect' => $collect,
    };
    print "Finished submitting phage assembly jobs.  Final output should reside in
the collection directory upon completion.\n";
    sleep(10);
    return($ret)
}

1;
