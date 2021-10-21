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

sub Riboseq {
    my ($class, %args) = @_;
    my $fastqc_job = Bio::Adventure::QA::Fastqc($class, %args);
    my $cutadapt_job = Bio::Adventure::Trim::Cutadapt($class, %args);
    $args{jdepends} = $cutadapt_job->{job_id};
    my $biopieces = Bio::Adventure::QA::Biopieces_Graph($class, %args);
    my $rrna_job = Bio::Adventure::Map::Bowtie_RRNA($class, %args);
    $args{jdepends} = $rrna_job->{job_id};
    my $bt_jobs = Bio::Adventure::Map::Bowtie($class, %args);
    my $ret = {
        fastqc => $fastqc_job,
        cutadapt => $cutadapt_job,
        biopieces => $biopieces,
        rrna => $rrna_job,
        bt => $bt_jobs};
    return($ret);
}

sub Bowtie {
    my ($class, %args) = @_;
    $args{aligner} = 'bowtie';
    my $rnaseq_jobs = $class->Pipeline_RNAseq(%args);
    return($rnaseq_jobs);
}

sub Bowtie2 {
    my ($class, %args) = @_;
    $args{aligner} = 'bowtie2';
    my $rnaseq_jobs = $class->Pipeline_RNAseq(%args);
    return($rnaseq_jobs);
}

sub BWA {
    my ($class, %args) = @_;
    $args{aligner} = 'bwa';
    my $rnaseq_jobs = $class->Pipeline_RNAseq(%args);
    return($rnaseq_jobs);
}

sub Hisat {
    my ($class, %args) = @_;
    $args{aligner} = 'hisat';
    my $rnaseq_jobs = $class->Pipeline_RNAseq(%args);
    return($rnaseq_jobs);
}

sub Kallisto {
    my ($class, %args) = @_;
    $args{aligner} = 'kallisto';
    my $rnaseq_jobs = $class->Pipeline_RNAseq(%args);
    return($rnaseq_jobs);
}

sub RNAseq {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        htseq_type => 'gene',
        htseq_id => 'ID',);
    my $prefix = sprintf("%02d", 1);
    my $final_locustag = basename(cwd());

    print "Starting fastqc.\n";
    my $fastqc = $class->Bio::Adventure::QA::Fastqc(
        jprefix => $prefix,);
    sleep(1);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nStarting trimmer.\n";
    my $trim = $class->Bio::Adventure::Trim::Trimomatic(
        input => $fastqc->{input},
        jprefix => $prefix,
        jname => 'trimomatic',);
    sleep(1);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nStarting mapper.\n";
    my $mapper;
    if ($args{aligner} eq 'bowtie') {
        $mapper = $class->Bio::Adventure::Map::Bowtie(
            jdepends => $trim->{job_id},
            input => $trim->{output},
            jprefix => $prefix,
            species => $options->{species},
            htseq_type => $options->{htseq_type},
            htseq_id => $options->{htseq_id},);
    } elsif ($args{aligner} eq 'bowtie2') {
        $mapper = $class->Bio::Adventure::Map::Bowtie2(
            jdepends => $trim->{job_id},
            input => $trim->{output},
            jprefix => $prefix,
            species => $options->{species},
            htseq_type => $options->{htseq_type},
            htseq_id => $options->{htseq_id},);
    } elsif ($args{aligner} eq 'bwa') {
        $mapper = $class->Bio::Adventure::Map::BWA(
            jdepends => $trim->{job_id},
            input => $trim->{output},
            jprefix => $prefix,
            species => $options->{species},
            htseq_type => $options->{htseq_type},
            htseq_id => $options->{htseq_id},);
    } elsif ($args{aligner} eq 'kallisto') {
        $mapper = $class->Bio::Adventure::Map::Kallisto(
            jdepends => $trim->{job_id},
            input => $trim->{output},
            jprefix => $prefix,
            species => $options->{species},);
    } elsif ($args{aligner} eq 'salmon') {
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

    return($mapper);
}

sub TNseq {
    my ($class, %args) = @_;
    my $fastqc_job = Bio::Adventure::QA::Fastqc($class, %args);
    $args{type} = 'tnseq';
    my $cutadapt_job = Bio::Adventure::Trim::Cutadapt($class, %args);
    $args{jdepends} = $cutadapt_job->{job_id};
    my $biopieces_job = Bio::Adventure::QA::Biopieces_Graph($class, %args);
    my $bt_jobs = Bio::Adventure::Map::Bowtie($class, %args);
    my $ret = {
        fastqc => $fastqc_job,
        cutadapt => $cutadapt_job,
        biopieces => $biopieces_job,
        bt => $bt_jobs,
    };
    return($ret);
}


sub Assemble {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        );
    my $final_locustag = basename(cwd());

    my $prefix = sprintf("%02d", 1);
    print "\nStarting trimmer.\n";
    my $trim = $class->Bio::Adventure::Trim::Trimomatic(
        input => $options->{input},
        jprefix => $prefix,
        jname => 'trimomatic',);
    sleep(1);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nStarting fastqc.\n";
    my $fastqc = $class->Bio::Adventure::QA::Fastqc(
        input => $trim->{input},
        jmem => 3,
        jnice => 100,
        jprefix => $prefix,);
    sleep(1);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nStarting read correction.\n";
    my $correct = $class->Bio::Adventure::Trim::Racer(
        jdepends => $trim->{job_id},
        input => $trim->{output},
        jprefix => $prefix,
        jname => 'racer',);
    sleep(1);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nClassifying sequences with Kraken2.\n";
    my $kraken = $class->Bio::Adventure::Annotation::Kraken(
        jdepends => $correct->{job_id},
        input => $correct->{output},
        jprefix => $prefix,
        jname => 'kraken',
        library => 'viral',);
    sleep(1);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nCreating initial assembly with Unicycler.\n";
    my $assemble = $class->Bio::Adventure::Assembly::Unicycler(
        jdepends => $correct->{job_id},
        input => $correct->{output},
        jprefix => $prefix,
        jname => 'unicycler',);
    sleep(1);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning phastaf.\n";
    my $phastaf = $class->Bio::Adventure::Phage::Phastaf(
        jdepends => $assemble->{job_id},
        input => $assemble->{output},
        jprefix => $prefix,
        jname => 'phastaf',);
    sleep(1);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nPerforming initial prokka annotation.\n";
    my $prokka = $class->Bio::Adventure::Annotation::Prokka(
        jdepends => $assemble->{job_id},
        input => $assemble->{output},
        jprefix => $prefix,
        jname => 'prokka',
        locus_tag => $final_locustag,);
    sleep(1);
}

=head2 C<PhageAssemble>

Time to collapse the set of tools I try out and make a simplified, but subject to change
pipeline of assembly tasks for phage assembly.

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
    sleep(1);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nStarting fastqc.\n";
    my $fastqc = $class->Bio::Adventure::QA::Fastqc(
        input => $trim->{input},
        jmem => 3,
        jnice => 100,
        jprefix => $prefix,);
    sleep(1);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nStarting read correction\n";
    my $correct = $class->Bio::Adventure::Trim::Racer(
        jdepends => $trim->{job_id},
        input => $trim->{output},
        jmem => 3,
        jprefix => $prefix,);
    sleep(1);

    ## Take a moment and get the host species for sequence filtering
    my $host_species = $options->{host_species};
    if (-r "host_species.txt") {
        my $host_file = FileHandle->new("<host_species.txt");
        while (my $line = <$host_file>) {
            chomp $line;
            $host_species = $line;
        }
        $host_file->close();
    }
    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nStarting filter with hisat2.\n";
    my $filter = $class->Bio::Adventure::Map::Hisat2(
        jdepends => $correct->{job_id},
        input => $correct->{output},
        jprefix => $prefix,
        do_htseq => 0,
        jname => 'hisatfilter',
        species => $host_species,);
    my $filtered_reads = $filter->{unaligned_comp};
    my $filter_id = $filter->{job_id};
    sleep(1);

    ## 05
    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nClassifying sequences with Kraken2 using the viral database.\n";
    my $kraken = $class->Bio::Adventure::Annotation::Kraken(
        jdepends => $filter->{job_id},
        input => $filtered_reads,
        jprefix => $prefix,
        jname => 'kraken',
        library => 'viral',);
    sleep(1);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nClassifying sequences with Kraken2 using the standard database.\n";
    my $kraken_std = $class->Bio::Adventure::Annotation::Kraken(
        jdepends => $filter->{job_id},
        input => $filtered_reads,
        jprefix => $prefix,
        jname => 'krakenstd',
        library => 'standard',);
    sleep(1);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nCreating initial assembly with Unicycler.\n";
    my $assemble = $class->Bio::Adventure::Assembly::Unicycler(
        jdepends => $filter->{job_id},
        input => $filtered_reads,
        jprefix => $prefix,
        jname => 'unicycler',);
    sleep(1);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nDepth filtering initial assembly.\n";
    my $depth_filtered = $class->Bio::Adventure::Assembly::Filter_Depth(
        jdepends => $assemble->{job_id},
        input => $assemble->{output},
        jprefix => $prefix,
        jname => 'depth_filter',);
    sleep(1);

    ## 08
    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning phastaf.\n";
    my $phastaf = $class->Bio::Adventure::Phage::Phastaf(
        jdepends => $depth_filtered->{job_id},
        input => $depth_filtered->{output},
        jprefix => $prefix,
        jname => 'phastaf',);
    sleep(1);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning virus ICTV classifier.\n";
    my $ictv = $class->Bio::Adventure::Phage::Classify_Phage(
        jdepends => $phastaf->{job_id},
        input => $depth_filtered->{output},
        jprefix => $prefix,
        jname => 'ictv',);
    sleep(1);

    ## 10
    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nSetting the Watson strand to the one with the most ORFs.\n";
    my $watsonplus = $class->Bio::Adventure::Annotation::Watson_Plus(
        jdepends => $ictv->{job_id},
        input => $depth_filtered->{output},
        jprefix => $prefix,
        jname => 'watsonplus',);
    sleep(1);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nUsing Phageterm to reorient the assembly to the terminii.\n";
    my $phageterm = $class->Bio::Adventure::Phage::Phageterm(
        jdepends => $watsonplus->{job_id},
        input => $filtered_reads,
        jprefix => $prefix,
        jname => 'phageterm',
        library => $watsonplus->{output},);
    sleep(1);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nPerforming reordering by putative terminase.\n";
    my $termreorder = $class->Bio::Adventure::Phage::Terminase_ORF_Reorder(
        jdepends => $phageterm->{job_id},
        fasta_tool => 'fastx36',
        input => $phageterm->{output},
        jprefix => $prefix,
        jname => 'reorder',
        library => 'terminase',
        species => 'phages',
        test_file => $phageterm->{test_file},);
    sleep(1);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nPerforming initial prokka annotation.\n";
    my $prokka = $class->Bio::Adventure::Annotation::Prokka(
        jdepends => $termreorder->{job_id},
        input => $termreorder->{output},
        jprefix => $prefix,
        jname => 'prokka',
        locus_tag => $final_locustag,);
    sleep(1);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning Jellyfish on the assembly.\n";
    my $jelly = $class->Bio::Adventure::Count::Jellyfish(
        jdepends => $prokka->{job_id},
        input => $prokka->{output_assembly},
        jprefix => $prefix,
        jname => 'jelly',);
    sleep(1);

    ## 15
    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nInvoking trinotate.\n";
    my $trinotate = $class->Bio::Adventure::Annotation::Trinotate(
        jdepends => $prokka->{job_id},
        input => $prokka->{output},
        jprefix => $prefix,
        jname => 'trinotate',
        config => 'phage.txt',);
    sleep(1);

    ## 16
    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nSearching for resistance genes with abricate.\n";
    my $abricate = $class->Bio::Adventure::Resistance::Abricate(
        jprefix => $prefix,
        jname => 'abricate',
        input => $prokka->{output},
        jdepends => $trinotate->{job_id},);
    sleep(1);

    ## 17
    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning interproscan.\n";
    my $interpro = $class->Bio::Adventure::Annotation::Interproscan(
        jprefix => $prefix,
        jname => 'interproscan',
        input => $prokka->{output_peptide},
        jdepends => $abricate->{job_id},);
    sleep(1);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning prodigal to get RBSes.\n";
    my $prodigal = $class->Bio::Adventure::Annotation::Prodigal(
        input => $prokka->{output_fsa},
        jdepends => $interpro->{job_id},
        jprefix => $prefix);
    sleep(1);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nMerging annotation files.\n";
    my $merge = $class->Bio::Adventure::Annotation::Merge_Annotations(
        input_abricate => $abricate->{output},
        input_classifier => $ictv->{output},
        input_fsa => $prokka->{output_fsa},
        input_genbank => $prokka->{output_genbank},
        input_interpro => $interpro->{output_tsv},
        input_phageterm => $phageterm->{test_file},
        input_prodigal => $prodigal->{output},
        input_prokka_tsv => $prokka->{output_tsv},
        input_trinotate => $trinotate->{output},
        jprefix => $prefix,
        jname => 'mergeannot',
        jdepends => $prodigal->{job_id},);
    sleep(1);

    ## 20
    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning cgview.\n";
    my $cgview = $class->Bio::Adventure::Visualization::CGView(
        input => $merge->{output_gbk},
        jdepends => $merge->{job_id},
        jprefix => $prefix);
    sleep(1);

    ## 21
    ## An extra invocation of merge_annotations which will not modify the final gbk file.
    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nMerging annotation files a second time.\n";
    my $merge2 = $class->Bio::Adventure::Annotation::Merge_Annotations(
        input_abricate => $abricate->{output},
        input_classifier => $ictv->{output},
        input_fsa => $prokka->{output_fsa},
        input_genbank => $prokka->{output_genbank},
        input_interpro => $interpro->{output_tsv},
        input_phageterm => $phageterm->{test_file},
        input_prodigal => $prodigal->{output},
        input_prokka_tsv => $prokka->{output_tsv},
        input_trinotate => $trinotate->{output},
        jprefix => $prefix,
        jnice => '1000',
        jname => 'mergeannot2',
        evalue => undef,
        jdepends => $prodigal->{job_id},);
    sleep(1);

    my $ret = {
        trim => $trim,
        fastqc => $fastqc,
        correction => $correct,
        filter => $filter,
        kraken_viral => $kraken,
        kraken_standard => $kraken_std,
        assembly => $assemble,
        depth_filter => $depth_filtered,
        phastaf => $phastaf,
        ictv => $ictv,
        watsonplus => $watsonplus,
        phageterm => $phageterm,
        terminase_reorder => $termreorder,
        prokka => $prokka,
        jellyfish => $jelly,
        trinotate => $trinotate,
        abricate => $abricate,
        interproscan => $interpro,
        prodigal => $prodigal,
        cgview => $cgview,
        merge_qualities => $merge,
        merge_unmodified => $merge2,
    };
    return($ret)
}

sub Annotate_Assembly {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input',],);
    my $prefix = sprintf("%02d", 10);
    my $final_locustag = basename(cwd());
    my $filtered_dir = dirname($options->{input});

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning phastaf.\n";
    my $phastaf = $class->Bio::Adventure::Phage::Phastaf(
        input => $options->{input},
        jprefix => $prefix,
        jname => 'phastaf',);
    sleep(1);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning virus ICTV classifier.\n";
    my $ictv = $class->Bio::Adventure::Phage::Blast_Classify(
        input => $options->{input},
        jprefix => $prefix,
        jname => 'ictv',);
    sleep(1);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nPerforming initial prokka annotation.\n";
    my $prokka = $class->Bio::Adventure::Annotation::Prokka(
        input => $options->{input},
        jprefix => $prefix,
        jname => 'prokka',
        locus_tag => $final_locustag,);
    sleep(1);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning Jellyfish on the assembly.\n";
    my $jelly = $class->Bio::Adventure::Count::Jellyfish(
        jdepends => $prokka->{job_id},
        input => $prokka->{output_assembly},
        jprefix => $prefix,
        jname => 'jelly',);
    sleep(1);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nInvoking trinotate.\n";
    my $trinotate = $class->Bio::Adventure::Annotation::Trinotate(
        jdepends => $prokka->{job_id},
        input => $prokka->{output},
        jprefix => $prefix,
        jname => 'trinotate',
        config => 'phage.txt',);
    sleep(1);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nSearching for resistance genes with abricate.\n";
    my $abricate = $class->Bio::Adventure::Resistance::Abricate(
        jprefix => $prefix,
        jname => 'abricate',
        input => $prokka->{output},
        jdepends => $prokka->{job_id},);
    sleep(1);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning interproscan.\n";
    my $interpro = $class->Bio::Adventure::Annotation::Interproscan(
        jprefix => $prefix,
        jname => 'interproscan',
        input => $prokka->{output_peptide},
        jdepends => $prokka->{job_id},);
    sleep(1);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning prodigal to get RBSes.\n";
    my $prodigal = $class->Bio::Adventure::Annotation::Prodigal(
        input => $prokka->{output_fsa},
        jdepends => $prokka->{job_id},
        jprefix => $prefix);
    sleep(1);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nMerging annotation files.\n";
    my $merge = $class->Bio::Adventure::Annotation::Merge_Annotations(
        input_abricate => $abricate->{output},
        input_fsa => $prokka->{output_fsa},
        input_genbank => $prokka->{output_genbank},
        input_interpro => $interpro->{output_tsv},
        input_prodigal => $prodigal->{output},
        input_prokka_tsv => $prokka->{output_tsv},
        input_trinotate => $trinotate->{output},
        jprefix => $prefix,
        jname => 'mergeannot',
        jdepends => $interpro->{job_id},);
    sleep(1);

    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nRunning cgview.\n";
    my $cgview = $class->Bio::Adventure::Visualization::CGView(
        input => $merge->{output_gbk},
        jdepends => $merge->{job_id},
        jprefix => $prefix);
    sleep(1);

    ## An extra invocation of merge_annotations which will not modify the final gbk file.
    $prefix = sprintf("%02d", ($prefix + 1));
    print "\nMerging annotation files a second time.\n";
    my $merge2 = $class->Bio::Adventure::Annotation::Merge_Annotations(
        input_abricate => $abricate->{output},
        input_fsa => $prokka->{output_fsa},
        input_genbank => $prokka->{output_genbank},
        input_interpro => $interpro->{output_tsv},
        input_prodigal => $prodigal->{output},
        input_prokka_tsv => $prokka->{output_tsv},
        input_trinotate => $trinotate->{output},
        jprefix => $prefix,
        jnice => '1000',
        jname => 'mergeannot2',
        evalue => undef,
        jdepends => $interpro->{job_id},);
    sleep(1);

    my $ret = {
        phastaf => $phastaf,
        ictv => $ictv,
        prokka => $prokka,
        jellyfish => $jelly,
        trinotate => $trinotate,
        abricate => $abricate,
        interproscan => $interpro,
        prodigal => $prodigal,
        cgview => $cgview,
        merge_qualities => $merge,
        merge_unmodified => $merge2,
    };
    return($ret)
}

1;
