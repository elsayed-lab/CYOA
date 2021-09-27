package Bio::Adventure::Annotation;
## LICENSE: gplv2
## ABSTRACT:  Kitty!
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use Bio::SeqIO;
use Bio::Seq;
use Bio::SeqFeature::Generic;
use Bio::Tools::GuessSeqFormat;
use Cwd qw"abs_path getcwd cwd";
use Data::Table;
use Data::Table::Excel qw"tables2xlsx";
use File::Basename;
use File::Copy qw"cp";
use File::Spec;
use File::Path qw"make_path";
use File::Which qw"which";
use File::ShareDir qw":ALL";
use List::MoreUtils qw"any";
use Template;
use Text::CSV_XS::TSV;

=head2 C<Aragorn>

Use aragorn to search for tRNA genes in a sequence database.  By
default this will explictly search for tmRNA as well.

=cut
sub Aragorn {
    my ($class, %args) = @_;
    my $check = which('aragorn');
    die("Could not find aragorn in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        modules => ['aragorn'],
        species => undef,
        arbitrary => ' -rp -fasta -w -m -t -mt ',
        );
    my $aragorn_args = $options->{arbitrary};
    my $job_name = $class->Get_Job_Name();
    my $output_dir = qq"outputs/aragorn";
    my $species_string = qq"";
    my $comment = qq!## This is a script to run aragorn.
!;
    my $jstring = qq!mkdir -p ${output_dir} && \\
  aragorn $options->{arbitrary} \\
    -o ${output_dir}/aragorn.txt \\
    $options->{input}
!;

    my $aragorn = $class->Submit(
        cpus => 6,
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => "aragorn_${job_name}",
        jprefix => "64",
        jstring => $jstring,
        jmem => 24,
        modules => $options->{modules},
        output => qq"${output_dir}/aragorn.txt",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => 'workstation',);
    return($aragorn);
}

=head2 C<Extend_Kraken_DB>

Add some more sequences to an existing kraken2 database.

=cut
sub Extend_Kraken_DB {
    my ($class, %args) = @_;
    my $check = which('kraken2');
    die("Could not find kraken2 in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        library => 'viral',
        modules => ['kraken'],);
    ## kraken2 --db ${DBNAME} --paired --classified-out cseqs#.fq seqs_1.fq seqs_2.fq
    my $job_name = $class->Get_Job_Name();
    my $output_dir = qq"outputs/extend_kraken";

    my $comment = qq!## This is a script to extend an existing kraken2 library with some new sequences.
!;
    my $jstring = qq!mkdir -p ${output_dir}
kraken2-build --download-taxonomy --db \${KRAKEN_DB_PATH}/$options->{library} \\
              2>${output_dir}/kraken2-build.out 1>&2
kraken2-build --add-to-library $options->{input} --db \${KRAKEN_DB_PATH}/$options->{library} \\
              2>>${output_dir}/kraken2-build.out 1>&2
kraken2-build --download-library $options->{library} \\
              --db \${KRAKEN_DB_PATH}/$options->{library} \\
              2>>${output_dir}/kraken2-build.out 1>&2
kraken2-build --build --db \${KRAKEN_DB_PATH}/$options->{library} \\
              2>>${output_dir}/kraken2-build.out 1>&2
!;
    my $kraken = $class->Submit(
        cpus => 6,
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => "kraken_${job_name}",
        jprefix => "99",
        jstring => $jstring,
        jmem => 96,
        modules => $options->{modules},
        output => qq"${output_dir}/kraken2-build.out",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => "large",);
    return($kraken);
}

=head2 C<Glimmer>

Use glimmer in two passes to search for ORFs in a sequence database.

=cut
sub Glimmer {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        modules => ['glimmer'],
        jprefix => '16',);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('glimmer3');
    die("Could not find glimmer in your PATH.") unless($check);
    my $job_name = $class->Get_Job_Name();
    my $output_dir = qq"outputs/$options->{jprefix}glimmer";

    my $comment = qq!## This is a script to run glimmer.
!;
##    my $jstring = qq!mkdir -p ${output_dir}
##long-orfs -n -t 1.15 $options->{input} ${output_dir}/first_run_longorfs.txt \\
##  2> ${output_dir}/first_run_longorfs.out 1>&2
##extract -t $options->{input} ${output_dir}/first_run_longorfs.txt > \\
##  ${output_dir}/first_run_training.txt
##build-icm -r ${output_dir}/first_run.icm < ${output_dir}/first_run_training.txt
##glimmer3 -o50 -g110 -t30 \\
##  $options->{input} \\
##  ${output_dir}/first_run.icm \\
##  ${output_dir}/first_run.out \\
##  2>${output_dir}first_run_glimmer3.out 1>&2
##
#### Use this first run output to go again
##
##extract -t $options->{input} train.coords > ${output_dir}/second_run_training.txt
##build-icm -r ${output_dir}/second_run.icm < ${output_dir}/second_run_training.txt
##upstream-coords.awk 25 0 train.coords | extract $options->{input} - > \\
##  ${output_dir}/second_run_upstream.txt
##elph ${output_dir}/second_run_upstream.txt LEN=6 | get-motif-counts.awk > \\
##  ${output_dir}/second_run_motif.txt
##startuse='start-codon-distrib -3 $options->{input} train.coords'
##glmmer3 -o50 -g110 -t30 -b ${output_dir}/second_run_motif.txt -P \${startuse} \\
##  $options->{input} ${output_dir}/second_run.icm \\
##  ${output_dir}/second_run.out
##!;
    my $jstring = qq!
cyoa_invoke_glimmer.pl --input $options->{input} --jprefix $options->{jprefix}
!;

    ## FIXME: There are a bunch of potentially useful glimmer outputs which should be put here.

    my $glimmer = $class->Submit(
        cpus => $options->{cpus},
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => "glimmer_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => 24,
        modules => $options->{modules},
        output => qq"${output_dir}/${job_name}_glimmer.out",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => "workstation",);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($glimmer);
}

=head2 C<Interproscan>

Use interproscan to look for existing annotations which are similar to
a set of provided ORFs.

=cut
sub Interproscan {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jprefix => '21',
        modules => ['interproscan'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('interproscan.sh');
    die("Could not find interproscan in your PATH.") unless($check);

    my $job_name = $class->Get_Job_Name();
    my $inputs = $class->Get_Paths($options->{input});
    my $cwd_name = basename(cwd());

    my $interproscan_exe_dir = dirname($check);
    ## Hey, don't forget abs_path requires a file which already exists.
    my $input_filename = basename($options->{input});
    my $output_filename = qq"${input_filename}.tsv";
    my $input_dir = dirname($options->{input});
    my $input_dirname = basename($input_dir);
    my $input_path = abs_path($input_dir);
    $input_path = qq"${input_path}/${input_filename}";
    my $output_dir = qq"outputs/$options->{jprefix}interproscan_${input_dirname}";
    my $comment = qq!## This is a interproscan submission script
!;
    my $jstring = qq!mkdir -p ${output_dir}
start=\$(pwd)
cd ${output_dir}
interproscan.sh -i ${input_path} 2>interproscan.err \\
  1>interproscan.out
ln -s ${output_filename} interproscan.tsv
cd \${start}
!;
    my $interproscan = $class->Submit(
        cpus => 6,
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => "interproscan_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => 80,
        modules => $options->{modules},
        output => qq"${output_dir}/interproscan.tsv",
        output_gff => qq"${output_dir}/${input_dirname}.faa.gff3",
        output_tsv => qq"${output_dir}/interproscan.tsv",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => "large",
        walltime => "144:00:00",);

    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($interproscan);
}

=head2 C<Kraken>

Use kraken2 to taxonomically classify reads.

=cut
sub Kraken {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        library => 'viral',
        jprefix => '11',
        modules => ['kraken'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('kraken2');
    die("Could not find kraken2 in your PATH.") unless($check);
    ## kraken2 --db ${DBNAME} --paired --classified-out cseqs#.fq seqs_1.fq seqs_2.fq
    my $job_name = $class->Get_Job_Name();
    my $input_directory = basename(cwd());
    my $output_dir = qq"outputs/$options->{jprefix}kraken_$options->{library}";
    make_path($output_dir);
    my $input_string = "";
    if ($options->{input} =~ /\:|\;|\,|\s+/) {
        my @in = split(/\:|\;|\,|\s+/, $options->{input});
        $input_string = qq" --paired <(less $in[0]) <(less $in[1]) ";
    } else {
        $input_string = qq"<(less $options->{input}) ";
    }
    my $comment = qq!## This is a kraken2 submission script
!;
    my $jstring = qq!kraken2 --db $ENV{KRAKEN2_DB_PATH}/$options->{library} \\
  --report ${output_dir}/kraken_report.txt --use-mpa-style \\
  --use-names ${input_string} \\
  --classified-out ${output_dir}/classified#.fastq.gz \\
  --unclassified-out ${output_dir}/unclassified#.fastq.gz \\
  2>${output_dir}/kraken.out 1>&2
!;
    my $kraken = $class->Submit(
        cpus => 6,
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => "kraken_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => 96,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => 'large',
        modules => $options->{modules},
        output => qq"${output_dir}/kraken_report.txt",);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($kraken);
}

sub Merge_Annotations {
   my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input_fsa', 'input_genbank', 'input_prokka_tsv'],
        input_abricate => '',
        input_interpro => '',
        input_classifier => '',
        input_phageterm => '',
        input_prodigal => '',
        input_trinotate => '',
        jprefix => '15',
        evalue => '1e-10',
        primary_key => 'locus_tag',
        product_columns => ['trinity_sprot_Top_BLASTX_hit', 'inter_Pfam', 'inter_TIGRFAM',],
        product_transmembrane => 'inter_TMHMM',
        product_signalp => 'trinity_SignalP',);

   my $output_name = basename($options->{input_fsa}, ('.fsa'));
   my $output_dir =  qq"outputs/$options->{jprefix}mergeannot";
   my $output_fsa = qq"${output_dir}/${output_name}.fsa";
   my $output_xlsx = qq"${output_dir}/${output_name}.xlsx";
   my $output_gbf = qq"${output_dir}/${output_name}.gbf";
   my $output_tbl = qq"${output_dir}/${output_name}.tbl";
   my $output_gbk = basename($output_gbf, ('.gbf'));
   $output_gbk = qq"${output_dir}/${output_gbk}";

   my $jstring = qq?
use Bio::Adventure::Annotation;
Bio::Adventure::Annotation::Merge_Annotations_Make_Gbk(\$h,
  input_abricate => '$options->{input_abricate}',
  input_classifier => '$options->{input_classifier}',
  input_fsa => '$options->{input_fsa}',
  input_genbank => '$options->{input_genbank}',
  input_phageterm => '$options->{input_phageterm}',
  input_interpro => '$options->{input_interpro}',
  input_prokka_tsv => '$options->{input_prokka_tsv}',
  input_prodigal => '$options->{input_prodigal}',
  input_trinotate => '$options->{input_trinotate}',
  jdepends => '$options->{jdepends}',
  jprefix => '$options->{jprefix}',
  jname => 'merge_annotations',
  );
?;
   my $merge_job = $class->Submit(
       input_abricate => $options->{input_abricate},
       input_classifier => $options->{input_classifier},
       input_fsa => $options->{input_fsa},
       input_genbank => $options->{input_genbank},
       input_interpro => $options->{input_interpro},
       input_phageterm => $options->{input_phageterm},
       input_prodigal => $options->{input_prodigal},
       input_prokka_tsv => $options->{input_prokka_tsv},
       input_trinotate => $options->{input_trinotate},
       jdepends => $options->{jdepends},
       jname => 'merge_annotations',
       jprefix => $options->{jprefix},
       jstring => $jstring,
       language => 'perl',
       library => $options->{library},
       output_dir => $output_dir,
       output_fsa => $output_fsa,
       output_gbf => $output_gbf,
       output_gbk =>  $output_gbk,
       output_tbl => $output_tbl,
       output_xlsx => $output_xlsx,
       primary_key => $options->{primary_key},
       shell => '/usr/bin/env perl',);
   $class->{language} = 'bash';
   $class->{shell} = '/usr/bin/env bash';
   return($merge_job);
}

sub Merge_Annotations_Make_Gbk {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input_fsa', 'input_genbank', 'input_prokka_tsv'],
        input_abricate => '',
        input_classifier => '',
        input_interpro => '',
        input_phageterm => '',
        input_prodigal => '',
        input_trinotate => '',
        evalue => '1e-10',
        jprefix => '15',
        primary_key => 'locus_tag',
        product_columns => ['trinity_sprot_Top_BLASTX_hit', 'inter_Pfam', 'inter_TIGRFAM',],
        product_transmembrane => 'inter_TMHMM',
        product_signalp => 'inter_signal',
        template_sbt => '/bio/reference/tbl2asn_template.sbt',
        );
    my $high_confidence = undef;
    my $likely_maximum = undef;
    my $possible_maximum = undef;
    if (defined($options->{evalue})) {
        $high_confidence = 1e-90;
        $likely_maximum = $options->{evalue} * $options->{evalue};
        $possible_maximum = $options->{evalue} * 100;
    }
    my $primary_key = $options->{primary_key};
    my $output_dir = qq"outputs/$options->{jprefix}mergeannot";
    make_path($output_dir);
    my $output_name = basename($options->{input_fsa}, ('.fsa'));
    my $output_fsa = qq"${output_dir}/${output_name}.fsa";
    my $output_xlsx = qq"${output_dir}/${output_name}.xlsx";
    my $output_gbf = qq"${output_dir}/${output_name}.gbf";
    my $output_tbl = qq"${output_dir}/${output_name}.tbl";
    my $log = qq"${output_dir}/${output_name}_runlog.txt";

    my $log_fh = FileHandle->new(">$log");
    print $log_fh "Merging annotations and writing new output files:
gbf: ${output_gbf}, tbl: ${output_tbl}, xlsx: ${output_xlsx}.\n";
    ## Here is a list of columns that will be created when merging these data sources.
    ## Make a hash of keys in the tsv file and names I want to use for adding notes/to the
    ## Bio::SeqIO annotations, keep in mind that the prokka ones are already there.
    my %note_keys = (
        ## Data types acquired from trinotate
        'trinity_gene_ontology_Pfam' => 'pfam_go',
        'trinity_Pfam' => 'pfam_hit',
        'trinity_phage_pep_BLASTX' => 'phage_blastx',
        'trinity_sprot_Top_BLASTX_hit' => 'swissprot_blastp',
        'trinity_gene_ontology_BLASTX' => 'blastx_go',
        'trinity_eggnog' => 'eggnog',
        'trinity_SignalP' => 'signalp',
        'trinity_RNAMMER' => 'rnammer',
        'trinity_terminase_BLASTX' => 'terminase_blast',
        'trinity_TmHMM' => 'tmhmm',
        ## The following hit types are from interpro
        'inter_MobiDBLite' => 'interpro_mobidb',
        'inter_SUPERFAMILY' => 'interpro_superfamily',
        'inter_Pfam' => 'interpro_pfam',
        'inter_Gene3D' => 'interpro_gene3d',
        'inter_SMART' => 'interpro_smart',
        'inter_PANTHER' => 'interpro_panther',
        'inter_CDD' => 'interpro_cdd',
        'inter_ProSitePatterns' => 'interpro_prositepatterns',
        'inter_Phobius' => 'interpro_phobius',
        'inter_PRINTS' => 'interpro_prints',
        'inter_signalp' => 'interpro_signalp',
        'inter_Coils' => 'interpro_coils',
        'inter_ProSiteProfiles' => 'interpro_prositeprofiles',
        'inter_TMHMM' => 'interpro_tmhmm',
        'inter_TIGRFAM' => 'interpro_tigrfam',
        ## And the abricate databases
        'abricate_argannot' => 'argannot',
        'abricate_card' => 'card',
        'abricate_ecoh' => 'ecoh',
        'abricate_ecoli_vf' => 'ecoli_vf',
        'abricate_megares' => 'megares',
        'abricate_ncbi' => 'ncbi_resistance',
        'abricate_plasmidfinder' => 'plasmidfinder',
        'abricate_resfinder' => 'resfinder',
        'abricate_vfdb' => 'vfdb',
        );

    ## A quick rundown of the various inputs and why they are being used:
    ## 1.   Data sources to create the merged table of annotations:
    ##  a.  input_prokka_tsv: The tsv output from prokka contains the base
    ##      annotation set used for everything else.
    ##  b.  input_trinotate: Outputs from trinotate are flexible and may
    ##      contain hits from arbitrary local databases.
    ##  c.  input_interpro: The interpro database provides a wide ranging
    ##      array of potentially interesting hits.
    ##  d.  input_abricate: Abricate provides a set of potentially helpful
    ##      resistance gene database queries.
    ##
    ## 2.   The merged table from #1 is used along with the following to
    ##      create a tbl file for tbl2asn's creation of the genbank output.
    ##  a.  input_genbank: This is the genbank output from prokka, it is
    ##      used to make writing the tbl input for tbl2asn easier because
    ##      it nicely sets the order of sequences.
    ##  b.  input_fsa: This is the entire sequence of the assembly, and is
    ##      used by tbl2asn to build the genbank file.

    ## This is #1 above and should contain all the interesting stuff from 1a-1d
    my $merged_data = {};

    print $log_fh "Reading prokka tsv data from $options->{input_prokka_tsv} to start.\n";
    ## 1a above, the starter set of annotations.
    unless (-r $options->{input_prokka_tsv}) {
        die("Unable to find the prokka tsv output, this is required.\n");
    }
    unless (-r $options->{input_genbank}) {
        die("Unable to find the prokka genbank output, this is required.\n");
    }

    ## The template sbt file was written with template toolkit variables which will need
    ## to be filled in, ideally with some information from the classifier.
    my $final_sbt = qq"${output_dir}/${output_name}.sbt";
    print $log_fh "Checking for ICTV classification data from $options->{input_classifier}.\n";
    my $taxonomy_information = {};
    ($merged_data, $taxonomy_information) = Merge_Classifier(
        input => $options->{input_classifier},
        output => $merged_data,
        primary_key => $options->{primary_key},
        template_sbt => $options->{template_sbt},
        template_out => $final_sbt,);
    print $log_fh "Wrote ${final_sbt} with variables filled in.\n";

    $merged_data = Merge_Prokka(
        primary_key => $options->{primary_key},
        output => $merged_data,
        input => $options->{input_prokka_tsv});

    if ($options->{input_trinotate} && -r $options->{input_trinotate}) {
        print $log_fh "Adding trinotate annotations from $options->{input_trinotate}.\n";
        ## 1b above, trinotate annotations.
        $merged_data = Merge_Trinotate(
            primary_key => $options->{primary_key},
            output => $merged_data,
            input => $options->{input_trinotate});
    } else {
        print $log_fh "Not including trinotate annotations.\n";
    }

    if ($options->{input_interpro} && -r $options->{input_interpro}) {
        print $log_fh "Adding interproscan annotations from $options->{input_interpro}.\n";
        ## 1c above, information from interproscan.
        ## The interpro output file is TSV, but it is super-annoying and not amendable to
        ## reading with Text::CSV
        $merged_data = Merge_Interpro(
            primary_key => $options->{primary_key},
            output => $merged_data,
            input => $options->{input_interpro});
    } else {
        print $log_fh "Not including interpro annotations.\n";
    }

    if ($options->{input_abricate} && -r $options->{input_abricate}) {
        print $log_fh "Adding abricate annotations from $options->{input_abricate}.\n";
        ## 1d above, resistance gene info provided by abricate.
        $merged_data = Merge_Abricate(
            primary_key => $options->{primary_key},
            output => $merged_data,
            input => $options->{input_abricate});
    } else {
        print $log_fh "Not including abricate annotations.\n";
    }

    ## print $log_fh "Reading prodigal genbank file from $options->{input_prodigal} to get RBS scores.\n";
    ## $merged_data = Merge_Prodigal(
    ##     primary_key => $options->{primary_key},
    ##     output => $merged_data,
    ##     input => $options->{input_prodigal});

    my %dtr_features = ();
    if ($options->{input_phageterm} && -r $options->{input_phageterm}) {
    ## Pull the direct-terminal-repeats from phageterm if they exist.
        %dtr_features = Bio::Adventure::Phage::Get_DTR(
            $class,
            input => $options->{input_phageterm});
    } else {
        print $log_fh "Not adding phageterm DTRs.\n";
    }

    ## This section uses the genbank input file to aid writing an output tbl file.

    ## The primary thing I need to recall when writing this is that tbl2asn, run with the arguments:
    ## tbl2asn -V -b -a r10k -l paired-ends -M n -N 1 -y 'modified from prokka' -Z testing.err -i EAb01.fsa -f EAb01.tbl
    ## Reads the following files: EAb01.fsa (the nucleotide assembly with a modified id for each contig),
    ## Eab01.tbl (the feature table, which prokka creates from the array of SeqIO objects).

    ## So, I need to create my own array of SeqIOs, presumably with the prokka genbank file.
    ## Then merge in the trinotate/etc annotations as notes and/or further inference tags.
    ## Finally, steal the tbl creator from prokka and invoke tbl2asn and prokka's gbf fixer.

    ## While I am reading the information from the existing genbank file, grab the sequences out
    ## and send them to the merged data so that I can write the starts/ends/strands/sequences.
    print $log_fh "Reading genbank file $options->{input_genbank} to acquire features.\n";
    my $in = FileHandle->new("less $options->{input_genbank} |");
    my $seqio = Bio::SeqIO->new(-format => 'genbank', -fh => $in);
    my $seq_count = 0;
    my $total_nt = 0;
    my $feature_count = 0;

    my @new_seq = ();
    my @new_features = ();
    while (my $seq = $seqio->next_seq) {
        my $seqid = $seq->id;
        push(@new_seq, $seqid);
        $seq_count++;
        my @feature_list = $seq->get_SeqFeatures();
        ## Check if we have the phageterm dtr feature, if so, put it at the beginning.
        foreach my $d (keys %dtr_features) {
            unshift(@feature_list, $dtr_features{$d});
        }

        for my $feat (@feature_list) {
            my $contig = $feat->seq_id(); ## The contig ID
            my $primary = $feat->primary_tag(); ## E.g. the type gene/cds/misc/etc
            my $annot = $feat->annotation();
            my $locus = 'failed_locustag';
            ## I want to mess with the CDS entries and leave the others alone.
            if ($feat->primary_tag eq 'CDS') {
                ## Get the information from our extra data source and add it as notes.
                my @loci = $feat->get_tag_values('locus_tag');
                $locus = $loci[0];
                my $new_info = $merged_data->{$locus};

                ## Pull out some sequence information to send back to $merged_data
                $merged_data->{$locus}->{start} = $feat->start;
                $merged_data->{$locus}->{end} = $feat->end;
                $merged_data->{$locus}->{strand} = $feat->strand;
                $merged_data->{$locus}->{cds} = $feat->seq->seq;
                $merged_data->{$locus}->{aaseq} = $feat->seq->translate->seq;
                ## Ok, back to our regularly scheduled programming of adding notes.
                ## to the Features from the annotation data.

                ## This is the likely place for mixing and matching where we want to add
                ## the annotation information, we could switch the add_tag_value() to inference,
                ## etc...
                my $note_string = '';
                ## my $match_number = 1;
                ROWDATA: foreach my $k (keys %{$new_info}) {
                    next ROWDATA if (!defined($new_info->{$k}));
                    next ROWDATA if ($new_info->{$k} eq '.');
                    if (exists($note_keys{$k})) {
                        my ($first_note, $stuff) = split(/\`/, $new_info->{$k});
                        my $note_string = qq"$note_keys{$k} hit: ${first_note}.";
                        $feat->add_tag_value('note', $note_string);

                        ## Make a simplified version of the note string as a match:
                        ## my @tag_arr = split(/\^/, $new_info->{$k});
                        ## my $tag_name = qq"match${match_number}";
                        ## my $tag_string = qq"$tag_arr[0] $tag_arr[3]";
                        ## print "ADDING TAG: $tag_name $tag_string\n";
                        ## $feat->add_tag_value($tag_name, $tag_string);
                        ## $match_number++;
                    }
                }

                ## Let us look through the array of product_columns
                ## ['trinity_sprot_Top_BLASTX_hit', 'inter_Pfam', 'inter_TIGRFAM',],
                ## Followed by the TmHMM and signalp colums
                ## product_transmembrane 'inter_TMHMM',
                ## product_signalp => 'inter_signal',
                my $tm_string = qq"Putative transmembrane domain containing protein";
                my $signal_string = "Putative signal peptide";
                my $product_string = undef;
                my $chosen_column = 0;
                my $tmhmm_column = $options->{product_transmembrane};
                if ($new_info->{$tmhmm_column}) {
                    my @tmhmm_data = split(/\^/, $new_info->{$tmhmm_column});
                    my $tm_region = $tmhmm_data[2];
                    $tm_region =~ s/^Q://g;
                    $product_string = $tm_string;
                }
                my $signalp_column = $options->{product_signalp};
                if ($new_info->{$signalp_column}) {
                    my @signal_data = split(/\^/, $new_info->{$signalp_column});
                    my $likelihood = $signal_data[2];
                    if ($likelihood > 0.9) {
                        $product_string = $signal_string;
                    }
                }
                my @product_columns = @{$options->{product_columns}};
                my $best_hit = 0;
                my $current_hit = 0;
                for my $col (@product_columns) {
                    ## Accession^Accession^Query,Hit^Identity^Evalue^RecName^Taxonomy
                    next unless(defined($new_info->{$col}));
                    my ($accession, $id, $query_hit, $identity, $evalue, $db_name, $taxonomy) = split(/\^/, $new_info->{$col});
                    next unless(defined($evalue));
                    $evalue =~ s/^(E|e)://g;
                    my $this_string;
                    if (defined($high_confidence)) {
                        if ($evalue <= $high_confidence) {
                            $this_string = qq"_High confidence ${db_name}";
                            $current_hit = 3;
                        } elsif ($evalue <= $likely_maximum) {
                            $this_string = qq"_Likely ${db_name}";
                            $current_hit = 2;
                        } elsif ($evalue <= $possible_maximum) {
                            $this_string = qq"_Potential ${db_name}";
                            $current_hit = 1;
                        }
                        ## If we are not attempting to discriminate among potential hits.
                    } else {
                        $this_string = qq"${db_name}";
                    }
                    if ($current_hit > $best_hit) {
                        $product_string = $this_string;
                        $chosen_column = $col;
                    }
                } ## Iterating over product columns looking for the most fun hit.

                ## See if $product_string has been filled
                if (defined($product_string)) {
                    $product_string =~ s/RecName: Full=//g;
                    my @current_values = $feat->remove_tag('product');
                    $product_string = qq"${product_string}";
                    my $new = $feat->add_tag_value('product', $product_string);
                    my $inf;
                    if ($product_string eq $signal_string) {
                        $inf = $feat->add_tag_value('inference', 'ab initio prediction:SignalP');
                    } elsif ($product_string eq $tm_string) {
                        $inf = $feat->add_tag_value('inference', 'ab initio prediction:TmHMM');
                    } else {
                        ## The %note_keys has more-readable versions of the column names, use one.
                        my $inference_string = $note_keys{$chosen_column};
                        $inference_string = qq"ab initio prediction:${inference_string}\n";
                        $inf = $feat->add_tag_value('inference', $inference_string);
                    }
                }

            } ## End looking for CDS entries
            push(@new_features, $feat);
        } ## End iterating over the feature list

        ## In theory, we now have a set of features with some new notes.
        ## So now let us steal the tbl writer from prokka and dump this new stuff...
        print $log_fh "Writing new tbl file to ${output_tbl}.\n";
        use Data::Dumper;
        print Dumper $taxonomy_information;
        my $tbl_written = Tbl_Writer(tbl_file => $output_tbl,
                                     taxonomy_information => $taxonomy_information,
                                     sequences => \@new_seq,
                                     features => \@new_features);
    } ## End Iterating over every sequence


    ## Remember that tbl2asn assumes the input files are all in the same directory.
    ## Before running tbl2asn, write a fresh copy of the fsa file containing the detected phage taxonomy.
    my $fsa_in = FileHandle->new("<$options->{input_fsa}");
    my $fsa_out = FileHandle->new(">$output_fsa");
    use Data::Dumper;
    print Dumper $taxonomy_information;
    while (my $line = <$fsa_in>) {
        ## chomp $line;  ## probably not needed for this.
        if ($line =~ /^\>/) {
            ## Replace the  [organism=Phage species] [strain=strain]
            ## in the prokka-derived fasta header, which looks like:
            ## >gnl|Prokka|test_1 [gcode=11] [organism=Phage species] [strain=strain]
            $line =~ s/(\[organism=.*\])/\[organism=Phage similar to $taxonomy_information->{taxon}\] \[strain=Similar to $taxonomy_information->{hit_accession}\]/g;
            ## $line =~ s/(\[organism=.*\])/\[organism=Phage similar to $taxonomy_information->{hit_description}\]/g;
            $line =~ s/\[strain=strain\]//g;
        }
        print $fsa_out $line;
    }
    $fsa_in->close();
    $fsa_out->close();

    ## This little section is stolen pretty much verbatim from prokka.
    my $tbl2asn_m_option = '-M n';
    if (scalar(@new_seq) > 10_000) {
        $tbl2asn_m_option = '-M b';
    }
    my $outdir = '.';
    my $accver = '1';
    my $EXE = 'slightly modified prokka';
    my $VERSION = '1.14.6';
    my $URL = 'https://github.com/tseemann/prokka';
    print $log_fh "tbl2asn uses the full assembly: $options->{input_fsa}, the tbl file ${output_tbl},
and modified template sbt file: ${final_sbt} to write a new gbf file: ${output_dir}/${output_name}.gbf.\n";
    my $tbl2asn_comment = qq"Annotated using $EXE $VERSION from $URL; Most similar taxon to this strain: $taxonomy_information->{taxon}";
    my $tbl_command = qq"tbl2asn -V b -c f -S F -a r10k -l paired-ends ${tbl2asn_m_option} -N ${accver} -y '${tbl2asn_comment}'" .
        " -Z ${output_dir}/${output_name}.err -t ${final_sbt} -i ${output_fsa} 1>${output_dir}/tbl2asn.log 2>${output_dir}/tbl2asn.err";
    print "TESTME: $tbl_command\n";
    print $log_fh "Running ${tbl_command}\n";
    my $tbl2asn_result = qx"${tbl_command}";
    my $sed_command = qq"sed 's/COORDINATES: profile/COORDINATES:profile/' ${output_dir}/${output_name}.gbf | sed 's/product=\"_/product=\"/g' >${output_dir}/${output_name}.gbk";
    my $sed_result = qx"${sed_command}";

    ## Now lets pull everything from the merged data and make a hopefully pretty xlsx file.
    print $log_fh "Writing final xlsx file of the annotations to $output_xlsx\n";

    my $written = Bio::Adventure::Annotation::Write_XLSX(
        $class,
        input => $merged_data,
        output => $output_xlsx,
        primary_key => $primary_key);

    ## Functions like this one should still return a job-like data structure so that I can
    ## putatively chain them, even though they are running primarily for their side-effects
    ## and the job that calls them is the one returning the fun information.
    my $output_gbk = basename($output_gbf, ('.gbf'));
    my $gbk_dir = dirname($output_gbf);
    $output_gbk = qq"${gbk_dir}/${output_gbk}.gbk";
    my $job = {
        output_dir => $output_dir,
        output_fsa => $output_fsa,
        output_gbf => $output_gbf,
        output_gbk =>  $output_gbk,
        output_tbl => $output_tbl,
        output_xlsx => $output_xlsx,};
    return($job);
}

sub Merge_Trinotate {
    my (%args) = @_;
    my $primary_key = $args{primary_key};
    my $merged_data = $args{output};
    my $trinotate_tsv = Text::CSV_XS::TSV->new({ binary => 1, });
    open(my $trinotate_fh, "<:encoding(utf8)", $args{input});
    my $trinotate_header = $trinotate_tsv->getline($trinotate_fh);
    $trinotate_tsv->column_names($trinotate_header);
    while (my $row = $trinotate_tsv->getline_hr($trinotate_fh)) {
        my $gene = $row->{'#gene_id'};
        foreach my $colname (@{$trinotate_header}) {
            next if ($colname eq '#gene_id');
            my $trin_col = qq"trinity_${colname}";
            $merged_data->{$gene}->{$trin_col} = $row->{$colname};
        }
    }
    close $trinotate_fh;
    ## Trinotate encodes as: CAF34166.1^CAF34166.1^Q:1135-1746,H:530-734^36.098%ID^E:3.33e-35^.^.
    ## TSP_BPKVM^TSP_BPKVM^Q:1000-1740,H:401-661^31.801%ID^E:7.61e-27^RecName: Full=Tail sheath protein;^Viruses; Duplodnaviria; Heunggongvirae; Uroviricota; Caudoviricetes; Caudovirales; Myoviridae; Tevenvirinae; Schizotequatrovirus
    ## Accession^Accession^Query,Hit^Identity^Evalue^RecName^Taxonomy
    return($merged_data);
}

sub Merge_Abricate {
    my (%args) = @_;
    my $primary_key = $args{primary_key};
    my $merged_data = $args{output};
    my $abricate_fh = FileHandle->new("<$args{input}");
    while (my $line = <$abricate_fh>) {
        chomp $line;
        my ($file, $sequence, $start, $end, $strand, $gene_id, $coverage, $covmap, $gaps, $pctcov, $pctid,
            $db, $accession, $resistance_name) = split(/\t/, $line);
        ## Look in the tsv_data hash for this annotation type for this gene.
        if (!defined($gene_id)) {
            next
        }
        if (defined($merged_data->{$gene_id}->{$db})) {
            ## If it already exists, for the moment just skip it
            next;
        } else {
            my $cell_data = qq"${accession}^${resistance_name}^Q:${start}-${end}^ID:${pctid}^Cov:${pctcov}";
            my $colname = qq"abricate_${db}";
            $merged_data->{$gene_id}->{$colname} = $cell_data;
        }
    }
    $abricate_fh->close();
    return($merged_data);
}

sub Merge_Classifier {
    my %args = @_;
    my $primary_key = $args{primary_key};
    my $template = $args{template_sbt};
    my $output = $args{template_out};
    my $merged_data = $args{output};
    my $default_values = {
        last_name => 'Margulieux',
        first_name => 'Katie',
        affiliation => 'Walter Reed Army Institute of Research',
        division => 'Bacterial Disease Branch/Wound Infections Department',
        city => 'Silver Spring',
        state => 'Maryland',
        country => 'United States of America',
        street => '503 Robert Grant Avenue',
        email => 'katie.r.margulieux.ctr@mail.mil',
        zip => '20910',
        bioproject => 'undefined bioproject',
        biosample => 'undefined biosample',
        second_email => 'abelew@umd.edu',
        ## The following entries will be changed by reading the input tsv file.
        taxon => 'Unknown taxonomy.',
        hit_length => 0,
        hit_accession => 'No matching accession.',
        hit_description => 'Unclassified phage genome.',
        hit_bit => 'No tblastx bit score found.',
        hit_sig => 'No significance score found.',
        hit_score => 'No tblastx score found.',
        user_comment => 'Phage genome with no detected similar taxonomy.',
    };

    my $input_tsv = Text::CSV_XS::TSV->new({ binary => 1, });
    print "TESTME: input file: $args{input}\n";
    open(my $classifier_fh, "<:encoding(utf8)", $args{input});
    my $input_header = $input_tsv->getline($classifier_fh);
    $input_tsv->column_names($input_header);
    my $rows_read = 0;
  ROWS: while (my $row = $input_tsv->getline_hr($classifier_fh)) {
      ## Just read the first row.
      if ($rows_read > 0) {
          last ROWS;
      }
      $rows_read++;

      my @wanted_columns = ('taxon', 'hit_length', 'hit_accession', 'hit_description',
                            'hit_bit', 'hit_sig', 'hit_score');
      foreach my $colname (@wanted_columns) {
          ## Check that we already filled this data point
          $default_values->{$colname} = $row->{$colname};
      }
  }
    close $classifier_fh;
    if ($default_values->{taxon} ne 'Unknown taxonomy.') {
        my $comment_string =  qq"tblastx derived taxonomy: $default_values->{taxon}, description: $default_values->{hit_description}, accession: $default_values->{hit_accession}, significance: $default_values->{hit_bit}, hit length: $default_values->{hit_length}, hit score: $default_values->{hit_score}";
        $default_values->{user_comment} = $comment_string;
    }
    my $tt = Template->new({
        ABSOLUTE => 1,});
    my $written = $tt->process($template, $default_values, $output) or print $tt->error(), "\n";
    ## Now write out template sbt file to the output directory with the various values filled in.
    return($merged_data, $default_values);
}

sub Merge_Prodigal {
    my (%args) = @_;
    my $primary_key = $args{primary_key};
    my $merged_data = $args{output};
    my $in = FileHandle->new("less $args{input} |");
    while (my $line = <$in>) {
        chomp $line;
        my ($contig, $prod, $cds, $start, $end, $score, $strand, $phase, $tags) = split(/\t/, $line);
        my %tag_hash = ();
        my @tag_pieces = split(/;/, $tags);
        for my $t (@tag_pieces) {
            my ($name, $value) = split(/=/, $t);
            $tag_hash{$name} = $value;
        }
    }  ## End iterating over the gff file.
    $in->close();
    return($merged_data);
}

sub Merge_Interpro {
    my (%args) = @_;
    my $primary_key = $args{primary_key};
    my $merged_data = $args{output};
    my $inter_fh = FileHandle->new("<$args{input}");
    while (my $line = <$inter_fh>) {
        chomp $line;
        my ($gene_id, $md5, $gene_length, $source, $hit_id, $hit_name,
            $hit_start, $hit_end, $hit_score, $hit_boolean, $hit_date,
            $interpro_id, $interpro_name) = split(/\t/, $line);
        if ($source =~ /^SignalP/) {
            $source = 'signalp';
        }
        $source = qq"inter_${source}";
        ## Look in the tsv_data hash for this annotation type for this gene.
        if (defined($merged_data->{$gene_id}->{$source})) {
            ## If it already exists, for the moment just skip it
            next;
        } else {
            ## Trinotate does:
            ## Accession^Accession^Query,Hit^Identity^Evalue^RecName^Taxonomy
            ## So let us do the same.
            my $cell_data = qq"${hit_id}^${interpro_id}^Q:${hit_start}-${hit_end}^^E:${hit_score}^${interpro_name}^.";
            $merged_data->{$gene_id}->{$source} = $cell_data;
        }
    }
    $inter_fh->close();
    return($merged_data);
}

sub Merge_Prokka {
    my (%args) = @_;
    my $primary_key = $args{primary_key};
    my $merged_data = $args{output};
    my $prokka_tsv = Text::CSV_XS::TSV->new({ binary => 1, });
    open(my $prokka_fh, "<:encoding(utf8)", $args{input});
    my $prokka_header = $prokka_tsv->getline($prokka_fh);
    $prokka_tsv->column_names($prokka_header);
    ROWS: while (my $row = $prokka_tsv->getline_hr($prokka_fh)) {
        next ROWS if ($row->{ftype} eq 'gene');
        my $key = $row->{$primary_key};
        $merged_data->{$key}->{$primary_key} = $key;
        foreach my $colname (@{$prokka_header}) {
            ## Check that we already filled this data point
            if ($merged_data->{$key}->{$colname}) {
                ## There is something here.
            } else {
                if ($row->{$colname}) {
                    $merged_data->{$key}->{$colname} = $row->{$colname};
                }
            }
        }  ## Run over the columns
    }
    close $prokka_fh;
    return($merged_data);
}

sub Write_XLSX {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'output'],
        primary_key => 'locus_tag',);

    my $tsv_data = $options->{input};
    my $output = $options->{output};
    my $primary_key = $options->{primary_key};

    ## Write the combined annotations to an xlsx output file.
    ## Write a quick xlsx file of what we have merged together.
    my @table;
    ## I am tired, so I will just do a two-pass over the keys to get all of their names.
    ##my @column_ids = ($primary_key);
  ##COLHUNT: foreach my $rowname (keys %{$tsv_data}) {
  ##    $tsv_data->{$rowname}->{$primary_key} = $rowname;
  ##    next COLHUNT unless (defined($tsv_data->{$rowname}));
  ##    my %internal = %{$tsv_data->{$rowname}};
  ##    foreach my $k (keys %internal) {
  ##        next if ($k eq $primary_key);
  ##        push(@column_ids, $k) unless any { $_ eq $k } @column_ids
  ##    }
    ##}
    my @column_ids = (
        'locus_tag',
        'product',
        'start',
        'end',
        'strand',
        'length_bp',
        'COG',
        'EC_number',
        'trinity_phage_pep_BLASTX',
        'trinity_terminase_BLASTX',
        'trinity_sprot_Top_BLASTX_hit',
        'inter_Pfam',
        'trinity_Pfam',
        'inter_TIGRFAM',
        'inter_CDD',
        'inter_TMHMM',
        'trinity_TmHMM',
        'inter_signalp',
        'trinity_SignalP',
        'inter_Coils',
        'inter_ProSitePatterns',
        'trinity_RNAMMER',
        'trinity_eggnog',
        'trinity_Kegg',
        'trinity_gene_ontology_BLASTX',
        'trinity_gene_ontology_Pfam',
        'abricate_argannot',
        'abricate_card',
        'abricate_ecoh',
        'abricate_ecoli_vf',
        'abricate_megares',
        'abricate_ncbi',
        'abricate_plasmidfinder',
        'abricate_resfinder',
        'abricate_vfdb',
        'aaseq',
        'cds',);

    foreach my $locus_tag (sort keys %{$tsv_data}) {
        my @row;
        for my $col (@column_ids) {
            my $cell;
            if (defined($tsv_data->{$locus_tag}->{$col})) {
                $cell = $tsv_data->{$locus_tag}->{$col};
            } else {
                $cell = '';
            }
            push(@row, $cell);
        }
        push(@table, \@row);
    }

    my $xlsx_data = Data::Table->new(\@table, \@column_ids, 0);
    ## tables2xlsx ($fileName, $tables, $names, $colors, $portrait, $columnHeaders)
    my $written = tables2xlsx($output, [$xlsx_data], ["Features"],
                              [["white","lightblue","blue"]], [1], [1]);
    my $out = basename($output, ('.xlsx'));
    my $dir = dirname($output);
    my $outtsv = FileHandle->new(">${dir}/${out}.tsv");
    my $tsv_string = $xlsx_data->tsv;
    print $outtsv $tsv_string;
    $outtsv->close();
}


sub Tbl_Writer {
    my %args = @_;
    my $file = $args{tbl_file};
    my $taxonomy_information = $args{taxonomy_information};
    my @seq = @{$args{sequences}};
    my @features = @{$args{features}};
    my $tbl_fh = FileHandle->new(">${file}");
    my $hypothetical_string = 'hypothetical protein';
    ## Ok, so it turns out that Torsten uses %seq, @seq, and $seq separately in prokka.
    ## That is more than a little evil.
    ## @seq is used in prokka as an array of IDs to keep the contig order correct.
    ## %seq is a hash of hashes where the first key is the ID and the second is a series
    ##   of data containers for sequences, features, etc.
    ## $seq is a series of Bio::SeqIO objects
    ## While I definitely like the idea of filling a hash with the various feature information,
    ## naming it %seq when you already have @seq and $seq is a bit insane and confusing.

    for my $sid (@seq) {
        print $tbl_fh ">Feature gnl|Prokka|${sid}\n";
        for my $f (@features) {
            if ($f->primary_tag eq 'source') {
                $f->add_tag_value('organism', "Phage similar to $taxonomy_information->{taxon}.");
                $f->add_tag_value('strain', "Similar accession: $taxonomy_information->{hit_accession}.");
                $f->seq_id("Phage species similar to $taxonomy_information->{taxon}.");
                use Data::Dumper;
                print Dumper $f;
            }
            if ($f->primary_tag eq 'CDS' and not $f->has_tag('product')) {
                $f->add_tag_value('product', $hypothetical_string);
            }
            if (my $name = TAG($f, 'gene')) {
                $f->add_tag_value('Name', $name);
            }
            # Make sure we have valid frames/phases (GFF column 8)
            $f->frame( $f->primary_tag eq 'CDS' ? 0 : '.' );
            my ($L, $R) = ($f->strand >= 0) ? ($f->start, $f->end) : ($f->end, $f->start);
            print $tbl_fh "$L\t$R\t", $f->primary_tag, "\n";
            for my $tag ($f->get_all_tags) {
                # remove GFF specific tags (start with uppercase letter)
                next if $tag =~ m/^[A-Z]/ and $tag !~ m/EC_number/i;
                for my $value ($f->get_tag_values($tag)) {
                    if (!defined($value)) {
                        print $tbl_fh "\t\t\t${tag}\t\n";
                    } else {
                        print $tbl_fh "\t\t\t${tag}\t${value}\n";
                    }
                }
            }
        }
    }
}


sub TAG {
  my($f, $tag) = @_;
  # very important to "return undef" here and not just "return"
  # otherwise it becomes non-existent/collapsed in a list
  return undef unless $f->has_tag($tag);
  return ($f->get_tag_values($tag))[0];
}



=head2 C<Prodigal>

Invoke prodigal on an assembly to search for ORFs.  This will by
default look for an existing training file provided by
Train_Prodigal()'.

=cut
sub Prodigal {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        species => undef,
        gcode => '11',
        output_dir => undef,
        jprefix => '17',
        modules => ['prodigal'],
        prodigal_outname => undef,);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('prodigal');
    die("Could not find prodigal in your PATH.") unless($check);

    my $inputs = $class->Get_Paths($options->{input});
    my $job_name = $class->Get_Job_Name();
    $job_name = basename($job_name, ('.fsa'));
    my $train_string = '';
    my $library_file;
    if ($options->{species}) {
        $library_file = qq"$options->{libdir}/hmm/$options->{species}_gc$options->{gcode}.training";
        if (!-r $library_file) {
            print "Could not find the training file for this species.\n";
            print "Sleeping for a moment so that you can do something.\n";
            sleep(5);
            $train_string = '';
        }
        $train_string = qq" -t ${library_file} ";
    }

    my $in_name = basename($options->{input}, ('.fasta'));
    my $output_dir;
    if (defined($options->{output_dir})) {
        if ($options->{output_dir}) {
            $output_dir = $options->{output_dir};
        }
    } else {
        $output_dir = qq"outputs/$options->{jprefix}prodigal_${in_name}";
    }

    my ($cds_file, $translated_file, $scores_file, $gff_file, $gbk_file);
    if ($options->{prodigal_outname}) {
        $cds_file = qq"${output_dir}/$options->{prodigal_outname}_cds.fasta";
        $translated_file = qq"${output_dir}/$options->{prodigal_outname}_translated.fasta";
        $scores_file = qq"${output_dir}/$options->{prodigal_outname}_scores.txt";
        $gff_file = qq"${output_dir}/$options->{prodigal_outname}.gff";
        $gbk_file = qq"${output_dir}/$options->{prodigal_outname}.gb";
    } else {
        $cds_file = qq"${output_dir}/predicted_cds.fasta";
        $translated_file = qq"${output_dir}/predicted_translated.fasta";
        $scores_file = qq"${output_dir}/predicted_scores.txt";
        $gff_file = qq"${output_dir}/predicted_cds.gff";
        $gbk_file = qq"${output_dir}/predicted_cds.gb";
    }

    my $comment = qq!## This is a script to run prodigal.
!;
    my $jstring = qq!mkdir -p ${output_dir}
prodigal ${train_string} \\
  -i $options->{input} \\
  -a ${translated_file} \\
  -d ${cds_file} \\
  -s ${scores_file} \\
  -f gff -o ${gff_file} \\
  2>${output_dir}/prodigal_gff.err \\
  1>${output_dir}/prodigal_gff.out
prodigal ${train_string} \\
  -i $options->{input} \\
  -f gbk -o ${gbk_file} \\
  2>${output_dir}/prodigal_gbk.err \\
  1>${output_dir}/prodigal_gbk.out
!;
    my $prodigal = $class->Submit(
        cpus => 1,
        comment => $comment,
        jdepends => $options->{jdepends},
        jmem => 24,
        jname => "prodigal_${job_name}",
        jprefix => $options->{jprefix},
        jqueue => "workstation",
        jstring => $jstring,
        modules => $options->{modules},
        output => $gbk_file,
        output_cds => $cds_file,
        output_gff => $gff_file,
        output_scores => $scores_file,
        output_translated => $translated_file,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        training_input => $library_file,);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($prodigal);
}

=head2 C<Train_Prodigal>

Some assemblies I have been performing are on sets of sequence which
are too small for prodigal to train itself sufficiently; so this
function was written to provide an opportunity for one to collate a
larger sequence database for training.

=cut
sub Train_Prodigal {
    my ($class, %args) = @_;
    my $check = which('prodigal');
    die("Could not find prodigal in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        gcode => '11',
        modules => ['prodigal'],);
    my $job_name = $class->Get_Job_Name();
    my $kingdom_string = '';
    my $output_dir = qq"$options->{libdir}/hmm";
    my $output = qq"${output_dir}/$options->{species}_gc$options->{gcode}.training";
    my $comment = qq!## This is a script to train prodigal.
!;
    my $jstring = qq!mkdir -p ${output_dir}
prodigal -i $options->{input} \\
  -t ${output} \\
  2>${output_dir}/prodigal_training.stderr \\
  1>${output_dir}/prodigal_training.stdout
!;
    my $prodigal = $class->Submit(
        cpus => 1,
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => "prodigal_training_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => 24,
        modules => $options->{modules},
        output => $output,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => "workstation",);
    return($prodigal);
}

=head2 C<Prokka>

Prokka is an automagic annotation tool for bacterial assemblies.  It
seems useful for other relatively small genomes.

=cut
sub Prokka {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        arbitrary => '',
        coverage => 30,
        evalue => '1e-05',
        gcode => '11',
        genus => 'phage',
        jprefix => '19',
        kingdom => 'bacteria',
        locus_tag => 'unknownphage',
        modules => ['prokka'],
        species => 'virus',);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('prokka');
    die("Could not find prokka in your PATH.") unless($check);

    my $job_name = $class->Get_Job_Name();
    my $kingdom_string = '';
    if ($options->{kingdom} ne '') {
        $kingdom_string = qq" --kingdom $options->{kingdom} --gcode $options->{gcode} ";
    }

    my $training_file = qq"$options->{libdir}/hmm/$options->{species}_gc$options->{gcode}.training";
    my $cwd_name = basename(cwd());
    my $input_name = basename(dirname($options->{input}));
    my $locus_tag;
    if ($options->{locus_tag}) {
        $locus_tag = $options->{locus_tag};
    } else {
        $locus_tag = basename($input_name, ('.fasta'));
    }
    my $output_dir = qq"outputs/$options->{jprefix}prokka_${input_name}";
    my $comment = qq!## This is a script to run prokka.
!;

##    prokka --addgenes -rfam --force --locustag EAb03 --genus Phage --species virus --cdsrnaolap --usegenus --kingdom Bacteria --gcode 11 --prodigaltf /home/trey/libraries/hmm/abaumannii_phage_samples_gc11.training --outdir outputs/19prokka_15termreorder_13unicycler --prefix EAb03 outputs/15termreorder_13unicycler/final_assembly_reordered.fasta --evalue '1e-05' --coverage 30 2>&1 | less

    my $jstring = qq!mkdir -p ${output_dir}
prokka --addgenes --rfam --force ${kingdom_string} \\
  --locustag ${locus_tag} --genus $options->{genus} \\
  --compliant --cdsrnaolap --usegenus \\
  --prodigaltf ${training_file} \\
  --outdir ${output_dir} \\
  --prefix ${cwd_name} \\
  $options->{input} \\
  2>${output_dir}/prokka.stderr \\
  1>${output_dir}/prokka.stdout
!;

    my $error_file =  qq"${output_dir}/${cwd_name}.err";
    my $peptide_file = qq"${output_dir}/${cwd_name}.faa";
    my $cds_file = qq"${output_dir}/${cwd_name}.ffn";
    my $assembly_copy = qq"${output_dir}/${cwd_name}.fna";
    my $assembly_renamed_contigs = qq"${output_dir}/${cwd_name}.fsa";
    my $genbank_file = qq"${output_dir}/${cwd_name}.gbk";
    my $gff_file = qq"${output_dir}/${cwd_name}.gff";
    my $sqn_file = qq"${output_dir}/${cwd_name}.sqn";
    my $tbl_file = qq"${output_dir}/${cwd_name}.tbl";
    my $tsv_file = qq"${output_dir}/${cwd_name}.tsv";

    my $prokka = $class->Submit(
        cpus => 6,
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => "prokka_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => 24,
        modules => $options->{modules},
        output => $cds_file,
        output_error => $error_file,
        output_peptide => $peptide_file,
        output_cds => $cds_file,
        output_assembly => $assembly_copy,
        output_assemblyrenamed => $assembly_renamed_contigs,
        output_fsa => $assembly_renamed_contigs,
        output_genbank => $genbank_file,
        output_gff => $gff_file,
        output_sqn => $sqn_file,
        output_tbl => $tbl_file,
        output_tsv => $tsv_file,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => "workstation",);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($prokka);
}


=head2 C<tRNAScan>

Alternative to aragorn.  Search for tRNAs!

=cut
sub tRNAScan {
    my ($class, %args) = @_;
    my $check = which('trnascan');
    die("Could not find trnascan in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        modules => ['trnascan'],
        species => undef,
        arbitrary => ' -G ',);
    my $trnascan_args = $options->{arbitrary};
    my $job_name = $class->Get_Job_Name();
    my $output_dir = qq"outputs/trnascan";
    my $species_string = qq"";
    my $comment = qq!## This is a script to run trnascan.
!;
    my $jstring = qq!mkdir -p ${output_dir} && \\
  tRNAscan-SE $options->{arbitrary} \\
    -o ${output_dir}/trnascan.txt \\
    $options->{input}
!;

    my $trnascan = $class->Submit(
        cpus => 6,
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => "trnascan_${job_name}",
        jprefix => "64",
        jstring => $jstring,
        jmem => 24,
        modules => $options->{modules},
        output => qq"${output_dir}/trnascan.txt",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => "workstation",);
    return($trnascan);
}

=head2 C<Extract_Annotations>

Used by Extract_Trinotate to extract annotations from the strangely encoded trinotate outputs.

=over

=item I<evalue> Evalue cutoff for trinotate output.

=item I<identity> Minimum percent identity cutoff for trinotate output.

=item I<ids> IDs to extract.

=item I<fh> Filehandle to parse.

=back

=cut
sub Extract_Annotations {
    my ($class, $datum, %args) = @_;
    my $min_identity = $args{identity};
    my $eval_max = $args{evalue};
    my $ids = $args{ids};
    my $fh = $args{fh};
    my $options = $class->Get_Vars();

    my $id = $datum->{prot_id};
    $id = $datum->{transcript_id} if ($id eq '.');
    if (defined($options->{species})) {
        my $species = $options->{species};
        $species =~ s/\s+/_/g;
        $id =~ s/TRINITY/${species}/g;
    }
    $id =~ s/\:\:/; /g;

    if (defined($ids->{$id})) {
        ## Then this record has been processed, just increment it and move on.
        $ids->{$id}++;
    } elsif (!defined($datum->{swissprot} && !defined($datum->{rnammer}))) {
        $ids->{$id} = 1;
    } elsif (defined($datum->{rnammer})) {
        my $seq = $datum->{sequence};
        $seq = join("\n", ($seq =~ m/.{1,80}/g));
        $ids->{$id} = 1;
        my $header_string = qq"${id}; ${id}; undef; undef; $datum->{rnammer}->[0]->{name}, $datum->{rnammer}->[0]->{region}; undef; undef";
        print $fh qq">${header_string}
${seq}
";
    } elsif ($datum->{swissprot}->[0]->{identity} >= $min_identity &&
                 $datum->{swissprot}->[0]->{e_value} <= $eval_max) {
        $ids->{$id} = 1;
        my $seq = $datum->{sequence};
        $seq = join("\n", ($seq =~ m/.{1,80}/g));
        my $header_string = qq"${id}; $datum->{swissprot}->[0]->{fullname}; $datum->{swissprot}->[0]->{species}; $datum->{swissprot}->[0]->{e_value}; $datum->{swissprot}->[0]->{identity}";
        print $fh qq">${header_string}
${seq}
";
    }
    return($ids);
}

=head2 C<Extract_Trinotate>

The trinotate output format is a bit... unwieldy.  This seeks to parse out the
useful information from it.

=over

=item I<input> * Input csv file from trinotate.

=item I<output> (interesting.fasta) Output for the parsed csv.

=item I<evalue> (1e-10) Evalue cutoff for trinotate output.

=item I<identity> (70) Minimum percent identity cutoff for trinotate output.

=back

=head3 C<Invocation>

> cyoa --task assembly --method extract --input trinotate_output.csv

=cut
sub Extract_Trinotate {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jname => "trin_rsem",
        output => 'interesting.fasta',
        evalue => 1e-10,
        identity => 70,);
    my $job_name = $class->Get_Job_Name();
    my $trinity_out_dir = qq"outputs/trinity_${job_name}";
    my $input = FileHandle->new("<$options->{input}");
    my $parser = Parse::CSV->new(
        handle => $input,
        sep_char => "\t",
        names => 1,
    );

    my $count = 0;
    my $all_annotations = [];
    my $out = FileHandle->new(">$options->{output}");
    while (my $object = $parser->fetch) {
        my $filled_in = Read_Write_Annotation($class, $object,
                                              db => $all_annotations,
                                              count => $count,
                                              evalue => $options->{evalue},
                                              identity => $options->{identity},
                                              fh => $out);
        $count = $count + 1;
    }
    $out->close();
    $input->close();
    return($count);
}


=head2 C<Read_Write_Annotation>

Called by Extract_Trinotate() to help write out the trinotate csv information
into an easier-to-read format.

=cut
sub Read_Write_Annotation {
    my ($class, $object, %args) = @_;
    my $annotations = $args{db};
    my $element_number = $args{count};
    my $fh = $args{fh};
    ## $object is a hash reference containing the various csv fields.
    ## Extract them by name and dump them to the array of material
    my $element = {};
    my $filled_in = 0;
    my $ids = {};

    foreach my $key (keys %{$object}) {
        my @result_list = ();
        my @field_elements = ();
        my $field_text = $object->{$key};
        if ($field_text =~ /\`/) {
            @field_elements = split(/\`/, $field_text);
        } else {
            $field_elements[0] = $field_text;
        }

        for my $t (0..$#field_elements) {
            my $ret = {};
            if ($key eq 'sprot_Top_BLASTX_hit') {
                if ($field_text eq '.') {
                    ## Null case: if a . then just drop out.
                    $element->{swissprot} = undef;
                } else {
                    ## The swissprot elements look something like:
                    ## Name        ## ARIA_ARATH^
                    ## Name again? ## ARIA_ARATH^
                    ## Query range ## Q:7-573,H:513-701^
                    ## %identity   ## 75.66%ID^
                    ## E-value     ## E:3e-102^
                    ## Full name   ## RecName: Full=ARM REPEAT PROTEIN INTERACTING WITH ABF2;^
                    ## Ontology    ## Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliophyta; eudicotyledons; Gunneridae; Pentapetalae; rosids; malvids; Brassicales; Brassicaceae; Camelineae; Arabidopsis',
                    ## Start by pulling apart the data using the chosen (somewhat odd) separator '^'.
                    my ($name, $name_again, $query_coords, $identity, $e_value, $fullname, $ontology) = split(/\^/, $field_text);
                    ## Clean up the fields a little
                    $identity =~ s/\%ID//g;
                    $e_value =~ s/E://g;
                    $fullname =~ s/RecName: //g;
                    $fullname =~ s/Full=//g;
                    $fullname =~ s/\;//g;
                    my ($query, $h) = split(/\,/, $query_coords);
                    $query =~ s/Q\://g;
                    $h =~ s/H\://g;
                    $name =~ s/Full=//g;
                    my @ontology_list = split(/; /, $ontology);
                    my $species = qq"$ontology_list[$#ontology_list - 1] $ontology_list[$#ontology_list]";
                    ## Put them back into the $ret
                    $ret->{name} = $name;
                    $ret->{coords} = $query_coords;
                    $ret->{identity} = $identity;
                    $ret->{e_value} = $e_value;
                    $ret->{fullname} = $fullname;
                    $ret->{species} = $species;
                    $filled_in = $filled_in + 6;
                    ## And refill field_elements with this information.
                    $field_elements[$t] = $ret;
                    ## Finally, fill in the swissprot entry with this information.
                    $element->{swissprot} = \@field_elements;
                }
            } elsif ($key eq 'gene_ontology_pfam') {
                if ($field_text eq '.') {
                    $element->{pfam_go} = undef;
                } else {
                    my ($id, $ontology, $name) = split(/\^/, $field_text);
                    $ret->{id} = $id;
                    $ret->{ontology} = $ontology;
                    $ret->{name} = $name;
                    $filled_in = $filled_in + 3;
                    $field_elements[$t] = $ret;
                    $element->{pfam_go} = \@field_elements;
                }
            } elsif ($key eq 'gene_ontology_blast') {
                if ($field_text eq '.') {
                    $element->{blast_go} = undef;
                } else {
                    my ($id, $ontology, $name) = split(/\^/, $field_text);
                    $ret->{id} = $id;
                    $ret->{ontology} = $ontology;
                    $ret->{name} = $name;
                    $filled_in = $filled_in + 3;
                    $field_elements[$t] = $ret;
                    $element->{blast_go} = \@field_elements;
                }
            } elsif ($key eq 'RNAMMER') {
                if ($field_text eq '.') {
                    $element->{rnammer} = undef;
                } else {
                    my ($name, $region) = split(/\^/, $field_text);
                    $ret->{name} = $name;
                    $ret->{region} = $region;
                    $filled_in = $filled_in + 2;
                    $field_elements[$t] = $ret;
                    $element->{rnammer} = \@field_elements;
                }
            } elsif ($key eq 'eggnog') {
                if ($field_text eq '.') {
                    $element->{eggnog} = undef;
                } else {
                    my ($id, $name) = split(/\^/, $field_text);
                    $ret->{id} = $id;
                    $ret->{name} = $name;
                    $filled_in = $filled_in + 2;
                    $field_elements[$t] = $ret;
                    $element->{eggnog} = \@field_elements;
                }
            } elsif ($key eq 'transcript_id') {
                $filled_in = $filled_in + 1;
                $element->{transcript_id} = $field_text;
            } elsif ($key eq 'Kegg') {
                if ($field_text eq '.') {
                    $element->{kegg} = undef;
                } else {
                    my ($id, $name) = split(/\^/, $field_text);
                    $ret->{id} = $id;
                    $ret->{name} = $name;
                    $filled_in = $filled_in + 2;
                    $field_elements[$t] = $ret;
                    $element->{kegg} = \@field_elements;
                }
            } elsif ($key eq 'prot_coords') {
                if ($field_text eq '.') {
                    $element->{prot_coords} = undef;
                } else {
                    $filled_in = $filled_in + 1;
                    $element->{prot_coords} = $field_text;
                }
            } elsif ($key eq 'prot_id') {
                $element->{prot_id} = $field_text;
            } elsif ($key eq 'transcript') {
                $element->{sequence} = $field_text;
            } elsif ($key eq 'TmHMM') {
                ## An example field
                ## ExpAA=124.31^PredHel=6^Topology=i59-78o83-105i126-148o152-174i187-209o229-251i
                if ($field_text eq '.') {
                    $element->{tmhmm} = undef;
                } else {
                    my ($expaa, $pred_helixes, $topology) = split(/\^/, $field_text);
                    $expaa =~ s/ExpAA=//g;
                    $pred_helixes =~ s/PredHel=//g;
                    $topology =~ s/Topology=//g;
                    $ret->{expaa} = $expaa;
                    $ret->{pred_helixes} = $pred_helixes;
                    $ret->{topology} = $topology;
                    $filled_in = $filled_in + 3;
                    $field_elements[$t] = $ret;
                    $element->{tmhmm} = \@field_elements;
                }
            } elsif ($key eq 'SignalP') {
                ## sigP:1^18^0.613^YES
                if ($field_text eq '.') {
                    $element->{signalp} = undef;
                } else {
                    my ($first, $second, $third, $boolean) = split(/\^/, $field_text);
                    $first =~ s/sigP\://g;
                    $ret->{first} = $first;
                    $ret->{second} = $second;
                    $ret->{third} = $third;
                    $ret->{boolean} = $boolean;
                    $filled_in = $filled_in + 4;
                    $field_elements[$t] = $ret;
                    $element->{tmhmm} = \@field_elements;
                }
            }
            ## Set the nth annotation to its annotation.
            $annotations->[$element_number] = $element;
        } ## End iterating over every element in a multi-element field
    } ## End the foreach every key in the object
    $ids = Extract_Annotations($class, $element,
                               ids => $ids,
                               fh => $fh,
                               evalue => $args{evalue},
                               identity => $args{identity},
                           );
    return($filled_in);
}

=head2 C<Transdecoder>

$hpgl->Transdecoder() submits a trinity denovo sequence assembly and runs its
default post-processing tools.

=over

=item I<input> * Output from trinity for post processing.

=back

=head3 C<Invocation>

> cyoa --task assembly --method transdecoder --input trinity.fasta

=cut
sub Transdecoder {
    my ($class, %args) = @_;
    my $check = which('TransDecoder.LongOrfs');
    die("Could not find transdecoder in your PATH.") unless($check);
    my $transdecoder_exe_dir = dirname($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        modules => ['transdecoder'],);
    my $transdecoder_input = File::Spec->rel2abs($options->{input});
    my $job_name = $class->Get_Job_Name();
    my $output_dir = qq"outputs/trinity_${job_name}";
    my $comment = qq!## This is a transdecoder submission script
!;
    my $jstring = qq!mkdir -p ${output_dir} && cd ${output_dir} && \\
  TransDecoder.LongOrfs -t ${transdecoder_input} \\
    2>${output_dir}/transdecoder_longorfs_${job_name}.err \\
    1>${output_dir}/transdecoder_longorfs_${job_name}.out
TransDecoder.Predict -t ${transdecoder_input} \\
  2>${output_dir}/transdecoder_predict_${job_name}.err \\
  1>${output_dir}/transdecoder_predict_${job_name}.out
${transdecoder_exe_dir}/util/cdna_alignment_orf_to_genome_orf.pl \\
  ${output_dir}/transcripts.fasta.transdecoder.gff3 \\
  ${output_dir}/transcripts.gff3 \\
  ${transdecoder_input} \\
  2>${output_dir}/cdna_alignment_${job_name}.err \\
  1>${output_dir}/transcripts.fasta.transdecoder.genome.gff3
!;
    my $transdecoder = $class->Submit(
        comment => $comment,
        jname => "transdecoder_${job_name}",
        jprefix => "47",
        jstring => $jstring,
        modules => $options->{modules},
        output => qq"${output_dir}/transcripts.fasta.transdecoder.genome.gff3",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => "workstation",);
    return($transdecoder);
}

=head2 C<Trinotate>

Submit a trinity denovo sequence assembly to trinotate.

=over

=item I<input> * Input fasta from trinity.

=back

=head3 C<Invocation>

> cyoa --task assembly --method trinotate --input trinity.fasta

=cut
sub Trinotate {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jprefix => '20',
        modules => ['divsufsort', 'transdecoder', 'blast', 'blastdb', 'signalp', 'hmmer',
                    'tmhmm', 'rnammer', 'trinotate', ],
        trinotate => 'autoTrinotate.pl',
        config => 'conf.txt',
        );

    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('Trinotate');
    die("Could not find Trinotate in your PATH:
 $ENV{PATH}.") unless($check);

    my $job_name = $class->Get_Job_Name();
    my $cwd_name = basename(cwd());
    my $trinotate_exe_dir = dirname($check);
    ## Once again, abs_path only works on stuff which already exists.
    ## So create the output directory, and use that.
    my $input_paths = $class->Get_Paths($options->{input});
    my $input_full = $input_paths->{fullpath};
    my $output_name = basename($input_full, ('.fasta', '.fa', '.fna', '.fsa'));
    $output_name = qq"${output_name}.tsv";
    my $output_dir = qq"outputs/$options->{jprefix}trinotate_$input_paths->{dirname}";
    my $comment = qq!## This is a trinotate submission script
!;
    my $jstring = qq!mkdir -p ${output_dir}
start=\$(pwd)
cd ${output_dir}
ln -sf $input_paths->{fullpath} .
if [ -f $input_paths->{filename}.gene_trans_map ]; then
  ln -s $input_paths->{filename}.gene_trans_map .
else
  ids=\$(grep "^>" $input_paths->{fullpath} | sed 's/>//g' | awk '{print \$1}')
  rm -f $input_paths->{filename}.gene_trans_map
  for i in \${ids}; do
    echo "\${i}	\${i}" >> $input_paths->{filename}.gene_trans_map
  done
fi

${trinotate_exe_dir}/auto/$options->{trinotate} \\
  --conf ${trinotate_exe_dir}/auto/$options->{config} \\
  --Trinotate_sqlite ${trinotate_exe_dir}/sample_data/Trinotate.boilerplate.sqlite \\
  --transcripts $input_paths->{filename} \\
  --gene_to_trans_map $input_paths->{filename}.gene_trans_map \\
  --CPU 6 \\
  2>trinotate_${job_name}.err \\
  1>trinotate_${job_name}.out
mv Trinotate.tsv ${output_name}
cd \${start}
!;
    my $trinotate = $class->Submit(
        cpus => 6,
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => qq"trinotate_$input_paths->{filename}_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => 80,
        modules => $options->{modules},
        output => qq"${output_dir}/${output_name}",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => 'large',
        jwalltime => '144:00:00',
        );
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($trinotate);
}


=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO


=cut
sub Watson_Plus {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jprefix => '44',
        );

    my $input_seq = $options->{input};
    my $job_name = 'watsonplus';
    my $output_dir = qq"outputs/$options->{jprefix}${job_name}";
    make_path($output_dir);
    my $watson_prodigal = $class->Bio::Adventure::Annotation::Prodigal(
        gcode => '11',
        input => $options->{input},
        jdepends => $options->{jdepends},
        jname => $job_name,
        jprefix => $options->{jprefix},
        modules => ['prodigal'],
        output_dir => $output_dir,
        species => 'phages',);
    $options->{jdepends} = $watson_prodigal->{job_id};
    my $input_gff = $watson_prodigal->{output_gff};
    my $output_file = basename($options->{input});
    $output_file = qq"${output_dir}/${output_file}";
    my $comment_string = qq"## This takes the prodigal output and checks to see which strand has
## more ORFs, if it is the crick strand, then the chromsome is reverse complemented.
";
    my $jstring = qq?
use Bio::Adventure::Annotation;
\$result = Bio::Adventure::Annotation::Watson_Rewrite(\$h,
  gff => '${input_gff}',
  input => '$options->{input}',
  jdepends => '$options->{jdepends}',
  jname => '$options->{jname}',
  jprefix => '$options->{jprefix}',
  output => '${output_file}',
  output_dir => '${output_dir}',);
?;
    my $rewrite = $class->Submit(
        comment => $comment_string,
        gff => $options->{gff},
        input => $options->{input},
        jdepends => $options->{jdepends},
        jname => $options->{jname},
        jprefix => $options->{jprefix},
        jstring => $jstring,
        language => 'perl',
        output => $output_file,
        output_dir => $output_dir,
        shell => '/usr/bin/env perl',);
    return($rewrite);
}

sub Watson_Rewrite {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],);

    my $input_gff = $options->{gff};
    my $input_seq = $options->{input};
    my $job_name = 'watsonplus';
    my $output_dir = qq"outputs/$options->{jprefix}${job_name}";
    make_path($output_dir);

    my $orf_count = {};
    my $annotation_in = FileHandle->new("<${input_gff}");
    my $write_file = basename($options->{input});
    $write_file = qq"${output_dir}/${write_file}";
    my $log_file = basename($write_file, ('.fasta'));
    $log_file = qq"${output_dir}/${log_file}.log";
    print "TESTME: $log_file\n";
    my $read_fasta = Bio::SeqIO->new(-file => qq"<$options->{input}", -format => 'Fasta');
    my $write_fasta = Bio::SeqIO->new(-file => qq">${write_file}", -format => 'Fasta');
    my $write_log = FileHandle->new(">${log_file}");
    print $write_log "Counting ORFs on the current watson and crick strand.
If the crick strand is larger, flipping them.
";
    ## A prodigal line looks like:
    ## EAb06   Prodigal_v2.6.3 CDS     3205    3414    8.7     -       0       ID=1_10;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.248;conf=88.17;score=8.74;cscore=-0.38;sscore=9.12;rscore=8.34;uscore=-0.93;tscore=2.86;
    while (my $line = <$annotation_in>) {
        next if ($line =~ /^#/);
        my ($name, $prod, $type, $start, $end, $score, $strand, $phase, $annot) = split(/\t/, $line);
        my @annot_list = split(/;/, $annot);
        my $id_annot = $annot_list[0];
        my ($contig, $orf);
        if ($id_annot =~ /ID=(.*)?_(\d+)$/) {
            $contig = $1;
            $orf = $2;
        } else {
            die("Could not get contig and orf.");
        }
        if (!defined($orf_count->{$contig})) {
            $orf_count->{$contig} = {
                plus => 0,
                minus => 0,
            };
        }
        if ($strand eq '+' || $strand eq '1') {
            $orf_count->{$contig}->{plus}++;
        } else {
            $orf_count->{$contig}->{minus}++;
        }
    } ## End iterating over every detected ORF
    $annotation_in->close();
    ## Now we should have some idea of which strand has the most ORFs for each contig.

  TESTLOOP: while (my $seq = $read_fasta->next_seq) {
      my $id = $seq->id;
      my $test = $orf_count->{$id};
      if ($test->{plus} >= $test->{minus}) {
          ## Then leave this contig alone.
          $write_fasta->write_seq($seq);
          print $write_log "$id was unchanged and has $test->{plus} plus and $test->{minus} minus ORFs.\n";
      } else {
          my $tmp = $seq->revcom();
          $write_fasta->write_seq($tmp);
          print $write_log "$id was reverse-complemented and now has $test->{minus} plus and $test->{plus} minus ORFs.\n";
      }
  } ## Done flipping sequence files.
    $write_log->close();
    return($orf_count);
}

1;
