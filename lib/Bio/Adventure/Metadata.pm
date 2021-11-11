package Bio::Adventure::Metadata;
## LICENSE: gplv2
## ABSTRACT:  Kitty!
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

=head1 C<Metadata>

This package is intended to collect and collate data from other tools into some
useful and readable formats.  Practically speaking, this mostly relates to
the annotation pipelines, which gather data from many sources into a single
genbank output file.

One other set of functions which I will move here include the various XYZ_Stats()
functions which gather the stdout/stderr from various tools in order to write a
csv summary.

=cut

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

=head2 C<Get_Aragorn>

=cut
sub Get_Aragorn {
    my ($class, %args) = @_;
    unless (-r $args{input}) {
        return(undef);
    }
    my $input_txt = FileHandle->new($args{input});
    my %aragorn_features = ();
    my $current_contig = '';
    my $current_description = '';
    my $current_sequence = '';
    my $last_line = '';
    my $counter = 0;
    my $number = '';
    my $type = '';
    my $coords = '';
    my $unknown = '';
    my $anticodon = '';
    my $annotation = '';
    my $trna_seq = '';
    my $start = '';
    my $end = '';
    my $strand = '+1';
    my $seq_portion = '';
    my $finished_annotation = 0;
    my @trna_annots = ();
    my %single_annot = ();
    my $orf_number = 0;
    my $formatted_number = sprintf("%04d", $orf_number);
  INPUT: while (my $line = <$input_txt>) {
      $counter++;
      chomp $line;
      if ($line =~ /^\d+ genes found/) {
          $current_contig = $last_line;
          $current_contig =~ s/^>//g;
          ## print "genes found line: $current_contig\n";
          $counter = 0;
      }
      ## if ($counter == 1) {
      if ($line =~ /^\d+\s+\S+\s+\S+\s+\d+\s+\(\w+\)$/) {
          $counter = 1;
          ($number, $type, $coords, $unknown, $anticodon) = split(/\s+/, $line);
          ($start, $end) = split(/\,/, $coords);
          ## print "TESTME PRE: $start\n";
          my $complement = $start;
          $complement =~ s/^(\w*).*$/$1/g;
          $strand = '+1';
          if ($complement eq 'c') {
              $strand = '-1';
          }
          ## print "TESTME AFTER COMP: $strand\n";

          $start =~ s/^\w*\[(\d+)$/$1/g;
          ## print "TESTME POST: $start\n";
          $end =~ s/\]//g;
          $anticodon =~ s/\(//g;
          $anticodon =~ s/\)//g;
          ## print "counter1: n:$number t:$type c:$coords u:$unknown a:$anticodon f:$start t:$end\n";
          $formatted_number = sprintf("%04d", $number);
      }
      if ($counter == 2) {
          my ($annot, $start_end) = split(/\s+c*\[/, $line);
          $annotation = $annot;
          $annotation =~ s/^>//g;
          ## print "counter2: $annotation\n";
      }
      if ($counter >= 3) {
          ## print "counter3+ $trna_seq\n";
          $seq_portion .= $line;
          $seq_portion =~ s/\s+//g;
          if ($seq_portion =~ /^\w{50}$/) {
              ## print "This is the max length sequence line.\n";
              $trna_seq = $seq_portion;
          } elsif ($seq_portion =~ /^\w+/) {
              $trna_seq .= $seq_portion;

              $counter = 1;
              my $orf_id = qq"${current_contig}_tRNA_${formatted_number}";
              $single_annot{contig} = $current_contig;
              $single_annot{formatted} = $orf_id;
              $single_annot{type} = $type;
              $single_annot{start} = $start;
              $single_annot{end} = $end;
              $single_annot{strand} = $strand;
              $single_annot{annotation} = $annotation;
              $single_annot{sequence} = $seq_portion;
              $single_annot{random} = $unknown;
              $single_annot{anticodon} = $anticodon;
              my %this_annot = %single_annot;
              push(@trna_annots, \%this_annot);
              $seq_portion = '';
          }
      }
      $last_line = $line;
  }
    $input_txt->close();

    ## We have an array of tRNA annotations, make them into SeqFeatures;
    my @aragorn_features = ();
    for my $a (@trna_annots) {
        my %trna = %{$a};
        my $gene_feature = Bio::SeqFeature::Generic->new(
            -primary => 'gene',
            -seq_id => $trna{contig},
            -start => $trna{start},
            -end => $trna{end},
            -strand => $trna{strand},
            -display_name => $trna{formatted},
            -tag => {
                'locus_tag' => $trna{formatted},
            });
        push(@aragorn_features, $gene_feature);
        my $trna_feature = Bio::SeqFeature::Generic->new(
            -primary => 'tRNA',
            -seq_id => $trna{contig},
            -source => 'aragorn',
            -start => $trna{start},
            -end => $trna{end},
            -strand => $trna{strand},
            -display_name => $trna{formatted},
            -tag => {
                'locus_tag' => $trna{formatted},
                'product' => $trna{type},
                'note' => $trna{annotation},
            },);
        push(@aragorn_features, $trna_feature);
    }
    return(\@aragorn_features);
}

=head2 C<Kraken_Accession>

Read the kraken_report.txt to extract the lowest common ancestor which has the most reads.
In order to get the taxonomy ID while doing this, we will need to read
the kraken output file which contains the status of each read; then
extract the third column and find out which taxon (except 0) which
has the most reads.  Once we have that taxonomy ID, we can go to


=cut
sub Kraken_Best_Hit {
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

=head2 C<Merge_Annotations>

Pull a series of annotation data sources into a single table of data.
The resulting table should be useful for tbl2asn.

This function basically just puts Merge_Annotations_Make_Gbk onto the
dependency chain.

=cut
sub Merge_Annotations {
   my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input_fsa', 'input_genbank', 'input_tsv'],
        input_abricate => '',
        input_classifier => '',
        input_glimmer => '',
        input_interpro => '',
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
\$h->Bio::Adventure::Metadata::Merge_Annotations_Make_Gbk(
  input_fsa => '$options->{input_fsa}',
  input_genbank => '$options->{input_genbank}',
  input_tsv => '$options->{input_tsv}',
  input_abricate => '$options->{input_abricate}',
  input_classifier => '$options->{input_classifier}',
  input_glimmer => '$options->{input_glimmer}',
  input_phageterm => '$options->{input_phageterm}',
  input_interpro => '$options->{input_interpro}',
  input_prodigal => '$options->{input_prodigal}',
  input_trinotate => '$options->{input_trinotate}',
  jdepends => '$options->{jdepends}',
  jprefix => '$options->{jprefix}',
  jname => 'merge_annotations',);
?;
   my $merge_job = $class->Submit(
       input_fsa => $options->{input_fsa},
       input_genbank => $options->{input_genbank},
       input_tsv => $options->{input_tsv},
       input_abricate => $options->{input_abricate},
       input_classifier => $options->{input_classifier},
       input_glimmer => $options->{input_glimmer},
       input_interpro => $options->{input_interpro},
       input_phageterm => $options->{input_phageterm},
       input_prodigal => $options->{input_prodigal},
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


=head2 C<Merge_Annotations_Make_Gbk>

This does the actual work of merging various annotation sources into a
new and more interesting genbank file.

=cut
sub Merge_Annotations_Make_Gbk {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input_fsa', 'input_genbank', 'input_tsv'],
        input_abricate => '',
        input_aragorn => '',
        input_classifier => '',
        input_glimmer => '',
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
        template_sbt => '/bio/reference/tbl2asn_template.sbt',);
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
        'abricate_vfdb' => 'vfdb',);

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

    print $log_fh "Reading tsv data from $options->{input_tsv} to start.\n";
    ## 1a above, the starter set of annotations.
    unless (-r $options->{input_tsv}) {
        die("Unable to find the tsv input: $options->{input_tsv}, this is required.\n");
    }
    unless (-r $options->{input_genbank}) {
        die("Unable to find the prokka genbank output: $options->{input_genbank}, this is required.\n");
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

    $merged_data = Merge_Start_Data(
        primary_key => $options->{primary_key},
        output => $merged_data,
        input => $options->{input_tsv});

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

    my %dtr_features = ();
    if ($options->{input_phageterm} && -r $options->{input_phageterm}) {
    ## Pull the direct-terminal-repeats from phageterm if they exist.
        %dtr_features = $class->Bio::Adventure::Phage::Get_DTR(
            input => $options->{input_phageterm});
        print $log_fh "Adding phageterm DTRs.\n";
    } else {
        print $log_fh "Not adding phageterm DTRs.\n";
    }

    ## An array reference;
    my $aragorn_features;
    if ($options->{input_aragorn} && -r $options->{input_aragorn}) {
    ## Pull the direct-terminal-repeats from phageterm if they exist.
        $aragorn_features = $class->Bio::Adventure::Metadata::Get_Aragorn(
            input => $options->{input_aragorn});
        print $log_fh "Adding aragorn tRNA annotations.\n";
    } else {
        print $log_fh "Not adding aragorn tRNA annotations.\n";
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
        foreach my $ara (@{$aragorn_features}) {
            unshift(@feature_list, $ara);
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
        my $tbl_written = Bio::Adventure::Annotation_Genbank::Write_Tbl_from_SeqFeatures(
            tbl_file => $output_tbl,
            taxonomy_information => $taxonomy_information,
            sequences => \@new_seq,
            features => \@new_features);
    } ## End Iterating over every sequence


    ## Remember that tbl2asn assumes the input files are all in the same directory.
    ## Before running tbl2asn, write a fresh copy of the fsa file containing the detected phage taxonomy.
    my $fsa_in = FileHandle->new("<$options->{input_fsa}");
    my $fsa_out = FileHandle->new(">$output_fsa");
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
    print $log_fh "Running ${tbl_command}\n";
    my $tbl2asn_result = qx"${tbl_command}";
    my $sed_command = qq"sed 's/COORDINATES: profile/COORDINATES:profile/' ${output_dir}/${output_name}.gbf | sed 's/product=\"_/product=\"/g' >${output_dir}/${output_name}.gbk";
    my $sed_result = qx"${sed_command}";

    ## Now lets pull everything from the merged data and make a hopefully pretty xlsx file.
    print $log_fh "Writing final xlsx file of the annotations to $output_xlsx\n";

    my $written = $class->Bio::Adventure::Metadata::Write_XLSX(
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
        output_gbk => $output_gbk,
        output_tbl => $output_tbl,
        output_xlsx => $output_xlsx,};
    return($job);
}

=head2 C<Merge_Trinotate>

Given a hash of annotation data, read a trinotate tsv file and pull
the information from it into that hash and pass it back to the caller.

=cut
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

=head2 C<Merge_Abricate>

Given a hash of annotation data, read the results from abricate and pull
that information into the input hash and pass it back to the caller.

=cut
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

=head2 C<Merge_Classifier>

Given a hash of annotation data, read ICTV classification information
and pass it back.

=cut
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

=head2 C<Merge_Prodigal>

Given a hash of annotation data, read some prodigal information and pass
the information from it into that hash and pass it back to the caller.

=cut
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

sub Merge_Glimmer {
    ## Rewrite this to read the predict file.
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

=head2 C<Merge_Interpro>

Given a hash of annotation data, read an interpro tsv file and pull
the information from it into that hash and pass it back to the caller.

=cut
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

=head2 C<Merge_Prokka>

Given a hash of annotation data, read a prokka log and pull
the information from it into that hash and pass it back to the caller.

=cut
sub Merge_Start_Data {
    my (%args) = @_;
    my $primary_key = $args{primary_key};
    my $merged_data = $args{output};
  ##  my $in = Bio::FeatureIO->new(-file => $args{input}, -format => 'GFF', -version => 3);
  ##FEATURES: while (my $f = $in->next_feature()) {
  ##    next FEATURES if ($f->primary_tag eq 'source');
  ##    next FEATURES if ($f->primary_tag eq 'gene');
  ##    my @tags = $f->get_all_tags();
  ##    my $id = $f->seq_id();
  ##    my $start = $f->start;
  ##    my $end = $f->end;
  ##    my $strand = $f->strand;
  ##    my $key = $f->display_name;
  ##    $merged_data->{$key}->{locus_tag} = $key;
  ##    $merged_data->{$key}->{start} = $start;
  ##    $merged_data->{$key}->{end} = $end;
  ##    $merged_data->{$key}->{strand} = $strand;
  ##    foreach my $t (@tags) {
  ##        my @values = $f->get_tag_values($t);
  ##        my $value = $values[0];
  ##        $merged_data->{$key}->{$t} = $value;
  ##    }
  ##}
    my $input_tsv = Text::CSV_XS::TSV->new({ binary => 1, });
    open(my $input_fh, "<:encoding(utf8)", $args{input});
    my $header = $input_tsv->getline($input_fh);
    $input_tsv->column_names($header);
  ROWS: while (my $row = $input_tsv->getline_hr($input_fh)) {
      ## next ROWS if ($row->{ftype} eq 'gene');
      my $key = $row->{$primary_key};
      $merged_data->{$key}->{$primary_key} = $key;
      foreach my $colname (@{$header}) {
          ## Check that we already filled this data point
          if ($merged_data->{$key}->{$colname}) {
              ## There is something here.
          } else {
              if ($row->{$colname}) {
                  $merged_data->{$key}->{$colname} = $row->{$colname};
              }
          }
      } ## Run over the columns
  }
    close $input_fh;
    return($merged_data);
}

=head2 C<Write_XLSX>

Given a pile of annotations in tabular format, send them to an xlsx
file.  I chose a set of annotation sources and columns that I like,
but if I were smarter I would make it mutable.

=cut
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

=head2 C<BT1_Stats>

Collect some alignment statistics from bowtie1.

This uses a little creative grepping in order to extract information out of the STDERR file
from bowtie, which happens to contain messages with the number of reads mapped etc.  The authors
of bowtie2/tophat/hisat2/etc all continued to write error logs with similar messages, so this
function is essentially the template for gathering statistics from all of those tools.

It requires the 'input' argument, which is the bowtie error file.

=cut
sub BT1_Stats {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input']);
    my $bt_input = $options->{input};
    my $bt_type = "";
    $bt_type = $options->{bt_type} if ($options->{bt_type});
    my $jname = "stats";
    $jname = $options->{jname} if ($options->{jname});
    my $jobid = qq"$options->{jbasename}_stats";
    my $count_table = "";
    $count_table = $options->{count_table} if ($options->{count_table});
    my $stat_output = qq"outputs/bowtie_stats.csv";
    my $comment = qq!## This is a stupidly simple job to collect alignment statistics.!;
    my $jstring = qq!
if [ \! -r "${stat_output}" ]; then
  echo "name,type,original_reads,reads,one_hits,failed,samples,rpm,count_table" > ${stat_output}
fi
original_reads=0
if [ -r "outputs/trimomatic_stats.csv" ]; then
  original_reads_tmp=\$(tail -n 1 outputs/trimomatic_stats.csv | awk -F, '{print \$2}')
  original_reads=\${original_reads_tmp:-0}
fi
reads_tmp=\$(grep "^# reads processed" ${bt_input} | awk -F: '{print \$2}' | sed 's/ //g')
reads=\${reads_tmp:-0}
one_align_tmp=\$(grep "^# reads with at least one reported" ${bt_input} | awk -F": " '{print \$2}' | sed 's/ .*//g')
one_align=\${one_align_tmp:-0}
failed_tmp=\$(grep "^# reads that failed to align" ${bt_input} | awk -F": " '{print \$2}' | sed 's/ .*//g')
failed=\${failed_tmp:-0}
sampled_tmp=\$(grep "^# reads with alignments sampled" ${bt_input} | awk -F": " '{print \$2}' | sed 's/ .*//g')
sampled=\${sampled_tmp:-0}
rpm_tmp=\$(perl -e "printf(1000000 / \${one_align})" 2>/dev/null)
rpm=\${rpm_tmp:-0}
stat_string=\$(printf "$options->{jbasename},${bt_type},%s,%s,%s,%s,%s,%s,${count_table}" "\${original_reads}" "\${reads}" "\${one_align}" "\${failed}" "\${sampled}" "\$rpm")
echo "\$stat_string" >> ${stat_output}!;
    my $stats = $class->Submit(
        comment => $comment,
        input => $bt_input,
        jdepends => $options->{jdepends},
        jname => $jname,
        jprefix => $options->{jprefix},
        jstring => $jstring,
        cpus => 1,
        jmem => 1,
        jqueue => 'throughput',);
    return($stats);
}

=head2 C<BT2_Stats>

Collects alignment statistics from bowtie2.  It is mostly a copy/paste from BT1_Stats().

=cut
sub BT2_Stats {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args,
        required => ['input']);
    my $bt_input = $options->{input};
    my $bt_type = $options->{bt_type};
    my $jname = "bt2_stats";
    $jname = $options->{jname} if ($options->{jname});
    my $jobid = qq"$options->{jbasename}_stats";
    my $count_table = "";
    $count_table = $options->{count_table} if ($options->{count_table});
    my $comment = qq!## This is a stupidly simple job to collect alignment statistics.!;
    my $output = "outputs/bowtie2_stats.csv";
    my $jstring = qq!
if [ \! -e "${output}" ]; then
    echo "original reads, single hits, failed reads, multi-hits, rpm" > ${output}
fi
original_reads_tmp=\$(grep " reads; of these" "${bt_input}" 2>/dev/null | awk '{print \$1}' | sed 's/ //g')
original_reads=\${original_reads_tmp:-0}
one_align_tmp=\$(grep " aligned exactly 1 time" "${bt_input}" | awk '{print \$1}' | sed 's/ .*//g')
one_align=\${one_align_tmp:-0}
failed_tmp=\$(grep " aligned 0 times" "${bt_input}" | awk '{print \$1}' | sed 's/ .*//g')
failed=\${failed_tmp:-0}
sampled_tmp=\$(grep " aligned >1 times" "${bt_input}" | awk '{print \$1}' | sed 's/ .*//g')
sampled=\${sampled_tmp:-0}
rpm_tmp=\$(perl -e "printf(1000000 / \$(( \${one_align} + \${sampled} )) ) " 2>/dev/null)
rpm=\${rpm_tmp:-0}
stat_string=\$(printf "$options->{jbasename},${bt_type},%s,%s,%s,%s,%s" "\${original_reads}" "\${one_align}" "\${failed}" "\${sampled}" "\${rpm}")
echo "\$stat_string" >> ${output}!;
    my $stats = $class->Submit(
        comment => $comment,
        input => $bt_input,
        jname => $jname,
        jdepends => $options->{jdepends},
        jprefix => $options->{jprefix},
        jstring => $jstring,
        cpus => 1,
        jmem => 1,
        jqueue => 'throughput',);
    return($stats);
}

=head2 C<BWA_Stats>

Collect some alignment statistics from bwa.

=cut
sub BWA_Stats {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $aln_input = $options->{aln_output};
    $aln_input = qq"${aln_input}.stats";
    my $mem_input = $options->{mem_output};
    $mem_input = qq"${mem_input}.stats";
    my $stat_output = qq"outputs/bwa_stats.csv";

    my $jname = "bwa_stats";
    $jname = $options->{jname} if ($options->{jname});
    my $jobid = qq"$options->{jbasename}_stats";
    my $count_table = "";
    $count_table = $options->{count_table} if ($options->{count_table});
    my $comment = qq!## This is a stupidly simple job to collect alignment statistics.!;
    my $jstring = qq!
if [ \! -r "${stat_output}" ]; then
    echo "# original reads, reads used, aln-aligned reads, mem-aligned reads, rpm" > ${stat_output}
fi
original_reads=0
if [ -r "outputs/trimomatic_stats.csv" ]; then
    original_reads_tmp=\$(tail -n 1 outputs/trimomatic_stats.csv | awk -F, '{print \$2}')
    original_reads=\${original_reads_tmp:-0}
fi
reads_tmp=\$(grep "^Total reads: " ${aln_input} | awk '{print \$3}' | sed 's/ //g')
reads=\${reads_tmp:-0}
aln_aligned_tmp=\$(grep "^Mapped reads" ${aln_input} | awk '{print \$3}' | sed 's/ .*//g')
aln_aligned=\${aln_aligned_tmp:-0}
mem_aligned_tmp=\$(grep "^Mapped reads" ${mem_input} | awk '{print \$3}' | sed 's/ .*//g')
mem_aligned=\${mem_aligned_tmp:-0}
rpm_tmp=\$(perl -e "printf(1000000 / \${aligned})" 2>/dev/null)
rpm=\${rpm_tmp:-0}
stat_string=\$(printf "$options->{jbasename},%s,%s,%s,%s,%s,${count_table}" "\${original_reads}" "\${reads}" "\${aln_aligned}" "\${mem_aligned}" "\$rpm")
echo "\${stat_string}" >> ${stat_output}!;
    my $stats = $class->Submit(
        comment => $comment,
        input => $aln_input,
        depends => $options->{jdepends},
        jname => $jname,
        jprefix => $options->{jprefix},
        jstring => $jstring,
        cpus => 1,
        jmem => 1,
        jqueue => "throughput",);
    return($stats);
}

=head2 C<Fastqc_Stats>

Collect some information from the fastqc output files and present them in a
simple-to-read csv file.

=cut
sub Fastqc_Stats {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jname => 'fqcst',
        paired => 1,);
    ## Dereferencing the options to keep a safe copy.
    my $jname = $options->{jname};
    ## Dereferencing the options to keep a safe copy.
    my $input_file = qq"$options->{input}/fastqc_data.txt";
    if ($options->{paired}) {
        $input_file = qq"$options->{input}/fastqc_data.txt";
    }
    my $stat_output = qq"outputs/fastqc_stats.csv";
    if ($options->{direction}) {
        $stat_output = qq"outputs/fastqc_$options->{direction}_stats.csv";
        $jname = qq"$options->{jname}_$options->{direction}";
    }
    my $comment = qq!## This is a stupidly simple job to collect alignment statistics.!;
    my $jstring = qq!
if [ \! -r "${stat_output}" ]; then
  echo "name,total_reads,poor_quality,per_quality,per_base_content,per_sequence_gc,per_base_n,per_seq_length,over_rep,adapter_content,kmer_content" > $stat_output
fi
total_reads_tmp=\$(grep "^Total Sequences" ${input_file} | awk -F '\\\\t' '{print \$2}')
total_reads=\${total_reads_tmp:-0}
poor_quality_tmp=\$(grep "^Sequences flagged as poor quality" ${input_file} | awk -F '\\\\t' '{print \$2}')
poor_quality=\${poor_quality_tmp:-0}
per_quality_tmp=\$(grep "Per base sequence quality" ${input_file} | awk -F '\\\\t' '{print \$2}')
per_quality=\${per_quality_tmp:-0}
per_base_content_tmp=\$(grep "Per base sequence content" ${input_file} | awk -F '\\\\t' '{print \$2}')
per_base_content=\${per_base_content_tmp:-0}
per_sequence_gc_tmp=\$(grep "Per sequence GC content" ${input_file} | awk -F '\\\\t' '{print \$2}')
per_sequence_gc=\${per_sequence_gc_tmp:-0}
per_base_n_tmp=\$(grep "Per base N content" ${input_file} | awk -F '\\\\t' '{print \$2}')
per_base_n=\${per_base_n_tmp:-0}
per_seq_length_tmp=\$(grep "Sequence Length Distribution" ${input_file} | awk -F '\\\\t' '{print \$2}')
per_seq_length=\${per_seq_length_tmp:-0}
over_rep_tmp=\$(grep "Overrepresented sequences" ${input_file} | awk -F '\\\\t' '{print \$2}')
over_rep=\${over_rep_tmp:-0}
adapter_content_tmp=\$(grep "Adapter Content" ${input_file} | awk -F '\\\\t' '{print \$2}')
adapter_content=\${adapter_content_tmp:-0}
kmer_content_tmp=\$(grep "Kmer Content" ${input_file} | awk -F '\\\\t' '{print \$2}')
kmer_content=\${kmer_content_tmp:-0}

stat_string=\$(printf "$options->{jname},%s,%s,%s,%s,%s,%s,%s,%s,%s,%s" "\${total_reads}" "\${poor_quality}" "\${per_quality}" "\${per_base_content}" "\${per_sequence_gc}" "\${per_base_n}" "\${per_seq_length}" "\${over_rep}" "\${adapter_content}" "\${kmer_content}")
echo "\$stat_string" >> ${stat_output}
!;
    my $stats = $class->Submit(
        comment => $comment,
        cpus => 1,
        input => $input_file,
        jdepends => $options->{jdepends},
        jmem => 1,
        jname => $jname,
        jprefix => $options->{jprefix},
        jqueue => 'throughput',
        jstring => $jstring,
        jwalltime => '00:1:00',
        output => qq"${stat_output}",);
    ## Added to return the state of the system to what it was
    ## before we messed with the options.
    return($stats);
}

=head2 C<HT2_Stats>

Collect alignment statistics from hisat 2.

=cut
sub HT2_Stats {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        output_dir => 'outputs',);
    my $ht_input = $options->{ht_input};
    my $jname = "ht2_stats";
    $jname = $options->{jname} if ($options->{jname});
    my $jobid = qq"$options->{jbasename}_stats";
    my $count_table = "";
    $count_table = $options->{count_table} if ($options->{count_table});
    my $comment = qq!## This is a stupidly simple job to collect alignment statistics.!;
    my $output = "$options->{output_dir}/hisat2_stats.csv";
    my $jstring = qq!
if [ \! -e "${output}" ]; then
    echo "id, original reads, single hits, failed reads, multi-hits, rpm" > ${output}
fi
original_reads_tmp=\$(grep " reads; of these" "${ht_input}" 2>/dev/null | awk '{print \$1}' | sed 's/ //g')
original_reads=\${original_reads_tmp:-0}
one_align_tmp=\$(grep " aligned exactly 1 time" "${ht_input}" | awk '{print \$1}' | sed 's/ .*//g')
one_align=\${one_align_tmp:-0}
failed_tmp=\$(grep " aligned 0 times" "${ht_input}" | tail -n 1 | awk '{print \$1}' | sed 's/ .*//g')
failed=\${failed_tmp:-0}
sampled_tmp=\$(grep " aligned >1 times" "${ht_input}" | awk '{print \$1}' | sed 's/ .*//g')
sampled=\${sampled_tmp:-0}
rpm_tmp=\$(perl -e "printf(1000000 / \$(( \${one_align} + \${sampled} )) ) " 2>/dev/null)
rpm=\${rpm_tmp:-0}
stat_string=\$(printf "$options->{jbasename},%s,%s,%s,%s,%s" "\${original_reads}" "\${one_align}" "\${failed}" "\${sampled}" "\${rpm}")
echo "\$stat_string" >> ${output}!;
    my $stats = $class->Submit(
        comment => $comment,
        input => $ht_input,
        jname => $jname,
        jdepends => $options->{jdepends},
        jprefix => $options->{jprefix},
        jstring => $jstring,
        output => $output,
        cpus => 1,
        jmem => 1,
        jqueue => 'throughput',);
    return($stats);
}

=head2 C<Salmon_Stats>

Collect some summary statistics from a salmon run.

=cut
sub Salmon_Stats {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $jname = "stats";
    $jname = $options->{jname} if ($options->{jname});
    my $jobid = qq"$options->{jbasename}_stats";
    my $outdir = dirname($options->{input});
    my $output = qq"${outdir}/salmon_stats.csv";
    my $comment = qq"## This is a stupidly simple job to collect salmon alignment statistics.";
    my $jstring = qq!
if [ \! -r "${output}" ]; then
  echo "basename,species,fragments,assigned,consistent,inconsistent,bias" > ${output}
fi
reads_tmp=\$(grep "num_compatible" $options->{input} | awk '{print \$2}' | sed 's/\,//g')
reads=\${reads_tmp:-0}
aligned_tmp=\$(grep "num_assigned" $options->{input} | awk '{print \$2}' | sed 's/\,//g')
aligned=\${aligned_tmp:-0}
consistent_tmp=\$(grep "concordant" $options->{input} | awk '{print \$2}' | sed 's/\,//g')
consistent=\${consistent_tmp:-0}
inconsistent_tmp=\$(grep "inconsistent" $options->{input} | awk '{print \$2}' | sed 's/\,//g')
inconsistent=\${inconsistent_tmp:-0}
bias_tmp=\$(grep "mapping_bias" $options->{input} | awk '{print \$2}' | sed 's/\,//g')
bias=\${bias_tmp:-0}
stat_string=\$(printf "$options->{jbasename},$options->{species},%s,%s,%s,%s,%s" "\${reads}" "\${aligned}" "\${consistent}" "\${inconsistent}" "\${bias}")
echo "\$stat_string" >> "${output}"!;
    my $stats = $class->Submit(
        comment => $comment,
        cpus => 1,
        input => $options->{input},
        jname => $jname,
        jdepends => $options->{jdepends},
        jprefix => $args{jprefix},
        jstring => $jstring,
        jmem => 1,
        jqueue => 'throughput',
        output => $output,);
    return($stats);
}

=head2 C<Tophat_Stats>

Collect alignment statistics from the accepted_hits.bam/unaligned.bam files
generated by a tophat run.

=cut
sub Tophat_Stats {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $accepted_input = $options->{accepted_input};
    my $accepted_output = qq"${accepted_input}.stats";
    my $unaccepted_input = $options->{unaccepted_input};
    my $unaccepted_output = qq"${unaccepted_input}.stats";
    my $read_info = $options->{prep_input};
    my $jname = "stats";
    $jname = $options->{jname} if ($options->{jname});
    my $jobid = qq"$options->{jbasename}_stats";
    my $count_table = "";
    $count_table = $options->{count_table} if ($options->{count_table});
    my $output = "outputs/tophat_stats.csv";
    my $comment = qq!## This is a stupidly simple job to collect tophat alignment statistics.!;
    my $jstring = qq!
if [ \! -r "${output}" ]; then
  echo "basename,species,original_reads,aligned_reads,failed_reads,rpm,count_table" > ${output}
fi
bamtools stats < "${accepted_input}" \\
    2>${accepted_output} 1>&2 && \\
  bamtools stats < "${unaccepted_input}" \\
    2>${unaccepted_output} 1>&2

original_reads=0
if [ -r "outputs/trimomatic_stats.csv" ]; then
  original_reads_tmp=\$(tail -n 1 outputs/trimomatic_stats.csv | awk -F, '{print \$2}')
  original_reads=\${original_reads_tmp:-0}
fi
reads_tmp=\$(grep "^reads_in " ${read_info} | awk -F= '{print \$2}' | sed 's/ //g')
reads=\${reads_tmp:-0}
aligned_tmp=\$(grep "^Total reads" ${accepted_output} | awk '{print \$3}' | sed 's/ .*//g')
aligned=\${aligned_tmp:-0}
failed_tmp=\$(grep "^Total reads" ${unaccepted_output} | awk '{print \$3}' | sed 's/ .*//g')
failed=\${failed_tmp:-0}
rpm_tmp=\$(perl -e "printf(1000000 / \${aligned})" 2>/dev/null)
rpm=\${rpm_tmp:-0}
stat_string=\$(printf "$options->{jbasename},$options->{species},%s,%s,%s,%s,%s,${count_table}" "\${original_reads}" "\${reads}" "\${aligned}" "\${failed}" "\$rpm")
echo "\$stat_string" >> "${output}"!;
    my $stats = $class->Submit(
        comment => $comment,
        cpus => 1,
        input => $accepted_input,
        jname => $jname,
        jdepends => $options->{jdepends},
        jprefix => $args{jprefix},
        jstring => $jstring,
        jmem => 1,
        jqueue => 'throughput',);
    return($stats);
}

=head2 C<Trimomatic_Stats>

Collect the trimming statistics from the output file 'trimomatic.out' and report
them in a .csv file by library.

=cut
sub Trimomatic_Stats {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        output_dir => 'outputs',
    );
    ## Dereferencing the options to keep a safe copy.
    my $basename = $options->{basename};
    my $input_file = "$options->{output_dir}/${basename}-trimomatic.out";
    my $jname = 'trimst';
    $jname = $options->{jname} if ($options->{jname});
    my $comment = qq!## This is a stupidly simple job to collect trimomatic statistics!;
    my $stat_output = qq"$options->{output_dir}/trimomatic_stats.csv";
    my $jstring = qq!
if [ \! -r ${stat_output} ]; then
  echo "total_reads,surviving_reads,dropped_reads" > ${stat_output}
fi
total_reads_tmp=\$(grep "^Input Reads" ${input_file} | awk '{print \$3}')
total_reads=\${total_reads_tmp:-0}
surviving_reads_tmp=\$(grep "^Input Reads" ${input_file} | awk '{print \$5}')
surviving_reads=\${surviving_reads_tmp:-0}
dropped_reads_tmp=\$(grep "^Input Reads" ${input_file} | awk '{print \$8}')
dropped_reads=\${dropped_reads_tmp:-0}

stat_string=\$(printf "${basename},%s,%s,%s" "\${total_reads}" "\${surviving_reads}" "\${dropped_reads}")
echo "\$stat_string" >> ${stat_output}
!;
    if ($options->{pairwise}) {
        ## The output looks a bit different for pairwise input:
        ## Input Read Pairs: 10000 Both Surviving: 9061 (90.61%) Forward Only Surviving: 457 (4.57%) Reverse Only Surviving: 194 (1.94%) Dropped: 288 (2.88%)
        $jstring = qq!
if [ \! -r ${stat_output} ]; then
  echo "total_reads,surviving_both,surviving_forward,surviving_reverse,dropped_reads" > ${stat_output}
fi
total_reads_tmp=\$(grep "^Input Read Pairs" ${input_file} | awk '{print \$4}')
total_reads=\${total_reads_tmp:-0}
surviving_both_tmp=\$(grep "^Input Read Pairs" ${input_file} | awk '{print \$7}')
surviving_both=\${surviving_both_tmp:-0}
surviving_forward_tmp=\$(grep "^Input Read Pairs" ${input_file} | awk '{print \$12}')
surviving_forward=\${surviving_forward_tmp:-0}
surviving_reverse_tmp=\$(grep "^Input Read Pairs" ${input_file} | awk '{print \$17}')
surviving_reverse=\${surviving_reverse_tmp:-0}
dropped_reads_tmp=\$(grep "^Input Read Pairs" ${input_file} | awk '{print \$20}')
dropped_reads=\${dropped_reads_tmp:-0}

stat_string=\$(printf "${basename},%s,%s,%s,%s,%s" "\${total_reads}" "\${surviving_both}" "\${surviving_forward}" "\${surviving_reverse}" "\${dropped_reads}")
echo "\$stat_string" >> ${stat_output}
!;
    }
    my $stats = $class->Submit(
        comment => $comment,
        cpus => 1,
        input => $input_file,
        jdepends => $options->{jdepends},
        jname => $jname,
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => 1,
        jqueue => 'throughput',
        jwalltime => '00:10:00',
        output => $stat_output,);
    return($stats);
}

=back

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

    L<bowtie> L<bowtie2> L<tophat> L<bwa> L<kallisto> L<samtools> L<htseq>

=cut

1;