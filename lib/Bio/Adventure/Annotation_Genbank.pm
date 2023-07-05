package Bio::Adventure::Annotation_Genbank;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';
use Data::Dumper;
use Data::Printer;

use Cwd qw"abs_path getcwd cwd";
use Data::Table;
use Data::Table::Excel qw"tables2xlsx";
use FileHandle;
use File::Basename;
use File::Copy qw"cp";
use File::Spec;
use File::Path qw"make_path";
use File::Which qw"which";
use File::ShareDir qw":ALL";
use List::MoreUtils qw"any uniq";
use Math::SigFigs qw"FormatSigFigs";
use Template;
use Text::CSV_XS::TSV;

use Bio::FeatureIO;
use Bio::Seq;
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::Lite;
use Bio::SeqIO;
use Bio::Tools::GuessSeqFormat;

=head2 C<Combine_CDS_Features>

 Merge features from multiple annotation sources.

 This is a generalized Feature-merger.  It attempts to keep a
 heirarchy of source priorities as well as maintain the annotations
 of the CDS predictions for the methods which provide useful or
 interesting information.

 It currently has some tortured logic to allow one to choose the
 preferred annotation source and also set the first/second source.
 I am pretty sure this logic is redundant.

=over

=item C<Arguments>

 first: First annotation source.
 first_name('prokka'): Set the source name from the first set.
 second: Guess!
 second_name('phanotate'): Set the source name for the second set.
 dominant('first'): Choose the annotation source which wins when there is a disagreement.
 notes: A hash of note annotations
 keep_both(0): Keep the union of the sources.

=back
=cut
sub Combine_CDS_Features {
    my %args = @_;
    my @first = @{$args{first}};
    my @second = @{$args{second}};
    my $keep_both = 0;
    $keep_both = $args{keep_both} if ($args{keep_both});
    my $assembly_type = 'phage';
    if (defined($args{assembly_type})) {
        $assembly_type = $args{assembly_type};
    }
    my $first_name = 'prokka';
    if (defined($args{first_name})) {
        $first_name = $args{first_name};
    }
    my $second_name = 'phanotate';
    if (defined($args{second_name})) {
        $second_name = $args{second_name};
    }
    my $dominant = 'first';
    if (defined($args{dominant})) {
        $dominant = $args{dominant};
    }
    my $notes = $dominant;
    if (defined($args{notes})) {
        $notes = $args{notes};
    }

    my @merged_cds = ();
    my $count = 0;
    my %finished_first = ();
    my %finished_second = ();
  FIRST: for my $first_cds (@first) {
      $count++;
      my $first_contig = $first_cds->seq_id();
      ## and the information about the location
      my $first_strand = $first_cds->strand();
      my $first_start = $first_cds->start();
      my $first_end = $first_cds->end();

      my ($first_fivep, $first_threep) = ($first_start, $first_end);
      if ($first_strand < 0) {
          ($first_fivep, $first_threep) = ($first_end, $first_start);
      }

      ## Get the set of tags and the locus tag specifically
      my @first_tags = $first_cds->get_all_tags();
      my @locus_tags = $first_cds->get_tag_values('locus_tag');
    SECOND: for my $second_cds (@second) {
        ## Another filtering opportunity presents itself here.
        my $second_contig = $second_cds->seq_id();
        ## Make sure we are dealing with the same contig.
        next SECOND unless ($second_contig eq $first_contig);
        ## We do not have to think about the type of these features.
        my $second_strand = $second_cds->strand();
        my $second_start = $second_cds->start;
        my $second_end = $second_cds->end;
        my ($second_fivep, $second_threep) = ($second_start, $second_end);
        if ($second_strand < 0) {
            ($second_fivep, $second_threep) = ($second_end, $second_start);
        }
        if (defined($finished_second{$second_threep})) {
            next SECOND;
        }

        my $first_length = scalar(@first);
        my $second_length = scalar(@second);
        if ($count < $first_length) {
            if ($second_end > $first_end) {
                last SECOND;
            }
        }

        if (($first_fivep eq $second_fivep) and ($first_threep eq $second_threep)) {
            if ($dominant eq 'first') {
                ## my $tag = $first_cds->add_tag_value('note', qq"cds_prediciton: ${first_name}");
                if ($notes eq 'both') {
                    my @second = $second_cds->get_tag_values('note');
                    $first_cds->add_tag_value('note', $second[0]);
                }
                push(@merged_cds, $first_cds);
            } else {
                if ($notes eq 'both') {
                    my @first = $first_cds->get_tag_values('note');
                    $second_cds->add_tag_values('note', $first[0]);
                }
                push(@merged_cds, $second_cds);
            }
            $finished_first{$first_threep} = 1;
            $finished_second{$first_threep} = 1;
        } elsif ($first_threep eq $second_threep) {
            if ($dominant eq 'first') {
                ## my $first_tag = $first_cds->add_tag_value('note', qq"cds_prediction: ${first_name}");
                push(@merged_cds, $first_cds);
                if ($keep_both) {
                    ## my $second_tag = $second_cds->add_tag_value('note', qq"cds_prediction: ${second_name}");
                    push(@merged_cds, $second_cds);
                }
            } else {
                ## my $second_tag = $second_cds->add_tag_value('note', qq"cds_prediction: ${second_name}");
                push(@merged_cds, $second_cds);
                if ($keep_both) {
                    ## my $first_tag = $first_cds->add_tag_value('note', qq"cds_prediction: ${first_name}");
                    push(@merged_cds, $first_cds);
                }
            }
            $finished_first{$first_threep} = 1;
            $finished_second{$first_threep} = 1;

            ## Important note here
            ## There is a missing piece of logic
            ## I have not taken into account the occasion when we have:
            ##        -------------------  phanotate
            ##             ----            glimmer
            ##   nor
            ##        ------------------- phanotate
            ##        -----               glimmer  (which should never happen)
            ##   nor
            ##        ------------------- phanotate
            ##                     ------ glimmer  (which is unlikely but possible)

        } elsif ($first_start < $second_start && $first_end >= $second_end) {
            print "The second feature is inside the first.  Skipping it.\n";
            print "1st start   2nd start         2nd end   1st end\n";
            print "$first_start  $second_start         $second_end  $first_end\n";
            $finished_second{$second_threep} = 1;
            next SECOND;
        } elsif ($first_start <= $second_start && $first_end > $second_end) {
            print "The second feature is inside the first.  Skipping it.\n";
            print "1st_start  2nd_start      2nd_end  1st_end\n";
            print "$first_start  $second_start         $second_end  $first_end\n";
            $finished_second{$second_threep} = 1;
            next SECOND;
        } else {
            $finished_second{$second_threep} = 1;
            ## my $second_tag = $second_cds->add_tag_value('note', qq"cds_prediction: ${second_name}");
            push(@merged_cds, $second_cds);
        }
    } ## End inner iteration (SECOND)
      if (!defined($finished_first{$first_threep})) {
          push(@merged_cds, $first_cds);
    }
      my $merged_cds_length = scalar(@merged_cds);
  } ## End outer iteration (FIRST)
    return(\@merged_cds);
}

=head2 C<Extract_Features>

  Extract keyword-based sequences from a genbank file.

=cut
sub Extract_Features {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        term => 'tail',
        type => 'CDS',
        id_tag => 'locus_tag',
        search_tag => 'note');
    my $type = $options->{type};
    unless (-r $options->{input}) {
        die("Unable to find input genbank file.");
    }
    my $output_dir = qq"outputs/extract_sequences";
    make_path($output_dir) unless (-d $output_dir);

    my $search = lc($options->{term});
    my $output_name = basename($options->{input}, ('.gz', '.xz'));
    $output_name = basename($output_name , ('.genbank', '.gbk', '.gb'));
    my $output_file = qq"${output_dir}/${output_name}_extracted.fasta";
    my $output = Bio::SeqIO->new(-file => qq">${output_file}",
                                 -format => 'Fasta');
    my $pep_file = qq"${output_dir}/${output_name}_extracted_pep.fasta";
    my $output_pep = Bio::SeqIO->new(-file => qq">${pep_file}",
                                     -format => 'Fasta');
    my $in = FileHandle->new("less $options->{input} |");
    my $seqio = Bio::SeqIO->new(-format => 'genbank', -fh => $in);
    my $num_hits = 0;
  SEQUENCES: while (my $seq = $seqio->next_seq) {
      my $seqid = $seq->id;
      my @feature_list = $seq->get_SeqFeatures();
    FEATURES: for my $feat (@feature_list) {
        my $primary = $feat->primary_tag(); ## E.g. the type gene/cds/misc/etc
        next FEATURES unless ($primary eq $options->{type});
        my $contig = $feat->seq_id(); ## The contig ID
        my $annot = $feat->annotation();
        my $type = $feat->primary_tag;
        my @names = $feat->get_tag_values($options->{id_tag});
        my $name = $names[0];
        ## use Data::Dumper;
        ## print Dumper $feat;
        my @tags = $feat->get_all_tags();
        my $chosen_tag = undef;
        if (defined($options->{id_tag}) && $options->{id_tag} ne '') {
            if (grep /^$options->{search_tag}$/, @tags) {
                $chosen_tag = $options->{search_tag};
            }
        }
        my @loci;
        if ($chosen_tag) {
            @loci = $feat->get_tag_values($chosen_tag);
            for my $l (@loci) {
                $l = lc($l);
                my $found = $l =~ m/${search}/;
                if ($found) {
                    $num_hits++;
                    print "Observed a potential hit, gene: ${l}\n";
                    my $cds_obj = Bio::Seq->new(-display_id => $name,
                                                -seq => $feat->seq->seq);
                    my $aa_obj = Bio::Seq->new(-display_id => $name,
                                               -seq => $feat->seq->translate->seq);

                    $output->write_seq($cds_obj);
                    $output_pep->write_seq($aa_obj);
                }
            }
        } else {
            for my $t (@tags) {
                @loci = $feat->get_tag_values($t);
                for my $l (@loci) {
                    $l = lc($l);
                    my $found = $l =~ m/${search}/;
                    if ($found) {
                        $num_hits++;
                        my $cds_obj = Bio::Seq->new(-display_id => $name,
                                                    -seq => $feat->seq->seq);
                        my $aa_obj = Bio::Seq->new(-display_id => $name,
                                                   -seq => $feat->seq->translate->seq);
                        $output->write_seq($cds_obj);
                        $output_pep->write_seq($aa_obj);
                    }
                }
            }
        }
    } ## End looking at features on this sequence
  } ## End looking at this sequence
    return($num_hits);
}

sub Extract_Notes {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        note_tag => 'note',
        gff_tag => 'locus_tag',
        required => ['input'],
        string => 'tail');

    my $string_name = $options->{string};
    $string_name =~ s/\s+/_/g;
    my $output_directory = 'outputs/extract_notes';
    make_path($output_directory, {verbose => 0}) if (!-d $output_directory);
    my $all_file = qq"${output_directory}/all_peptides.tsv";
    my $filtered_file = qq"${output_directory}/${string_name}_peptides.tsv";

    my $in = FileHandle->new("less $options->{input} |");
    my $all_out = FileHandle->new(">${all_file}");
    my $filtered_out = FileHandle->new(">${filtered_file}");
    my $seqio = Bio::SeqIO->new(-format => 'genbank', -fh => $in);

  SEQUENCES: while (my $seq = $seqio->next_seq) {
      my $seqid = $seq->id;
      my @feature_list = $seq->get_SeqFeatures();
    FEATURES: for my $feat (@feature_list) {
        my $noted = $feat->has_tag($options->{note_tag});
        my $locus = $feat->has_tag($options->{gff_tag});
        next FEATURES unless (defined($locus));
        my @name_array = $feat->get_tag_values($options->{gff_tag});
        my @notes;
        if (!defined($noted)) {
            @notes = @name_array;
        } else {
            my @notes = $feat->get_tag_values($options->{note_tag});
        }
        my $name = $name_array[0];
        my $uc_string = ucfirst($options->{tag_string});
        my $found = grep(/$options->{tag_string}|$uc_string/, @notes);
        if ($found) {
            print $filtered_out qq"${name}\t@notes\n";
        }
        print $all_out qq"${name}\t@notes\n";
    }
  }
    $filtered_out->close();
    $all_out->close();
    $in->close();
}

sub Filter_Edge_Features {
    my %args = @_;
    my @features = @{$args{features}};
    my $source = $args{source};
    my @new_features;
    my $count = 0;
    for my $f (@features) {
        $count++;
        my $contig_id = $f->seq_id;
        my $contig = $source->{$contig_id};
        my $good_feature = Query_Edge_Feature(feature => $f, contig => $contig);
        if ($good_feature) {
            ## print "Keeping this feature!\n";
            push(@new_features, $f);
        } else {
            print "Filter_Edge_Feature() Dropping bad feature, number ${count}.\n";
        }
    }
    return(\@new_features);
}

=head2 C<Merge_CDS_Predictions>

 Master function used to merge multiple CDS annotation sources.

 This should set up a job which is able start a job capable of merging
 CDS/ORF predictions from prokka/prodigal and glimmer.  Note that
 prokka uses prodigal, but by having a separate merge for prodigal, it
 becomes possible to pull in the confidence estimate provided by
 prodigal as well as use a separate trainer or different parameters.

 This requires input genbank file from prokka, the gff features from
 prodigal, and the prediction file from glimmer3.

=over

=item C<Arguments>

 input(required): Prokka output annotations (e.g. untrained prodigal).
 input_glimmer(''): Output file from a glimmer invocation.
 input_phanotate(''): Output file from a phanotate invocation.
 input_prodigal(''): Output file from a (trained) prodigal invocation.
 primary_key('locus_tag'): Set the primary key for the hash of annotations.
 jmem(8): Expected memory usage.
 jprefix('19'): Prefix for the job name and output directories.


=cut
sub Merge_CDS_Predictions {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        input_glimmer => '',
        input_phanotate => '',
        input_prodigal => '',
        primary_key => 'locus_tag',
        jmem => 8,
        jprefix => '19',
        modules => ['ncbi_tools/6.1']);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    ## Even though it is a bit redundant, record all the output filenames now
    ## so that they are easily accessible for the various downstream search tools.
    my $output_dir =  qq"outputs/$options->{jprefix}merge_cds_predictions";
    my $output_basename = basename($options->{input}, ('.gbk'));
    ## The 'raw' genome
    my $output_genome = qq"${output_dir}/${output_basename}.fna";
    ## The full assembly
    my $output_fsa = qq"${output_dir}/${output_basename}.fsa";
    ## CDS sequences
    my $output_cds = qq"${output_dir}/${output_basename}.ffn";
    ## Amino acid sequences
    my $output_faa = qq"${output_dir}/${output_basename}.faa";
    ## Feature file
    my $output_gff = qq"${output_dir}/${output_basename}.gff";
    ## tbl file (e.g. input for tbl2asn
    my $output_tbl = qq"${output_dir}/${output_basename}.tbl";
    ## initial genbank file
    my $output_gbf = qq"${output_dir}/${output_basename}.gbf";
    ## final genbank file
    my $output_gbk = qq"${output_dir}/${output_basename}.gbk";
    ## A tsv containing arbitrarily chosen annotation information.
    my $output_tsv = qq"${output_dir}/${output_basename}.tsv";
    ## After we do the various annotation searches, we will add
    ## the metadata template file (sbt), the tsv, and xlsx outputs.
    my $output_log = qq"${output_dir}/${output_basename}_runlog.txt";
    my $comment = '## This will hopefully merge CDS predictions from prodigal/glimmer/prokka.';

    my $jstring = qq?
use Bio::Adventure::Annotation;
my \$result = \$h->Bio::Adventure::Annotation_Genbank::Merge_CDS_Predictions_Worker(
  input => '$options->{input}',
  input_glimmer => '$options->{input_glimmer}',
  input_phanotate => '$options->{input_phanotate}',
  input_prodigal => '$options->{input_prodigal}',
  jprefix => '$options->{jprefix}',
  jname => 'merge_orfs',
  output_dir => '${output_dir}',
  output_basename => '${output_basename}',
  output_genome => '${output_genome}', ## .fna, initial genome
  output_fsa => '${output_fsa}', ## .fsa, assembly
  output_cds => '${output_cds}', ## .ffn, CDS nucleotide sequences by locus tag.
  output_faa => '${output_faa}', ## .faa, CDS amino acids by locus tag.
  output_gff => '${output_gff}', ## .gff, merged CDS predictions.
  output_tbl => '${output_tbl}', ## .tbl, input for tbl2asn
  output_gbf => '${output_gbf}', ## .gbf, initial genbank file.
  output_gbk => '${output_gbk}', ## .gbk, cleaned genbank file.
  output_tsv => '${output_tsv}', ## .tsv
  output_log => '${output_log}', ## Run log
  primary_key => '$options->{primary_key}',);
?;
    my $merge_orfs = $class->Submit(
        comment => $comment,
        input => $options->{input},
        input_glimmer => $options->{input_glimmer},
        input_phanotate => $options->{input_phanotate},
        input_prodigal => $options->{input_prodigal},
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
        jname => 'merge_orfs',
        jprefix => $options->{jprefix},
        jstring => $jstring,
        language => 'perl',
        output_dir => $output_dir,
        output_basename => $output_basename,
        output_genome => $output_genome,
        output_fsa => $output_fsa,
        output_cds => $output_cds,
        output_faa => $output_faa,
        output_gff => $output_gff,
        output_tbl => $output_tbl,
        output_gbf => $output_gbf,
        output_gbk => $output_gbk,
        output_tsv => $output_tsv,
        output_log => $output_log,
        primary_key => $options->{primary_key},);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($merge_orfs);
}

=back

=head2 C<Merge_CDS_Predictions_Worker>

 This does the actual work spawned by Merge_CDS_Predictions().

 In its current form it does the following:
 1.  Reads the phanotate tsv file into an array of SeqFeatures.
 2.  Reads the prokka gbk file into another array and splits it up
     by type (source/other/cds/etc)
 3.  Reads the prodigal GFF file into an array of SeqFeatures.
     (This may seem redundant with prokka, but I promise it is not,
     partially because it pulls in the neat Shine-Dalgarno information,
     and partially because this separate run is trained a little differently).
 4.  Reads the Glimmer predict file into yet another array of Features.
 5.  Runs the function Combine_CDS_Features() multiple times, each time
     attempting to add a few more features and/or merge in the interesting
     stuff.  This is a place where some logic should be added to allow one to
     choose arbitrary sets of features rather than the currently hard-coded logic
     to state that glimmer<prokka<prodigal<phanotate.
 6.  Pull the non-CDS and assembly features from prokka into the array.
 7.  Cleanup the names of the features so they are consistent across methods.
 8.  Write out the various output files of interest, including:  .fsa of the
     assembly, .ffn of the nucleotide coding sequences, .faa of the amino
     acid sequences, .gff of the features, .tbl of the features in the format
     expected by tbl2asn(1).
 9.  Actually run tbl2asn in order to generate genbank flat files.

 Currently it does not return anything useful, which is dumb.

=over

=item C<Arguments>

 I don't feel like typing this right now, the important ones are in Merge_CDS_Predictions().

=cut
sub Merge_CDS_Predictions_Worker {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        input_glimmer => '',
        input_phanotate => '',
        input_prodigal => '',
        output_dir => '.',
        output_basename => 'assembly',
        output_genome => 'assembly.fna',
        output_fsa => 'assembly.fsa',
        output_cds => 'assembly.ffn',
        output_faa => 'assembly.faa',
        output_gff => 'assembly.gff',
        output_tbl => 'assembly.tbl',
        output_gbf => 'assembly.gbf',
        output_gbk => 'assembly.gbk',
        output_tsv => 'assembly.tsv',
        jprefix => '19',);

    ## Start gathering features from each of the annotation methods:
    ## Phanotate first.
    print "Reading phanotate seqfeatures.\n";
    my $phanotate_features_ref = Read_Phanotate_to_SeqFeatures(input => $options->{input_phanotate});
    my @phanotate_features = @{$phanotate_features_ref};
    ## Follow that up with prokka, keep in mind that prokka provides
    ## features of multiple types, the genes, cds, source, and others.
    print "Reading prokka seqfeatures.\n";
    my ($pfeat, $pseq, $seqids) = Read_Prokka_Gbk_to_SeqFeatures(input => $options->{input});
    my @prokka_sequences = @{$pseq};
    my @sequence_ids = @{$seqids};
    my $prokka_features = Separate_Prokka_Features(input => $pfeat);
    my @prokka_cds = @{$prokka_features->{cds}};
    my %source_features = %{$prokka_features->{source}};
    my %gene_features = %{$prokka_features->{gene}};
    my %other_features = %{$prokka_features->{other}};
    ## Add the prodigal features next.
    my $prodfeat = Read_Prodigal_GFF_to_SeqFeatures(input => $options->{input_prodigal});
    my @prodigal_features = @{$prodfeat};
    ## Read the glimmer predict file, transform it into a feature list, and return that.
    my $glimmer_features_ref = Predict_to_Features(in => $options->{input_glimmer});
    my @glimmer_features = @{$glimmer_features_ref};

    ## Add a filter for bad Edge features
    my $phanotate_filtered = Filter_Edge_Features(features => \@phanotate_features,
                                                  source => \%source_features);
    ##my $prokka_cds_filtered = Filter_Edge_Features(features => \@prokka_cds,
    ##                                               source => \%source_features);
    my $glimmer_filtered = Filter_Edge_Features(features => \@glimmer_features,
                                                source => \%source_features);
    my $prodigal_filtered = Filter_Edge_Features(features => \@prodigal_features,
                                                 source => \%source_features);

    my $first_merged = Combine_CDS_Features(first => $phanotate_filtered, second => \@prokka_cds);
    my $second_merged = Combine_CDS_Features(first => $first_merged, second => $glimmer_filtered,);
    my $final_merged_ref = Combine_CDS_Features(first => $second_merged, second => $prodigal_filtered,
                                                notes => 'both');
    my @final_merged = @{$final_merged_ref};

    ## Now we should have a data structure containing:
    ##  1.  A hash of the sources (e.g. the contigs)
    ##  2.  A hash of genes from prokka, which I may discard and rewrite.
    ##  3.  A hash of other features from prokka, which will be good to keep.
    ##  4.  A hash of cds from the merge of prokka and glimmer.

    ## Add the non CDS features
    my $final_length = scalar(@final_merged);
    for my $other (keys %other_features) {
        my $other_feature = $other_features{$other};
        push(@final_merged, $other_feature);
    }
    $final_length = scalar(@final_merged);

    ## Pull the source features into an array.
    my @assembly_features = ();
    for my $contig (keys %source_features) {
        my $source_feature = $source_features{$contig};
        ## Put the source at the front.
        push(@assembly_features, $source_feature);
    }

    ## Lets make the CDS names just be the basename of the sample
    ## followed by an incrementer.

    ## There is one important corner case that needs to be addressed.
    ## Sometimes glimmer (and maybe prodigal?) picks up a potential ORF
    ## which crosses the 0 mark of the assembly.  Since this is used
    ## primarily with circular genomes, that is perfectly reasonable.
    ## But I am not certain how to properly handle it.
    ## It seems like this is the best place to handle that case though,
    ## since at this point we have the full catalog of features in one place.
    my $prefix = basename($options->{basedir});
    my @renamed = Rename_Features(
        features => \@final_merged, assembly => \@assembly_features,
        prefix => $prefix, add_genes => 1);
    ## Write the various output files.

    my $fsa_written = Write_Fsa_from_SeqFeature_Generic(
        input_seq => \@assembly_features,
        output_fsa => $options->{output_fsa});
    my $cds_written = Write_CDS_from_SeqFeatures(
        input_features => \@renamed,
        input_seq => $fsa_written,
        output_cds => $options->{output_cds},
        output_faa => $options->{output_faa},
        output_tsv => $options->{output_tsv});
    my $gff_writen = Write_Gff_from_SeqFeatures(
        input_features => \@renamed,
        output_gff => $options->{output_gff});
    my $tbl_written = Write_Tbl_from_SeqFeatures(
        tbl_file => $options->{output_tbl}, features => \@renamed,
        sequences => \@sequence_ids,
        taxonomy_information => undef,);

    my $gbk_written = Write_Gbk(
        input_seq => \@prokka_sequences,
        input_tbl => $options->{output_tbl},
        output_gbk => $options->{output_gbk},
        output_fsa => $options->{output_fsa},);
    ## TODO: Write the other output files.

    ## This should return something useful -- maybe parse the error file from tbl2asn?
}

=back

=head2 C<Predict_to_Features>

 Read the nasty glimmer output format into something sensible.

 This reads the .predict file producted by glimmer and sends the data
 to an array of SeqFeature::Generic objects.

=cut
sub Predict_to_Features {
    my %args = @_;
    my $glimmer = $args{in};
    my $glimmer_in = FileHandle->new("less ${glimmer} |");
    my $contig_id;
    my $fixed_id;
    my @glimmer_features = ();
  ENTRIES: while (my $line = <$glimmer_in>) {
      chomp $line;
      if ($line =~ /^>/) {
          ## Then this tells us the contig
          $contig_id = $line;
          $fixed_id = Remove_Contig_Cruft($contig_id);
          next ENTRIES;
      }
      ## Example prediction
      ## orf00001   109351       25  +2    12.06
      my ($orf_id, $start, $end, $strand_frame, $score) = split(/\s+/, $line);
      my ($str, $frame) = split(//, $strand_frame);
      my $strand;
      if ($str eq '+') {
          $strand = '1';
      } elsif ($str eq '-') {
          $strand = '-1';
      } else {
          $strand = '0';
      }

      my $start_phase = $frame;
      $start_phase =~ s/^.{1}(\d)$/$1/g;
      my $phase = $start_phase - 1;
      my $fivep = $start;
      my $threep = $end;
      if ($strand eq '-1') {
          $fivep = $end;
          $threep = $start;
      }
      my $feature = Bio::SeqFeature::Generic->new(
          -primary_tag => 'CDS',
          -display_name => $orf_id,
          -seq_id => $fixed_id,
          -start => $fivep,
          -end => $threep,
          -strand => $strand,
          -score => $score,
          -frame => $phase,
          -tag => {
              note => 'cds_prediction: glimmer',
              locus_tag => $orf_id,
              transl_table => 11, });
      $feature->strand($strand);
      push(@glimmer_features, $feature);
  } ## End iterating over every line of the glimmer3 output.
    $glimmer_in->close();
    return(\@glimmer_features);
}

sub Query_Edge_Feature {
    my %args = @_;
    my $feature = $args{feature};
    my $contig = $args{contig};
    ## This is a test for weirdo features from phanotate/prodigal.
    ## I think the answer lies in a piece of the scoring provided by prodigal:
    ## Here is an annotation line from progial:
    ## >gnl|Prokka|EPr2_1_1 # 1 # 117 # 1 # ID=1_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.504
    ## I think therefore the easiest solution is to query if the start is 1 and the amino acid is not 'M'.
    my $start = $feature->start;
    my $end = $feature->end;
    my $strand = $feature->strand;
    my $contig_sequence = $contig->seq;
    my $contig_length = length($contig_sequence);
    ## Another peculiar edge case:
    ## Sometimes a feature is picked up which is a little before the origin and continues after the origin.
    ## e.g. the strand may be +1, the start is ~ 100nt before 0 and the end is ~ 100nt after 0.
    my $straddles = 0;
    $straddles = 1 if ($start > $end);
    if ($straddles) {
        print "This feature straddles the origin east and is problematic for genbank.\n";
        return(0);
    }

    my $feature_seq = $contig_sequence->trunc($start, $end);
    $feature_seq = $feature_seq->revcom if ($strand < 0);
    my $amino_acids = $feature_seq->translate->seq;
    my $first_amino = substr($amino_acids, 0, 1);
    my $edge = 0;
    $edge = 1 if ($start == 1 && $strand > 0);
    $edge = 1 if ($end == $contig_length && $strand < 0);

    my $ret = 1;
    my $not_start_codon = 1;
    $not_start_codon = 0 if ($first_amino eq 'M' or $first_amino eq 'L');
    ## print "The start is $start and the first amino acid is $first_amino\n";
    if ($edge && $not_start_codon) {
        print "This feature starts at 1 and not with 'M|L'! and should be skipped\n";
        $ret = 0;
    }
    return($ret);
}

=head2 C<Read_Phanotate_to_SeqFeatures>

 Read Phanotate output into a set of features.

 Either this or Predict_to_Features() should probably be renamed,
 because they do much the same; except of course this reads
 the phanotate tsv output file and sends it to an array of
 SeqFeatures.

=cut
sub Read_Phanotate_to_SeqFeatures {
    my %args = @_;
    my $phanotate = $args{input};
    my $phanotate_in = FileHandle->new("less ${phanotate} |");
    my $contig_id;
    my $fixed_id;
    my @phanotate_features = ();
    my $count = 0;
  ENTRIES: while (my $line = <$phanotate_in>) {
      chomp $line;
      next ENTRIES if ($line =~ /^#/);
      $count++;
      my $orf_number = sprintf("%04d", $count);
      my ($start, $end, $frame, $contig, $score) = split(/\t/, $line);
      my $fixed_id = Remove_Contig_Cruft($contig);
      my $orf_id = qq"${fixed_id}_${orf_number}";
      my $strand = '+';
      if ($frame eq '+') {
          $strand = '1';
      } elsif ($frame eq '-') {
          $strand = '-1';
      } else {
          $strand = '0';
      }
      my $fivep = $start;
      my $threep = $end;
      if ($end < $start) {
          $fivep = $end;
          $threep = $start;
      }
      $score = FormatSigFigs($score, 3);

      my $feature = Bio::SeqFeature::Generic->new(
          -primary_tag => 'CDS',
          -display_name => $orf_id,
          -seq_id => $fixed_id,
          -start => $fivep,
          -end => $threep,
          -strand => $strand,
          -score => $score,
          -tag => {
              note => qq"cds_prediction: phanotate, score: ${score}",
              locus_tag => $orf_id,
              transl_table => 11, });
      $feature->strand($strand);
      push(@phanotate_features, $feature);
  } ## End iterating over every line of the glimmer3 output.
    $phanotate_in->close();
    return(\@phanotate_features);
}

=head2 C<Read_Prokka_Gbk_to_SeqFeatures>

 Dump prokka genbank output into a seqfeature array.

 Read the genbank output from prokka into a set of SeqFeature's.
 This will be the basis for quite a few of the downstream methods
 along the way to creating a merged set of annotations.

=cut
sub Read_Prokka_Gbk_to_SeqFeatures {
    my %args = @_;
    my $prokka_in = FileHandle->new("less $args{input} |");
    my $seqio = Bio::SeqIO->new(-format => 'genbank', -fh => $prokka_in);
    my (@prokka_sequences, @prokka_features, @sequence_ids);
    while (my $seq = $seqio->next_seq) {
        push(@prokka_sequences, $seq);
        my $seq_id = $seq->id;
        push(@sequence_ids, $seq_id);
        my @feature_list = $seq->get_SeqFeatures();
        for my $f (@feature_list) {
            push(@prokka_features, $f);
        }
    } ## End grabbing prokka sequences
    $prokka_in->close();
    return(\@prokka_features, \@prokka_sequences, \@sequence_ids);
}

=head2 C<Read_Prodigal_GFF_to_SeqFeatures>

 Dump prodigal gff output into a seqfeature array.

 Read the gff output from prodigal into an array of SeqFeature's.
 Collect the Shine-Dalgarno information and scores while at it.

=cut
sub Read_Prodigal_GFF_to_SeqFeatures {
    my %args = @_;
    my $prodigal_in = FileHandle->new("less $args{input} |");
    my $prodigal_io = Bio::FeatureIO->new(-format => 'gff', -fh => $prodigal_in);
    my @prodigal_features = ();
    while (my $feature = $prodigal_io->next_feature()) {
        my $contig_id = $feature->seq_id;
        my $fixed_id = Remove_Contig_Cruft($contig_id);
        my $reset_contig = $feature->seq_id($fixed_id);

        ## Pull out the scoring information
        my $spacer = '';
        my $rbs = '';
        my $uscore = '';
        my $tscore = '';
        my $sscore = '';
        my $rscore = '';
        my $cscore = '';
        my $start_type = '';
        my @spacers = $feature->get_tag_values('rbs_spacer');
        $feature->remove_tag('rbs_spacer');
        $spacer = $spacers[0] if ($spacers[0]);
        my @rbses = $feature->get_tag_values('rbs_motif');
        $feature->remove_tag('rbs_motif');
        $rbs = $rbses[0] if ($rbses[0]);
        my @uscores = $feature->get_tag_values('uscore');
        $feature->remove_tag('uscore');
        $uscore = $uscores[0] if ($uscores[0]);
        my @tscores = $feature->get_tag_values('tscore');
        $feature->remove_tag('tscore');
        $tscore = $tscores[0] if ($tscores[0]);
        my @sscores = $feature->get_tag_values('sscore');
        $feature->remove_tag('sscore');
        $sscore = $sscores[0] if ($sscores[0]);
        my @rscores = $feature->get_tag_values('rscore');
        $feature->remove_tag('rscore');
        $rscore = $rscores[0] if ($rscores[0]);
        my @cscores = $feature->get_tag_values('cscore');
        $feature->remove_tag('cscore');
        $cscore = $cscores[0] if ($cscores[0]);
        my @start_types = $feature->get_tag_values('start_type');
        $feature->remove_tag('start_type');
        $start_type = $start_types[0] if ($start_types[0]);

        my $score_string = qq"cds_prediction: trained prodigal(rbs:${rbs} spacer:${spacer} start:${start_type} c:${cscore} r:${rscore} s:${sscore} t:${tscore} u:${uscore})";
        $feature->add_tag_value('note', $score_string);
        push(@prodigal_features, $feature);
    }
    $prodigal_in->close();
    return(\@prodigal_features);
}

=head2 C<Remove_Contig_Cruft>

 Standardize contig IDs upon merging them.

 This takes the various contig IDs producted by prokka and friends and
 attempts to simplify/standardize them so that downstream reader methods
 will not get confused.

=cut
sub Remove_Contig_Cruft {
    my $contig_id = shift;
    ## Take into account the possible contig lines
    ## >gnl|Prokka|EPr2_1 [gcode=11] [organism=Phage species] [strain=strain]
    ## vs.
    ## >gnl|Prokka|EPr2_1
    ## vs
    ## >EPr2_1
    ## Get rid of the ^>
    $contig_id =~ s/^>//g;
    if ($contig_id =~ /\[/) {  ## Then it is the long annoying one.
        ## And the stuff after the whitespace
        $contig_id =~ s/^(\S+)?\s+.*$/$1/g;
        ## Then pull out anything left after the last |
        $contig_id =~ s/^(.*\|)(\w+)\s*$/$2/g;
    } elsif ($contig_id =~ /\|/) {
        ## grab anything after the last |
        $contig_id =~ s/^(.*\|)?(\S+)\s*$/$2/g;
    } else {
        ## Just drop the >
        $contig_id =~ s/^>//g;
    }
    return($contig_id);
}

=head2 C<Rename_Features>

 Make a consistent set of feature names after merging multiple data sources.

 This function takes a big pile of features from different data sources
 (prokka/prodigal/glimmer/phanotate/aragorn/tnrascan) and attempts to give
 them a single set of canonical names.

 One detail worth noting, I decided to _not_ begin numbering again when
 going from one contig to another.  In addition, tRNA features have
 a separate number from CDS.

=cut
sub Rename_Features {
    my %args = @_;
    my $prefix = $args{prefix};
    my @features = @{$args{features}};
    my @assembly = @{$args{assembly}};
    my $count = 0;
    my $contig_id;

    ## Set aside a feature in case we need to move one to the end.
    my $terminal_feature = undef;

    my %keyed_assembly = ();
    for my $a (@assembly) {
        my $key = $a->seq_id;
        $keyed_assembly{$key} = $a;
    }

    my %orfs_by_contig = ();
    my @renamed = ();
  RENAME: for my $n (@features) {
      $count++;
      if ($n->primary_tag eq 'source') {
          next RENAME;
      }
      my $display_number = sprintf("%04d", $count);
      $contig_id = $n->seq_id();
      my $display_name = qq"${prefix}_${display_number}";
      my $protein_name = qq"Phokka:${display_name}";
      my $start = $n->start;
      my $end = $n->end;
      my $strand = $n->strand;
      $orfs_by_contig{$display_name} = $contig_id;
      ## Handle the weird corner case of features which bridge the origin here:
      ## I think these are happening only when phageterm reorients the genome to a
      ## DTR and prodigal/glimmer picks up an ORF which ends at the beginning of
      ## the DTR.
      ## With that in mind, I am going to copy the logic (maybe move it) from
      ## Write_CDS_from_SeqFeatures() which handles this case, but end the feature
      ## at the origin.
      if ($end < $start) {
          ## I think any features which fall in this category should end very close to position 1.
          ## In theory, probably significantly less than 100.
          $count--; ## Reset the ORF counter.
          if ($end < 1000) {
              my $source_feature = $keyed_assembly{$n->seq_id};
              ## Assume this is caused by phageterm-reorganization.
              my $moved_end = $n->end($source_feature->end);
              my $note = qq"This CDS is actually ${start} nt. longer, but was truncated because it crossed the origin.  The primary hypothesis for this is due to phageterm-reorganization of the genome, resulting in a DTR which begins at the end of this CDS.";
              $n->add_tag_value('note', $note);
              $terminal_feature = $n;
          } else {
              print "I am not sure how to handle this feature, start:${start} end:${end}\n";
          }
          next RENAME;
      }

      my $id = $n->display_name($display_name);
      $n->remove_tag('locus_tag');
      $n->add_tag_value('locus_tag', $display_name);
      if ($n->has_tag('protein_id')) {
          $n->remove_tag('protein_id');
      }
      $n->add_tag_value('protein_id', $protein_name);
      push(@renamed, $n);
  } ## End of the rename loop.

    $count++;
    ## Final step in addressing phageterm-annoying features,
    ## Pick up the modified feature and stick it on the end of the list.
    if (defined($terminal_feature)) {
        my $display_number = sprintf("%04d", $count);
        my $contig_id = $terminal_feature->seq_id();
        my $display_name = qq"${prefix}_${display_number}";
        my $protein_name = qq"Phokka:${display_name}";
        my $start = $terminal_feature->start;
        my $end = $terminal_feature->end;
        my $strand = $terminal_feature->strand;
        my $id = $terminal_feature->display_name($display_name);
        $terminal_feature->remove_tag('locus_tag');
        $terminal_feature->add_tag_value('locus_tag', $display_name);
        if ($terminal_feature->has_tag('protein_id')) {
            $terminal_feature->remove_tag('protein_id');
        }
        $terminal_feature->add_tag_value('protein_id', $protein_name);
        push(@renamed, $terminal_feature);
    }

    ## Finally, re-iterate over the features and add genes if requested.
    my @final_features = ();
    if ($args{add_genes}) {
      FEATLOOP: for my $f (@renamed) {
          if ($f->primary_tag() eq 'source') {
              push(@final_features, $f);
              next FEATLOOP;
          }

          my $this_id = $f->display_name;
          my $this_contig_id = $orfs_by_contig{$this_id};
          my $gene = Bio::SeqFeature::Generic->new(
              -primary_tag => 'gene',
              -seq_id => $this_contig_id,
              -display_name => $f->display_name,
              -start => $f->start,
              -strand => $f->strand,
              -end => $f->end,
              -tag => { locus_tag => $f->display_name },);
          push(@final_features, $gene);
          push(@final_features, $f);
      } ## Finished iterating over the features and adding genes.
    }

    my $len = scalar(@final_features);
    return(@final_features);
}

=head2 C<Separate_Prokka_Features>

 Separate prokka feature types into multiple data structures.

 Given a pile of features extracted from a prokka annotation, pull them apart
 so that the assembly, genes, cds, etc may be addressed separately.

=cut
sub Separate_Prokka_Features {
    my %args = @_;
    my @prokka = @{$args{input}};
    my $source_remove = 'gnl\|Prokka\|';
    if (defined($args{source_remove})) {
        $source_remove = $args{source_remove};
    }
    ## new_features will be the result of the merge.
    my %source_features = ();
    my %gene_features = ();
    my %other_features = ();
    my @cds_features = ();
  PROKKA: for my $prokka_f (@prokka) {
      my $prokka_tag = $prokka_f->primary_tag();
      my $prokka_name = $prokka_f->display_name();
      my $new_name = $prokka_name;
      $new_name =~ s/$source_remove//g;
      my $display_renamed = $prokka_f->display_name($new_name);
      my $prokka_source_id = $prokka_f->seq_id;
      if ($prokka_tag eq 'source') {
          $source_features{$prokka_source_id} = $prokka_f;
      } elsif ($prokka_tag eq 'gene') {
          $gene_features{$prokka_name} = $prokka_f;
      } elsif ($prokka_tag eq 'CDS') {
          $prokka_f->add_tag_value('note', 'cds_prediction: prodigal via prokka');
          push(@cds_features, $prokka_f);
      } else {
          $other_features{$prokka_name} = $prokka_f;
      }
  }  ## End initial iteration over prokka features.

    my %ret = (
        source => \%source_features,
        cds => \@cds_features,
        gene => \%gene_features,
        other => \%other_features,);
    return(\%ret);
}

=head2 C<Write_CDS_from_SeqFeature>

 What it says on the tin!

 Given a pile of SeqFeatures, write out the individual sequences
 from it.  This will ideally write out a .faa of amino acids,
 .ffn of the CDS nucleotides, and a .tsv summarizing them.

=cut
sub Write_CDS_from_SeqFeatures {
    my %args = @_;
    ## The full sequences
    my @sequences = @{$args{input_seq}};
    my %named_seq = ();
    my $translated = 0;
    $translated = $args{translated} if (defined($args{translated}));
    my @features = @{$args{input_features}};
    my $out_cds = Bio::SeqIO->new(-file => qq">$args{output_cds}",
                                  -format => 'Fasta');
    my $out_faa = Bio::SeqIO->new(-file => qq">$args{output_faa}",
                                  -format => 'Fasta');
    my $out_tsv = FileHandle->new(">$args{output_tsv}");
    my $tsv_header = qq"locus_tag\tcontig\ttype\tsource\tstart\tend\tstrand\tcds_prediciton\taa_sequence\n";
    print $out_tsv $tsv_header;

    ## Make a hash of the contigs by name
  CONTIGS: for my $s (@sequences) {
      my $id = $s->display_name;
      $named_seq{$id} = $s;
  }

    my $written = 0;
    my $count = 0;
    my $full_sequence;
    my $maximum_length = 0;
  FEATURES: for my $in (@features) {
      $count++;
      my $type = $in->primary_tag;
      next FEATURES if ($type eq 'gene');
      my $source = $in->source_tag;
      my $contig = $in->seq_id;
      my $name = $in->display_name;
      my $start = $in->start;
      my $end = $in->end;
      my $strand = $in->strand;
      my $full_sequence = $named_seq{$contig};
      ## This is hopefully only happening when a phageterm-reorganized
      ## genome is putting an ORF at the 'top of the clock'
      ## My solution therefore will be to check that it is bound at the
      ## new beginning.
      my $cds_obj = '';
      if ($end < $start) {
          ## I am going to assume that any real orf in this context must
          ## be no more than 5000 before and/or 5000 after the 0 point, thus
          ## a 6k ORF ending at position 100 will fail this test.
          my $post_dist = $end;
          my $pre_dist = $full_sequence->length - $start;
          if ($pre_dist < 5000 && $post_dist < 5000) {
              print "Found an overlap with 12 on the clock.\n";
              my $pre_sequence = $full_sequence->trunc($end, $pre_dist);
              my $post_sequence = $full_sequence->trunc(1, $start);
              if ($strand < 0) {
                  $pre_sequence = $pre_sequence->revcom;
                  $post_sequence = $post_sequence->revcom;
              }
              my $tmp_sequence_string = $pre_sequence->seq . $post_sequence->seq;
              $cds_obj = Bio::Seq->new(-display_id => $name, -seq => $tmp_sequence_string);
          } else {
              print "Cannot deal with this sequence right now, start:${start}, end:${end}\n";
              next FEATURES;
          }
      } else { ## Normal sequence where start < end
          $cds_obj = $full_sequence->trunc($start, $end);
          if ($strand < 0) {
              $cds_obj = $cds_obj->revcom;
          }
      }
      my $aa_obj = $cds_obj->translate;
      my $cds_sequence_string = $cds_obj->seq;
      my $aa_sequence_string = $aa_obj->seq;
      $aa_sequence_string =~ s/\*$//g;

      my $cds_seq_obj = Bio::Seq->new(-display_id => $name, -seq => $cds_sequence_string);
      $out_cds->write_seq($cds_seq_obj);
      my $aa_seq_obj = Bio::Seq->new(-display_id => $name, -seq => $aa_sequence_string);
      $out_faa->write_seq($aa_seq_obj);
      $written++;

      ## Now figure out what information to write to the tsv file.
      my @notes;
      my $note;
      if ($in->has_tag('note')) {
          @notes = $in->get_tag_values('note');
          $note = $notes[0];
          $note =~ s/^cds_prediction:\s+//g;
      }

      my $tsv_line = qq"${name}\t${contig}\t${type}\t${source}\t${start}\t${end}\t${strand}\t${note}\t${aa_sequence_string}\n";
      print $out_tsv $tsv_line;
  }
    $out_tsv->close();
    return($written);
}

=head2 C<Write_Fsa_from_SeqFeature_Generic>

 Also what it says on the tin!

 Given an array of SeqFeature::Generic's, write out a .fsa
 assembly file.

=cut
sub Write_Fsa_from_SeqFeature_Generic {
    my %args = @_;
    my $id_prefix = 'gnl|Phokka|';
    if (defined($args{id_prefix})) {
        $id_prefix = $args{id_prefix};
    }
    my @in = @{$args{input_seq}};
    my $out = Bio::SeqIO->new(-file => qq">$args{output_fsa}",
                              -format => 'Fasta');
    my $written = 0;
    my @seq_objects = ();
    for my $s (@in) {
        $written++;
        my $id = $s->seq_id;
        my $entire_sequence = $s->entire_seq();
        my $sequence_string = $entire_sequence->seq();
        my $seq_obj = Bio::Seq->new(-display_id => $id,
                                    -seq => $sequence_string);
        $out->write_seq($seq_obj);
        push(@seq_objects, $seq_obj);
    }
    return(\@seq_objects);
}

=head2 C<Write_Gbk>

 This mostly runs tbl2asn.

 This function is a little bit of a misnomer, it does not actually
 write a Gbk file, but instead takes a .tbl file from Write_Tbl()
 and runs tbl2asn(1) on it; invoking a little logic on the way
 to hopefully make sure that the final result is sane.

=cut
sub Write_Gbk {
    my %args = @_;
    $args{contigs} = 1 if (!defined($args{contigs}));
    my $run_sed = 1;
    $run_sed = $args{sed} if (defined($args{sed}));

    my $tbl2asn_m_option = '-M n';
    if ($args{contigs} > 10_000) {
        $tbl2asn_m_option = '-M b';
    }
    my $log_fh = *STDOUT;
    $log_fh = $args{log_fh} if (defined($args{log_fh}));
    my $outdir = '.';
    my $accver = '1';
    my $EXE = 'slightly modified prokka';
    my $VERSION = '1.14.6';
    my $URL = 'https://github.com/tseemann/prokka';
    print $log_fh "Inputs for tbl2asn: $args{output_fsa}, $args{input_tbl}";
    my $sbt_arg = '';
    if (defined($args{sbt})) {
        print $log_fh " $args{sbt}\n";
        $sbt_arg = qq"-t $args{sbt} ";
    } else {
        print $log_fh "\n";
    }
    print $log_fh qq"Outputs from tbl2asn: $args{output_gbk}.\n";
    my $tbl2asn_comment = qq"Annotated using ${EXE} ${VERSION} from ${URL}";
    if (defined($args{taxonomy})) {
        my $taxonomy_information = $args{taxonomy};
        $tbl2asn_comment .= qq"; Most similar taxon to this strain: $taxonomy_information->{taxon}";
    }
    my $default_args = '-V b -c f -S F -a r10k -l paired-ends ';
    if (defined($args{tbl2asn_args})) {
        $default_args = $args{tbl2asn_args};
    }
    my $out_basedir = dirname($args{output_gbk});
    my $out_basefile = basename($args{output_gbk}, ('.gbk'));
    my $stderr_file = qq"${out_basedir}/${out_basefile}.stderr";
    my $error_file = qq"${out_basedir}/${out_basefile}.err";
    my $stdout_file = qq"${out_basedir}/${out_basefile}.stdout";
    my $tbl_command = qq"tbl2asn ${default_args} ${tbl2asn_m_option} \\
  -N ${accver} -y '${tbl2asn_comment}' \\
  -Z ${error_file} ${sbt_arg} -i $args{output_fsa} \\
  1>${stdout_file} 2>${stderr_file}
";
    print $log_fh qq"Running ${tbl_command}\n";
    my $tbl2asn_result = qx"${tbl_command}";
    my $sed_result = undef;
    if ($run_sed && -r qq"${out_basedir}/${out_basefile}.gbf") {
        my $sed_command = qq"sed 's/COORDINATES: profile/COORDINATES:profile/' \\
  ${out_basedir}/${out_basefile}.gbf | \\
  sed 's/product=\"_/product=\"/g' > $args{output_gbk}";
        print $log_fh qq"Running ${sed_command}\n";
        $sed_result = qx"${sed_command}";
    } else {
        print $log_fh qq"Not running sed, copying file.\n";
        $sed_result = cp(qq"${out_basedir}/${out_basefile}.gbf", $args{output_gbk});
    }
    my $results = [$tbl2asn_result, $sed_result];
    return($results);
}

=head2 C<Write_Gff_from_SeqFeatures>

 Given an array of SeqFeatures, write out a hopefully useful .gff file.

=cut
sub Write_Gff_from_SeqFeatures {
    my %args = @_;
    my @in = @{$args{input_features}};
    use Bio::Tools::GFF;
    my $gff_out = FileHandle->new(">$args{output_gff}");
    my $gffio = Bio::Tools::GFF->new(-noparse => 1, -gff_version => 3);
    my $written = 0;
    for my $s (@in) {
        ## Standardize the feature set
        my $product = 'hypothetical protein';
        if ($s->has_tag('product')) {
            my @products = $s->get_tag_values('product');
            $product = $products[0];
        }
        my $inference = 'ab initio prediction';
        if ($s->has_tag('inference')) {
            my @inferences = $s->get_tag_values('inference');
            $inference = $inferences[0];
        }

        my $standardized = Bio::SeqFeature::Generic->new(
            -primary_tag => $s->primary_tag,
            -display_name => $s->display_name,
            -seq_id => $s->seq_id,
            -start => $s->start,
            -end => $s->end,
            -strand => $s->strand,
            -score => $s->score,
            -tag => {
                locus_tag => $s->display_name,
                transl_table => 11,
                inference => $inference,
                product => $product,
            });
        my $gff_string = $gffio->gff_string($standardized);
        print $gff_out "${gff_string}\n";
        $written++;
    }
    $gff_out->close();
    return($written);
}

=head2 C<Write_Tbl_from_SeqFeatures>

 Write a .tbl file from some seqfeatures.

 This is basically just stealing a small piece of code from Torsten
 Seeman's prokka tool which is capable of writing a tbl file in
 preparation for running NCBI's tbl2sqn and generating a genbank file.

=cut
sub Write_Tbl_from_SeqFeatures {
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

    my %contig_to_orf = ();
  FEATUREIDS: for my $j (@features) {
      my $orf_id = $j->display_name;
      my $contig_id = $j->seq_id;
      next FEATUREIDS if (!defined($orf_id));
      next FEATUREIDS if ($orf_id eq '');
      $contig_to_orf{$orf_id} = $contig_id;
  }

  OUTERSEQ: for my $sid (@seq) {
      print STDOUT "Working on feature: ${sid}\n";
      print $tbl_fh ">Feature ${sid}\n";
    INNERFEATURE: for my $f (@features) {
        my $feature_name = $f->display_name;
        ## use Data::Dumper;
        ## print Dumper $f;
        next INNERFEATURE unless ($contig_to_orf{$feature_name} eq $sid);
        print "In ${sid} working on $feature_name\n";
        if ($f->primary_tag eq 'source') {
            $f->strand(1);
            if (defined($taxonomy_information)) {
                my $org_set = $f->add_tag_value('organism', "Phage similar to $taxonomy_information->{taxon}.");
                my $strain_set = $f->add_tag_value('strain', "Similar accession: $taxonomy_information->{hit_accession}.");
                my $id_set = $f->seq_id("Phage species similar to $taxonomy_information->{taxon}.");
            }
        }
        if ($f->primary_tag eq 'CDS' and not $f->has_tag('product')) {
            my $product_set = $f->add_tag_value('product', $hypothetical_string);
        }
        if (my $name = TAG($f, 'gene')) {
            my $name_set = $f->add_tag_value('Name', $name);
        }
        # Make sure we have valid frames/phases (GFF column 8)
        $f->frame($f->primary_tag eq 'CDS' ? 0 : '.');
        if (!defined($f->strand) || !defined($f->start) || !defined($f->end)) {
            my $name = $f->display_name;
            my $strand = $f->strand;
            my $start = $f->start;
            my $end = $f->end;
            print "THERE IS A PROBLEM WITH: ${name}\n";
            print "ONE OF STRAND:${strand} START:${start} OR END:${end} IS UNDEFINED.\n";
        }
        my ($L, $R) = ($f->strand >= 0) ? ($f->start, $f->end) : ($f->end, $f->start);
        print $tbl_fh "${L}\t${R}\t", $f->primary_tag, "\n";
      WRITE_TAGS: for my $tag ($f->get_all_tags) {
          # remove GFF specific tags (start with uppercase letter)
          next WRITE_TAGS if $tag =~ m/^[A-Z]/ and $tag !~ m/EC_number/i;
          my @unwanted = ('gc_cont', 'sscore', 'phase', 'conf', 'tscore', 'partial', 'rbs_spacer', 'rbs_motif', 'cscore', 'start_type', 'uscore', 'seq_id', 'rscore', 'source', 'frame', 'score', 'type');
          next WRITE_TAGS if (grep $_ eq $tag, @unwanted);
          for my $value ($f->get_tag_values($tag)) {
              if (!defined($value)) {
                  print $tbl_fh "\t\t\t${tag}\t\n";
              } else {
                  print $tbl_fh "\t\t\t${tag}\t${value}\n";
              }
          } ## End iterating over the tag values.
      } ## End getting all tags
    } ## End iterating over features
    } ## End iterating over sequence objects.
    return(%contig_to_orf);
}

=head2 C<TAG>

 I just copy/pasted this outright from prokka so that the table writer
 will work properly.

=cut
sub TAG {
    my($f, $tag) = @_;
    # very important to "return undef" here and not just "return"
    # otherwise it becomes non-existent/collapsed in a list
    return undef unless $f->has_tag($tag);
    return ($f->get_tag_values($tag))[0];
}

1;

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

tbl2asn(1)  Bio::SeqFeature::Generic  Bio::Tools::GFF

=cut
