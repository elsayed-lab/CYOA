package Bio::Adventure::Annotation_Genbank;
## LICENSE: gplv2
## ABSTRACT:  Kitty!
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

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
use Template;
use Text::CSV_XS::TSV;

use Bio::FeatureIO;
use Bio::Seq;
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::Lite;
use Bio::SeqIO;
use Bio::Tools::GuessSeqFormat;

=head2 C<Merge_CDS_Predictions>

This should set up a job which is able start a job capable of merging
CDS/ORF predictions from prokka/prodigal and glimmer.  Note that
prokka uses prodigal, but by having a separate merge for prodigal, it
becomes possible to pull in the confidence estimate provided by
prodigal as well as use a separate trainer or different parameters.

This requires input genbank file from prokka, the gff features from
prodigal, and the prediction file from glimmer3.

=cut
sub Merge_CDS_Predictions {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        input_glimmer => '',
        input_prodigal => '',
        jprefix => '19',
        primary_key => 'locus_tag',);

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
    ## After we do the various annotation searches, we will add
    ## the metadata template file (sbt), the tsv, and xlsx outputs.
    my $comment = qq"## This will hopefully merge CDS predictions from prodigal/glimmer/prokka.";

    my $jstring = qq?
use Bio::Adventure::Annotation;
\$h->Bio::Adventure::Annotation_Genbank::Merge_CDS_Predictions_Worker(
  input => '$options->{input}',
  input_glimmer => '$options->{input_glimmer}',
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
  primary_key => '$options->{primary_key}',);
?;
    my $merge_orfs = $class->Submit(
        comment => $comment,
        input => $options->{input},
        input_glimmer => $options->{input_glimmer},
        input_prodigal => $options->{input_prodigal},
        jdepends => $options->{jdepends},
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
        primary_key => $options->{primary_key},
        shell => '/usr/bin/env perl',);
    $class->{language} = 'bash';
    $class->{shell} = '/usr/bin/env bash';
    return($merge_orfs);
}

sub Merge_CDS_Predictions_Worker {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        input_glimmer => '',
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
        jprefix => '19',);

    ## First read in the prokka data as a set of sequences and features
    ## Keep in mind that the features will have the full sequence attached.
    my @prokka_features = ();
    my @prokka_sequences = ();
    my @sequence_ids = ();

    my $prokka_in = FileHandle->new("less $options->{input} |");
    my $seqio = Bio::SeqIO->new(-format => 'genbank', -fh => $prokka_in);
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

    ## Read the prodigal gff as a feature list.
    my @prodigal_features = ();
    my $prodigal_in = FileHandle->new("less $options->{input_prodigal} |");
    my $prodigal_io = Bio::FeatureIO->new(-format => 'gff', -fh => $prodigal_in);
    while (my $feature = $prodigal_io->next_feature()) {
        push(@prodigal_features, $feature);
    }
    $prodigal_in->close();

    ## Read the glimmer predict file, transform it into a feature list, and return that.
    my @glimmer_features = Predict_to_Features(in => $options->{input_glimmer});
    ## This order of operations is important, I am effectively saying that I trust the
    ## run of prodigal that I manually performed more than either prokka or glimmer
    ## because it will override either of them.  Though to be fair, most of the time
    ## the prokka-derived prodigal data is expected to agree (except I did train them
    ## slightly differently and set a couple of options.
    my @prokka_glimmer = Gather_Glimmer_CDS(start => \@prokka_features,
                                            add => \@glimmer_features);
    my @prok_glim_prod = Gather_Prodigal_CDS(start => \@prokka_glimmer,
                                             add => \@prodigal_features);
    my $num = scalar(@prok_glim_prod);
    my @renamed = Rename_Features(features => \@prok_glim_prod, add_genes => 1);

    ## Write the various output files.
    my $fsa_written = Write_Fsa(input_seq => \@prokka_sequences, output_fsa => $options->{output_fsa});
    my $cds_written = Write_CDS(input_seq => \@renamed, output_cds => $options->{output_cds});
    my $faa_written = Write_Faa(input_seq => \@renamed, output_faa => $options->{output_faa});
    my $gff_writen = Write_Gff(input_seq => \@renamed, output_gff => $options->{output_gff});

    my $tbl_written = Write_Tbl(tbl_file => $options->{output_tbl}, features => \@renamed,
                                sequences => \@sequence_ids,
                                taxonomy_information => undef,);
    my $gbk_written = Write_Gbk(
        input_seq => \@prokka_sequences,
        input_tbl => $options->{output_tbl},
        output_gbk => $options->{output_gbk},
        output_fsa => $options->{output_fsa},);

    ## TODO: Write the other output files.
}

## Something in my logic for merging the locus_tag names is incorrect here.
sub Gather_Glimmer_CDS {
    my %args = @_;
    my @prokka = @{$args{start}};
    my @glimmer = @{$args{add}};

    ## new_features will be the result of the merge.
    my @new_features = ();
    my $current_source;
    my $pcount = 0;
    my $orf_number = 1;
  PROKKA: for my $prokka_f (@prokka) {
      $pcount++;
      my $prokka_tag = $prokka_f->primary_tag();
      if ($prokka_tag eq 'source') {
          push(@new_features, $prokka_f);
          $current_source = $prokka_f;
          next PROKKA;
      }
      if ($prokka_tag eq 'gene') {
          next PROKKA;
      }
      ## The contig ID
      my $prokka_contig = $prokka_f->seq_id();
      ## prediction location.
      my $prokka_strand = $prokka_f->strand();
      my $prokka_start = $prokka_f->start();
      my $prokka_end = $prokka_f->end();
      ## I want to make sure display_name is set.
      my @display_names = $prokka_f->get_tag_values('locus_tag');
      ## With the caveat that I am going to turn around and change it almost immediately.
      my $named = $prokka_f->display_name($display_names[0]);
      my $continue_scanning = 1;
    SCANNER: while ($continue_scanning) {
        my $found_glimmer = Scan_Glimmer_Features(
            features => \@glimmer, contig => $prokka_contig,
            strand => $prokka_strand, start => $prokka_start, end => $prokka_end,
            orf_number => $orf_number, source => 'glimmer', source_seq => $current_source);
        @glimmer = @{$found_glimmer->{remaining}};
        if ($found_glimmer->{operation} eq 'prepend') {
            my $operation = $found_glimmer->{operation};
            my $new_data = $found_glimmer->{feature};
            my $new_number = $found_glimmer->{number};
            ## Replace the glimmer array with an array which is 1 element smaller.
            $orf_number = $new_number;
            push(@new_features, $new_data);
        } elsif ($found_glimmer->{operation} eq 'none') {
            my $formatted_number = sprintf("%04d", $orf_number);
            my $formatted_orf = qq"${prokka_contig}_${formatted_number}";
            if ($prokka_f->has_tag('protein_id')) {
                $prokka_f->remove_tag('protein_id');
            }
            $prokka_f->add_tag_value('protein_id', qq"prokka:${formatted_orf}");
            $prokka_f->remove_tag('locus_tag');
            $prokka_f->add_tag_value('locus_tag', $formatted_orf);
            push(@new_features, $prokka_f);
            $continue_scanning = 0;
            next PROKKA;
        }
    }
  }
    return(@new_features);
}

sub Gather_Prodigal_CDS {
    my %args = @_;
    my @prokka = @{$args{start}};
    my @prodigal = @{$args{add}};
    my $annotation_tag = 'inference';
    if (defined($args{tag})) {
        $annotation_tag = $args{tag};
    }
    ## Do it again for prodigal
    my @final_features = ();
    my $current_source;
    my $pcount = 0;
    my $orf_number = 1;
  PROKKA: for my $prokka_f (@prokka) {
      $pcount++;
      my $prokka_tag = $prokka_f->primary_tag();
      if ($prokka_tag eq 'source') {
          push(@final_features, $prokka_f);
          $current_source = $prokka_f;
          next PROKKA;
      }
      if ($prokka_tag eq 'gene') {
          next PROKKA;
      }
      ## The contig ID
      my $prokka_contig = $prokka_f->seq_id();
      ## prediction location.
      my $prokka_strand = $prokka_f->strand();
      my $prokka_start = $prokka_f->start();
      my $prokka_end = $prokka_f->end();
      my $continue_scanning = 1;
      while ($continue_scanning) {
          my $found_prodigal = Scan_Prodigal_Features(
              features => \@prodigal, contig => $prokka_contig,
              strand => $prokka_strand, start => $prokka_start, end => $prokka_end,
              orf_number => $orf_number, source => 'prodigal', source_seq => $current_source);
          @prodigal = @{$found_prodigal->{remaining}};
          if ($found_prodigal->{operation} eq 'prepend') {
              my $operation = $found_prodigal->{operation};
              my $new_data = $found_prodigal->{feature};
              my $new_number = $found_prodigal->{number};
              ## print "Setting orf_number from $orf_number to $new_number\n";
              ## Replace the glimmer array with an array which is 1 element smaller.
              $orf_number = $new_number;
              push(@final_features, $new_data);
          } elsif ($found_prodigal->{operation} eq 'merge') {
              my $new_data = $found_prodigal->{feature};
              my $new_number = $found_prodigal->{number};
              my $formatted_number = sprintf("%04d", $new_number);
              $orf_number = $new_number;
              my $formatted_orf = qq"${prokka_contig}_${formatted_number}";
              $prokka_f->remove_tag('protein_id');
              $prokka_f->add_tag_value('protein_id', qq"prodigal:${formatted_orf}");
              $prokka_f->remove_tag('locus_tag');
              $prokka_f->add_tag_value('locus_tag', $formatted_orf);
              my @new_scores = $new_data->get_tag_values($annotation_tag);
              $prokka_f->add_tag_value($annotation_tag, $new_scores[0]);
              $prokka_f->display_name($formatted_orf);
              push(@final_features, $prokka_f);
              $continue_scanning = 0;
              next PROKKA;
          } elsif ($found_prodigal->{operation} eq 'none') {
              my $formatted_number = sprintf("%04d", $orf_number);
              my $formatted_orf = qq"${prokka_contig}_${formatted_number}";
              $prokka_f->remove_tag('protein_id');
              $prokka_f->add_tag_value('protein_id', qq"prokka:${formatted_orf}");
              $prokka_f->remove_tag('locus_tag');
              $prokka_f->add_tag_value('locus_tag', $formatted_orf);
              $prokka_f->display_name($formatted_orf);
              push(@final_features, $prokka_f);
              $continue_scanning = 0;
              $orf_number++;
              next PROKKA;
          }
      }
  }
    return(@final_features);
}

sub Rename_Features {
    my %args = @_;
    my @features = @{$args{features}};
    my $count = 0;
  RENAME: for my $n (@features) {
      next if ($n->primary_tag() eq 'source');
      $count++;
      my $display_number = sprintf("%04d", $count);
      my $primary_tag = $n->seq_id();
      my $display_name = qq"${primary_tag}_${display_number}";
      my $protein_name = qq"Phokka:${display_name}";
      my $end = $n->end;
      my $start = $n->start;
      my $strand = $n->strand;
      my $id = $n->display_name($display_name);
      $n->remove_tag('locus_tag');
      $n->add_tag_value('locus_tag', $display_name);
      $n->remove_tag('protein_id');
      $n->add_tag_value('protein_id', $protein_name);
  }
    if ($args{add_genes}) {
        my @new_features = ();
      FEATLOOP: for my $f (@features) {
          if ($f->primary_tag() eq 'source') {
              push(@new_features, $f);
              next FEATLOOP;
          }
          my $gene = Bio::SeqFeature::Generic->new(
              -primary_tag => 'gene',
              -start => $f->start,
              -strand => $f->strand,
              -end => $f->end,
              -tag => { locus_tag => $f->display_name },);
          push(@new_features, $gene);
          push(@new_features, $f);
      } ## Finished iterating over the features and adding genes.
        @features = @new_features;
  }
    my $len = scalar(@features);
    return(@features);
}

sub Scan_Glimmer_Features {
    my %args = @_;
    my @feat = @{$args{features}};
    my $contig = $args{contig};
    my $strand = $args{strand};
    my $start = $args{start};
    my $end = $args{end};
    my $number = $args{orf_number};
    my $source_seq = $args{source_seq};
    my $seq = $source_seq->seq();
    my $formatted_number = sprintf("%04d", $number);
    my $formatted_orf = qq"${contig}_${formatted_number}";
    $number++;

    my $fcount = 0;
    my $returned_feature = undef; ## The feature to send back
    my $operation = 'none';
  THESE: for my $f (@feat) {
      $fcount++;
      my $this_contig = $f->seq_id();
      $this_contig =~ s/^gnl\|Prokka\|//g;
      my $this_name = $f->name();
      my $this_strand = $f->strand();
      my $this_start = $f->start();
      my $this_end = $f->end();
      my $this_phase = $f->phase();
      if ($contig eq $this_contig) {
          if ($this_end < $end) {
              ## Remove the feature to be returned from the feature array.
              $operation = 'prepend';
              my $features_pre = scalar(@feat);
              my $tmp_feature = splice(@feat, $fcount - 1, 1);
              my $features_post = scalar(@feat);
              my $subsequence = $seq->trunc($this_start, $this_end);
              if ($this_strand eq '-1') {
                  $subsequence = $subsequence->revcom();
              }
              my $translation = $subsequence->translate->seq();
              $returned_feature = Bio::SeqFeature::Generic->new(
                  -start => $this_start,
                  -end => $this_end,
                  -strand => $this_strand,
                  -primary => 'CDS',
                  -phase => $this_phase,
                  -source_tag => $args{source},
                  -display_name => $formatted_orf,
                  -score => $f->score,
                  -tag => {
                      locus_tag => $formatted_orf,
                      protein_id => qq"$args{source}:${formatted_orf}",
                      product => 'hypothetical protein',
                      translation => $translation,
                      inference => "ab initio prediction:glimmer",
                      transl_table => 11,
                  },);
              $returned_feature->attach_seq($seq);
              $returned_feature->seq_id($source_seq->seq_id);
              ## And get out of this loop.
              last THESE;
          } elsif (($this_end == $end && $this_strand eq $strand && $this_start eq $start) ||
                   ($this_end == $end && $this_strand eq $strand)) {
              ## This glimmer feature is identical in location to an existing prokka feature.
              ## So pull it out of the array and discard it.
              my $tmp_feature = splice(@feat, $fcount - 1, 1);
              last THESE;
          }
      } else {
          next THESE;
      }
  } ## End iterating over the array to destroy.

    ## We should get here only through a last or if we ran out of features.
    my $ret;
    if (defined($returned_feature)) {
        $ret = {
            operation => $operation,
            number => $number,
            feature => $returned_feature,
            remaining => \@feat,
        };
    } else {
        $ret = {
            operation => 'none',
            remaining => \@feat,
        };
    }
    return($ret);
}

sub Scan_Prodigal_Features {
    my %args = @_;
    my @feat = @{$args{features}};
    my $contig = $args{contig};
    my $strand = $args{strand};
    my $start = $args{start};
    my $end = $args{end};
    my $number = $args{orf_number};
    my $source_seq = $args{source_seq};
    my $annotation_tag = 'inference';
    if (defined($args{tag})) {
        $annotation_tag = $args{tag};
    }
    my $seq = $source_seq->seq();
    my $formatted_number = sprintf("%04d", $number);
    my $formatted_orf = qq"${contig}_${formatted_number}";
    $number++;

    my $fcount = 0;
    my $returned_feature = undef; ## The feature to send back
    my $operation = 'none';
  THESE: for my $f (@feat) {
      $fcount++;
      my $this_contig = $f->seq_id();
      $this_contig =~ s/^gnl\|Prokka\|//g;
      my $this_name = $f->name();
      my $this_strand = $f->strand();
      my $this_start = $f->start();
      my $this_end = $f->end();
      my $this_phase = $f->phase();
      ## print "PRODIGAL THESE: $this_start $this_end $this_strand\n";
      if ($contig eq $this_contig) {
          if ($this_end < $end) {
              $operation = 'prepend';
              my $tmp_feature = splice(@feat, $fcount - 1, 1);
              my $subsequence = $seq->trunc($this_start, $this_end);
              if ($this_strand eq '-1') {
                  $subsequence = $subsequence->revcom();
              }
              my $translation = $subsequence->translate->seq();

              my @spacer = $f->get_tag_values('rbs_spacer');
              my @rbs = $f->get_tag_values('rbs_motif');
              my @uscore = $f->get_tag_values('uscore');
              my @tscore = $f->get_tag_values('tscore');
              my @sscore = $f->get_tag_values('sscore');
              my @rscore = $f->get_tag_values('rscore');
              my @cscore = $f->get_tag_values('cscore');
              my @start_type = $f->get_tag_values('start_type');
              my $score_string = qq"rbs:$rbs[0] spacer:$spacer[0] start:$start_type[0] c:$cscore[0] r:$rscore[0] s:$sscore[0] t:$tscore[0] u:$uscore[0]";
              $returned_feature = Bio::SeqFeature::Generic->new(
                  -start => $this_start,
                  -end => $this_end,
                  -strand => $this_strand,
                  -primary => 'CDS',
                  -phase => $this_phase,
                  -source_tag => $args{source},
                  -display_name => $formatted_orf,
                  -score => $f->score,
                  -tag => { locus_tag => $formatted_orf,
                            protein_id => qq"$args{source}:${formatted_orf}",
                            product => 'hypothetical protein',
                            translation => $translation,
                            inference => 'ab initio prediction:prodigal',
                            transl_table => 11,
                            $annotation_tag => $score_string,
                  },);
              $returned_feature->attach_seq($seq);
              $returned_feature->seq_id($source_seq->seq_id);
              ## And get out of this loop.
              last THESE;
          } elsif (($this_end == $end && $this_strand eq $strand && $this_start eq $start) ||
                   ($this_end == $end && $this_strand eq $strand)) {
              ## Remove the feature to be returned from the feature array.
              $operation = 'merge';
              my $tmp_feature = splice(@feat, $fcount - 1, 1);
              my $subsequence = $seq->trunc($this_start, $this_end);
              if ($this_strand eq '-1') {
                  $subsequence = $subsequence->revcom();
              }
              my $translation = $subsequence->translate->seq();

              my @spacer = $f->get_tag_values('rbs_spacer');
              my @rbs = $f->get_tag_values('rbs_motif');
              my @uscore = $f->get_tag_values('uscore');
              my @tscore = $f->get_tag_values('tscore');
              my @sscore = $f->get_tag_values('sscore');
              my @rscore = $f->get_tag_values('rscore');
              my @cscore = $f->get_tag_values('cscore');
              my @start_type = $f->get_tag_values('start_type');
              my $score_string = qq"rbs:$rbs[0] spacer:$spacer[0] start:$start_type[0] c:$cscore[0] r:$rscore[0] s:$sscore[0] t:$tscore[0] u:$uscore[0]";
              $returned_feature = Bio::SeqFeature::Generic->new(
                  -start => $this_start,
                  -end => $this_end,
                  -strand => $this_strand,
                  -primary => 'CDS',
                  -phase => $this_phase,
                  -source_tag => $args{source},
                  -display_name => $formatted_orf,
                  -score => $f->score,
                  -tag => { locus_tag => $formatted_orf,
                            protein_id => qq"$args{source}:${formatted_orf}",
                            product => 'hypothetical protein',
                            translation => $translation,
                            inference => 'ab initio prediction:prodigal',
                            transl_table => 11,
                            $annotation_tag => $score_string,
                  },);
              $returned_feature->attach_seq($seq);
              $returned_feature->seq_id($source_seq->seq_id);
              ## And get out of this loop.
              last THESE;
          }
      } else {
          next THESE;
      }
  } ## End iterating over the array to destroy.

    ## We should get here only through a last or if we ran out of features.
    my $ret;
    if (defined($returned_feature)) {
        $ret = {
            operation => $operation,
            number => $number,
            feature => $returned_feature,
            remaining => \@feat,
        };
    } else {
        $ret = {
            operation => 'none',
            remaining => \@feat,
        };
    }
    return($ret);
}

sub Predict_to_Features {
    my %args = @_;
    my $glimmer = $args{in};
    my $glimmer_in = FileHandle->new("less ${glimmer} |");
    my $contig_id;
    my @glimmer_features = ();
  ENTRIES: while (my $line = <$glimmer_in>) {
      chomp $line;
      if ($line =~ /^>/) {
          ## Then this tells us the contig
          $contig_id = $line;
          $contig_id =~ s/^>(\w+)?\s+(.*)$/$1/g;
          if ($contig_id =~ /^>/) {
              $contig_id =~ s/^>//g;
          }
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
      my $feature = Bio::SeqFeature::Lite->new(
          -seq_id => $contig_id,
          -name => $orf_id,
          -start => $start,
          -end => $end,
          -score => $score,
          -phase => $phase,
          -strand => $strand,
          -type => 'glimmer',);
      ## Something is weird when constructing one of these, it is always setting strand to +1.
      $feature->strand($strand);
      push(@glimmer_features, $feature);
  } ## End iterating over every line of the glimmer3 output.
    $glimmer_in->close();
    return(@glimmer_features);
}

sub Write_Fsa {
    my %args = @_;
    my $id_prefix = 'gnl|Phokka|';
    if (defined($args{id_prefix})) {
        $id_prefix = $args{id_prefix};
    }
    my @in = @{$args{input_seq}};
    my $out = Bio::SeqIO->new(-file => qq">$args{output_fsa}",
                              -format => 'Fasta');
    my $written = 0;
    for my $s (@in) {
        $written++;
        my $start_id = $s->id;
        my $new_id = $s->id(qq"gnl|Phokka|${start_id}");
        $out->write_seq($s);
    }
    return($written);
}

sub Write_CDS {
    my %args = @_;
    my @in = @{$args{input_seq}};
    my $out = Bio::SeqIO->new(-file => qq">$args{output_cds}",
                              -format => 'Fasta');
    my $written = 0;
    my $count = 0;
    my $full_sequence;
    for my $s (@in) {
        $count++;
        if ($count == 1) {
            $full_sequence = $s->seq;
            next;
        }
        next if ($s->primary_tag eq 'source');
        next if ($s->primary_tag eq 'gene');
        my $id = $s->display_name;
        my $start = $s->start;
        my $end = $s->end;
        my $strand = $s->strand;
        my $sequence = $full_sequence->trunc($start, $end);
        if ($strand < 0) {
            $sequence = $sequence->revcom;
        }
        my $sequence_string = $sequence->seq;
        my $seq = Bio::Seq->new(-display_id => $id, -seq => $sequence_string);
        $out->write_seq($seq);
        $written++;
    }
    return($written);
}

sub Write_Faa {
    my %args = @_;
    my @in = @{$args{input_seq}};
    my $out = Bio::SeqIO->new(-file => qq">$args{output_faa}",
                              -format => 'Fasta');
    my $written = 0;
    my $count = 0;
    my $full_sequence;
    for my $s (@in) {
        $count++;
        if ($count == 1) {
            $full_sequence = $s->seq;
            next;
        }
        next if ($s->primary_tag eq 'source');
        next if ($s->primary_tag eq 'gene');
        my $id = $s->display_name;
        my $start = $s->start;
        my $end = $s->end;
        my $strand = $s->strand;
        my $sequence = $full_sequence->trunc($start, $end);
        if ($strand < 0) {
            $sequence = $sequence->revcom;
        }
        $sequence = $sequence->translate;
        my $sequence_string = $sequence->seq;
        my $seq = Bio::Seq->new(-display_id => $id, -seq => $sequence_string);
        $out->write_seq($seq);
        $written++;
    }
    return($written);
}

=head2 C<Write_Tbl>

This is basically just stealing a small piece of code from Torsten
Seeman's prokka tool which is capable of writing a tbl file in
preparation for running NCBI's tbl2sqn and generating a genbank file.

=cut
sub Write_Tbl {
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
    my $tbl2asn_comment = qq"Annotated using $EXE $VERSION from $URL";
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
    if ($run_sed) {
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

sub Write_Gff {
    my %args = @_;
    my @in = @{$args{input_seq}};
    use Bio::Tools::GFF;
    my $gff_out = FileHandle->new(">$args{output_gff}");
    my $gffio = Bio::Tools::GFF->new(-noparse => 1, -gff_version => 3);
    my $written = 0;
    for my $s (@in) {
        my $gff_string = $gffio->gff_string($s);
        print $gff_out $gff_string;
        $written++;
    }
    $gff_out->close();
    return($written);
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
