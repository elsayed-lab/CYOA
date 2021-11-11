package Bio::Adventure::Annotation_Genbank;
## LICENSE: gplv2
## ABSTRACT:  Kitty!
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
    ## A tsv containing arbitrarily chosen annotation information.
    my $output_tsv = qq"${output_dir}/${output_basename}.tsv";
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
  output_tsv => '${output_tsv}', ## .tsv
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
        output_tsv => $output_tsv,
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
        output_tsv => 'assembly.tsv',
        jprefix => '19',);

    ## First read in the prokka data as a set of sequences and features
    ## Keep in mind that the features will have the full sequence attached.
    my ($pfeat, $pseq, $seqids) = Read_Prokka_Gbk_to_SeqFeatures(
        input => $options->{input});
    my @prokka_features = @{$pfeat};
    my @prokka_sequences = @{$pseq};
    my @sequence_ids = @{$seqids};

    my $prodfeat = Read_Prodigal_GFF_to_SeqFeatures(
        input => $options->{input_prodigal});
    my @prodigal_features = @{$prodfeat};

    ## Read the glimmer predict file, transform it into a feature list, and return that.
    my @glimmer_features = Predict_to_Features(in => $options->{input_glimmer});
    ## This order of operations is important, I am effectively saying that I trust the
    ## run of prodigal that I manually performed more than either prokka or glimmer
    ## because it will override either of them.  Though to be fair, most of the time
    ## the prokka-derived prodigal data is expected to agree (except I did train them
    ## slightly differently and set a couple of options.

    my $prokka_glimmer = Gather_Glimmer_CDSv3(prokka => \@prokka_features,
                                              glimmer => \@glimmer_features);
    ## Now we should have a data structure containing:
    ##  1.  A hash of the sources (e.g. the contigs)
    ##  2.  A hash of genes from prokka, which I may discard and rewrite.
    ##  3.  A hash of other features from prokka, which will be good to keep.
    ##  4.  A hash of cds from the merge of prokka and glimmer.
    my @merged_features = @{$prokka_glimmer->{cds}};
    my $final_features = Gather_Prodigal_CDSv3(merged => \@merged_features,
                                               prodigal => \@prodigal_features);
    my @final_set = @{$final_features};
    ## Add the non CDS features
    my $final_length = scalar(@final_set);
    for my $other (keys %{$prokka_glimmer->{other}}) {
        my $other_feature = $prokka_glimmer->{other}->{$other};
        push(@final_set, $other_feature);
    }
    $final_length = scalar(@final_set);
    ##my $prokka_glimmer_len = scalar(@prokka_glimmer);
    ##print "TESTME: $prokka_glimmer_len\n";
    ##my @prok_glim_prod = Gather_Prodigal_CDS(prokka => \@prokka_glimmer,
    ##                                         prodigal => \@prodigal_features);
    ##my $all_len = scalar(@prok_glim_prod);
    ##print "TESTME: $all_len\n";
    ##my $num = scalar(@prok_glim_prod);

    ## Lets make the CDS names just be the basename of the sample
    ## followed by an incrementer.
    my $prefix = basename($options->{basedir});
    my @renamed = Rename_Features(features => \@final_set,
                                  prefix => $prefix,
                                  add_genes => 1);
    ## Write the various output files.
    my @assembly_features = ();
    for my $contig (keys %{$prokka_glimmer->{source}}) {
        my $source_feature = $prokka_glimmer->{source}->{$contig};
        ## Put the source at the front.
        push(@assembly_features, $source_feature);
    }

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
}

## Something in my logic for merging the locus_tag names is incorrect here.
sub Gather_Glimmer_CDSv3 {
    my %args = @_;
    my @prokka = @{$args{prokka}};
    my @glimmer = @{$args{glimmer}};
    my $source_remove = 'gnl\|Prokka\|' unless defined($args{source_remove});
    my $keep_both = 0;
    $keep_both = $args{keep_both} if ($args{keep_both});

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
          push(@cds_features, $prokka_f);
      } else {
          $other_features{$prokka_name} = $prokka_f;
      }
  }

    ## These three arrays will hold the final information
    ## I may drop the genes.
    my @final_cds = ();
    my $prokka_count = 0;
    my %finished_glimmer = ();
    my %finished_prokka = ();
  CDS: for my $p_cds (@cds_features) {
      $prokka_count++;
      ## Definitely need the contig:
      my $p_contig = $p_cds->seq_id();
      ## and the information about the location
      my $p_strand = $p_cds->strand();
      my $p_start = $p_cds->start();
      my $p_end = $p_cds->end();
      my ($p_fivep, $p_threep) = ($p_start, $p_end);
      if ($p_strand < 0) {
          ($p_fivep, $p_threep) = ($p_end, $p_start);
      }
      ## Get the set of tags and the locus tag specifically
      my @tags = $p_cds->get_all_tags();
      my @locus_tags = $p_cds->get_tag_values('locus_tag');
    GLIMMER: for my $g_cds (@glimmer) {
        ## We do not have to think about the type of these features.
        my $g_strand = $g_cds->strand();
        my $g_start = $g_cds->start;
        my $g_end = $g_cds->end;
        my ($g_fivep, $g_threep) = ($g_start, $g_end);
        if ($g_strand < 0) {
            ($g_fivep, $g_threep) = ($g_end, $g_start);
        }
        if (defined($finished_glimmer{$g_threep})) {
            next GLIMMER;
        }
        ## We are scanning through the prokka features in lexical order (start->end),
        ## therefore scan the glimmer data _only_ until we get to the current end.
        ## Thus, if we pass the current prokka end, get out of this inner loop
        ## With one important caveat, what if there are glimmer annotations
        ## after the final prokka end?
        if ($prokka_count < scalar(@cds_features)) {
            ## As long as we are not on the last feature
            ## test to see if we should break out to the next.
            my $cds_featlen = scalar(@cds_features);
            if ($g_end > $p_end) {
                ## print "$prokka_count is less than $cds_featlen and glimmer: $g_end is bigger than prokka end: $p_end\n";
                last GLIMMER;
            }
        }

        ## Now think through the possibilities:
        ##  1.  They are the same, if so, take the prokka feature
        ##  2.  They agree on the stop codon, but not the start, decide what to do.
        ##  3.  They disagree, then in this inner loop add the glimmer feature
        ##  4.  They disagree, and the prokka feature is never found, then add the prokka feature after the inner loop.
        if (($p_fivep eq $g_fivep) and ($p_threep eq $g_threep)) {
            ## print "prokka: $p_fivep and $p_threep are the same as glimmer.\n";
            ## Then use the prokka annotation.
            $p_cds->add_tag_value('note', 'cds_prediction: prodigal via prokka.');
            push(@final_cds, $p_cds);
            $finished_glimmer{$g_threep} = 1;
            $finished_prokka{$p_threep} = 1;
        } elsif ($p_threep eq $g_threep) {
            ## print "The threep annotations are the same: $p_threep\n";
            ## When they do not agree on the start codon, I say go with prokka.
            $g_cds->add_tag_value('note', 'cds_prediction: untrained glimmer.');
            push(@final_cds, $p_cds);
            $finished_glimmer{$g_threep} = 1;
            $finished_prokka{$p_threep} = 1;
            if ($keep_both) {
                push(@final_cds, $g_cds);
            }
        } else {
            ## print "The glimmer is different: $g_threep and $g_fivep\n";
            ## The final possibility is that the glimmer entry is different.
            $finished_glimmer{$g_threep} = 1;
            $g_cds->add_tag_value('note', 'cds_prediction: untrained glimmer.');
            push(@final_cds, $g_cds);
        }
    } ## End inner iteration over glimmer features.
      if (!defined($finished_prokka{$p_threep})) {
          ## Then this prokka entry was never found.
          ## print "The entry $p_start $p_end is only in prokka.\n";
          $p_cds->add_tag_value('note', 'cds_prediction: prodigal via prokka.');
          push(@final_cds, $p_cds);
    }
      my $final_cds_length = scalar(@final_cds);
      ## print "End outer loop, cds array is: $final_cds_length entries.\n";
  } ## End the outer iteration over prokka features.

    my %ret = (
        source => \%source_features,
        cds => \@final_cds,
        gene => \%gene_features,
        other => \%other_features,);
    return(\%ret);
}


sub Gather_Prodigal_CDSv3 {
    my %args = @_;
    my @merged = @{$args{merged}};
    my @prodigal = @{$args{prodigal}};
    my $keep_both = 0;
    $keep_both = $args{keep_both} if ($args{keep_both});

    ## These three arrays will hold the final information
    ## I may drop the genes.
    my @final_cds = ();
    my $merged_count = 0;
    my %finished_merged = ();
    my %finished_prodigal = ();
  MERGED: for my $m_cds (@merged) {
      $merged_count++;
      my $m_contig = $m_cds->seq_id();
      ## and the information about the location
      my $m_strand = $m_cds->strand();
      my $m_start = $m_cds->start();
      my $m_end = $m_cds->end();
      my ($m_fivep, $m_threep) = ($m_start, $m_end);
      if ($m_strand < 0) {
          ($m_fivep, $m_threep) = ($m_end, $m_start);
      }
      ## Get the set of tags and the locus tag specifically
      my @tags = $m_cds->get_all_tags();
      my @locus_tags = $m_cds->get_tag_values('locus_tag');
    PRODIGAL: for my $p_cds (@prodigal) {
        ## We do not have to think about the type of these features.
        my $p_strand = $p_cds->strand();
        my $p_start = $p_cds->start;
        my $p_end = $p_cds->end;
        my ($p_fivep, $p_threep) = ($p_start, $p_end);

        ## Do a little work to standardize the prodigal feature:
        my $spacer = '';
        my $rbs = '';
        my $uscore = '';
        my $tscore = '';
        my $sscore = '';
        my $rscore = '';
        my $cscore = '';
        my $start_type = '';
        my @spacers = $p_cds->get_tag_values('rbs_spacer');
        $p_cds->remove_tag('rbs_spacer');
        $spacer = $spacers[0] if ($spacers[0]);
        my @rbses = $p_cds->get_tag_values('rbs_motif');
        $p_cds->remove_tag('rbs_motif');
        $rbs = $rbses[0] if ($rbses[0]);
        my @uscores = $p_cds->get_tag_values('uscore');
        $p_cds->remove_tag('uscore');
        $uscore = $uscores[0] if ($uscores[0]);
        my @tscores = $p_cds->get_tag_values('tscore');
        $p_cds->remove_tag('tscore');
        $tscore = $tscores[0] if ($tscores[0]);
        my @sscores = $p_cds->get_tag_values('sscore');
        $p_cds->remove_tag('sscore');
        $sscore = $sscores[0] if ($sscores[0]);
        my @rscores = $p_cds->get_tag_values('rscore');
        $p_cds->remove_tag('rscore');
        $rscore = $rscores[0] if ($rscores[0]);
        my @cscores = $p_cds->get_tag_values('cscore');
        $p_cds->remove_tag('cscore');
        $cscore = $cscores[0] if ($cscores[0]);
        my @start_types = $p_cds->get_tag_values('start_type');
        $p_cds->remove_tag('start_type');
        $start_type = $start_types[0] if ($start_types[0]);

        my $score_string = qq"cds_prediction: trained prodigal(rbs:${rbs} spacer:${spacer} start:${start_type} c:${cscore} r:${rscore} s:${sscore} t:${tscore} u:${uscore})";
        if (!$p_cds->has_tag('note')) {
            $p_cds->add_tag_value('note', $score_string);
        }

        if ($p_strand < 0) {
            ($p_fivep, $p_threep) = ($p_end, $p_start);
        }
        if (defined($finished_prodigal{$p_threep})) {
            next PRODIGAL;
        }
        ## We are scanning through the merged features in lexical order (start->end),
        ## therefore scan the prodigal data _only_ until we get to the current end.
        ## Thus, if we pass the current merged end, get out of this inner loop
        ## With one important caveat, what if there are prodigal annotations
        ## after the final merged end?
        if ($merged_count < scalar(@merged)) {
            ## As long as we are not on the last feature
            ## test to see if we should break out to the next.
            my $cds_featlen = scalar(@merged);
            if ($p_end > $m_end) {
                last PRODIGAL;
            }
        }

        ## Now think through the possibilities:
        ##  1.  They are the same, in this case, we take the prodigal feature
        ##  2.  They agree on the stop codon, but not the start, decide what to do.
        ##  3.  They disagree completely, then add the prodigal feature.
        ##  4.  They disagree and the merged feature is never found, then add the merged after.
        if (($m_fivep eq $p_fivep) and ($m_threep eq $p_threep)) {
            ## print "merged: $m_fivep and $m_threep are the same as prodigal.\n";
            ## Then use the prodigal annotation.
            push(@final_cds, $p_cds);
            $finished_merged{$m_threep} = 1;
            $finished_prodigal{$m_threep} = 1;
        } elsif ($m_threep eq $p_threep) {
            ## print "The threep annotations are the same: $p_threep\n";
            ## When they do not agree on the start codon, I say go with prodigal.
            push(@final_cds, $p_cds);
            $finished_prodigal{$p_threep} = 1;
            $finished_merged{$p_threep} = 1;
            if ($keep_both) {
                push(@final_cds, $m_cds);
            }
        } else {
            ## print "The prodigal is different: $p_threep and $p_fivep\n";
            ## The final possibility is that the prodigal entry is different.
            $finished_prodigal{$p_threep} = 1;
            push(@final_cds, $p_cds);
        }
    } ## End inner iteration over glimmer features.
      if (!defined($finished_merged{$m_threep})) {
          ## Then this merged entry was not found.
          ## print "The entry $m_start $m_end was only in the merged data.\n";
          push(@final_cds, $m_cds);
    }
      my $final_cds_length = scalar(@final_cds);
      ## print "End outer loop, cds array is: $final_cds_length entries.\n";
  } ## End the outer iteration over prokka features.
    return(\@final_cds);
}

sub Gather_Glimmer_CDS {
    my %args = @_;
    my @prokka = @{$args{prokka}};
    my @glimmer = @{$args{glimmer}};

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
    my @prokka = @{$args{prokka}};
    my @prodigal = @{$args{prodigal}};
    my $annotation_tag = 'note';
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
    my $prefix = $args{prefix};
    my @features = @{$args{features}};
    my $count = 0;
    my $contig_id;
  RENAME: for my $n (@features) {
      next if ($n->primary_tag() eq 'source');
      $count++;
      my $display_number = sprintf("%04d", $count);
      $contig_id = $n->seq_id();
      my $display_name = qq"${prefix}_${display_number}";
      my $protein_name = qq"Phokka:${display_name}";
      my $end = $n->end;
      my $start = $n->start;
      my $strand = $n->strand;
      my $id = $n->display_name($display_name);
      $n->remove_tag('locus_tag');
      $n->add_tag_value('locus_tag', $display_name);
      if ($n->has_tag('protein_id')) {
          $n->remove_tag('protein_id');
      }
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
              -seq_id => $contig_id,
              -display_name => $f->display_name,
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
    my @gfeat = @{$args{features}};
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
  THESE: for my $glimmer_f (@gfeat) {
      $fcount++;
      my $this_contig = $glimmer_f->seq_id();
      $this_contig =~ s/^gnl\|Prokka\|//g;
      my $this_name = $glimmer_f->name();
      my $this_strand = $glimmer_f->strand();
      my $this_start = $glimmer_f->start();
      my $this_end = $glimmer_f->end();
      my $this_phase = $glimmer_f->phase();
      if ($contig eq $this_contig) {
          if ($this_end < $end) {
              ## Remove the feature to be returned from the feature array.
              $operation = 'prepend';
              my $features_pre = scalar(@gfeat);
              my $tmp_feature = splice(@gfeat, $fcount - 1, 1);
              my $features_post = scalar(@gfeat);
              my $subsequence = $seq->trunc($this_start, $this_end);
              if ($this_strand eq '-1') {
                  $subsequence = $subsequence->revcom();
              }
              my $translation = $subsequence->translate->seq();
              $returned_feature = Bio::SeqFeature::Generic->new(
                  -display_name => $this_contig,
                  -start => $this_start,
                  -end => $this_end,
                  -strand => $this_strand,
                  -primary => 'CDS',
                  -phase => $this_phase,
                  -source_tag => $args{source},
                  -score => $glimmer_f->score,
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
              next THESE;
          } elsif (($this_end == $end && $this_strand eq $strand && $this_start eq $start) ||
                   ($this_end == $end && $this_strand eq $strand)) {
              ## This glimmer feature is identical in location to an existing prokka feature.
              ## So pull it out of the array and discard it.
              my $tmp_feature = splice(@gfeat, $fcount - 1, 1);
              next THESE;
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
            remaining => \@gfeat,
        };
    } else {
        $ret = {
            operation => 'none',
            remaining => \@gfeat,
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
    my $annotation_tag = 'note';
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
      ##my $feature = Bio::SeqFeature::Lite->new(
      ##    -seq_id => $fixed_id,
      ##    -name => $orf_id,
      ##    -start => $start,
      ##    -end => $end,
      ##    -score => $score,
      ##    -phase => $phase,
      ##    -strand => $strand,
      ##    -type => 'CDS',);
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
              locus_tag => $orf_id,
              transl_table => 11, });
      $feature->strand($strand);
      push(@glimmer_features, $feature);
  } ## End iterating over every line of the glimmer3 output.
    $glimmer_in->close();
    return(@glimmer_features);
}

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
        ## my $new_id = qq"${id_prefix}${start_id}";
        my $entire_sequence = $s->entire_seq();
        my $sequence_string = $entire_sequence->seq();
        ## my $seq_obj = Bio::Seq->new(-display_id => $new_id,
        my $seq_obj = Bio::Seq->new(-display_id => $id,
                                    -seq => $sequence_string);
        $out->write_seq($seq_obj);
        push(@seq_objects, $seq_obj);
    }
    return(\@seq_objects);
}

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
      ##for my $k (keys %{$s}) {
      ##    print "TESTME NAMED: $k\n";
      ##}
      my $id = $s->display_name;
      ##print "TESTME: THE ID IS: $id\n";
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
     my $cds_sequence = $full_sequence->trunc($start, $end);
     if ($strand < 0) {
         $cds_sequence = $cds_sequence->revcom;
     }
     my $aa_sequence = $cds_sequence->translate;
     my $cds_sequence_string = $cds_sequence->seq;
     my $aa_sequence_string = $aa_sequence->seq;
     $aa_sequence_string =~ s/\*$//g;

     my $cds_seq_obj = Bio::Seq->new(-display_id => $name, -seq => $cds_sequence_string);
     $out_cds->write_seq($cds_seq_obj);
     my $aa_seq_obj = Bio::Seq->new(-display_id => $name, -seq => $aa_sequence_string);
     $out_faa->write_seq($aa_seq_obj);
     $written++;

     ## Now figure out what information to write to the tsv file.

     ## $aa_sequence_string
     ##my @all_tags = $in->get_all_tags;
     ##print Dumper @all_tags;
     my @notes = $in->get_tag_values('note');
     my $note = $notes[0];
     $note =~ s/^cds_prediction:\s+//g;
     my $tsv_line = qq"${name}\t${contig}\t${type}\t${source}\t${start}\t${end}\t${strand}\t${note}\t$aa_sequence_string\n";
     print $out_tsv $tsv_line;
 }
    $out_tsv->close();
    return($written);
}

=head2 C<Write_Tbl>

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

    for my $sid (@seq) {
        print $tbl_fh ">Feature ${sid}\n";
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
            if (!defined($f->strand) || !defined($f->start) || !defined($f->end)) {
                my $name = $f->display_name;
                my $strand = $f->strand;
                my $start = $f->start;
                my $end = $f->end;
                print "THERE IS A PROBLEM WITH: $name\n";
                print "ONE OF STRAND:$strand START:$start OR END:$end IS UNDEFINED.\n";
            }
            my ($L, $R) = ($f->strand >= 0) ? ($f->start, $f->end) : ($f->end, $f->start);
            print $tbl_fh "$L\t$R\t", $f->primary_tag, "\n";
            WRITE_TAGS: for my $tag ($f->get_all_tags) {
                # remove GFF specific tags (start with uppercase letter)
                next WRITE_TAGS if $tag =~ m/^[A-Z]/ and $tag !~ m/EC_number/i;
                my @unwanted = ('gc_cont', 'sscore', 'phase', 'conf', 'tscore', 'partial', 'rbs_spacer', 'rbs_motif', 'cscore', 'start_type', 'uscore', 'seq_id', 'rscore', 'source', 'frame', 'score', 'type');
                ## next WRITE_TAGS if ($tag ~~ @unwanted);
                ## next WRITE_TAGS if any { $_ eq $tag } @unwanted;
                next WRITE_TAGS if (grep $_ eq $tag, @unwanted);
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

sub Read_Prodigal_GFF_to_SeqFeatures {
    my %args = @_;
    my $prodigal_in = FileHandle->new("less $args{input} |");
    my $prodigal_io = Bio::FeatureIO->new(-format => 'gff', -fh => $prodigal_in);
    my @prodigal_features = ();
    while (my $feature = $prodigal_io->next_feature()) {
        my $contig_id = $feature->seq_id;
        my $fixed_id = Remove_Contig_Cruft($contig_id);
        my $reset_contig = $feature->seq_id($fixed_id);
        push(@prodigal_features, $feature);
    }
    $prodigal_in->close();
    return(\@prodigal_features);
}

1;
