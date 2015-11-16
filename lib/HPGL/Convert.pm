package HPGL;
use common::sense;
use autodie;

=head1 NAME

    HPGL::Convert - Perform conversions between various formats:
    sam->bam, gff->fasta, genbank->fasta, etc.

=head1 SYNOPSIS

    use HPGL;
    my $hpgl = new HPGL;
    $hpgl->Gff2Fasta(genome => 'mmusculus.fasta', gff => 'mmusculus.gff');

=head2 Methods

=over 4

=item C<Gff2Fasta>

    $hpgl->Gff2Fasta(genome => 'something.fasta', gff => 'something.gff')
    will read the genome's fasta and annotation files and return a set
    of fasta files which include the coding sequence regions which are
    members of args{feature_type} in the gff file.

    It writes one file of amino acid sequence and one of nucleotides.
    Upon completion, it returns the number of entries written.

=cut
sub Gff2Fasta {
    my $me = shift;
    my %args = @_;
    my $genome = $args{genome};
    my $gff = $args{gff};
    my $genome_basename = basename($genome, ('.fasta'));
    my $chromosomes = $me->Read_Genome(genome => $genome);
    open(GFF, "<${gff}");
    my $out_fasta_amino = new FileHandle();
    my $out_fasta_nt = new FileHandle();
    my $aa_out_name = qq"${genome_basename}_cds_aa.fasta";
    my $nt_out_name = qq"${genome_basename}_cds_nt.fasta";
    $out_fasta_amino->open(">${aa_out_name}");
    $out_fasta_nt->open(">${nt_out_name}");
    my $tag = 'ID';
    $tag = $me->{tag} if ($me->{tag});
    my $feature_type = 'CDS';
    $feature_type = $me->{feature_type} if ($me->{feature_type});
    my $annotation_in = new Bio::Tools::GFF(-fh => \*GFF, -gff_version => 3);
    my $features_written = 0;
  LOOP: while(my $feature = $annotation_in->next_feature()) {
      ##print "TAGS: $feature->{_primary_tag} vs $me->{feature_type}\n" if ($me->{debug} == 1);
      next LOOP unless ($feature->{_primary_tag} eq $feature_type);
      my $location = $feature->{_location};
      my $start = $location->start();
      my $end = $location->end();
      my $strand = $location->strand();
      my @something = $feature->each_tag_value();
      ##print Dumper $feature if ($me->{debug} == 1);
      my @ids = $feature->each_tag_value($tag);
      my $gff_chr = $feature->{_gsf_seq_id};
      my $gff_string = $annotation_in->{$gff_chr};
      if (!defined($chromosomes->{$gff_chr})) {
          print STDERR "Something is wrong with $gff_chr\n";
          next LOOP;
      }
      my $id = "";
      foreach my $i (@ids) {
          $id .= "$i ";
      }
      $id =~ s/(\s+$)//g;
      my $cds = $chromosomes->{$gff_chr}->subseq($start, $end);
      if ($strand == -1) {
          $cds = reverse($cds);
          $cds =~ tr/ATGCatgc/TACGtacg/;
      }
      my $seq_obj = new Bio::Seq();
      $seq_obj->seq($cds);
      my $aa_cds = $seq_obj->translate->seq();
      print $out_fasta_amino ">${gff_chr} ${id}
${aa_cds}
";
      print $out_fasta_nt ">${gff_chr} ${id}
${cds}
";
      $features_written++;
  } ## End LOOP
    close(GFF);
    $out_fasta_amino->close();
    $out_fasta_nt->close();
    return($features_written);
}

=item C<Gff2Gtf>

    $hpgl->Gff2Gtf(gff => 'mmusculus.gff')
    reads a given gff file and writes a gtf file from the features
    found therein.  It returns the number of features written.

=cut
sub Gff2Gtf {
    my $me = shift;
    my %args = @_;
    my $input = $args{gff};
    use File::Basename;
    use Bio::FeatureIO;

    my ($name, $path, $suffix) = fileparse($input, qr/\.gff/);
    my $out_file = $path . $name . ".gtf";

    my $in_gff = new Bio::FeatureIO('-file' => "$input",
                                    '-format' => 'GFF',
                                    '-version' => 3);
    my $out_gtf = new FileHandle();
    $out_gtf->open(">$out_file");

    my $features_written = 0;
    while (my $feature = $in_gff->next_feature()) {
        # the output handle is reset for every file
        ##According to UCSC, GTF file contains 9 column:
        ##<seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes]
        ## An example working gtf line: (except the double spaces are tabs...)
        ##TcChr20-S  TriTrypDB  CDS  14765  15403  .  -  0    gene_id "cds_TcCLB.397937.10-1"; transcript_id "cds_TcCLB.397937.10-1";

        my @tags = $feature->get_all_tags();
        my $seqid = $feature->seq_id();
        my $location = $feature->location();
        my $start = $feature->start();
        my $end = $feature->end();
        my $strand = $feature->strand();
        if ($strand == 1) {
            $strand = '+';
        } else {
            $strand = '-';
        }
        my $source = $feature->source_tag();
        my $primary_id = $feature->primary_tag();
        my $string = qq"$seqid\t$source\tCDS\t$start\t$end\t.\t$strand\t0\t";

        my $last_column = "";
        my $annot = $feature->{annotation};
        my $stringified;
        my $key;
        foreach $key ($annot->get_all_annotation_keys()) {
            my @values = $annot->get_Annotations($key);
            foreach my $value (@values) {
                $stringified = $value->{value};
                $stringified = $value->{term}->{name} unless($stringified);
##                print "TSTME: $value $key $stringified\n";
                $last_column .= qq!${key} "${stringified}"; !;
            }
        }
        $string .= "$last_column\n";
        print $out_gtf $string;
    }
    $out_gtf->close();
    return($features_written);
}

=head2
    Read_Genome()
=cut
sub Read_Genome {
    my $me = shift;
    my %args = @_;
    my $genome = $args{genome};
    my $chromosomes = {};
    print "TESTME: REading $genome\n";
    my $input = new Bio::SeqIO(-file => $genome, -format => 'Fasta');
    while (my $genome_seq = $input->next_seq()) {
        next unless(defined($genome_seq->id));
        my $id = $genome_seq->id;
        $chromosomes->{$id} = $genome_seq;
    }
    return($chromosomes);
}

=head2
    Sam2Bam()
=cut
sub Sam2Bam {
    my $me = shift;
    $me->Check_Options(["species"]);
    my %args = @_;
    ## A little logic to see if the invocation of this is for an extent .sam file
    ## Or calling it on an existing .fastq(.gz), in which case one must assume
    ## It is being called on one or more files in the bowtie_out/ directory
    ## And which start with $basename, include something like -trimmed-v0M1.sam
    my $basename = $me->{basename};
    my @input_list = ();
    if ($args{input}) {
        push(@input_list, $args{input});
    } elsif (-r $me->{input} and $me->{input} =~ /\.sam$/) {
        push(@input_list, $me->{input});
        my $sam = $me->Samtools(sam => \@input_list);
    } elsif (-r $me->{input} and $me->{input} =~ /\.fastq$/) {
        if (-d "bowtie_out") {
            find({ wanted => sub { push(@input_list, "bowtie_out/$_") if ($_ =~ /\.sam$/); }, follow => 1 }, 'bowtie_out/');
            my $sam = $me->Samtools(sam => \@input_list);
        } else {
            foreach my $k (%{$me->{bt_args}}) {
                my $output_string = "bowtie_out/${basename}-${k}.sam";
                push(@input_list, $output_string);
            }
            my $bt = $me->Bowtie();
            my $sam = $me->Samtools(depends => $bt->{pbs_id}, sam => \@input_list);
        }
    } else {
        die("I don't know what to do without a .fastq file nor a .sam file.\n");
    }
}


=head2
    Samtools()
=cut
sub Samtools {
    my $me = shift;
    my %args = @_;
    my $basename = $me->{basename};
    my $depends = "";
    if ($args{depends}) {
        $depends = $args{depends};
    }
    my $input = $args{input};
    my $output = $input;
    $output =~ s/\.sam$/\.bam/g;
    my $sorted = $input;
    $sorted =~ s/\.sam$//g;
    $sorted = qq"${sorted}-sorted";
    print "Converting to a compressed/sorted bam file.\n";
    my $job_string = qq!samtools view -u -t $me->{libdir}/genome/$me->{species}.fasta -S $input 1>$output
samtools sort -l 9 $output $sorted
rm $output && rm $input && mv ${sorted}.bam $output && samtools index $output
bamtools stats -in $output 2>${output}.stats 1>&2
!;
    my $comment = qq!## Converting the text sam to a compressed, sorted, indexed bamfile.
## Also printing alignment statistics to ${output}.stats!;
    my $samtools = $me->Qsub(job_name => "sam",
                             depends => $depends,
                             job_string => $job_string,
                             input => $input,
                             output => $output,
                             comment => $comment,
			     prescript => $args{prescript},
			     postscript => $args{postscript},
        );
    return($samtools);
}

sub Gb2Gff {
    my $me = shift;
    my %args = @_;
    my $input = $me->{input};

    if (!defined($input)) {
        my $term = new Term::ReadLine('>');
        my $attribs = $term->Attribs;
        $attribs->{completion_suppress_append} = 1;
        my $OUT = $term->OUT || \*STDOUT;
        $input = $term->readline("Please provide an input genbank file: ");
        $input =~ s/\s+$//g;
    }
    my @suffix = (".gb", ".genbank");
    my $base = basename($input, @suffix);

    my $seqio = new Bio::SeqIO(-format => 'genbank', -file => $input);
    my $seq_count = 0;
    my $total_nt = 0;
    my $feature_count = 0;
    my $fasta = new Bio::SeqIO(-file => qq">${base}.fasta", -format => 'fasta', -flush => 0);
    my $gffout = new Bio::Tools::GFF(-file => ">${base}_all.gff", -gff_version => 3);
    my $gene_gff = new Bio::Tools::GFF(-file => ">${base}.gff", -gff_version => 3);
    my $cds_gff = new Bio::Tools::GFF(-file => ">${base}_cds.gff", -gff_version => 3);
    my $inter_gffout = new Bio::Tools::GFF(-file => ">${base}_interCDS.gff", -gff_version => 3);
    my $cds_fasta = new Bio::SeqIO(-file => qq">${base}_cds.fasta", -format => 'Fasta');
    my $pep_fasta = new Bio::SeqIO(-file => qq">${base}_pep.fasta", -format => 'Fasta');
    while(my $seq = $seqio->next_seq) {
        $seq_count++;
        $total_nt = $total_nt + $seq->length();
        $fasta->write_seq($seq);
        my @feature_list = ();
        # defined a default name

        foreach my $feature ($seq->top_SeqFeatures()) {
            $gffout->write_feature($feature);
            my $feat_count = 0;
        }
        for my $feat_object ($seq->get_SeqFeatures) {
            $feature_count++;
            my $feat_count = 0;
            if ($feat_object->primary_tag eq "CDS") {
                $cds_gff->write_feature($feat_object);
                $feat_count++;
                my $start = $feat_object->start;
                my $end = $feat_object->end;
                my $len = $feat_object->length;
                my $id = $feat_object->{_gsf_tag_hash}->{locus_tag}->[0];
                my $seq = $feat_object->spliced_seq->seq;
                my $pep = $feat_object->spliced_seq->translate->seq;
                my $size = {
                    id => $id,
                    start => $start,
                    end => $end,
                    len => $len,
                    feat => $feat_object,
                };
                push(@feature_list, $size);
                my $seq_object = new Bio::PrimarySeq(-id => $id, -seq => $seq);
                my $pep_object = new Bio::PrimarySeq(-id => $id, -seq => $pep);
                $cds_fasta->write_seq($seq_object);
                ##print "TESTME: $seq_object\n";
                my $ttseq = $seq_object->seq;
                ##print "TESTME: $id $ttseq\n";
                $pep_fasta->write_seq($pep_object);
            } elsif ($feat_object->primary_tag eq "gene") {
                $gene_gff->write_feature($feat_object);
            }
        } ## End looking at every feature and putting them into the @feature_list
        ## Now make a hash from it and fill in the inter-cds data

        my %inter_features = ();
        my ($p_st, $p_en, $c_st, $c_en, $inter_start, $inter_end) = 0;
      INTER: for my $c (0 .. $#feature_list) {
          my $p_size = {};
          if ($c == 0) {
              $p_size = $feature_list[$#feature_list];
          } else {
              $p_size = $feature_list[$c - 1];
          }

          my $c_size = $feature_list[$c];
          next INTER if (!defined($p_size->{start}));
          next INTER if (!defined($p_size->{end}));
          next INTER if (!defined($c_size->{start}));
          next INTER if (!defined($c_size->{end}));

          my $p_start = $p_size->{start};
          my $p_end = $p_size->{end};
          my $c_start = $c_size->{start};
          my $c_end = $c_size->{end};
          my $interstart = 0;
          my $interend = 0;

          if ($p_start > $p_end) {
              $interstart = $p_start;
          } else {
              $interstart = $p_end;
          }
          if ($c_start < $c_end) {
              $interend = $c_start;
          } else {
              $interend = $c_end;
          }
          my $tmp_start = $interstart;
          my $tmp_end = $interend;
          if ($tmp_start > $tmp_end) {
              $interstart = $tmp_end;
              $interend = $tmp_start;
          }
          ## I need a little logic to catch the first feature
          ## As it is currently written, it grabs the position of the last feature
          ## and uses that number rather than the beginning of the chromosome.
          ## Therefore we end up with a first feature which contains the entire chromosome.
          if ($c <= 2) {  ## Arbitrarily set it to check the first 2 features
              if ($interend >= 100000) {  ## If the first feature take >= 1/2 the chromosome
                  $interend = $interstart;
                  $interstart = 1;
              }
          }

          my $new_feature = $c_size->{feat};
          my $inter_location = new Bio::Location::Atomic(-start => $interstart, -end => $interend, -strand => 1);
          $new_feature->{_location} = $inter_location;
          $inter_gffout->write_feature($new_feature);
      }
    } ## End while every sequence
    my $ret_stats = {
        num_sequences => $seq_count,
        total_nt => $total_nt,
        num_features => $feature_count,
    };
    return($ret_stats);
}

sub TriTryp2Text {
    my $me = shift;
    my %args = @_;
    use List::MoreUtils qw"uniq";

    my $input = $me->{input};
    if (!defined($input)) {
        my $term = new Term::ReadLine('>');
        my $attribs = $term->Attribs;
        $attribs->{completion_suppress_append} = 1;
        my $OUT = $term->OUT || \*STDOUT;
        $input = $term->readline("Please provide an input tritryp text file: ");
        $input =~ s/\s+$//g;
    }
    my $species = $me->{species};
    if (!defined($species)) {
        my $term = new Term::ReadLine('>');
        my $attribs = $term->Attribs;
        $attribs->{completion_suppress_append} = 1;
        my $OUT = $term->OUT || \*STDOUT;
        $species = $term->readline("Please provide a basename for the output files (I am just using species): ");
        $species =~ s/\s+$//g;
    }
    my $base = $species;

    my $master_table = $me->TriTryp_Master(input => $input, base => $species);
    my ($ortho_table, $go_table);
    if ($args{ortho}) {
        $ortho_table = $me->TriTryp_Ortho(input => $input, base => $species);
    }
    if ($args{go}) {
        $go_table = $me->TriTryp_GO(input => $input, base => $species);
    }

    my $ret = {
        master_lines => $master_table->{lines_written},
        ortho_lines => $ortho_table->{lines_written},
        go_lines => $go_table->{lines_written},
    };
    return($ret);
}

sub TriTryp_DL_Text {
    my $me = shift;
    my %args = @_;
    my $species = $args{species};

    my $ua = new LWP::UserAgent;
    $ua->agent("HPGL/downloader ");

    my $final_request;
    # Create a request
    my $species = 'LmajorFriedlin';
    my $req = HTTP::Request->new(GET => qq"http://tritrypdb.org/common/downloads/Current_Release/${species}/txt/");
    $req->content_type('text/html');
    my $res = $ua->request($req);
    my $text_listing = "";
    if ($res->is_success) {
        $text_listing = $res->content;
    } else {
        print $res->status_line, "\n";
    }
    my $line = "";
    if ($text_listing) {
        my @lines = split(/\n/, $text_listing);
        my $filename = "";
      LISTING: foreach $line (@lines) {
          chomp $line;
          next LISTING unless ($line =~ m/\<a href="TriTrypDB.+Gene\.txt">/);
          $line =~ s/^.*\<a href=.*"\>(TriTryp.+Gene\.txt).*$/\1/g;
          if ($line =~ /^TriTrypDB/) {
              $final_request = $line;
          }
      }
    }
    if ($text_listing) {
        my $output = new FileHandle;
        my $url = qq"http://tritrypdb.org/common/downloads/Current_Release/${species}/txt/${final_request}";
        print "Going to write to $final_request with the output from: $url\n";
        unless (-r $final_request) { ## Don't redownload if we already did.
            $output->open(">$final_request");
            ##open(OUT, ">$text_listing");
            my $req = HTTP::Request->new(GET => qq"$url");
            $req->content_type('text/html');
            my $res = $ua->request($req);
            if ($res->is_success) {
                ##print OUT $res->content;
                print $output $res->content;
            } else {
                print $res->status_line, "\n";
            }
            $output->close();
            ##close(OUT);
        } ## End unless
    } ## End checking to see that we got a downloadable filename
    return($final_request);
}


sub TriTryp_Ortho {
    my $me = shift;
    my %args = @_;
    my $input = $args{input};
    my $base = $args{base};
    my $orth = {};
    my $current_id = "";
    my $in = new FileHandle;
    my $out = new FileHandle;
    $in->open("<$input");
    $out->open(">${base}_ortho.txt");
    my $reading_ortho = 0;
    my $ortho_count = 0;
    my $lines_written = 0;
    LINES: while(my $line = <$in>) {
	chomp $line;
	next if ($line =~ /\[/);
	next if ($line =~ /^\s*$/);
	if ($line =~ /:/) {
	    my ($key, $datum) = split(/:\s+/, $line);
	    if ($key eq 'Gene ID') {
		$current_id = $datum;
		$reading_ortho = 0;
	    } elsif ($key eq 'TABLE' and $datum eq 'Orthologs and Paralogs within TriTrypDB') {
		$reading_ortho = 1;
		$ortho_count = 0;
		next LINES;
	    } elsif ($key eq 'TABLE') {
		$reading_ortho = 0;
	    }
	}
	if ($reading_ortho == 1) {
	    $ortho_count++;
	    my($gene, $organism,$product,$syntenic,$comments) = split(/\t+/, $line);
#	    $orth->{$current_id}->{$ortho_count}->{gene} = $gene;
#	    $orth->{$current_id}->{$ortho_count}->{organism} = $organism;
#	    $orth->{$current_id}->{$ortho_count}->{product} = $product;
#	    $orth->{$current_id}->{$ortho_count}->{syntenic} = $syntenic;
#	    $orth->{$current_id}->{$ortho_count}->{comments} = $comments;
#	    print "Adding $current_id $ortho_count $gene\n";
	    print $out "$current_id\t$ortho_count\t$gene\t$organism\t$product\t$syntenic\t$comments\n";
            $lines_written++;
	}
    }
    $in->close();
    $out->close();
    $orth->{lines_written} = $lines_written;
    return($orth);
}

sub TriTryp_Master {
    my $me = shift;
    my %args = @_;
    my $input = $args{input};
    my $base = $args{base};
    my $master_table = {};
    my $current_id = "";
    my $current_key = "";
    my @key_list = ();
    my $in = new FileHandle;
    my $out = new FileHandle;
    print "TESTME: reading $input\n";
    $in->open("<$input");
    my $count = 0;
    my $waiting_for_next_gene = 0;
    my $lines_written = 0;
    LINES: while(my $line = <$in>) {
	chomp $line;
        if ($line =~ /^TABLE/) {
            $waiting_for_next_gene = 1;
        }
        if ($line =~ /^Gene ID/) {
            $waiting_for_next_gene = 0;
        }
        next LINES if ($waiting_for_next_gene);
	next unless($line =~ /:/);
	my ($key, $datum) = split(/:\s+/, $line);
	if ($key eq 'Gene ID') {
	    $current_id = $datum;
	} else {
	    $current_key = $key;
	    $current_key =~ s/\s+/_/g;
            push(@key_list, $current_key);
	    $master_table->{$current_id}->{$current_key} = $datum;
	}
    }
    $in->close();
##    print "New test: Starting uniq\n";
    @key_list = uniq(@key_list);
    my $num_keys = scalar(@key_list);
##    print "End test: Ended uniq $num_keys @key_list\n";
    my $header_string = "ID\t";
    foreach my $key (@key_list) {
        $header_string .= "$key\t";
    }
    $header_string =~ s/\s+$//g;
    $out->open(">${base}_master.tab");
    print $out "${header_string}\n";
    foreach my $id (keys %{$master_table}) {
        my %datum = %{$master_table->{$id}};
        my $line_string = "$id\t";
        foreach my $key (@key_list) {
            my $dat = "undef";
            if (defined($datum{$key})) {
                $dat = $datum{$key};
            }
            $line_string .= "$dat\t";
        }
        $line_string =~ s/\s+$//g;
        print $out "${line_string}\n";
        $lines_written++;
    }
    $out->close();
    $master_table->{lines_written} = $lines_written;
    return($master_table);
}


sub TriTryp_GO {
    my $me = shift;
    my %args = @_;
    my $input = $args{input};
    my $base = $args{base};
    my $go = {};
    my $current_id = "";
    my $in = new FileHandle;
    my $out = new FileHandle;
    my $gaf = new FileHandle;
    $in->open("<$input");
    $out->open(">${base}_go.txt");
    $gaf->open(">${base}_go.gaf");
    my $reading_go = 0;
    my $go_count = 0;
    my $lines_written = 0;
    LINES: while(my $line = <$in>) {
	chomp $line;
	next if ($line =~ /\[/);
	next if ($line =~ /^\s*$/);
	if ($line =~ /:/) {
	    my ($key, $datum) = split(/:\s+/, $line);
	    if ($key eq 'Gene ID') {
		$current_id = $datum;
		$reading_go = 0;
	    } elsif ($key eq 'TABLE' and $datum eq 'GO Terms') {
		$reading_go = 1;
		$go_count = 0;
		next LINES;
	    } elsif ($key eq 'TABLE') {
		$reading_go = 0;
	    }
	}
	if ($reading_go == 1) {
	    $go_count++;
	    my($gene, $ontology, $go_term_name, $source, $evidence_code, $low_evidence_code) = split(/\t+/, $line);
#	    $orth->{$current_id}->{$ortho_count}->{gene} = $gene;
#	    $orth->{$current_id}->{$ortho_count}->{organism} = $organism;
#	    $orth->{$current_id}->{$ortho_count}->{product} = $product;
#	    $orth->{$current_id}->{$ortho_count}->{syntenic} = $syntenic;
#	    $orth->{$current_id}->{$ortho_count}->{comments} = $comments;
#	    print "Adding $current_id $ortho_count $gene\n";
            print $out "$current_id\t$gene\t$ontology\t$go_term_name\t$source\t$evidence_code\n";
####                   1      2           3     4 5       6      7                8            9       10        11     12     13                      14       15
####                   source gene   shortname qual|go  dbref  evidencecode     withfrom       P    goname      gosym  type   taxon                   date      assignedby
            print $gaf "TRI\t$current_id\tundef\t\t$gene\tundef\t$evidence_code\t$go_term_name\tP\t$go_term_name\tundef\tgene\ttaxon:$me->{taxid}}\t20050101\tTriTryp\n";
            $lines_written++;
	}
    }
    $in->close();
    $out->close();
    $gaf->close();
    $go->{lines_written} = $lines_written;
    return($go);
}

=back

=head1 AUTHOR - atb

Email abelew@gmail.com

=head1 SEE ALSO

    L<samtools> L<Bio::FeatureIO>

=cut

1;
