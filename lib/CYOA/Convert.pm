package CYOA;
use common::sense;
use autodie;
use File::Basename;
use Bio::FeatureIO;

=head1 NAME

    CYOA::Convert - Perform conversions between various formats:
    sam->bam, gff->fasta, genbank->fasta, etc.

=head1 SYNOPSIS

    use CYOA;
    my $hpgl = new CYOA;
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
    $me->Check_Options(["gff","genome"]);
    my $genome = $me->{genome};
    $genome = $args{genome} if (!defined($genome));
    my $gff = $me->{gff};
    $gff = $args{gff} if (!defined($gff));
    my $genome_basename = basename($genome, ('.fasta'));
    my $chromosomes = $me->Read_Genome(genome => $genome);

    my $gff_handle = new FileHandle;
    $gff_handle->open("less ${gff} |");
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
    $feature_type = $args{feature_type} if ($args{feature_type});
    my $annotation_in = new Bio::Tools::GFF(-fh => $gff_handle, -gff_version => 3);
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
      print $out_fasta_amino ">${gff_chr}_${id}
${aa_cds}
";
      print $out_fasta_nt ">${gff_chr}_${id}
${cds}
";
      $features_written++;
  } ## End LOOP
    $gff_handle->close();
    $out_fasta_amino->close();
    $out_fasta_nt->close();
    return($features_written);
}

sub Read_GFF {
    my $me = shift;
    my %args = @_;
    my $chromosomes = $args{chromosomes};
    open(GFF, "<$args{gff}");
    use Bio::Tools::GFF;
    my $annotation_in = new Bio::Tools::GFF(-fh => \*GFF, -gff_version => 3);
    my $gff_out = {};
    print "Starting to read gff: $args{gff}\n";
  LOOP: while(my $feature = $annotation_in->next_feature()) {
      next LOOP unless ($feature->{_primary_tag} eq $args{gff_type});
      my $location = $feature->{_location};
      my $start = $location->start();
      my $end = $location->end();
      my $strand = $location->strand();
      my @ids = $feature->each_tag_value("ID");
      my $id = "";
      my $gff_chr = $feature->{_gsf_seq_id};
      my $gff_string = $annotation_in->gff_string($feature);
      if (!defined($chromosomes->{$gff_chr})) {
          print STDERR "Something is wrong with $gff_chr\n";
          next LOOP;
      }
      foreach my $i (@ids) {
          $i =~ s/^cds_//g;
          $i =~ s/\-\d+$//g;
          $id .= "$i ";
      }
      $id =~ s/\s+$//g;
      my @gff_information = split(/\t+/, $gff_string);
      my $description_string = $gff_information[8];
      my $orf_chromosome = $gff_chr;
      my $annot = {
          id => $id,
          start => $start,  ## Genomic coordinate of the start codon
          end => $end,      ## And stop codon
          strand => $strand,
          description_string => $description_string,
          chromosome => $gff_chr,
      };
      $gff_out->{$gff_chr}->{$id} = $annot;
  } ## End looking at every gene in the gff file
    close(GFF);
    return($gff_out);
}


=item C<Gff2Gtf>

    $hpgl->Gff2Gtf(gff => 'mmusculus.gff')
    reads a given gff file and writes a gtf file from the features
    found therein.  It returns the number of features written.

    note that I only use gtf files for tophat, thus they must have a tag 'transcript_id'!
    This is woefully untrue for the tritrypdb gff files.  Thus I need to have a regex in this
    to make sure that web_id or somesuch is changed to transcript_id.

=cut
sub Gff2Gtf {
    my $me = shift;
    my %args = @_;
    my $input = $args{gff};

    my ($name, $path, $suffix) = fileparse($input, qr/\.gff/);
    my $out_file = $path . $name . ".gtf";

    my $in_gff = new Bio::FeatureIO('-file' => "$input",
                                    '-format' => 'GFF',
                                    '-version' => 3);
    open(OUT_GTF, ">${out_file}");
    select(OUT_GTF);
     $| = 1;
    ##my $out_gtf = new FileHandle();
    ##$out_gtf->open(">$out_file");

    my $features_written = 0;
    my $in_features = 0;
    FEATURES: while (my $feature = $in_gff->next_feature()) {
        $in_features++;
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
        my $phase = '.';
        if ($primary_id eq 'CDS') {
            $phase = $feature->phase();
        }
        next FEATURES if ($primary_id eq 'supercontig' or $primary_id eq 'region' or $primary_id eq 'chromosome');
        ##        my $string = qq"${seqid}\t${source}\tCDS\t$start\t$end\t.\t$strand\t0\t";
        my $string = qq"${seqid}\t${source}\t${primary_id}\t${start}\t${end}\t.\t${strand}\t${phase}\t";

        my $last_column = "";
        my $annot = $feature->{annotation};
        my $stringified;
        my $key;
        my $transcript_string = "";
        my $geneid_string = "";
        foreach $key ($annot->get_all_annotation_keys()) {
            my @values = $annot->get_Annotations($key);
            foreach my $value (@values) {
                $stringified = $value->{value};
                $stringified = $value->{term}->{name} unless($stringified);
                $last_column .= qq!${key} "${stringified}"; !;
            }
            if ($key eq 'ID') {
                $transcript_string = qq!transcript_id "${stringified}"; !;
            }
            if ($key eq 'seq_id') {
                $geneid_string = qq!gene_id "${stringified}"; !;
            }
        }
        my $seq_string = $string . ${transcript_string} . ${geneid_string} . $last_column;
        $seq_string =~ s/\s+$//g;
        ##print STDOUT "$seq_string\n";
        print "$seq_string\n";
        ##print OUT_GTF $seq_string;
        $features_written++;
    } ## End iterating over every FEATURES
    ##    $out_gtf->close();
    close(OUT_GTF);
    select(STDOUT);
    print STDERR "Finished writing gtf file.\n";
    return($features_written);
}

=item C<Read_Genome>

    $hpgl->Read_Genome(genome => 'mmusculus.fasta');
    read the given genome file and return a hash of the chromosomes.

=cut
sub Read_Genome {
    my $me = shift;
    my %args = @_;
    my $genome = $args{genome};
    my $chromosomes = {};
    my $fh = new FileHandle;
    $fh->open("less $genome |");
    my $input = new Bio::SeqIO(-fh => $fh, -format => 'Fasta');
    while (my $genome_seq = $input->next_seq()) {
        next unless(defined($genome_seq->id));
        my $id = $genome_seq->id;
        $chromosomes->{$id} = $genome_seq;
    }
    $fh->close();
    return($chromosomes);
}

=item C<Sam2Bam>

    $hpgl->Sam2Bam();
    Used to invoke samtools to take the sam output from bowtie/bwa and
    convert it to an compressed-sorted-indexed bam alignment.

    This function just calls $me->Samtools(), but has a little logic
    to see if the invocation of this is for an extent .sam file or
    calling it on an existing .fastq(.gz), in which case one must
    assume it is being called on one or more files in the bowtie_out/
    directory and which start with $basename, include something like
    -trimmed-v0M1.sam.

=cut
sub Sam2Bam {
    my $me = shift;
    $me->Check_Options(["species"]);
    my %args = @_;
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


=item C<Samtools>

    $hpgl->Samtools() calls (in order): samtools view, samtools sort,
    and samtools index.  Upon completion, it invokes bamtools stats to
    see what the alignments looked like.

    It explicitly does not pipe one samtools invocation into the next,
    not for any real reason but because when I first wrote it, it
    seemed like the sorting was taking too long if I did not already
    have the alignments in a bam file.

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
## Also printing alignment statistics to ${output}.stats
## This job depended on: $depends!;
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

=item C<Gb2Gff>

    $hpgl->Gb2Gff() takes a genbank genome file and splits it into:
    A genomic fasta file, CDS fasta, peptide fasta, gff file of all
    entries, CDS, and interCDS regions.

=cut
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

    my $in = new FileHandle;
    $in->open("less ${input} |");
    my $seqio = new Bio::SeqIO(-format => 'genbank', -fh => $in);
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

        ## $feature is an object of type: Bio::SeqFeatureI
        ## Some of the things you can call on it are:
        ## display_name(), primary_tag(), has_tag(something), get_tag_values(), get_all_tags(), attach_seq(Bio::Seq),
        ## seq(), entire_seq(), seq_id(), gff_string(), spliced_seq(), primary_id(), phase(), 
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
                my $id_string = "";
                my $id_hash = {};
                foreach my $thing (keys %{$feat_object->{_gsf_tag_hash}}) {
                    my @arr = @{$feat_object->{_gsf_tag_hash}->{$thing}};
                    ##print Dumper $feat_object->{_gsf_tag_hash}->{$thing};
                    ##print "WTF: $arr[0]\n";
                    $id_hash->{$thing} = $arr[0];
                }
                my $id = qq"$id_hash->{protein_id}";
                my $desc = qq"$id_hash->{gene} ; $id_hash->{db_xref} ; $id_hash->{product}";
                ## print "TESTME: $id\n";
                my $start = $feat_object->start;
                my $end = $feat_object->end;
                my $len = $feat_object->length;
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
                my $seq_object = new Bio::PrimarySeq(-id => $id, -seq => $seq, description => $desc);
                my $pep_object = new Bio::PrimarySeq(-id => $id, -seq => $pep, description => $desc);
                $cds_fasta->write_seq($seq_object);
                my $ttseq = $seq_object->seq;
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
    close($in);
    return($ret_stats);
}

=item C<TriTryp2Text>

    $hpgl->TriTryp2Text() generates some simple text tables from the
    much more elaborate text files provided by the TriTrypDB.

=cut
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

=item C<TriTryp_DL_Text>

    $hpgl->TriTryp_DL_Text() downloads the newest TriTrypDb text,
    genomic fasta, gff, and annotated CDS files for $args{species}.
    The species must be in their format, so unlike other invocations
    of 'species' where 'lmajor' would be fine, this would require:
    'LmajorFriedlin' or 'TcruziCLBrener'

=cut
sub TriTryp_Download {
    my $me = shift;
    my %args = @_;
    my $species = $args{species};

    my $ua = new LWP::UserAgent;
    $ua->agent("CYOA/downloader ");

    my $final_text_request;
    my $final_fasta_request;
    my $final_cds_request;
    my $final_gff_request;
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
              $final_text_request = $line;
              $final_gff_request = $final_text_request;
              $final_gff_request =~ s/Gene\.txt/\.gff/g;
              $final_fasta_request = $final_gff_request;
              $final_fasta_request =~ s/\.gff/_Genome\.fasta/g;
              $final_cds_request = $final_fasta_request;
              $final_cds_request =~ s/_Genome\.fasta/_AnnotatedCDSs\.fasta/g;
                        ##TriTrypDB-26_LmajorFriedlin_AnnotatedCDSs.fasta
          }
      }
    }
    if ($text_listing) {
        my $text_output = new FileHandle;
        my $gff_output = new FileHandle;
        my $fasta_output = new FileHandle;
        my $cds_output = new FileHandle;
        my $text_url = qq"http://tritrypdb.org/common/downloads/Current_Release/${species}/txt/${final_text_request}";
        print "Going to write to $final_text_request with $text_url\n";
        unless (-r $final_text_request) { ## Don't redownload if we already did.
            $text_output->open(">$final_text_request");
            my $req = HTTP::Request->new(GET => qq"$text_url");
            $req->content_type('text/html');
            my $res = $ua->request($req);
            if ($res->is_success) {
                print $text_output $res->content;
            } else {
                print $res->status_line, "\n";
            }
            $text_output->close();
        } ## End unless

        my $gff_url = qq"http://tritrypdb.org/common/downloads/Current_Release/${species}/gff/data/${final_gff_request}";
        print "Going to write to $final_gff_request with $gff_url\n";
        unless (-r $final_gff_request) { ## Don't redownload if we already did.
            $gff_output->open(">$final_gff_request");
            my $req = HTTP::Request->new(GET => qq"$gff_url");
            $req->content_type('text/html');
            my $res = $ua->request($req);
            if ($res->is_success) {
                print $gff_output $res->content;
            } else {
                print $res->status_line, "\n";
            }
            $gff_output->close();
        } ## End unless

        my $fasta_url = qq"http://tritrypdb.org/common/downloads/Current_Release/${species}/fasta/data/${final_fasta_request}";
        print "Going to write to $final_fasta_request with $fasta_url\n";
        unless (-r $final_fasta_request) { ## Don't redownload if we already did.
            $fasta_output->open(">$final_fasta_request");
            my $req = HTTP::Request->new(GET => qq"$fasta_url");
            $req->content_type('text/html');
            my $res = $ua->request($req);
            if ($res->is_success) {
                print $fasta_output $res->content;
            } else {
                print $res->status_line, "\n";
            }
            $fasta_output->close();
        } ## End unless

        my $cds_url = qq"http://tritrypdb.org/common/downloads/Current_Release/LmajorFriedlin/fasta/data/${final_cds_request}";

        print "Going to write to $final_cds_request with $cds_url\n";
        unless (-r $final_cds_request) { ## Don't redownload if we already did.
            $cds_output->open(">$final_cds_request");
            my $req = HTTP::Request->new(GET => qq"$cds_url");
            $req->content_type('text/html');
            my $res = $ua->request($req);
            if ($res->is_success) {
                print $cds_output $res->content;
            } else {
                print $res->status_line, "\n";
            }
            $cds_output->close();
        } ## End unless

    } ## End checking to see that we got a downloadable filename
    return($final_fasta_request);
}

=item C<TriTryp_Ortho>

    $hpgl->TriTryp_Ortho() generates a text table of
    orthologs/paralogs from the large TriTrypDB text files.

=cut
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

=item C<TriTryp_Master>

    $hpgl->TriTryp_Master() generates a text table of
    the ~50 columns of data from the large TriTrypDB text files.

=cut
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

=item C<TriTryp_GO>

    $hpgl->TriTryp_GO() generates a text table of gene ontology data
    from the large TriTrypDB text files.

=cut
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
