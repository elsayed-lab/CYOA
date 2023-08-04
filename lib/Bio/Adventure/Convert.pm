package Bio::Adventure::Convert;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use feature 'try';
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';
use Bio::FeatureIO;
use Bio::Tools::GFF;
use Bio::Root::Exception;
use File::Basename;
use File::Which qw"which";
use List::MoreUtils qw"uniq";
use Text::CSV_XS::TSV;

no warnings qw"experimental::try";

=head1 NAME

 Bio::Adventure::Convert - Perform conversions between various formats.


=head1 SYNOPSIS

 The functions here handle the various likely conversions one may perform.
 sam to bam, gff to fasta, genbank to fasta, etc.

=head1 METHODS

=head2 C<Gb2Gff>

 A poorly named genbank flat file converter.

 Gb2Gff() takes a genbank genome file and splits it into:
 A genomic fasta file, CDS fasta, peptide fasta, gff file of all
 entries, CDS, and interCDS regions.

=cut
sub Gb2Gff {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jprefix => '78');
    my $base = basename($options->{input}, ('.xz', '.gz', '.bz2'));
    $base = basename($base, ('.gb', '.gba', '.gbk', '.genbank'));
    my $dir = dirname($options->{input});
    if (defined($options->{output_dir})) {
        $dir = $options->{output_dir};
    }
    my $output_fasta = qq"${dir}/${base}.fsa";
    my $output_all_gff = qq"${dir}/${base}_all.gff";
    my $output_gene_gff = qq"${dir}/${base}.gff";
    my $output_cds_gff = qq"${dir}/${base}_cds.gff";
    my $output_inter_gff = qq"${dir}/${base}_interCDS.gff";
    my $output_cds_fasta = qq"${dir}/${base}.ffn";
    my $output_pep_fasta = qq"${dir}/${base}.faa";
    my $output_rrna_gff = qq"${dir}/${base}_rrna.gff";
    my $output_rrna_fasta = qq"${dir}/${base}_rrna.fasta";

    my $stderr = qq"${base}_convert.stderr";
    my $stdout = qq"${base}_convert.stdout";
    my $comment = '## Convert from genbank to other formats.';
    my $jname = qq"${base}_convert";
    my $jstring = qq!use Bio::Adventure::Convert;
my \$result = \$h->Bio::Adventure::Convert::Gb2Gff_Worker(
  input => '$options->{input}',
  jdepends => '$options->{jdepends}',
  jname => '${jname}',
  jprefix => '$options->{jprefix}',
  output => '${output_fasta}',
  output_fasta => '${output_fasta}',
  output_all_gff => '${output_all_gff}',
  output_gene_gff => '${output_gene_gff}',
  output_cds_gff => '${output_cds_gff}',
  output_inter_gff => '${output_inter_gff}',
  output_cds_fasta => '${output_cds_fasta}',
  output_pep_fasta => '${output_pep_fasta}',
  output_rrna_gff => '${output_rrna_gff}',
  output_rrna_fasta => '${output_rrna_fasta}',
  stdout => '${stdout}',
  stderr => '${stderr}',);
!;
    my $convert = $class->Submit(
        comment => $comment,
        input => $options->{input},
        jcpu => 1,
        jdepends => $options->{jdepends},
        jname => $jname,
        jprefix => $options->{jprefix},
        jstring => $jstring,
        output => $output_fasta,
        output_fasta => $output_fasta,
        output_all_gff => $output_all_gff,
        output_gene_gff => $output_gene_gff,
        output_cds_gff => $output_cds_gff,
        output_inter_gff => $output_inter_gff,
        output_cds_fasta => $output_cds_fasta,
        output_pep_fasta => $output_pep_fasta,
        output_rrna_gff => $output_rrna_gff,
        output_rrna_fasta => $output_rrna_fasta,
        stdout => $stdout,
        stderr => $stderr,
        language => 'perl',);
    return($convert);
}

sub Gb2Gff_Worker {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        output_fasta => 'output.fsa',
        output_all_gff => 'all.gff',
        output_gene_gff => 'gene_gff',
        required => ['input']);

    my $in = FileHandle->new("less $options->{input} |");
    my $seqio = Bio::SeqIO->new(-format => 'genbank', -fh => $in);
    my $seq_count = 0;
    my $total_nt = 0;
    my $feature_count = 0;

    my $fasta = Bio::SeqIO->new(-file => qq">$options->{output_fasta}",
                                -format => 'fasta', -flush => 0);
    my $gffout = Bio::Tools::GFF->new(-file => ">$options->{output_all_gff}",
                                      -gff_version => 3);
    my $gene_gff = Bio::Tools::GFF->new(-file => ">$options->{output_gene_gff}",
                                        -gff_version => 3);
    my $cds_gff = Bio::Tools::GFF->new(-file => ">$options->{output_cds_gff}",
                                       -gff_version => 3);
    my $inter_gffout = Bio::Tools::GFF->new(-file => ">$options->{output_inter_gff}",
                                            -gff_version => 3);
    my $cds_fasta = Bio::SeqIO->new(-file => qq">$options->{output_cds_fasta}",
                                    -format => 'Fasta');
    my $pep_fasta = Bio::SeqIO->new(-file => qq">$options->{output_pep_fasta}",
                                    -format => 'Fasta');
    my $rrna_gffout = Bio::Tools::GFF->new(-file => ">$options->{output_rrna_gff}", -gff_version => 3);
    my $rrna_fasta = Bio::SeqIO->new(-file => qq">$options->{output_rrna_fasta}", -format => 'Fasta');
    while (my $seq = $seqio->next_seq) {
        $seq_count++;
        $total_nt = $total_nt + $seq->length();
        $fasta->write_seq($seq);
        print "Wrote ${seq_count} features.\n";
        my @feature_list = ();
        ## defined a default name

        ## $feature is an object of type: Bio::SeqFeatureI
        ## Some of the things you can call on it are:
        ## display_name(), primary_tag(), has_tag(something), get_tag_values(), get_all_tags(), attach_seq(Bio::Seq),
        ## seq(), entire_seq(), seq_id(), gff_string(), spliced_seq(), primary_id(), phase(),
        foreach my $feature ($seq->top_SeqFeatures()) {
            $gffout->write_feature($feature);
            my $feat_count = 0;
        }
      FEAT: for my $feat_object ($seq->get_SeqFeatures) {
          $feature_count++;
          my $feat_count = 0;
          if ($feat_object->primary_tag eq 'rRNA') {
              $rrna_gffout->write_feature($feat_object);
              my $id_string = '';
              my $desc = '';
              my $id_hash = {};
              foreach my $thing (keys %{$feat_object->{_gsf_tag_hash}}) {
                  my @arr = @{$feat_object->{_gsf_tag_hash}->{$thing}};
                  $id_hash->{$thing} = $arr[0];
              }
              my $id = '';
              if (defined($id_hash->{protein_id})) {
                  $id = qq"$id_hash->{protein_id}";
              }
              if (defined($id_hash->{gene})) {
                  $desc .= "$id_hash->{gene} ; ";
              }
              if (defined($id_hash->{db_xref})) {
                  $desc .= "$id_hash->{db_xref} ; ";
              }
              if (defined($id_hash->{product})) {
                  $desc .= "$id_hash->{product}";
              }
              my $start = $feat_object->start;
              my $end = $feat_object->end;
              my $len = $feat_object->length;
              my $seq;
              try {
                  $seq = $feat_object->spliced_seq->seq;
              }
              catch ($e) {
                  print qq"Something went wrong getting the sequence.\n";
                  next FEAT;
              }
              my $size = {
                  id => $id_string,
                  start => $start,
                  end => $end,
                  len => $len,
                  feat => $feat_object,
              };
              my $rrna_object = Bio::PrimarySeq->new(-id => $id_string, -seq => $seq,
                                                     description => $desc);
              $rrna_fasta->write_seq($rrna_object);
          } elsif ($feat_object->primary_tag eq 'CDS') {
              $cds_gff->write_feature($feat_object);
              $feat_count++;
              my $id_string = '';
              my $id_hash = {};
              foreach my $thing (keys %{$feat_object->{_gsf_tag_hash}}) {
                  my @arr = @{$feat_object->{_gsf_tag_hash}->{$thing}};
                  $id_hash->{$thing} = $arr[0];
              }
              my $id = '';
              if (defined($id_hash->{protein_id})) {
                  $id = qq"$id_hash->{protein_id}";
              }
              my $desc = '';
              if (defined($id_hash->{gene})) {
                  $desc .= "$id_hash->{gene} ; ";
              }
              if (defined($id_hash->{db_xref})) {
                  $desc .= "$id_hash->{db_xref} ; ";
              }
              if (defined($id_hash->{product})) {
                  $desc .= "$id_hash->{product}";
              }
              $desc .= "\n";
              my $start = $feat_object->start;
              my $end = $feat_object->end;
              my $len = $feat_object->length;
              ## This is in response to the puzzling error:
              ## "Error::throw("Bio::Root::Exception", "Location end (601574) exceeds length (0) of called sequence C"...) called at /sw/local/perl/5.28.1/perl5/Bio/Root/Root.pm line 449"
              my $seq = '';
              my $pep = '';

              my $e;
              try {
                  $seq = $feat_object->spliced_seq->seq;
              }
              catch ($e) {
                  print "Something went wrong getting the sequence.\n";
                  next FEAT;
              }
              try {
                  $pep = $feat_object->spliced_seq->translate->seq;
              }
              catch ($e) {
                  print "Something went wrong getting the peptide sequence.\n";
                  next FEAT;
              }
              my $size = {
                  id => $id,
                  start => $start,
                  end => $end,
                  len => $len,
                  feat => $feat_object,
              };
              push(@feature_list, $size);
              my $seq_object = Bio::PrimarySeq->new(-id => $id, -seq => $seq, description => $desc);
              my $pep_object = Bio::PrimarySeq->new(-id => $id, -seq => $pep, description => $desc);
              $cds_fasta->write_seq($seq_object);
              my $ttseq = $seq_object->seq;
              $pep_fasta->write_seq($pep_object);
          } elsif ($feat_object->primary_tag eq 'gene') {
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
          if ($c <= 2) {   ## Arbitrarily set it to check the first 2 features
              if ($interend >= 100000) { ## If the first feature take >= 1/2 the chromosome
                  $interend = $interstart;
                  $interstart = 1;
              }
          }

          my $new_feature = $c_size->{feat};
          my $inter_location = Bio::Location::Atomic->new(-start => $interstart,
                                                          -end => $interend,
                                                          -strand => 1);
          $new_feature->{_location} = $inter_location;
          $inter_gffout->write_feature($new_feature);
      }
    }                           ## End while every sequence
    my $ret_stats = {
        num_sequences => $seq_count,
        total_nt => $total_nt,
        num_features => $feature_count,
        fasta_out => $options->{output_fasta},
        all_gff_out => $options->{output_all_gff},
        gene_gff => $options->{output_gene_gff},
        inter_cds => $options->{output_inter_gff},
        cds_fasta => $options->{output_cds_fasta},
        cds_gff => $options->{output_cds_gff},
        pep_fasta => $options->{output_pep_fasta},
        rrna_gff => $options->{output_rrna_gff},
        rrna_fasta => $options->{output_rrna_fasta},
    };
    close($in);
    $fasta->close();
    $gffout->close();
    $gene_gff->close();
    $cds_gff->close();
    $inter_gffout->close();
    $cds_fasta->close();
    $pep_fasta->close();
    $rrna_gffout->close();
    $rrna_fasta->close();
    return($ret_stats);
}

=head2 C<Gff2Fasta>

 Extract feature types from a gff/fasta pair.

 Gff2Fasta(genome => 'something.fasta', gff => 'something.gff')
 will read the genome's fasta and annotation files and return a set
 of fasta files which include the coding sequence regions which are
 members of args{feature_type} in the gff file.

 It writes one file of amino acid sequence and one of nucleotides.
 Upon completion, it returns the number of entries written.

=over

=item C<Arguments>

 species: Add the .gff and .fasta suffixes to this and look in $libdir.
 gff(required): Gff file to read.
 input(required): Genome fasta file to read.
 gff_tag('gene_id'): ID tag to use when extracting features.
 gff_type(undef): gff type to extract, left undefined it will count
  them up and choose the most highly represented.

=item C<Invocation>

> cyoa --method gff2 --input lmajor.fasta \
   --gff lmajor.gff --gff_type CDS --gff_tag ID

=cut
sub Gff2Fasta {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        species => '',
        input => '',
        gff => '',
        gff_tag => 'gene_id',
        gff_type => 'cds',);
    unless ($options->{species} || ($options->{input} && $options->{gff})) {
        die("This requires either a species or an input gff/fasta.")
    }
    my $genome;
    my $gff;
    my $species;
    if ($options->{species}) {
        my $input_base = qq"$options->{libpath}/genome/$options->{species}";
        $genome = qq"${input_base}.fasta";
        $gff = qq"${input_base}.gff";
        $species = $options->{species};
    } elsif ($options->{input} && $options->{gff}) {
        $genome = $options->{input};
        $gff = $options->{gff};
        $species = basename($genome, split(/,/, $class->{suffixes}));
        $species = basename($genome, ('.fasta', '.fna', '.faa'));
    }

    my $wanted_tag = $options->{gff_tag};
    my $wanted_type = $options->{gff_type};
    my $chromosomes = $class->Read_Genome_Fasta(genome => $genome);
    my $gff_handle = FileHandle->new("less ${gff} |");
    my $nt_out = Bio::SeqIO->new(
        -format => 'Fasta',
        -file => qq">${species}_${wanted_type}_${wanted_tag}_nt.fasta");
    my $aa_out = Bio::SeqIO->new(
        -format => 'Fasta',
        -file => qq">${species}_${wanted_type}_${wanted_tag}_aa.fasta");
    ## Note that this and the next line might not be a good idea,
    ## HT_Types only looks at the first n (40,000) records and uses that as a heuristic
    ## to see that the wanted type is actually in the gff file.
    my $feature_type = $class->Bio::Adventure::Count::HT_Types(
        annotation => $gff,
        feature_type => $wanted_type,);
    $feature_type = $feature_type->[0];
    my $annotation_in = Bio::Tools::GFF->new(-fh => $gff_handle, -gff_version => 3);
    my $features_written = 0;
    my $features_read = 0;
  LOOP: while (my $feature = $annotation_in->next_feature()) {
      $features_read++;
      my $primary_type = $feature->primary_tag();
      unless ($primary_type eq $wanted_type) {
          next LOOP;
      }
      my $start = $feature->start();
      my $end = $feature->end();
      my $strand = $feature->strand();
      my $gff_chr = $feature->seq_id();
      my @ids;
      my $e;
      try {
          @ids = $feature->each_tag_value($wanted_tag);
      }
      catch ($e) {
          print "Did not find the tag: ${wanted_tag}, perhaps check the gff file.\n";
          next LOOP;
      }
      my $id = $ids[0];
      my $genome_obj = $chromosomes->{$gff_chr}->{obj};
      my $cds = $genome_obj->trunc($start, $end);
      if ($strand == -1) {
          $cds = $cds->revcom;
      }
      my $id_string = qq"chr_${gff_chr}_id_${id}_start_${start}_end_${end}";
      my $nt_sequence = Bio::Seq->new(
          -id => $id_string,
          -seq => $cds->seq());
      $nt_out->write_seq($nt_sequence);
      my $aa_sequence = Bio::Seq->new(
          -id => $id_string,
          -seq => $cds->translate->seq());
      $aa_out->write_seq($aa_sequence);
      $features_written++;
  } ## End LOOP
    $gff_handle->close();
    print "Wrote ${features_written} features to the aa/nt files.\n";
    return($features_written);
}

=back

=head2 C<Gff2Gtf>

 Pretty much unused gff to gtf converter.

 Reads a given gff file and writes a gtf file from the features
 found therein.  It returns the number of features written.

 note that I only use gtf files for tophat, thus they must have a tag 'transcript_id'!
 This is woefully untrue for the tritrypdb gff files.  Thus I need to have a regex in this
 to make sure that web_id or somesuch is changed to transcript_id.

=cut
sub Gff2Gtf {
    my ($class, %args) = @_;
    my $input = $args{gff};

    my ($name, $path, $suffix) = fileparse($input, qr/\.gff/);
    my $out_file = $path . $name . '.gtf';

    my $in_gff = Bio::FeatureIO->new(-file => "${input}",
                                     -format => 'GFF',
                                     -version => 3);
    my $out_gtf = FileHandle->new(">${out_file}");
    local $| = 1;
    ##my $out_gtf = new FileHandle();
    ##$out_gtf->open(">$out_file");

    my $features_written = 0;
    my $in_features = 0;
  FEATURES: while (my $feature = $in_gff->next_feature()) {
      $in_features++;
      ## the output handle is reset for every file
      ## According to UCSC, GTF file contains 9 column:
      ## <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes]
      ## An example working gtf line: (except the double spaces are tabs...)
      ## TcChr20-S  TriTrypDB  CDS  14765  15403  .  -  0    gene_id "cds_TcCLB.397937.10-1"; transcript_id "cds_TcCLB.397937.10-1";

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
      next FEATURES if ($primary_id eq 'supercontig' or
                        $primary_id eq 'region' or
                        $primary_id eq 'chromosome');
      my $string = qq"${seqid}\t${source}\t${primary_id}\t${start}\t${end}\t.\t${strand}\t${phase}\t";

      my $last_column = '';
      my $annot = $feature->{annotation};
      my $stringified;
      my $transcript_string = '';
      my $geneid_string = '';
      foreach my $key ($annot->get_all_annotation_keys()) {
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
      print $out_gtf $seq_string;
      $features_written++;
  } ## End iterating over every FEATURES
    $out_gtf->close();
    return($features_written);
}

=head2 C<Read_GFF>

 use Bio::Tools::GFF to gather information from a gff file.

 This pulls information of potential interest from a gff file.  In
 theory it could be trivially modified to be robust to different
 formats and other shenanigans.

=cut
sub Read_GFF {
    my ($class, %args) = @_;
    my $chromosomes = $args{chromosomes};
    my $gff = FileHandle->new("<$args{gff}");

    my $annotation_in = Bio::Tools::GFF->new(-fh => \$gff, -gff_version => 3);
    my $gff_out = {};
    print "Starting to read gff: $args{gff}\n";
  LOOP: while(my $feature = $annotation_in->next_feature()) {
      next LOOP unless ($feature->{_primary_tag} eq $args{gff_type});
      my $location = $feature->{_location};
      my $start = $location->start();
      my $end = $location->end();
      my $strand = $location->strand();
      my @ids = $feature->each_tag_value('ID');
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
          $id .= "${i} ";
      }
      $id =~ s/\s+$//g;
      my @gff_information = split(/\t+/, $gff_string);
      my $description_string = $gff_information[8];
      my $orf_chromosome = $gff_chr;
      my $annot = {
          id => $id,
          start => $start, ## Genomic coordinate of the start codon
          end => $end, ## And stop codon
          strand => $strand,
          description_string => $description_string,
          chromosome => $gff_chr,
      };
      $gff_out->{$gff_chr}->{$id} = $annot;
  } ## End looking at every gene in the gff file
    $gff->close();
    return($gff_out);
}

=head2 C<Sam2Bam>

 Sort, compress, index a sam file into a usable bam file.

 Used to invoke samtools to take the sam output from bowtie/bwa and
 convert it to an compressed-sorted-indexed bam alignment. This
 function just calls $class->Samtools(), but has a little logic
 to see if the invocation of this is for an extent .sam file or
 calling it on an existing .fastq(.gz), in which case one must
 assume it is being called on one or more files in the bowtie_out/
 directory and which start with $basename, include something like
 -trimmed-v0M1.sam.

=over

=item C<Arguments>

 input(required): Input sam file.
 species(required): Species name to hunt down the appropriate
  fasta/gff files.
 modules('samtools', 'bamtools'): Environment modules to load.

=item C<Invocation>

> cyoa --task convert --method sam2bam --input bob.sam --species hg38_91

=cut
sub Sam2Bam {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input'],);
    my @input_list = ();
    my $paths = $class->Get_Paths($options->{input});
    if ($options->{input}) {
        push(@input_list, $options->{input});
    } elsif (-r $options->{input} and $options->{input} =~ /\.sam$/) {
        push(@input_list, $options->{input});
    } elsif (-r $options->{input} and $options->{input} =~ /\.fastq$/) {
        if (-d "bowtie_out") {
            find({ wanted => sub { push(@input_list, "bowtie_out/$_") if ($_ =~ /\.sam$/); }, follow => 1 }, 'bowtie_out/');

        } else {
            foreach my $k (%{$options->{bt_args}}) {
                my $output_string = qq"bowtie_out/$paths->{jbasename}-${k}.sam";
                push(@input_list, $output_string);
            }
            my $bt = $class->Bio::Adventure::Map::Bowtie(%args);
        }
    } else {
        die("I don't know what to do without a .fastq file nor a .sam file.\n");
    }
    my $loaded = $class->Module_Load(modules => $options->{modules});
    my $sam = $class->Bio::Adventure::Convert::Samtools(%args,
        jdepends => $options->{jdepends},
        sam => \@input_list);
    return($sam);
}

=back

=head2 C<Samtools>

 Does the actual work of converting a sam to a bam file.

 $hpgl->Samtools() calls (in order): samtools view, samtools sort,
 and samtools index.  Upon completion, it invokes bamtools stats to
 see what the alignments looked like.

 It explicitly does not pipe one samtools invocation into the next,
 not for any real reason but because when I first wrote it, it
 seemed like the sorting was taking too long if I did not already
 have the alignments in a bam file.

=over

=item C<Arguments>

 input(required): Input sam file.
 species(required): Input species prefix for finding a genome.

=cut
sub Samtools {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        jmem => 30,
        jname => 'sam',
        jprefix => '',
        paired => 1,
        mismatch => 1,);
    my $input = $options->{input};
    my $output = $input;
    $output =~ s/\.sam$/\.bam/g;
    my $sorted_name = $input;
    $sorted_name =~ s/\.sam$//g;
    $sorted_name = qq"${sorted_name}-sorted";
    my $paired_name = $sorted_name;
    $paired_name =~ s/\-sorted/\-paired/g;
    my $workdir = dirname($input);
    ## Add a samtools version check because *sigh*
    my $samtools_version = qx"samtools 2>&1 | grep 'Version'";
    ## Start out assuming we will use the new samtools syntax.

    my $stderr = qq"${output}_samtools.stderr";
    my $stdout = qq"${output}_samtools.stdout";
    my $samtools_first = qq!
## If a previous sort file exists due to running out of memory,
## then we need to get rid of them first.
## hg38_100_genome-sorted.bam.tmp.0000.bam
if [[ -f "${output}.tmp.000.bam" ]]; then
  rm -f ${output}.tmp.*.bam
fi
samtools view -u -t $options->{libdir}/genome/$options->{species}.fasta \\
  -S ${input} -o ${output}  \\
  2>${stderr} \\
  1>${stdout}
!;

    my $samtools_second = qq"samtools sort -l 9 ${output} \\
  -o ${sorted_name}.bam \\
  2>>${stderr} \\
  1>>${stdout}";
    ## If there is a 0.1 in the version string, then use the old syntax.
    if ($samtools_version =~ /0\.1/) {
        $samtools_first = qq"samtools view -u \\
  -t $options->{libdir}/genome/$options->{species}.fasta \\
  -S ${input} 1>${output}";
        $samtools_second = qq"samtools sort -l 9 ${output} \\
  ${sorted_name} \\
  2>>${stderr} \\
  1>>${stdout}";
    }
    my $jstring = qq!
echo "Starting samtools"
if [[ -f "${output}" && -f "${input}" ]]; then
  echo "Both the bam and sam files exist, rerunning."
elif [[ -f "${output}" ]]; then
  echo "The output file exists, quitting."
  exit 0
elif [[ \! -f "${input}" ]]; then
  echo "Could not find the samtools input file."
  exit 1
fi
${samtools_first}
echo "First samtools command finished with \$?"
${samtools_second}
rm ${output}
rm ${input}
mv ${sorted_name}.bam ${output}
samtools index ${output} \\
  2>>${stderr} \\
  1>>${stdout}
echo "Second samtools command finished with \$?"
bamtools stats -in ${output} \\
  2>>${output}_samtools.stats 1>&2
echo "Bamtools finished with \$?"
!;
    if ($options->{paired}) {
        $jstring .= qq!
## The following will fail if this is single-ended.
samtools view -b -f 2 \\
  -o ${paired_name}.bam \\
  ${output} \\
  2>>${stderr} \\
  1>>${stdout}
samtools index ${paired_name}.bam \\
  2>>${stderr} \\
  1>>${stdout}
bamtools stats -in ${paired_name}.bam \\
  2>>${output}_samtools.stats 1>&2
!;
    } else {
        ## If it is not paired, just set the paired output name to the output name so that
        ## jobs which run this as part of a chain do not need to go looking for the paired output.
        ## I am considering this primarily in the case of gatk deduplication which demands
        ## properly paired reads or single ended; therefore if I set this here I do not need to
        ## check later for different inputs.
        $paired_name = basename($output, ('.bam'));
    }

    unless ($options->{mismatch}) {
        $jstring .= qq!
bamtools filter -tag XM:0 \\
  -in ${output} \\
  -out ${sorted_name}_nomismatch.bam \\
  2>>${output}_samtools.stats 1>&2
echo "bamtools filter finished with: \$?"
samtools index \\
  ${sorted_name}_nomismatch.bam \\
  2>>${stderr} \\
  1>>${stdout}
echo "final samtools index finished with: \$?"
!;
    }

    if ($options->{samtools_mapped}) {
        my $mapped = $input;
        $mapped =~ s/\.sam/_samtools_mapped\.fastq/g;
        $jstring .= qq"
samtools fastq ${output} -F 4 \\
  1>${mapped} \\
  2>${mapped}.stderr
xz -9e -f ${mapped}
";
    }

    if ($options->{samtools_unmapped}) {
        my $unmapped = $input;
        $unmapped =~ s/\.sam/_samtools_unmapped\.fastq/g;
        $jstring .= qq"
samtools fastq ${output} -f 4 \\
  1>${unmapped} \\
  2>${unmapped}.stderr
xz -9e -f ${unmapped}
";
    }

    my $comment = qq!## Converting the text sam to a compressed, sorted, indexed bamfile.
## Also printing alignment statistics to ${output}.stats
!;
    if ($options->{jdepends}) {
        $comment .= qq!
## This job depended on: $options->{jdepends}!;
    }
    my $jobname = qq"$options->{jname}_$options->{species}";
    my $samtools = $class->Submit(
        comment => $comment,
        depends => $options->{jdepends},
        input => $input,
        output => $output,
        paired => $options->{paired},
        paired_output => qq"${paired_name}.bam",
        stderr => $stderr,
        stdout => $stdout,
        postscript => $options->{postscript},
        prescript => $options->{prescript},
        jcpu => 4,
        jmem => $options->{jmem},
        jname => $jobname,
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jwalltime => '18:00:00',);
    return($samtools);
}

=back

=head1 AUTHOR - atb

Email abelew@gmail.com

=head1 SEE ALSO

L<samtools> L<Bio::FeatureIO> L<Bio::Tools::GFF> L<Bio::Seq> L<Bio::SeqIO>

=cut

1;
