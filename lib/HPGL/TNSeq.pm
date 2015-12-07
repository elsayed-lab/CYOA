package HPGL;
use common::sense;
use autodie qw":all";
use diagnostics;
use PerlIO;
use PerlIO::gzip;
use File::Path qw"make_path";

use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Bio::DB::Sam;
use Cwd qw(abs_path getcwd);

=head1 Name
    HPGL::TNSeq - some functions for working with TNSeq data

=head1 SYNOPSIS

    The functions here were copy/pasted from some other code pieces I
    wrote while working on mgas 5448 tnseq data.  They have not been
    incorporated into the general HPGL.pm scheme yet.

=head2 Methods

=over 4

=item C<Sort_Indexes>

    $hpgl->Sort_Indexes(inputdir => '.');
    will assume there is a $(pwd)/indexes.txt with a list of
    sequencing index to sample names, it will then read all the fastq
    files in the inputdir and sort them into piles by index.  Most
    tasks/tools of course do not/should not need this, but the tnseq
    libraries used by the McIver lab are not quite standard.

=cut
sub Sort_Indexes {
    my $me = shift;
    my %args = @_;
    my $index_file = $me->{index_file};
    $index_file = $args{index_file} if (defined($args{index_file}));
    my $input = $me->{input};
    $input = $args{input} if (defined($args{input}));
    $input = '.' unless (defined($input));
    my $outdir = $me->{output_dir};
    $outdir = $args{output_dir} if (defined($args{output_dir}));
    $outdir = 'sorted' unless (defined($outdir));
    $me->{outdir} = $outdir;
    my $index_hash = {
    };

    if ($me->{help}) {
        print qq"This script is intended to split TNSeq reads as indexed by Yoann Le Breton.
It therefore makes some assumptions about the beginnings of reads and requires two arguments:
 --inputdir or -d <dir>:  A directory containing the raw sequences (compressed or not)
 --indexes or -i <file>:  A file containing the sequence barcodes and sample names, one per line.
";
        exit(0);
    }

    ## Make sure we have a directory to read sequences from
    if (!defined($input) or !-r $input) {
        my $term = new Term::ReadLine('>');
        my $attribs = $term->Attribs;
        $attribs->{completion_suppress_append} = 1;
        my $prompt = "Please provide a directory of fastq files or an individual file: ";
        my $OUT = $term->OUT or \*STDOUT;
        $input = $term->readline("Please provide a directory of fastq files or an individual file: ");
    }

    ## Make sure we have a file to read indexes from, and gather them.
    if (!defined($index_file) or !-r $index_file) {
        my $term = new Term::ReadLine('>');
        my $attribs = $term->Attribs;
        $attribs->{completion_suppress_append} = 1;
        my $prompt = "Please provide an index file: ";
        my $OUT = $term->OUT || \*STDOUT;
        $index_file = $term->readline("Please provide a file with delimited index->samples: ");
    }

    $me->{indexes} = $me->Read_Indexes(index_file => $index_file);
    if (-f $input) {
        $me->Sort_TNSeq_File(file => $input);
    } elsif (-d $input) {
        $me->Sort_TNSeq_Dir(dir => $input);
    } else {
        die("I need a directory containing some sequence.") unless($input);
    }

    my $log = new FileHandle;
    $log->open(">$me->{outdir}/tnseq_sorting_stats.txt");
    my $message = qq"
The total number of reads read in $input: $me->{indexes}->{total}->{read}.
The total number of reads sorted in $input: $me->{indexes}->{total}->{written}.
";
    print STDOUT $message;
    print $log $message;
    my %final = %{$me->{indexes}};
    foreach my $k (keys %final) {
        next if ($k eq 'total');
        $message = qq"Index $k had $me->{indexes}->{$k}->{written} reads written from $input to $me->{indexes}->{$k}->{filename}.\n";
        print STDOUT $message;
        print $log $message;
        $me->{indexes}->{$k}->{handle}->close();
    }
    $log->close();
    return($me->{indexes}->{total});
}

=item C<Sort_TNSeq_File>

  $hpgl->Sort_TNSeq_File(file => 'filename.fastq.gz', indexes =>*Index Hash*)
  is called by Sort_Indexes() when a single file is to be demultiplexed.
  It assumes the input file is gzipped and will use the PerlIO::gzip
  filter to uncompress the data.  It will then extract the sequences
  per the rules set up in the index hash.

=cut
sub Sort_TNSeq_File {
    my $me = shift;
    my %args = @_;
    my $input = $args{file};
    my $indexes = $me->{indexes};

    use String::Approx;
    return(undef) unless ($input =~ /\.fastq/);
    print STDOUT "Starting to read: $input\n";
    open(INPUT, "lesspipe $input |");  ## PerlIO::gzip is failing on some files.
    ##open(INPUT, "<:gzip", "$input");
    my $in = new Bio::SeqIO(-fh => \*INPUT, -format => 'Fastq');
    while (my $in_seq = $in->next_dataset()) {
        $me->{indexes}->{total}->{read}++;
        print "Read $me->{indexes}->{total}->{read} entries.\n" if ($me->{debug});
        my $id = $in_seq->{'-descriptor'};
        my $sequence = $in_seq->{'-seq'};
        my $qual = $in_seq->{'-raw_quality'};
        my $found = 0;
        foreach my $index (keys %{$indexes}) {
            my $match = "^\\w{0,5}$index\\w+";
            if ($sequence =~ m/$match/) {
                $found++;
                my $new_seq = $sequence;
                ## Remember that I truncated the index in case of problems
                ## Leaving two any-characters off of the beginning of the line attempts to ensure that
                ## These are not picked up erroneously.
                $new_seq =~ s/^\w{0,5}${index}..(\w+$)/$1/g;
                ## The next 5 lines are an overly verbose way of extracting the positions which matched the (\w+) above
                ## and using those positions to pull out the same substring from the quality scores
                my @match_start = @-;
                my @match_end = @+;
                my $match_start_pos = $match_start[1];
                my $match_end_pos = $match_end[1];
                my $new_qual = substr($qual, $match_start_pos, $match_end_pos);
                ## The following search was added due to some weird PCR artifacts in some datasets.
                ## This may complicate submitting this data to SRA.
                if ($me->{tnseq_trim}) {
                    my $pcr_match = "ACAGGTTGGATGATA";
                    if ($sequence =~ m/$pcr_match/) {
                        my $newer_seq = $new_seq;
                        $newer_seq =~ s/^(\w+)${pcr_match}(\w+)$/$1/g;
                        @match_start = @-;
                        @match_end = @+;
                        my $pcr_match_start = $match_start[1];
                        my $pcr_match_end = $match_end[1];
                        my $newer_qual = substr($qual, $pcr_match_start, $pcr_match_end);
                        $new_seq = $newer_seq;
                        $new_qual = $newer_qual;
                    }
                }

                $me->{indexes}->{total}->{written}++;
                if (($me->{indexes}->{total}->{written} % 10000) == 0) {
                    print "Wrote $me->{indexes}->{total}->{written} reads from $input\n";
                    ## dump_sizes();
                }
                my $fastq_string = qq"\@$id
$new_seq
+
$new_qual
";
                my $handle = $indexes->{$index}->{handle};
                print $handle $fastq_string;
                $indexes->{$index}->{written}++;
            } ## End looking for an individual match from %indexes
        }  ## End recursing over %indexes
        if ($found == 0) {
            my $handle = $indexes->{unknown}->{handle};
            my $fastq_string = qq"\@$id
$sequence
+
$qual
";
            print $handle $fastq_string;
            $indexes->{unknown}->{written}++;
        }
    } ## End reading the fastq file
    $me->{indexes} = $indexes;
    close(INPUT);
}

sub Sort_TNSeq_File_Approx {
    my $me = shift;
    my %args = @_;
    my $input = $args{file};
    my $indexes = $me->{indexes};

    use String::Approx qw"amatch";
    return(undef) unless ($input =~ /\.fastq/);
    print STDOUT "Starting to read: $input\n";
    open(INPUT, "lesspipe $input |");  ## PerlIO::gzip is failing on some files.
    ##open(INPUT, "<:gzip", "$input");
    my $in = new Bio::SeqIO(-fh => \*INPUT, -format => 'Fastq');
    my $count = 0;
    READS: while (my $in_seq = $in->next_dataset()) {
        $count++;
        $me->{indexes}->{total}->{read}++;
        print "Read $me->{indexes}->{total}->{read} entries.\n" if ($me->{debug});
        my $id = $in_seq->{'-descriptor'};
        my $sequence = $in_seq->{'-seq'};
        my $qual = $in_seq->{'-raw_quality'};
        my $found = 0;
        my @index_list = ();
        foreach my $index (keys %{$indexes}) {
            push(@index_list, $index) unless ($index eq 'total' or $index eq 'unknown' or $index eq 'ambiguous');
        }
        my $found = 0;
        my $found_id = '';
        my $match_substring = substr($sequence, 0, 9);
        foreach my $index (@index_list) {
            ##print "TESTME2: $match_substring vs $index\n";
            my @matches = amatch($index, [ 'S1' ], ($match_substring));
            if (scalar(@matches) == 1) {
                $found++;
                $found_id = $index;
            }
        } ## End each index
        if ($found == 1) {
            my $new_sequence = substr($sequence, 9,);
            my $new_quality = substr($qual, 9,);
            $me->{indexes}->{total}->{written}++;
            if (($me->{indexes}->{total}->{written} % 10000) == 0) {
                print "Wrote $me->{indexes}->{total}->{written} reads from ${input}\n";
            }
            my $fastq_string = qq"\@$id
$new_sequence
+
$new_quality
";
            my $handle = $indexes->{$found_id}->{handle};
            print $handle $fastq_string;
            $indexes->{$found_id}->{written}++;
        } elsif ($found > 1) {
            my $handle = $indexes->{ambiguous}->{handle};
            my $fastq_string = qq"\@$id
$sequence
+
$qual
";
            print $handle $fastq_string;
            $indexes->{ambiguous}->{written}++;
        } elsif ($found < 1) {
            my $handle = $indexes->{unknown}->{handle};
            my $fastq_string = qq"\@$id
$sequence
+
$qual
";
            print $handle $fastq_string;
            $indexes->{unknown}->{written}++;
        }
    } ## End reading the fastq file
    $me->{indexes} = $indexes;
    close(INPUT);
}

=item C<Sort_TNSeq_Dir>

  $hpgl->Sort_TNSeq_File(dir => '.', indexes =>*Index Hash*)
  is called by Sort_Indexes() when a directory of fastq files are to
  be demultiplexed. It uses File::Find to find all the .fastq.gz files
  in the target directory (recursively) which do not have the name of
  the output directory in the filename.
  It then passes each file to Sort_TNSeq_File()

=cut
sub Sort_TNSeq_Dir {
    my $me = shift;
    my %args = @_;

    my $cwd_dir = getcwd();
    my $searchdir = qq"$args{dir}";
    unless ($searchdir =~ /^\//) {
        $searchdir = qq"${cwd_dir}/${searchdir}";
    }
    print "Searching: $searchdir  for files to read.\n";
    my @directory = ($searchdir);
    my @file_list = ();
    find(sub { push(@file_list, $File::Find::name) if ($File::Find::name =~ /\.fastq\.gz/ and $file::Find::name !=~ /$me->{outdir}/); }, @directory);
    foreach my $file (@file_list) {
        next if ($file =~ /$me->{outdir}/);
        ## This might be incorrect, CHECKME!
        $file =~ s/\/\.\//\//;
        print "Processing $file\n";
        $me->Sort_TNSeq_File_Approx(file => $file);
    }
}

=item C<Read_Indexes>

  $h->Read_Indexes()
  reads a file containing the indexes and sets up the data structure
  used to keep count of the reads demultiplexed along with the set
  of file handles to which to write the data.

=cut
sub Read_Indexes {
    my $me = shift;
    my %args = @_;

    my $unknown = new FileHandle;
    my $ambiguous = new FileHandle;
    make_path($me->{outdir}) if (!-d $me->{outdir});
    unlink("$me->{outdir}/unknown.fastq.gz") if (-r "$me->{outdir}/unknown.fastq.gz");
    unlink("$me->{outdir}/ambiguous.fastq.gz") if (-r "$me->{outdir}/ambiguous.fastq.gz");
    $unknown->open("|gzip >> $me->{outdir}/unknown.fastq.gz");
    $ambiguous->open("|gzip >> $me->{outdir}/ambiguous.fastq.gz");
    my $indexes = {
        total => {
            read => 0,
            written => 0,
        },
        unknown => {
            name => 'unknown',
            filename => 'uknown.fastq',
            handle => $unknown,
            written => 0,
        },
        ambiguous => {
            name => 'ambiguous',
            filename => 'ambiguous.fastq',
            handle => $ambiguous,
            written => 0,
        },
    };

    my $index_file = new FileHandle;
    $index_file->open("<$args{index_file}");
    while (my $line = <$index_file>) {
        chomp $line;
        next unless ($line =~ /^A|T|G|C|a|t|g|c/);
        my ($ind, $sample) = split(/\s+|\,|;/, $line);
        my $output_filename = qq"$me->{outdir}/${sample}.fastq.gz";
        unlink($output_filename) if (-r $output_filename);
        my $handle = new FileHandle;
        $handle->open("|gzip >> $output_filename");
        print "Opened for writing: $output_filename\n" if ($me->{debug});
        $indexes->{$ind} = {
            name => $sample,
            filename => $output_filename,
            handle => $handle,
            written => 0,
        };
    }
    $index_file->close();
    return($indexes);
}



=item C<Essentiality_TAs>

    $hpgl->Essentiality_TAs(bam => 'lib1.bam');
    will read the $args{bam} bam file and look for reads which
    start/end with TA in the orientation expected for a successful
    insertion by the mariner retrotransposable element.

    essentiality_tas
    essentiality_tas -f genome.fasta -o output.txt -n genome_name -g genome.gff -b bamfile.bam

        Other options:  The following are possible
        --name  |-n : The name to print in the text file (required)
        --bam   |-b : The bam alignment (required)
        --fasta |-f : The fasta genome  (required)
        --output|-o : The output file   (required)
        --debug |-d : Print debugging information?
        --gff   |-g : The gff annotation (required)
        --kept_tas  : Keep the intermediate tas file (useful for coverage data)
        --type  |-t : Type of entry in the gff file to use (gene by default)
        --help      : Print this help information


=cut
sub Essentiality_TAs {
    my $me = shift;
    my %args = @_;

    my $output_directory = qq"outputs/essentiality";
    make_path($output_directory);
    ## A smarter way of handling counting inter-cds regions would be to use the Read_GFF() functions to pull the intercds regions rather than the hackish @inter_starts @inter_ends.
    $me->Check_Options(['species','input']);
    my $input = $me->{input};
    $input = $args{input} if (defined($args{input}));

    print STDERR "Hey, debug is on! $me->{debug}\n" if ($me->{debug});
    my $input_name = basename($input, ('.bam'));
    ## This is a bit of a mess, I need to clean up how these arguments are passed.
    my $output = $me->{output};
    $output = qq"${input}_out" unless (defined($args{output}));

    my $genome_fasta = qq"$me->{libdir}/$me->{libtype}/$me->{species}.fasta";
    my $genome_gff = $genome_fasta;
    $genome_gff =~ s/\.fasta/\.gff/;
    my $genome = $me->Allocate_Genome(fasta => $genome_fasta);
    my @data = @{$genome->{data}};
    ## These 5 arrays will hold the lists of start/ends of each region of interest
    ## This is a somewhat dirty method for this, but I just wanted something easy
    my $cds_gff = $me->Read_Genome_GFF(gff => $genome_gff, feature_type => 'gene');
    my $chr_name = $cds_gff->{stats}->{chromosomes}->[0];  ## Grab the name of the first chromosome
    my @names = @{$cds_gff->{stats}->{feature_names}};
    my @cds_starts = @{$cds_gff->{stats}->{cds_starts}};
    my @cds_ends = @{$cds_gff->{stats}->{cds_ends}};
    my @inter_starts = @{$cds_gff->{stats}->{inter_starts}};
    my @inter_ends = @{$cds_gff->{stats}->{inter_ends}};

    print STDERR "Counting TAs in ${input}, this takes a while.\n";
    my $datum = $me->Count_TAs(data => \@data);
    @data = @{$datum};

    my $tas_file = qq"${output_directory}/${input_name}_tas.txt";
    print STDERR "Printing list of #TAs observed by position to ${tas_file}.\n";
    my $ta_count = new FileHandle;
    $ta_count->open(">${tas_file}");
    foreach my $c (0 .. $#data) {
        if ($data[$c]->{fwd} or $data[$c]->{rev} or $data[$c]->{mismatch} or $data[$c]->{ta}) {
            my $fwd = $data[$c]->{fwd};
            my $rev = $data[$c]->{rev};
            my $mis = $data[$c]->{mismatch};
            print $ta_count "$c\t$fwd\t$rev\t$mis\n";
        }
    }
    $ta_count->close();

    my $inter_tas = new FileHandle;
    my $cds_tas = new FileHandle;
    $inter_tas->open(">${output_directory}/${input_name}_interCDS_tas.txt");
    $cds_tas->open(">${output_directory}/${input_name}_gene_tas.txt");
    my $file_start = qq"#type=COPY_NUMBER
Chromosome\tStart\tEnd\tFeature reads\t TAs
";
    print $inter_tas $file_start;
    print $cds_tas $file_start;

    print STDERR "Printing #TAs observed by (inter)CDS region to ${input_name}_(gene|interCDS)_tas.txt\n";
    foreach my $c (0 .. $#data) {
        if ($data[$c]->{ta}) {
            my $fwd = $data[$c]->{fwd};
            my $rev = $data[$c]->{rev};
            my $mis = $data[$c]->{mismatch};
            my $total = $fwd + $rev;
            my $d = $c + 2;
            my $feature = "unfound";
            ## $c and $d outline the position of the TA
            ## So I want to figure out if the TA is inside an ORF or intercds, and its name...
          NAMES: for my $i (0 .. $#names) {
              ## print "TESTME: inter_start: $inter_starts[$i] inter_end: $inter_ends[$i] vs. $c CDS $cds_starts[$i] $cds_ends[$i]\n";
              if ($c >= $inter_starts[$i] and $c <= $inter_ends[$i]) {
                  $feature = "inter_$names[$i]";
                  last NAMES;
              } if ($c >= $cds_starts[$i] and $c <= $cds_ends[$i]) {
                  $feature = "cds_$names[$i]";
                  last NAMES;
              }
          } ## End searching for feature names.
            ## Now look for inter vs cds features
            if ($feature =~ /^inter/) {
                print $inter_tas "${chr_name}\t$c\t$d\t$feature\t$total\t1\n";
            } else {
                print $cds_tas "${chr_name}\t$c\t$d\t$feature\t$total\t1\n";
            }
        }
    }
    $inter_tas->close();
    $cds_tas->close();
}

sub Count_TAs {
    my $me = shift;
    my %args = @_;
    my $length = $args{length};
    my $input = $me->{input};
    my $data = $args{data};
    my $sam = new Bio::DB::Sam(-bam => ${input},
                               -fasta => $me->{fasta},);
    my @targets = $sam->seq_ids;
    my $num = scalar(@targets);
    my $bam = Bio::DB::Bam->open(${input});
    my $header = $bam->header;
    my $target_count = $header->n_targets;
    my $target_names = $header->target_name;
    my $align_count = 0;
    my $million_aligns = 0;
    my $alignstats = qx"samtools idxstats $input";
    my @alignfun = split(/\n/, $alignstats);
    my @aligns = split(/\t/, $alignfun[0]);
    my @unaligns = split(/\t/, $alignfun[1]);
    my $number_reads = $aligns[2] + $unaligns[3];
    print STDERR "There are $number_reads alignments in ${input} made of $aligns[2] aligned reads and $unaligns[3] unaligned reads.\n";
  BAMLOOP: while (my $align = $bam->read1) {
      $align_count++;
      ##if ($me->{debug}) {  ## Stop after a relatively small number of reads when debugging.
      ##    last BAMLOOP if ($align_count > 200);
      ##}
      if (($align_count % 1000000) == 0) {
          $million_aligns++;
          print STDERR "Finished $million_aligns million alignments out of ${number_reads}.\n";
      }
      my $seqid = $target_names->[$align->tid];
      ## my $start = $align->pos + 1;
      ## my $end = $align->calend;
      my $start = $align->pos;
      my $end = $align->calend - 1;
      my $cigar = $align->cigar_str;
      my $strand = $align->strand;
      my $seq = $align->query->dna;
      my $qual= $align->qual;

      my @seq_array = split(//, $seq);
      ## At this point there are a few possibilities:
      ## 1.  We have only a forward hit and it is correct
      ## 2.  We have only a reverse hit and it is correct
      ## 3.  We have both and they are correct
      ###  The rest of the possibilities come up because we allow 1 mismatch in the alignment
      ## 4.  It is a forward hit but off by one base due to mismatch
      ## 5.  It is a reverse hit but off by one base due to mismatch
      ## 6.  It looks like both, but only one is correct due to mismatch
      ## 7.  No hit
      my $forward_hit = 0;
      my $reverse_hit = 0;
      if ($seq_array[0] eq 'T' and $seq_array[1] eq 'A') {
          $forward_hit = 1;
      }
      if ($seq_array[-1] eq 'A' and $seq_array[-2] eq 'T') {
          $reverse_hit = 1;
      }
      my $new_end = $end - 1;
      if ($cigar eq '') {
          ## print "This did not align.\n";
          next BAMLOOP;
      }
      if ($new_end < 0) {
          ##print "The alignment appears to end before the start of the genome.\n"; ## if ($me->{debug});
          ##print "start:$start end:$end cigar:$cigar seq:$seq $strand $qual\n"; ## if ($me->{debug});
          next BAMLOOP;
      }
      ## Now do the same check on the reference sequence (genome)
      my $ref_forward_hit = 0;
      my $ref_reverse_hit = 0;
      ## $data was passed to this function and
      if (!defined($data->[$start]->{nt})) {
          print "WTF UNDEF $start\n";
          $data->[$start]->{nt} = 'N';
      }
      if (!defined($data->[$start + 1]->{nt})) {
          print "WTF UNDEF $start\n";
          $data->[$start + 1]->{nt} = 'N';
      }
      if ($data->[$start]->{nt} eq 'T' and $data->[$start + 1]->{nt} eq 'A') {
          $ref_forward_hit = 1;
      }
      if ($data->[$new_end]->{nt} eq 'T' and $data->[$end]->{nt} eq 'A') {
          $ref_reverse_hit = 1;
      }

      if ($forward_hit == 0 and $reverse_hit == 0) {
          if ($ref_forward_hit == 1) {
              print "$start Mismatch reference forward hit $qual $cigar.\n" if ($me->{debug});
              $data->[$start]->{fwd}++;
              $data->[$start]->{mismatch}++;
          } elsif ($ref_reverse_hit == 1) {
              print "$new_end Mismatch reference reverse hit $qual $cigar.\n" if ($me->{debug});
              $data->[$new_end]->{rev}++;
              $data->[$new_end]->{mismatch}++;
          } else {
              print "$start $end No easily observable hit $qual $cigar.\n" if ($me->{debug});
              $data->[$start]->{mismatch}++;
              $data->[$new_end]->{mismatch}++;
          }
      } elsif ($forward_hit == 1 and $reverse_hit == 0) {
          if ($ref_forward_hit == 1) {
              print "$start fwd_only Query: ${seq} Ref:$data->[$start]->{nt}$data->[$start + 1]->{nt} $qual $cigar\n" if ($me->{debug});
              $data->[$start]->{fwd}++;
          } else {
              print "$start Mismatch forward hit Query: ${seq} Ref:$data->[$start]->{nt}$data->[$start + 1]->{nt} $qual $cigar\n" if ($me->{debug});
              $data->[$start]->{mismatch}++;
          }
      } elsif ($forward_hit == 0 and $reverse_hit == 1) {
          if ($ref_reverse_hit == 1) {
              print "$new_end rev_only Query: ${seq} Ref:$data->[$new_end]->{nt}$data->[$end]->{nt} $qual $cigar\n" if ($me->{debug});
              $data->[$new_end]->{rev}++;
          } else {
              print "$new_end Mismatch reverse hit Query: ${seq} Ref:$data->[$start]->{nt}$data->[$start + 1]->{nt} $qual $cigar\n" if ($me->{debug});
              $data->[$new_end]->{mismatch}++;
          }
      } elsif ($forward_hit == 1 and $reverse_hit == 1) {
          	     if ($ref_forward_hit == 1 and $ref_reverse_hit == 1) {
		 print "$start $new_end fwd and reverse hits: Query: ${seq} Ref:$data->[$start]->{nt}$data->[$start + 1]->{nt} $data->[$new_end]->{nt}$data->[$end]->{nt} $qual $cigar\n" if ($me->{debug});
		 $data->[$start]->{fwd}++;
##		 $data->[$new_end]->{rev}++  ## Removing this because double count
	     } elsif ($ref_forward_hit == 1) {
		 print "$start fwd hit confirmed: Query: ${seq} Ref:$data->[$start]->{nt}$data->[$start + 1]->{nt} $data->[$new_end]->{nt}$data->[$end]->{nt} $qual $cigar\n" if ($me->{debug});
		 $data->[$start]->{fwd}++;
	     } elsif ($ref_reverse_hit == 1) {
		 print "$new_end rev hit confirmed: Query: ${seq} Ref:$data->[$start]->{nt}$data->[$start + 1]->{nt} $data->[$new_end]->{nt}$data->[$end]->{nt} $qual $cigar\n" if ($me->{debug});
		 $data->[$new_end]->{rev}++;
	     } else {
		 print "$start $new_end Neither: Query: ${seq} Ref:$data->[$start]->{nt}$data->[$start + 1]->{nt} $data->[$new_end]->{nt}$data->[$end]->{nt} $qual $cigar\n" if ($me->{debug});
		 $data->[$start]->{mismatch}++;
##		 $data->[$new_end]->{mismatch}++;  ## No double count here either
	     }
	 } else {
	     die("This else statement should never happen.\n");
	 }
  }  ## End reading each bam entry
    return($data);
} ## End Count_Tas

sub Allocate_Genome {
    my $me = shift;
    my %args = @_;
    my $genome = {
        length => 0,
        data => [],
        sequences => 0,
    };
    my $data = [];
    my $infile = new Bio::SeqIO(-file => $args{fasta}, -format => 'Fasta');
    while (my $in_seq = $infile->next_seq()) {
        $genome->{sequences}++;
        $genome->{length} = $genome->{length} + $in_seq->length;
        my $seq = $in_seq->seq;
        my @seq_array = split(//, $seq);
        my $count = 0;
        print "Allocating array of length $genome->{length}\n";
      LOOP: foreach my $base (@seq_array) {
          $data->[$count]->{nt} = $base;
          my $next;

          if (defined($seq_array[$count + 1])) {
              $next = $seq_array[$count + 1];
          } else { ## This should only happen at the end of the genome.
              print "Setting $count + 1 to 'N', end of genome.\n";
              $seq_array[$count + 1] = 'N';
              $next = 'N';
              last LOOP;
          }
          if ($base eq 'T' and $next eq 'A') {
              $data->[$count]->{ta} = 1;
          } else {
              $data->[$count]->{ta} = 0;
          }
          $data->[$count]->{fwd} = 0;
          $data->[$count]->{rev} = 0;
          $data->[$count]->{mismatch} = 0;
          $count++;
      }
        print "Finished allocating array.\n";
    }
    $genome->{data} = $data;
    print "Found $genome->{sequences} sequence(s) in $me->{fasta} with $genome->{length} characters.\n";
    return($genome);
}

1;

## EOF;
