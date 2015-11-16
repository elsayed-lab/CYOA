package HPGL;

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

    our %options = (
        debug => undef,
        inputdir => undef,
        indexes => 'indexes.txt',
        index_hash => {},
        handles => {},
        help => undef,
        reads_read => 0,
        reads_written => 0,
        written_by_index => {unsorted => 0,},
        );
    our $opt = GetOptions(
        "debug" => \$options{debug},
        "inputdir|d:s" => \$options{inputdir},
        "indexes|i:s" => \$options{indexes},
        "help|h" => \$options{help},
        );

    if ($options{help}) {
        print qq"This script is intended to split TNSeq reads as indexed by Yoann Le Breton.
It therefore makes some assumptions about the beginnings of reads and requires two arguments:
 --inputdir or -d <dir>:  A directory containing the raw sequences (compressed or not)
 --indexes or -i <file>:  A file containing the sequence barcodes and sample names, one per line.
";
        exit(0);
    }

    ## Make sure we have a directory to read sequences from
    if (!defined($options{inputdir}) or !-r $options{inputdir}) {
        my $term = new Term::ReadLine('>');
        my $attribs = $term->Attribs;
        $attribs->{completion_suppress_append} = 1;
        my $prompt = "Please provide a directory: ";
        my $OUT = $term->OUT || \*STDOUT;
        $options{inputdir} = $term->readline("Please provide a directory from which to read unsorted sequence: ");
    }

    ## Make sure we have a file to read indexes from, and gather them.
    if (!defined($options{indexes}) or !-r $options{indexes}) {
        my $term = new Term::ReadLine('>');
        my $attribs = $term->Attribs;
        $attribs->{completion_suppress_append} = 1;
        my $prompt = "Please provide an index file: ";
        my $OUT = $term->OUT || \*STDOUT;
        $options{indexes} = $term->readline("Please provide a file with the space delimited index->samples: ");
    }

    Read_Indexes();
    die("I need a directory containing some sequence.") unless($options{inputdir});
    Sort_Indexes();

    sub Read_Indexes {
        my %indexes = ();
        open(INDEX, "<$options{indexes}");
        while (my $line = <INDEX>) {
            chomp $line;
            next unless ($line =~ /^A|T|G|C|a|t|g|c/);
            my ($ind, $sample) = split(/\s+/, $line);
            $indexes{$ind} = $sample;
            print "Index $ind goes with sample $sample.\n" if ($options{debug});
        }
        close(INDEX);
        print "End Read_Indexes() \%indexes\n";
        print Dumper %indexes;
        $options{index_hash} = \%indexes;
    }

    sub Sort_Indexes {
        my $indexes_names = $options{index_hash};
        print Dumper $indexes_names if ($options{debug});
        foreach my $k (keys %{$indexes_names}) {
            my $filename = "$indexes_names->{$k}.fasta";
            open(my $tmp, ">$filename");
            print "Opening: $filename\n" if ($options{debug});
            $options{handles}->{$k} = $tmp;
        }
        my $unsorted_filename = "unknown.fasta";
        open(my $unsorted, ">$unsorted_filename");
        print "Opening unsorted: $unsorted_filename\n";
        $options{handles}->{unsorted} = $unsorted;

        my $cwd_dir = getcwd();
        my $searchdir = qq"$options{inputdir}";
        unless ($searchdir =~ /^\//) {
            $searchdir = qq"${cwd_dir}/${searchdir}";
        }
        print "TESTME: $searchdir\n";
        my @directory = ($searchdir);
        find(\&wanted, @directory);

        sub wanted {
            my $input = $File::Find::name;
            my $indexes_names = $options{index_hash};
            return(undef) unless ($input =~ /\.fastq/);
            print STDOUT "Starting to read: $input\n";
            open(INPUT, "<:gzip", "$input");
            my $in = new Bio::SeqIO(-fh => \*INPUT, -format => 'Fastq');
            while (my $in_seq = $in->next_seq()) {
                $options{reads_read}++;
                my $id = $in_seq->id();
                my $seq = $in_seq->seq();
                my $found = 0;
                foreach my $index (keys %{$indexes_names}) {
                    my $match = "^\\w{0,5}$index\\w+";
                    if ($seq =~ /$match/) {
                        $found++;
                        my $new_seq = $seq;
                        ##print STDOUT "Original Sequence: $new_seq  ";
                        if (!defined($options{written_by_index}->{$index})) {
                            $options{written_by_index}->{$index} = 1;
                        } else {
                            $options{written_by_index}->{$index}++;
                        }
                        ## Note the .. in the following line
                        ## Remember that I truncated the index in case of problems
                        ## Leaving two any-characters off of the beginning of the line attempts to ensure that
                        ## These are not picked up erroneously.
                        $new_seq =~ s/^\w{0,5}${index}..(\w+$)/$1/g;
                        $new_seq =~ s/^(\w+)AGGTTGGATGATA/$1/g;
                        ##print STDOUT "New Sequence: $new_seq  goes to $index $options{index_hash}->{$index}\n";
                        my $hand = $options{handles}{$index};
                        select($hand);
                        $options{reads_written}++;
                        if (($options{reads_written} % 10000) == 0) {
                            print STDOUT "Wrote $options{reads_written} from $input\n";
                            dump_sizes();
                        }
                        print $hand ">$id
$new_seq
";
                    }
                }
                if ($found == 0) {
                    my $hand = $options{handles}->{unsorted};
                    $options{written_by_index}->{unsorted}++;
                    select($hand);
                    print $hand ">$id
$seq
";
                }
            }
            print STDOUT "The total number of reads read in $input: $options{reads_read}.\n";
            print STDOUT "The total number of reads sorted in $input: $options{reads_written}.\n";
            foreach my $k (keys %{$options{written_by_index}}) {
                print STDOUT "Index $k had $options{written_by_index}->{$k} reads from $input.\n";
            }
            close(INPUT);
        }

        foreach my $k (keys %{$options{handles}}) {
            close($options{handles}{$k});
        }
    }
}


=item C<Essentiality_TAs>

    $hpgl->Essentiality_TAs(bam => 'lib1.bam');
    will read the $args{bam} bam file and look for reads which
    start/end with TA in the orientation expected for a successful
    insertion by the mariner retrotransposable element.

=cut
sub Essentiality_TAs {

    use strict;
    use warnings;
    use autodie qw":all";
    use Bio::SeqIO;
    use Bio::Tools::GFF;
    use Data::Dumper;
    use Bio::DB::Sam;
    use Getopt::Long;
    use File::Basename;
    use Cwd qw(abs_path getcwd);
    use Term::ReadLine;
    use Pod::Usage;

    =head1 NAME

        essentiality_tas

        =head1 SYNOPSIS

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


        ## An important note, I recently trimmed off the leading/ending TAs from these reads, as a result I need to check if the next/previous 2 genomic bases are TA rather than the first/last 2 bases.  In addition, I have to make some assumption regarding the stranded-ness of each hit...
        ## I am going to leave the current logic here and put an if around it

        our $cwd = getcwd;
    our %options = (
        name => undef,
        bam => undef,
        fasta => undef,
        output => undef,
        debug => undef,
        gff => undef,
        kept_tas => 1,
        type => "gene",
        help => undef,
        );
    our $data = [];

    Check_Options();

    if ($options{help}) {
        pod2usage();
    }

    my $genome_length = Check_Fasta();
    ## These 5 arrays will hold the lists of start/ends of each region of interest
    ## This is a somewhat dirty method for this, but I just wanted something easy
    my @inter_starts;
    my @inter_ends;
    my @cds_starts;
    my @cds_ends;
    my @names;
    Read_GFF();

    if ($options{kept_tas}) {
        Read_Sam_TAs();
    } else {
        Read_Sam_NoTAs();
    }

    foreach my $c (0 .. $#$data) {
        if ($data->[$c]->{fwd} or $data->[$c]->{rev} or $data->[$c]->{mismatch} or $data->[$c]->{ta}) {
            my $fwd = $data->[$c]->{fwd};
            my $rev = $data->[$c]->{rev};
            my $mis = $data->[$c]->{mismatch};
            print "$c\t$fwd\t$rev\t$mis\n";
        }
    }


    open(INTER_TAS_ONLY, ">$options{output}-inter");
    open(CDS_TAS_ONLY, ">$options{output}");
    my $file_start = qq"#type=COPY_NUMBER
Chromosome\tStart\tEnd\tFeature reads\t TAs
";
    print INTER_TAS_ONLY $file_start;
    print CDS_TAS_ONLY $file_start;

    foreach my $c (0 .. $#$data) {
        if ($data->[$c]->{ta}) {
            my $fwd = $data->[$c]->{fwd};
            my $rev = $data->[$c]->{rev};
            my $mis = $data->[$c]->{mismatch};
            my $total = $fwd + $rev;
            my $d = $c + 2;
            my $feature = "unfound";
            ## $c and $d outline the position of the TA
            ## So I want to figure out if the TA is inside an ORF or intercds, and its name...
          NAMES: for my $i (0 .. $#names) {
              if ($c >= $inter_starts[$i] and $c <= $inter_ends[$i]) {
                  $feature = "inter_$names[$i]";
                  last NAMES;
              } if ($c >= $cds_starts[$i] and $c <= $cds_ends[$i]) {
                  $feature = "cds_$names[$i]";
                  last NAMES;
              }
          } ## End searching for feature names.

            if ($feature =~ /^inter/) {
                print INTER_TAS_ONLY "$options{name}\t$c\t$d\t$feature\t$total\t1\n";
            } else {
                print CDS_TAS_ONLY "$options{name}\t$c\t$d\t$feature\t$total\t1\n";
            }
        }
    }
    close(INTER_TAS_ONLY);
    close(CDS_TAS_ONLY);

    sub Check_Options{
        my $opt = GetOptions(
            "fasta|f:s" => \$options{fasta},
            "output|o:s" => \$options{output},
            "name|n:s" => \$options{name},
            "debug|d:s" => \$options{debug},
            "gff|g:s" => \$options{gff},
            "bam|b:s" => \$options{bam},
            "type|t:s" => \$options{type},
            "kept|k:s" => \$options{kept_tas},
            "help" => \$options{help},
            );
        my $term = new Term::ReadLine('>');
        my $attribs = $term->Attribs;
        $attribs->{completion_suppress_append} = 1;
        my $OUT = $term->OUT || \*STDOUT;
        if (!defined($options{fasta})) {
            $options{fasta} = $term->readline("Please provide a fasta for this genome: ");
            $options{fasta} =~ s/\s+$//g;
            $options{fasta} =~ s/\~/${ENV{HOME}}/g;
            $options{fasta} = abs_path(qq"$options{fasta}");
        }
        if (!defined($options{bam})) {
            $options{bam} = $term->readline("Please provide an input bamfile: ");
            $options{bam} =~ s/\s+$//g;
            $options{bam} =~ s/\~/${ENV{HOME}}/g;
            $options{bam} = abs_path(qq"$options{bam}");
        }
        my $input_sam = abs_path(qq"$options{bam}");
        my $output_txt;
        if (!$options{output}) {
            $output_txt = $input_sam;
            my @suffix = ('.bam');
            my $output_base = basename($output_txt, @suffix);
            $options{output} = $output_base . '-count_by_base.txt';
        } else {
            $options{output} = abs_path(qq"$options{output}");
        }
    }

    sub Read_GFF {
        open(GFF, "<$options{gff}");
        my $annotations = new Bio::Tools::GFF(-fh => \*GFF, -gff_version => 3);
        my $start = 1;
        my $old_start = 1;
      LOOP: while(my $feature = $annotations->next_feature()) {
          next LOOP unless ($feature->{_primary_tag} eq $options{type});
          my $location = $feature->{_location};
          $old_start = $start;
          $start = $location->start();
          push(@inter_starts, $old_start);
          push(@inter_ends, $start);
          my $end = $location->end();
          push(@cds_starts, $start);
          push(@cds_ends, $end);
          my $strand = $location->strand();
          my @ids = $feature->each_tag_value("locus_tag");
          ##      print "TESTME: $start $end @ids\n";
          my $id = "";
          foreach my $i (@ids) {
              $id .= "$i ";
          }
          $id =~ s/\s+$//g;
          push(@names, $id);
      } ## End looking at every feature of the gff
        close(GFF);
    } ## End Read_GFF

    sub Check_Fasta {
        my $length = 0;
        my $sequences = 0;
        my $infile = new Bio::SeqIO(-file => $options{fasta}, -format => 'Fasta');
        while (my $in_seq = $infile->next_seq()) {
            $sequences++;
            $length += $in_seq->length;
            my $seq = $in_seq->seq;
            my @seq_array = split(//, $seq);
            my $count = 0;
            print "Allocating array of length $length\n";
          LOOP: foreach my $base (@seq_array) {
              my $next;
              if (defined($seq_array[$count + 1])) {
                  $next = $seq_array[$count + 1];
              } else {
                  print "Setting $count + 1 to 'N'\n";
                  $seq_array[$count + 1] = 'N';
                  $next = 'N';
                  last LOOP;
              }
              if ($base eq 'T' and $next eq 'A') {
                  $data->[$count]->{ta} = 1;
              } else {
                  $data->[$count]->{ta} = 0;
              }
              $data->[$count]->{nt} = $base;
              $data->[$count]->{fwd} = 0;
              $data->[$count]->{rev} = 0;
              $data->[$count]->{mismatch} = 0;
              $count++;
          }
            print "Finished allocating array.\n";
        }
        print "Found $sequences sequence(s) in $options{fasta} with $length characters.\n";
        return($length);
    }

    sub Read_Sam_NoTAs {
        my %args = @_;
        my $length = $args{length};
        my $input = $options{input};
        my $sam = new Bio::DB::Sam(-bam => $options{bam},
                                   -fasta => $options{fasta},);
        my @targets = $sam->seq_ids;
        my $num = scalar(@targets);
        my $bam = Bio::DB::Bam->open($options{bam});
        my $header = $bam->header;
        my $target_count = $header->n_targets;
        my $target_names = $header->target_name;
        my $align_count = 0;
        my $million_aligns = 0;
      BAMLOOP: while (my $align = $bam->read1) {
          $align_count++;
          ##	 last BAMLOOP if ($align_count > 20000);
          if (($align_count % 1000000) == 0) {
              $million_aligns++;
              print STDERR "Finished $million_aligns million alignments out of 52,259,285.\n";
          }
          my $seqid = $target_names->[$align->tid];
          ##	 my $start = $align->pos + 1;
          ##	 my $end = $align->calend;
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
          my $query_forward_hit = 0;
          my $query_reverse_hit = 0;
          ## Note that in the TAs version of this function, I check for TA in the beginning/end of the sequence, this is no longer possible.
          if ($strand == 1) {
              $query_forward_hit = 1;
          } elsif ($strand == -1) {
              $query_reverse_hit = 1;
          } else {
              print "Neither forward nor reverse!?\n";
              next BAMLOOP;
          }

          my $new_end = $end - 1;
          if ($new_end < 0) {
              print "The alignment appears to end before the start of the genome.\n";
              print "$start $end $cigar $seq $strand $qual\n";
              next BAMLOOP;
          }
          ## Now do the same check on the reference sequence (genome)
          my $ref_start_hit = 0;
          my $ref_end_hit = 0;
          if (!defined($data->[$start]->{nt})) {
              print "WTF UNDEF $start\n";
              $data->[$start]->{nt} = 'N';
          }
          if (!defined($data->[$start + 1]->{nt})) {
              print "WTF UNDEF $start\n";
              $data->[$start + 1]->{nt} = 'N';
          }
          if ($data->[$start-2]->{nt} eq 'T' and
              $data->[$start-1]->{nt} eq 'A') {
              $ref_start_hit = 1;
              $ref_end_hit = 0;
          }
          elsif ($data->[$end+1]->{nt} eq 'T'
                 and $data->[$end+2]->{nt} eq 'A') {
              $ref_start_hit = 0;
              $ref_end_hit = 1;
          }
          else {
              $ref_start_hit = 0;
              $ref_end_hit = 0;
          }

          if ($query_forward_hit == 1 and $query_reverse_hit == 0) {
              print "q_fwd " if ($options{debug});
              if ($ref_start_hit == 1) {
                  print "r_st $start Query: ${seq} Ref:$data->[$start-2]->{nt}$data->[$start-1]->{nt} $qual $cigar\n" if ($options{debug});
                  $data->[$start]->{fwd}++;
              } elsif ($ref_end_hit == 1) {
                  print "r_en $start Query: ${seq} Ref:$data->[$end+1]->{nt}$data->[$end+2]->{nt} $qual $cigar\n" if ($options{debug});
                  $data->[$start]->{rev}++;
              } else {
                  print "r_no $start Query: ${seq} Ref:$data->[$start-2]->{nt}$data->[$start-1]->{nt}|$data->[$end+1]->{nt}$data->[$end+2]->{nt} $qual $cigar\n" if ($options{debug});
                  $data->[$start]->{mismatch}++;
              }

          } elsif ($query_forward_hit == 0 and $query_reverse_hit == 1) {
              print "q_rev " if ($options{debug});
              if ($ref_start_hit == 1) {
                  print "r_st $new_end Query: ${seq} Ref:$data->[$start-2]->{nt}$data->[$start-1]->{nt} $qual $cigar\n" if ($options{debug});
                  $data->[$new_end]->{fwd}++;
              } elsif ($ref_end_hit == 1) {
                  print "r_en $new_end Query: ${seq} Ref:$data->[$end-1]->{nt}$data->[$end+2]->{nt} $qual $cigar\n" if ($options{debug});
                  $data->[$new_end]->{rev}++;
              } else {
                  print "r_no $new_end Query: ${seq} Ref:$data->[$new_end]->{nt}$data->[$end]->{nt}|$data->[$start-2]->{nt}$data->[$start-1]->{nt} $qual $cigar\n" if ($options{debug});
                  $data->[$new_end]->{mismatch}++;
              }
          } else {  ## something weird happened.
              print "Something weird happened...\n";
              next BAMLOOP;
          }
      }
    }

    sub Read_Sam_TAs {
        my %args = @_;
        my $length = $args{length};
        my $input = $options{input};
        my $sam = new Bio::DB::Sam(-bam => $options{bam},
                                   -fasta => $options{fasta},);
        my @targets = $sam->seq_ids;
        my $num = scalar(@targets);
        my $bam = Bio::DB::Bam->open($options{bam});
        my $header = $bam->header;
        my $target_count = $header->n_targets;
        my $target_names = $header->target_name;
        my $align_count = 0;
        my $million_aligns = 0;
      BAMLOOP: while (my $align = $bam->read1) {
          $align_count++;
          ##	 last BAMLOOP if ($align_count > 20000);
          if (($align_count % 1000000) == 0) {
              $million_aligns++;
              print STDERR "Finished $million_aligns million alignments out of 52,259,285.\n";
          }
          my $seqid = $target_names->[$align->tid];
          ##	 my $start = $align->pos + 1;
          ##	 my $end = $align->calend;
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
          if ($new_end < 0) {
              print "The alignment appears to end before the start of the genome.\n";
              print "$start $end $cigar $seq $strand $qual\n";
              next BAMLOOP;
          }
          ## Now do the same check on the reference sequence (genome)
          my $ref_forward_hit = 0;
          my $ref_reverse_hit = 0;
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
                  print "$start Mismatch reference forward hit $qual $cigar.\n" if ($options{debug});
                  $data->[$start]->{fwd}++;
                  $data->[$start]->{mismatch}++;
              } elsif ($ref_reverse_hit == 1) {
                  print "$new_end Mismatch reference reverse hit $qual $cigar.\n" if ($options{debug});
                  $data->[$new_end]->{rev}++;
                  $data->[$new_end]->{mismatch}++;
              } else {
                  print "$start $end No easily observable hit $qual $cigar.\n" if ($options{debug});
                  $data->[$start]->{mismatch}++;
                  $data->[$new_end]->{mismatch}++;
              }
          } elsif ($forward_hit == 1 and $reverse_hit == 0) {
              if ($ref_forward_hit == 1) {
                  print "$start fwd_only Query: ${seq} Ref:$data->[$start]->{nt}$data->[$start + 1]->{nt} $qual $cigar\n" if ($options{debug});
                  $data->[$start]->{fwd}++;
              } else {
                  print "$start Mismatch forward hit Query: ${seq} Ref:$data->[$start]->{nt}$data->[$start + 1]->{nt} $qual $cigar\n" if ($options{debug});
                  $data->[$start]->{mismatch}++;
              }
          } elsif ($forward_hit == 0 and $reverse_hit == 1) {
              if ($ref_reverse_hit == 1) {
                  print "$new_end rev_only Query: ${seq} Ref:$data->[$new_end]->{nt}$data->[$end]->{nt} $qual $cigar\n" if ($options{debug});
                  $data->[$new_end]->{rev}++;
              } else {
                  print "$new_end Mismatch reverse hit Query: ${seq} Ref:$data->[$start]->{nt}$data->[$start + 1]->{nt} $qual $cigar\n" if ($options{debug});
                  $data->[$new_end]->{mismatch}++;
              }
          } elsif ($forward_hit == 1 and $reverse_hit == 1) {
              if ($ref_forward_hit == 1 and $ref_reverse_hit == 1) {
                  print "$start $new_end fwd and reverse hits: Query: ${seq} Ref:$data->[$start]->{nt}$data->[$start + 1]->{nt} $data->[$new_end]->{nt}$data->[$end]->{nt} $qual $cigar\n" if ($options{debug});
                  $data->[$start]->{fwd}++;
                  ##		 $data->[$new_end]->{rev}++  ## Removing this because double count
              } elsif ($ref_forward_hit == 1) {
                  print "$start fwd hit confirmed: Query: ${seq} Ref:$data->[$start]->{nt}$data->[$start + 1]->{nt} $data->[$new_end]->{nt}$data->[$end]->{nt} $qual $cigar\n" if ($options{debug});
                  $data->[$start]->{fwd}++;
              } elsif ($ref_reverse_hit == 1) {
                  print "$new_end rev hit confirmed: Query: ${seq} Ref:$data->[$start]->{nt}$data->[$start + 1]->{nt} $data->[$new_end]->{nt}$data->[$end]->{nt} $qual $cigar\n" if ($options{debug});
                  $data->[$new_end]->{rev}++;
              } else {
                  print "$start $new_end Neither: Query: ${seq} Ref:$data->[$start]->{nt}$data->[$start + 1]->{nt} $data->[$new_end]->{nt}$data->[$end]->{nt} $qual $cigar\n" if ($options{debug});
                  $data->[$start]->{mismatch}++;
                  ##		 $data->[$new_end]->{mismatch}++;  ## No double count here either
              }
          } else {
              die("This else statement should never happen.\n");
          }
      }
    }
}

1;
