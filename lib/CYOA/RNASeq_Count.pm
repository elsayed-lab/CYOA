package CYOA;
use Bio::Tools::GFF;
use String::Approx qw"amatch";

=head1 NAME
    CYOA::RNASeq_Count - Perform Sequence alignments counting with HTSeq

=head1 SYNOPSIS

    use CYOA;
    my $hpgl = new CYOA;
    $hpgl->HT_Multi();

=head2 Methods

=over 2

=item C<HT_Multi>

    some stuff here

=cut
sub HT_Multi {
    my $me = shift;
    my %args = @_;
    my @jobs = ();
    $me->Check_Options(args => \%args, needed => ["species", "htseq_type"]);
    my $basename = $me->{basename};
    my $species = $me->{species};
    my $script_suffix = qq"";
    if ($args{suffix}) {
        $script_suffix = $args{suffix};
    }
    my %gff_types = (
        misc => "transcript",
        rintron => "exon",
        nmd => "exon",
        linc => "gene",
        pseudo => "exon",
        antisense => "exon",
        exon => "exon",
        fiveputr => "five_prime_UTR",
        threeputr => "three_prime_UTR",
        sn => "exon",
        mi => "exon",
        sno => "exon",
        rrna => "gene",
        interCDS => "CDS",
        operons => "gene",
        all => $me->{htseq_type},
    );
    my $htseq_runs = 0;
    my $input = $args{input} || $me->{input};
##    my $aligntype = $args{aligntype};
    ## my $out_prefix = $input;
    my $out_basename = basename($input, (".bam", ".sam"));
    my $out_suffixdir = dirname($input);  ## The directory name we will change from bowtie/tophat/whatever to htseq
    my $out_prefixdir = dirname($out_suffixdir);  ## This we will keep the same
    my $new_outdir = qq"${out_prefixdir}/htseq";
    make_path($new_outdir, {verbose => 0}) unless (-r $new_outdir);
    foreach my $type (keys %gff_types) {
        my $gff = qq"$me->{libdir}/genome/${species}_${type}.gff";
        my $gtf = $gff;
        $gtf =~ s/\.gff/\.gtf/g;
        my $output = qq"${new_outdir}/${out_basename}_${type}.count";
        if (-r "$gff") {
            print "Found $gff, performing htseq with it.\n";
            my $ht = $me->HTSeq(htseq_type => $gff_types{$type},
                                htseq_gff => $gff,
                                jobname => "ht${type}${script_suffix}_${species}",
                                depends => $args{depends},
                                job_prefix => $args{job_prefix},
                                qsub_queue => 'throughput',
                                input => $input,
                                output => $output,
                                prescript => $args{prescript},
                                postscript => $args{postscript},
                                suffix => $args{suffix},
                );
            push(@jobs, $ht);
            $htseq_runs++;
        } elsif (-r "$gtf") {
            print "Found $gtf, performing htseq with it.\n";
            my $ht = $me->HTSeq(htseq_type => $gff_types{$type},
                                htseq_gff => $gtf,
                                jobname => "ht${type}${script_suffix}_${species}",
                                depends => $args{depends},
                                qsub_queue => 'throughput',
                                job_prefix => $args{job_prefix},
                                input => $input,
                                output => $output,
                                prescript => $args{prescript},
                                postscript => $args{postscript},
                                suffix => $args{suffix},
                );
            push(@jobs, $ht);
            $htseq_runs++;
        } else {
            print "Did not find ${gff} nor ${gtf}, skipping htseq_${type}.\n";
        }
    } ## End foreach type

    ## Also perform a whole genome count
    my $gff = qq"$me->{libdir}/genome/${species}.gff";
    my $gtf = $gff;
    $gtf =~ s/\.gff/\.gtf/g;
    my $output = $input;
    $output =~ s/\.bam$/\.count/g;
    if (-r "$gff") {
        print "Found $gff, performing htseq_all with it.\n";
        my $ht = $me->HTSeq(htseq_type => $gff_types{all},
                            htseq_gff => $gff,
                            jobname => "htall${script_suffix}_${species}",
                            depends => $args{depends},
                            input => $input,
                            output => $output,
                            suffix => $args{suffix},
                            prescript => $args{prescript},
                            postscript => $args{postscript},
            );
        push(@jobs, $ht);
    } elsif (-r "$gtf") {
        print "Found $gtf, performing htseq_all with it.\n";
        my $ht = $me->HTSeq(htseq_type => $gff_types{all},
                            htseq_gff => $gff,
                            jobname => "htall${script_suffix}_${species}",
                            depends => $args{depends},
                            input => $input,
                            output => $output,
                            suffix => $args{suffix},
                            prescript => $args{prescript},
                            postscript => $args{postscript},
            );
        push(@jobs, $ht);
    } else {
        print "Did not find ${gff} nor ${gtf}, not running htseq_all.\n";
    }
    return(\@jobs);
}

sub HT_Types {
    my $me = shift;
    my %args = @_;
    my $reader = new FileHandle;
    my $my_type = $args{type};
    $my_type = "" unless($my_type);
    $reader->open("<$args{annotation}");
    my $annotation_in = new Bio::Tools::GFF(-fh => $reader);
    my $gff_out = {};
    print "Checking types of $args{annotation}\n";
    my $max_lines = 100000;
    my $count = 0;
    my %found_types = ();
  LOOP: while(my $feature = $annotation_in->next_feature()) {
      $count++;
      my $type = $feature->{_primary_tag};
      if ($found_types{$type}) {
          $found_types{$type}++;
      } else {
          $found_types{$type} = 1;
      }
      if ($count > $max_lines) {
          last LOOP;
      }
  } ## End loop
    $reader->close();
    my $found_my_type = 0;
    my $max_type = "";
    my $max = 0;
    foreach my $type (keys %found_types) {

        if ($found_types{$type} > $max) {
            $max_type = $type;
            $max = $found_types{$type};
        }
        if ($found_types{$my_type}) {
            print "The specified type: ${my_type} is in the gff file, comprising $found_types{$my_type} of the first 40,000.\n";
            $found_my_type = 1;
        }
    } ## End the loop

    if ($found_my_type == 1) {
        return($my_type);
    } else {
        print "Did not find your specified type.  Changing it to: ${max_type} which had ${max} entries of the first 40,000.\n";
        return($max_type);
    }
}

=head2
    HTSeq()
=cut
sub HTSeq {
    my $me = shift;
    my %args = @_;
    $me->Check_Options(args => \%args, needed => ["species", "htseq_stranded", "htseq_args",]);
    my $basename = $me->{basename};
    my $species = $me->{species};
    my $stranded = $me->{htseq_stranded};
    my $htseq_type = $args{htseq_type};
    my $htseq_jobname = qq"hts_${species}";
    $htseq_jobname = $args{jobname} if ($args{jobname});
    my $depends = "";
    $depends = $args{depends} if ($args{depends});
    my $libtype = "genome";
    $libtype = $args{libtype} if ($args{libtype});
    my $gff = qq"$me->{libdir}/${libtype}/${species}.gff";
    $gff = $args{htseq_gff} if ($args{htseq_gff});
    my $gtf = $gff;
    $gtf =~ s/\.gff/\.gtf/g;
    my $htseq_args = "";
    if (defined($me->{htseq_args}->{$species})) {
        $htseq_args = $me->{htseq_args}->{$species};
    } else {
        $htseq_args = $me->{htseq_args}->{all};
    }
    my $input = "${basename}.bam";
    $input = $args{input} if ($args{input});
    my $output = $input;
    if ($args{suffix}) {
        my $suffix = $args{suffix};
        $output =~ s/\.bam/_${suffix}\.count/g;
    } else {
        $output =~ s/\.bam/\.count/g;
    }
    $output = $args{output} if ($args{output});
    my $error = $input;
    $error =~ s/\.bam/_htseq\.err/g;
    if (!-r "${gff}" and !-r "${gtf}") {
        die("Unable to read ${gff} nor ${gtf}, please fix this and try again.\n");
    }
    my $annotation = $gtf;
    if (!-r "${gtf}") {
        $annotation = $gff;
    }
    ## Set the '-t FEATURETYPE --type' argument used by htseq-count
    ## This may be provided by a series of defaults in %HT_Multi::gff_types, overridden by an argument
    ## Or finally auto-detected by HT_Types().
    ## This is imperfect to say the nicest thing possible, I need to consider more appropriate ways of handling this.
    if (!defined($htseq_type)) {
        $htseq_type = $me->{htseq_type};
    }
    if (!defined($htseq_type) or $htseq_type eq '' or $htseq_type eq 'auto') {
        $htseq_type = $me->HT_Types(annotation => $annotation, type => $htseq_type);
    }

    my $job_string = qq!htseq-count -q -f bam -s ${stranded} ${htseq_args} -t ${htseq_type} \\
  ${input} ${annotation} \\
  1>${output} 2>${error} && \\
    xz -f -9e ${output}
!;
    my $comment = qq!## Counting the number of hits in ${input} for each feature found in ${annotation}
## Is this stranded? ${stranded}.  The defaults of htseq are:
## $me->{htseq_args}->{default}
!;
    my $htseq = $me->Qsub(job_name => $htseq_jobname,
                          qsub_mem => 6,
                          depends => $depends,
                          job_string => $job_string,
                          comment => $comment,
                          job_prefix => $args{job_prefix},
                          qsub_queue => "throughput",
                          input => $input,
                          output => $output,
                          prescript => $args{prescript},
                          postscript => $args{postscript},
        );
    return($htseq);
}

sub Mi_Map {
    my $me = shift;
    my %args = @_;

    eval "use Bio::DB::Sam; 1;";
    $me->Check_Options(args => \%args, needed => ["mirbase_data", "mature_fasta",
                                                  "mi_genome", "bamfile",]);

    my $bam_base = basename($me->{bamfile}, (".bam", ".sam"));
    my $pre_map = qq"${bam_base}_mirnadb.txt";
    my $final_map = qq"${bam_base}_mature.count";

    ## Step 1:  Read the mature miRNAs to get the global IDs and sequences of the longer transcripts.
    ## This provides the IDs as 'mm_miR_xxx_3p' and short sequences ~21 nt.
    print "Starting to read miRNA sequences.\n";
    my $sequence_db = Read_Mi(seqfile => $me->{mature_fasta});
    ## print Dumper $sequence_db;  ## Works

    ## Step 2:  Collate those IDs against the set of miRNA_transcripts->miRNA_mature IDs using
    ## the tab delimited file downloaded from mirbase.org
    print "Starting to read miRNA mappings.\n";
    my $sequence_mappings = Read_Mappings(mappings => $me->{mirbase_data},
                                          output => $pre_map,
                                          seqdb => $sequence_db);
    ## print Dumper $sequence_mappings;
    ## At this point, we have brought together mature sequence/mature IDs to parent transcript IDs
    ## When we read the provided bam alignment, we will simultaneously read the immature miRNA database
    ## and use them to make a bridge from the mature sequence/IDs to immature miRNA IDs.

    my $final_hits = Read_Bam(mappings => $sequence_mappings,
                              mi_genome => $me->{mi_genome},
                              bamfile => $me->{bamfile});
    Final_Print(data => $final_hits, output => $final_map);

    ## Quick and dirty sequence reader
    sub Read_Mi {
        my %args = @_;
        my $fasta = new Bio::SeqIO(-file => $args{seqfile}, -format => "Fasta");
        my %sequences = ();
        while (my $mi_seq = $fasta->next_seq()) {
            next unless(defined($mi_seq->id));
            my $id = $mi_seq->id;
            my $length = $mi_seq->length;
            my $seq = $mi_seq->seq;
            $sequences{$id}->{sequence} = $seq;
        }
        ## print Dumper %sequences;
        return(\%sequences);
    }

    sub Read_Mappings {
        my %args = @_;
        my $output = $args{output};
        my $seqdb = $args{seqdb};
        my $inmap = new FileHandle;
        $inmap->open($args{mappings});
        my $mimap = {};
        while (my $line = <$inmap>) {
            chomp $line;
            $line =~ s/"//g;
            my ($hit_id, $ensembl, $mirbase_id, $fivep_mature, $fivep_id, $threep_mature, $threep_id) = split(/\s+/, $line);
            if (defined($fivep_id) and $fivep_id ne "") {
                ## print "TESTME: $mirbase_id : $fivep_id\n";
                $mimap->{$fivep_id}->{mirbase_id} = $mirbase_id;
                $mimap->{$fivep_id}->{hit_id} = $hit_id;
                $mimap->{$fivep_id}->{ensembl} = $ensembl;
                $mimap->{$fivep_id}->{mimat} = $fivep_id;
            }
            if (defined($threep_id) and $threep_id ne "") {
                ## print "TESTME: $mirbase_id : $threep_id\n";
                $mimap->{$threep_id}->{mirbase_id} = $mirbase_id;
                $mimap->{$threep_id}->{hit_id} = $hit_id;
                $mimap->{$threep_id}->{ensembl} = $ensembl;
                $mimap->{$threep_id}->{mimat} = $threep_id;
            }
        }
        $inmap->close();

        my $out = new FileHandle;
        $out->open(">${output}");
        ## Now re-key the database of miRNAs so that there are (potentially) multiple elements to each mirbase ID
        ## This is a bit redundant, but hopefully clearer.  The idea is to first gather all miR data, then make a list of specific
        ## ID/sequences for each mature miRNA child of the mirbase_id parent entry.  Eg:
        ## mirbase_id parent ->  [ mature_id1, mature_id2, ... ]
        ## where each mature_id is a { sequence, count, id, ensembl, ... }
        my $newdb = {};
      LOOP: foreach my $id (keys %{$seqdb}) {
          next LOOP unless ($id =~ /^mmu/);
          my @hit_list = ();
          if (defined($mimap->{$id})) {
              ##print "Found $id across the mappings.\n";
              my $new_id = $id;
              my $sequence = $seqdb->{$id}->{sequence};
              if ($mimap->{$id}->{hit_id}) {
                  $new_id = $mimap->{$id}->{hit_id};
              }
              my $element = {
                  sequence => $sequence,
                  count => 0,
                  mirbase_id => $mimap->{$id}->{mirbase_id},
                  id => $id,
                  ensembl => $mimap->{$id}->{ensembl},
                  hit_id => $mimap->{$id}->{hit_id},
                  mimat => $mimap->{$id}->{mimat},
              };
              my $ensembl_id = $element->{ensembl};
              if (defined($newdb->{$ensembl_id})) {
                  ## Then there should be an existing element, append to it.
                  @hit_list = @{$newdb->{$ensembl_id}};
                  push(@hit_list, $element);
                  ## $newdb->{$new_id} = \@hit_list;
                  $newdb->{$ensembl_id} = \@hit_list;
                  print $out "$id; $element->{mirbase_id}; $element->{ensembl}; $element->{hit_id}; $element->{mimat}; $element->{sequence}\n";

              } else {
                  @hit_list = ();
                  push(@hit_list, $element);
                  ## $newdb->{$new_id} = \@hit_list;
                  $newdb->{$ensembl_id} = \@hit_list;
              }

          } else {
              ## print "Did not find $id\n";
          }
      }
        $out->close();
        return($newdb);
    }

    sub Read_Bam {
        my %args = @_;
        my $mappings = $args{mappings};
        my $sam = new Bio::DB::Sam(-bam => $args{bamfile}, -fasta=> $args{mi_genome},);
        my $bam = $sam->bam;
        my $header = $bam->header;
        my $target_count = $header->n_targets;
        my $target_names = $header->target_name;
        my $align_count = 0;
        print "Started reading bamfile.\n";
        $options{debug} = 0;
        ## I probably don't need to acquire most of this information, as I am only really taking
        ## the read_seq and read_seqid
      BAM: while (my $align = $bam->read1) {
          if ($options{debug}) {
              last BAM if ($align_count > 4000);
          }
          my $read_seqid = $target_names->[$align->tid];
          ## my $read_start = $align->pos + 1;
          ## my $read_end = $align->calend;
          ## my $read_strand = $align->strand;
          ## my $read_cigar = $align->cigar_str;
          ## my @read_scores = $align->qscore;        # per-base quality scores
          ## my $read_match_qual= $align->qual;       # quality of the match
          ## my $read_length = $align->query->end;
          my $read_seq = $align->query->seq->seq; ## maybe?
          $align_count++;

          ## Check each element in the mappings to see if they are contained in the read's sequence and vice versa.
          ## If so, increment the count for that element and move on.
          $read_seqid =~ s/^chr[\d+]_//g;
          ## print "TESTME MAPPINGS: <$read_seqid> <$mappings->{$read_seqid}>\n";

          if ($mappings->{$read_seqid}) {
              ##print Dumper $mappings->{$read_seqid};
              my @element_list = @{$mappings->{$read_seqid}};
              my $found_element = 0;
              my $element_length = scalar(@element_list);
              ## print "Found map, checking against ${element_length} mature RNAs.\n";

              foreach my $c (0 .. $#element_list) {
                  my $element_datum = $element_list[$c];
                  my $element_seq = $element_list[$c]->{sequence};
                  ## print "Comparing: $read_seq vs. $element_seq\n";

                  my @read_vs_element = amatch($read_seq, [ 'S1' ], ($element_seq));
                  my @element_vs_read = amatch($element_seq, [ 'S1' ], ($read_seq));
                  if (scalar(@read_vs_element) > 0 or scalar(@element_vs_read) > 0) {
                      $mappings->{$read_seqid}->[$c]->{count}++;
                      ##print "Found hit with amatch against $element_seq: ";
                      ##use Data::Dumper;
                      ##print Dumper $element_datum;
                      ##print "We need to get the slot in position: {ensembl_id} -> [c] -> {id}\n";
                      ##print "is read_seqid ensembl? $read_seqid  \n";
                      ##print "Incremented $mappings->{$read_seqid}->[$c]->{id} to $mappings->{$read_seqid}->[$c]->{count}\n";
                      ##print Dumper $mappings;
                      $found_element = 1;
                  }
                  ## if ($read_seq =~ /$element_seq/ or $element_seq =~ /$read_seq/) {
                  ##     $mappings->{$read_seqid}->[$c]->{count}++;
                  ##     print "Using an exact match: ";
                  ##     print "Incremented $mappings->{$read_seqid}->[$c]->{mirbase_id} to $mappings->{$read_seqid}->[$c]->{count}\n";
                  ##     $found_element = 1;
                  ## }
              }
              ##if ($found_element == 0) {
              ##    print "No map for $read_seq\n";
              ##}
          }
      } ## Finish iterating through every read
        return($mappings);
    }

    sub Final_Print {
        my %args = @_;
        my $final = $args{data};
        my $output = new FileHandle;
        $output->open(">$args{output}");
        foreach my $immature (keys %{$final}) {
            foreach my $mature (@{$final->{$immature}}) {
                print $output "$mature->{mimat} $mature->{count}\n";
            }
        }
        $output->close();
    } ## End of Final_Print

} ## End of Mi_Map

=back

=head1 AUTHOR - atb

Email <abelew@gmail.com>

=head1 SEE ALSO

    L<htseq-count>

=cut

1;
