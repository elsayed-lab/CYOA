package HPGL;
use common::sense;
use autodie qw":all";

use Bio::DB::Sam;
use Bio::Tools::GFF;
use PerlIO;
use Data::Dumper;
use Math::BigFloat ':constant';
use Bio::SeqIO;
use POSIX;

sub Calibrate {

    print "This script attempts to calibrate a set of aligned reads using
my observations about the relationship between the observed position
of the read and the most likely position of the ribosome on an mRNA.
These observations are not likely to be the same for other organisms
nor even different states of the same organism.
The function Calibrate_Reads() below may therefore be replaced with
new observations.\n";

    my $all_state_prob =  {
        total => 0,
        last => 0,
        '0' => 0,
    };
    my $all_transition_prob = {
        total => 0,
        '00' => 0,
    };

    print "Starting to read the genome.\n";
    my $chromosome_number = 0;
    my $chromosomes = {};
    Read_Fasta();

    print "Finished Reading the genome, read $chromosome_number chromosomes.\n";
    Read_Bam();

    sub Read_Fasta {
        my $input_genome = new Bio::SeqIO(-file => $options{fasta}, -format => 'Fasta');
        while (my $genome_seq = $input_genome->next_seq()) {
            next unless(defined($genome_seq->id));
            $chromosome_number++;
            my $id = $genome_seq->id;
            print "Reading chromosome: $id\n";
            my $length = $genome_seq->length;
            my $empty_forward = [];
            my $empty_reverse = [];
            for my $c (0 .. $length) {
                $empty_forward->[$c] = 0;
                $empty_reverse->[$c] = 0;
            }
            $chromosomes->{$id}->{forward} = $empty_forward;
            $chromosomes->{$id}->{reverse} = $empty_reverse;
        }
    }

    sub Read_Bam {
        my %args = @_;
        system("mkdir -p $options{outdir}") if (!-d $options{outdir});
        my $sam = new Bio::DB::Sam(-bam => $options{input}, -fasta=> $options{fasta},);
        my $bam = $sam->bam;
        my $header = $bam->header;
        my $target_count = $header->n_targets;
        my $target_names = $header->target_name;
        my $align_count = 0;
        print "Started reading bamfile.\n";
      BAM: while (my $align = $bam->read1) {
          if ($options{debug}) {
              last BAM if ($align_count > 200000);
          }
          my $read_seqid = $target_names->[$align->tid];
          my $read_start = $align->pos + 1;
          my $read_end = $align->calend;
          my $read_strand = $align->strand;
          my $read_cigar = $align->cigar_str;
          my @read_scores = $align->qscore;        # per-base quality scores
          my $read_match_qual= $align->qual;       # quality of the match
          my $read_length = $align->query->end;
          next BAM unless(Good_Size($read_length));
          $align_count++;
          my $new_position;
          if ($read_strand == 1) {
              $new_position = Calibrate(start => $read_start, strand => 1, length => $read_length);
              $chromosomes->{$read_seqid}->{forward}->[$new_position]++ if ($options{asite});
              $chromosomes->{$read_seqid}->{forward}->[$new_position + 3]++ if ($options{psite});
              $chromosomes->{$read_seqid}->{forward}->[$new_position + 6]++ if ($options{esite});
          } else {
              $new_position = Calibrate(start => $read_end, strand => -1, length => $read_length);
              $chromosomes->{$read_seqid}->{reverse}->[$new_position]++ if ($options{asite});
              $chromosomes->{$read_seqid}->{reverse}->[$new_position + 3]++ if ($options{psite});
              $chromosomes->{$read_seqid}->{reverse}->[$new_position + 6]++ if ($options{esite});
          }

          if (($align_count % 10000) == 0) {
              my $a = $align_count / 1000;
              print "${a}k chr: $read_seqid $read_strand $new_position\n";
          }
      } ## End reading the bam file

        my $gff_entries;
        if ($options{gff}) {
            $gff_entries = Read_GFF();
        }

        print "Finished reading bam file.\n";
        for my $chr (1 .. 36) {
            $chr = sprintf("%02d", $chr);
            $chr = "LmjF.${chr}";
            print "Counting and printing hits on: $chr\n";

            my $fwd = $chromosomes->{$chr}->{forward};
            my $rev = $chromosomes->{$chr}->{reverse};

            my $forward_name = qq"$options{outdir}/f_$chr.tab";
            my $reverse_name = qq"$options{outdir}/r_$chr.tab";
            Print_Chromosome($forward_name, $fwd);
            Print_Chromosome($reverse_name, $rev);

            my ($fwd_states, $fwd_state_prob, $fwd_transition_prob) = Count_States($fwd);
            my ($rev_states, $rev_state_prob, $rev_transition_prob) = Count_States($rev);

            my $fwd_state_file = qq"$options{outdir}/f_${chr}.states";
            my $rev_state_file = qq"$options{outdir}/r_${chr}.states";
            Print_States($fwd_state_file, $fwd_states);
            Print_States($rev_state_file, $rev_states);
            my $fwd_stateprob_file = qq"$options{outdir}/f_${chr}_state.prob";
            my $rev_stateprob_file = qq"$options{outdir}/r_${chr}_state.prob";
            Print_State_Prob($fwd_stateprob_file, $fwd_state_prob);
            Print_State_Prob($rev_stateprob_file, $rev_state_prob);
            my $fwd_transprob_file = qq"$options{outdir}/f_${chr}_trans.prob";
            my $rev_transprob_file = qq"$options{outdir}/r_${chr}_trans.prob";
            Print_Trans_Prob($fwd_transprob_file, $fwd_transition_prob);
            Print_Trans_Prob($rev_transprob_file, $rev_transition_prob);

            my $all_stateprob_file = qq"$options{outdir}/all_state.prob";
            my $all_transprob_file = qq"$options{outdir}/all_trans.prob";
            Print_State_Prob($all_stateprob_file, $all_state_prob);
            Print_Trans_Prob($all_transprob_file, $all_transition_prob);

            if ($options{gff}) {
                Count_Tables($fwd, $gff_entries->{$chr}, 1);
                Count_Tables($rev, $gff_entries->{$chr}, -1);
            }

        }  ## End printing every chromosome
    }

    sub Count_Tables {
        my $data = shift;
        my $entries = shift;
        my $strand = shift;
        my $multiplier = $options{asite} + $options{psite} + $options{esite};
        my $count = new FileHandle;
        $count->open(">>$options{outdir}/options{gff_type}.count");
      CANDIDATES: foreach my $candidate (keys %{$entries}) {
          my $orf_start = $entries->{$candidate}->{start} - 1;
          my $orf_end = $entries->{$candidate}->{end} - 1;

          my $candidate_count = 0;
        LOOP: for my $c ($orf_start .. $orf_end) {
            next LOOP unless ($data->[$c] > 0);
            $candidate_count = $candidate_count + $data->[$c];
        }
          $candidate_count = floor($candidate_count / $multiplier);
          my $string = qq"$entries->{$candidate}->{id}\t$candidate_count\n";
          print $count $string;
          print "TESTME: $string";
      }
        $count->close();
    }

    sub Count_States {
        my $datum = shift;
        my @states = ();
        my $state_prob =  {
            total => 0,
            last => 0,
            '0' => 0,
        };
        my $transition_prob = {
            total => 0,
            '00' => 0,
        };

        my @data = @{$datum};
      DATA: while (@data) {
          last DATA if (scalar(@data) < 3);
          my $first = shift @data;
          my $second = shift @data;
          my $third = shift @data;
          my $current = 0;

          $all_state_prob->{total}++;
          $state_prob->{total}++;
          $all_transition_prob->{total}++;
          $transition_prob->{total}++;

          if ($first == 0 and $second == 0 and $third == 0) {
              $current = 0;
          } elsif ($first > $second and $first > $third) {
              $current = 1;
          } elsif ($second > $first and $second > $third) {
              $current = 2;
          } elsif ($third > $first and $third > $second) {
              $current = 3;
          } elsif ($first == $second and $first > $third) {
              $current = 4;
          } elsif ($second == $third and $second > $first) {
              $current = 5;
          } elsif ($first == $third and $first > $second) {
              $current = 6;
          } else {
              $current = 7;
          }
          push(@states, $current);

          ## Now increment the relevant total state information
          if ($all_state_prob->{$current}) {
              $all_state_prob->{$current}++;
          } else {
              $all_state_prob->{$current} = 1;
          }
          my $transition = "$all_state_prob->{last}$current";
          if ($all_transition_prob->{$transition}) {
              $all_transition_prob->{$transition}++;
          } else {
              $all_transition_prob->{$transition} = 1;
          }

          ## Now do the likelihoods for the current chromosome
          if ($state_prob->{$current}) {
              $state_prob->{$current}++;
          } else {
              $state_prob->{$current} = 1;
          }
          if ($transition_prob->{$transition}) {
              $transition_prob->{$transition}++;
          } else {
              $transition_prob->{$transition} = 1;
          }

          ## End the calculation and flip last->current
          $all_state_prob->{last} = $current;
      }
        return(\@states, $state_prob, $transition_prob);
    }

    sub Good_Size {
        my $size = shift;
        my $size_string = $options{sizes};
        my @allowed_sizes = split(/\,/, $size_string);
        foreach my $allowed (@allowed_sizes) {
            return(1) if ($size == $allowed);
        }
        return(0);
    }

    sub Calibrate {
        my %args = @_;
        my $rel_pos = $args{start};
        my $strand = $args{strand};
        my $length = $args{length};

        ## Beginning of reads to the start codon
        if ($options{orientation} eq 'start' and $strand eq '1') {
            if ($length == 24) {
                $rel_pos = $rel_pos + 9;
            } elsif ($length == 25) {
                $rel_pos = $rel_pos + 10;
            } elsif ($length == 26) {
                $rel_pos = $rel_pos + 10;
            } elsif ($length == 27) {
                $rel_pos = $rel_pos + 12;
            } else {
                $rel_pos = $rel_pos + 13;
            }
        } else {
            if ($length == 25) {
                $rel_pos = $rel_pos - 10;
            } elsif ($length == 27) {
                $rel_pos = $rel_pos - 12;
            } elsif ($length == 29) {
                $rel_pos = $rel_pos - 13;
            } elsif ($length == 30) {
                $rel_pos = $rel_pos - 13;
            } elsif ($length == 31) {
                $rel_pos = $rel_pos - 13;
            } elsif ($length == 32) {
                $rel_pos = $rel_pos - 13;
            } elsif ($length == 33) {
                $rel_pos = $rel_pos - 13;
            } elsif ($length == 34) {
                $rel_pos = $rel_pos - 13;
            } else {
                $rel_pos = $rel_pos - 13;
            }
        }
        return($rel_pos);
    }

    sub Print_Chromosome {
        my $name = shift;
        my $datum = shift;
        my $out = new FileHandle;
        $out->open(">${name}");
        my @data = @{$datum};
        foreach my $base (@data) {
            print $out "$base\n";
        }
        $out->close();
    }

    sub Print_States {
        my $file = shift;
        my $datum = shift;
        my $out = new FileHandle;
        $out->open(">${file}");
        foreach my $s (@{$datum}) {
            print $out $s;
        }
        print $out "\n";
        $out->close();
    }

    sub Print_State_Prob {
        my $file = shift;
        my $states = shift;
        my $out = new FileHandle;
        $out->open(">${file}");
        foreach my $k (sort keys %{$states}) {
            my $likelihood;
            if ($states->{total} == 0) {
                $likelihood = 0;
            } else {
                $likelihood = $states->{$k} / $states->{total};
                my $n = new Math::BigFloat($likelihood);
                $likelihood = $n->bround(8);
            }
            print $out "$k $states->{$k} $likelihood\n";
        }
        $out->close();
    }

    sub Print_Trans_Prob {
        my $file = shift;
        my $trans = shift;
        my $out = new FileHandle;
        $out->open(">${file}");
        foreach my $k (sort keys %{$trans}) {
            my $likelihood;
            if ($trans->{total} == 0) {
                $likelihood = 0;
            } else {
                $likelihood = $trans->{$k} / $trans->{total};
                my $n = new Math::BigFloat($likelihood);
                $likelihood = $n->bround(8);
            }
            print $out "$k $trans->{$k} $likelihood\n";
        }
        $out->close();
    }

    sub Read_GFF {
        open(GFF, "<$options{gff}");
        use Bio::Tools::GFF;
        my $annotation_in = new Bio::Tools::GFF(-fh => \*GFF, -gff_version => 3);
        my $gff_out = {};
        print "Starting to read gff: $options{gff}\n";
      LOOP: while(my $feature = $annotation_in->next_feature()) {
          next LOOP unless ($feature->{_primary_tag} eq $options{gff_type});
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
}



sub Graph_Reads {

    if ($options{help}) {
        print "This script attempts to simplify the graphing of read distributions for ribosome profiling experiments.
Likely options include:
--output|-o: directory to put files into
--input|-i: input text file containing lengths and read positions
--minsize: minimum size to consider
--maxsize: maximum size to consider
--maxpos: maximum distance from the position of interest to the right
--minpos: maximum distance from the position of interest to the left
--correction: include some simple attempted corrections from previous observations
--orientation: count to the stop codon or start codon
--anchor: count from the beginning or end of reads
--help: print this.
";
        exit(0);
    }

    if (!-d $options{output}) {  system("mkdir $options{output}") };
    #my $counted = Read_Counts(min => 24, max => 36);
    my $counted = Read_Counts();

    sub Read_Counts {
        my %args = @_;
        my $minsize = $options{minsize};
        my $maxsize = $options{maxsize};
        my $minpos = $options{minpos};
        my $maxpos = $options{maxpos};

        foreach my $size ($minsize .. $maxsize) {
            #	open($size, ">${options{output}}/${size}.tab");
            open($size, ">:gzip", "${options{output}}/${size}.tab.gz");
        }
        # open(TOTAL, ">$options{output}/total.tab");
        open(TOTAL, ">:gzip", "${options{output}}/total.tab.gz");
        open(IN, "zcat $options{input} |");
      LOOP: while (my $line = <IN>) {
          chomp $line;
          ##4902 LmjF.01.0010 -1 ORF: 3704 4702 REL: 3717 3705 985 13
          my ($something, $abs_pos, $id, $strand, $orf, $st, $en, $rel, $rel_start, $rel_end, $rel_pos, $length)  = split(/\s+/, $line);
          next LOOP if (!defined($id))
              next LOOP if ($id =~ /RNA/);
          next unless defined($length);
          if ($length >= $minsize and $length <= $maxsize and $rel_pos <= $maxpos and $rel_pos >= $minpos) {
              if ($options{correction}) {
                  $rel_pos = Correct(rel => $rel_pos, length => $length);
              }
              select($length);
              print TOTAL "$id $rel_pos\n";
              print "$id $rel_pos\n";
          } else {
              next LOOP;
          }
      } ## End foreach line of the input file
        foreach my $size ($minsize .. $maxsize) {
            close($size);
        }
        close(IN);
        close(TOTAL);
    }

    sub Correct {
        my %args = @_;
        my $rel_pos = $args{rel};
        my $length = $args{length};

        ## Beginning of reads to the start codon
        if ($options{orientation} eq 'start' and $options{anchor} eq 'start') {
            if ($length == 24) {
                $rel_pos = $rel_pos + 9;
            } elsif ($length == 25) {
                $rel_pos = $rel_pos + 10;
            } elsif ($length == 26) {
                $rel_pos = $rel_pos + 10;
            } elsif ($length == 27) {
                $rel_pos = $rel_pos + 12;
            } else {
                $rel_pos = $rel_pos + 13;
            }
        }

        elsif ($options{orientation} eq 'start' and $options{anchor} eq 'end') {
            if ($length == 24) {
                $rel_pos = $rel_pos - 11;
            } elsif ($length == 25) {
                $rel_pos = $rel_pos - 11;
            } elsif ($length == 26) {
                $rel_pos = $rel_pos - 12;
            } elsif ($length == 27) {
                $rel_pos = $rel_pos - 14;
            } elsif ($length == 28) {
                $rel_pos = $rel_pos - 14;
            } elsif ($length == 29) {
                $rel_pos = $rel_pos - 15;
            } elsif ($length == 30) {
                $rel_pos = $rel_pos - 16;
            } elsif ($length == 31) {
                $rel_pos = $rel_pos - 17;
            } elsif ($length == 32) {
                $rel_pos = $rel_pos - 18;
            } elsif ($length == 33) {
                $rel_pos = $rel_pos - 19;
            } elsif ($length == 34) {
                $rel_pos = $rel_pos - 20;
            } elsif ($length == 35) {
                $rel_pos = $rel_pos - 21;
            } else {
                $rel_pos = $rel_pos - ($length - 14);
            }
        }

        ## Beginning of reads to the stop codon
        elsif ($options{orientation} eq 'end' and $options{anchor} eq 'start') {
            if ($length == 24) {
                $rel_pos = $rel_pos + 14;
            } elsif ($length == 25) {
                $rel_pos = $rel_pos + 15;
            } elsif ($length == 26) {
                $rel_pos = $rel_pos + 15;
            } elsif ($length == 27) {
                $rel_pos = $rel_pos + 17;
            } else {
                $rel_pos = $rel_pos + 18;
            }
        }

        elsif ($options{orientation} eq 'end' and $options{anchor} eq 'end') {
            if ($length == 24) {
                $rel_pos = $rel_pos - 6;
            } elsif ($length == 25) {
                $rel_pos = $rel_pos - 6;
            } elsif ($length == 26) {
                $rel_pos = $rel_pos - 11;
            } elsif ($length == 27) {
                $rel_pos = $rel_pos - 12;
            } elsif ($length == 28) {
                $rel_pos = $rel_pos - 12;
            } elsif ($length == 29) {
                $rel_pos = $rel_pos - 11;
            } elsif ($length == 30) {
                $rel_pos = $rel_pos - 11;
            } elsif ($length == 31) {
                $rel_pos = $rel_pos - 12;
            } elsif ($length == 32) {
                $rel_pos = $rel_pos - 10;
            } elsif ($length == 33) {
                $rel_pos = $rel_pos - 11;
            } elsif ($length == 34) {
                $rel_pos = $rel_pos - 11;
            } else {
                $rel_pos = $rel_pos - 10;
            }
        }
        return($rel_pos);
    }

}

1;

__END__

=head1 NAME

count_reads.pl - Describe the usage of script briefly

=head1 SYNOPSIS

count_reads.pl [options] args

      -opt --long      Option description

=head1 DESCRIPTION

Stub documentation for count_reads.pl,

=head1 AUTHOR

Ashton Trey Belew, E<lt>abelew@gmail.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014 by Ashton Trey Belew

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut

## EOF
