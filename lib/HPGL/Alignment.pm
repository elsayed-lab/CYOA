package HPGL;
use common::sense;
use autodie qw":all";

use Cwd;
use Bio::Seq;
use Bio::SearchIO::blast;
use Bio::SearchIO::fasta;
use Bio::Tools::Run::StandAloneBlast;
use POSIX qw"ceil";

=head1 NAME
    HPGL::Alignment - Perform Sequence alignments with tools like Blast/Fasta

=head1 SYNOPSIS

    use HPGL;
    my $hpgl = new HPGL;
    $hpgl->Split_Align(query => 'amino_query.fasta', library => 'amino_library.fasta');

=head1 SYNOPSISold

  The following is the original pod for this:

  split_align_blast.pl -i some_multi_fasta_query.fasta -l some_multi_fasta_library.fasta
  split_align_blast --help

 Other options:  The following are some ways to pass options to this script and their results:
  --best_only : tells this to print a summary table of only the best hits
  --blast|-b : which blast program is desired? (blastn|blastp|tblastn|tblastx)
  --blast_params|-bp ' -e 100 ' : Sets the minimum e-score to 100 -- set any arbitrary blast options here
  --clean : Clean up the mess this script makes
  --formatdb|-f : write out formatdb options (really shouldn't be needed)
                  ('formatdb -p \$peptide -o T -n blastdb/\$lib -s -i ' by default)
  --help|h : print out the help
  --input|-i bob.fasta  : Sets the query library to bob.fasta
  --lib|-l clbrener.fasta : Sets the database library to clbrener.fasta
  --number|-n 10 : Sets the number of jobs to run to 10 (200 by default)
  --output|-o : A file to which to write the concatenated last outputs (split_align_output.txt.gz)
  --parse|-p : Add this argument to parse the blast output (1 by default, set to 0 to skip parsing)
  --peptide|-p : T for a peptide database, F for a nucleotide for formatdb -- this will be overridden depending on the blast prorgam used
  --queue|q : what pbs queue (workstation)
  --skip : Skip the alignments and just try parsing the outputs

  A couple of example invocations which I used in my own work:
  split_align_blast.pl -i esmer_annotated_cds.fasta -l TriTrypDB-9.0_LmajorFriedlin_Genome.fasta -q workstation
    That split the esmer_annotated_cds file into ~200 pieces and aligned them against the leishmania major genome on the workstation queue.

  split_align_blast.pl -i esmer_annotated_cds.fasta -l TriTrypDB-9.0_LmajorFriedlin_Genome.fasta -b tblastx --parse 0
    Same deal, but use tblastx and don't parse the output

=head1 DESCRIPTIONold

    This script seeks to make running large blast jobs easier.  It will split the input query file
    into many smaller files, format a blastdb appropriate for the alignment, submit them to the cluster,
    collect the alignments, parse the outputs, and cleanup the mess.

=head2 Methods

=over 4

=item C<Blast_Parse>

    $hpgl->Blast_Parse(query => 'nt_query.fasta', library => 'nt_library.fasta');
    will create a blast library from 'library' and perform individual
    blast searches for every sequence in nt_query.fasta against it.
    Finally, it writes a summary of the results in a set of tables.

=cut
sub Parse_Blast {
    my $me = shift;
    my %args = @_;
    print STDERR qq"Please note that this function calls blastall
    separately for every sequence in the database.  As a result it is
    not recommended fo use with large sequence libraries.  Instead use
    the separate functions 'Run_Blast()' or 'Split_Align_Blast()'
    followed by 'Parse_Blast()' which does these steps separately.";
    my $blast_program = 'blastp';
## Search libraries are the various combinations of esmer, nonesmer, and unassigned which are of potential interest
## thus all 7 of the possible states
## These libraries are generated in make_blast.sh
    my @search_libraries = ('nr',);
    #my @search_libraries = ('esmer',
    #                       'non',
    #                       'unas',
    #                       'esmer-non',
    #                       'esmer-unas',
    #                       'non-unas',
    #                       'all',);
    ## States to save: singles, doubles, few (3-10), many (10+)
    my $blast_output = new Bio::SearchIO(-format => 'blast',);
    my $query = $args{query};
    for my $library (@search_libraries) {
        my $query_library = new Bio::SeqIO(-file => ${query}, -format => 'Fasta');
        my $output_directory = qq"$me->{basedir}/outputs/${query}_vs_${library}";
        make_path("${output_directory}") unless (-d $output_directory);
        open(COUNTS, ">${output_directory}/counts.txt");
        open(ZEROS, ">${output_directory}/zeros.fasta");
        open(SINGLES, ">${output_directory}/singles.fasta");
        open(DOUBLES, ">${output_directory}/doubles.fasta");
        open(FEW, ">${output_directory}/few.fasta");
        open(MANY, ">${output_directory}/many.fasta");
        my $seq_count = 0;

        my @params = (-program => $blast_program,
                      ## Show GI's in definition lines?
                      -F => 'T',  ## Filter query sequence
                      -G => '-1', ## Cost to open a gap
                      -E => '-1', ## Cost to extend a gap
                      -X => '0',  ## X dropoff value for gapped alignment
                      -I => 'T',  ## Default in blast is F
##                      -q => '-3', ## Penalty for nucleotide mismatch (blastn only)
                      ## And many more
                      -database => $library,
                      -outfile => qq"blast_out",
            );

        while (my $query_seq = $query_library->next_seq()) {
            ## I think this while loop may not be needed.
            $seq_count++;
            my $id = $query_seq->id;
            my $desc = $query_seq->desc;
            my $seq = $query_seq->seq;
            my $search = new Bio::Tools::Run::StandAloneBlast(@params);
            print Dumper $search;
            my $blast_output;
            eval {
                ##    $blast_output = $search->blastall($query_library);
                $blast_output = $search->blastall($query_seq);
            };
            if ($@) {
                print "Error? $@\n";
            }
            my $result_count = 0;
            while (my $result = $blast_output->next_result()) {
                $result_count++;
                my $query_name = $result->query_name();
                my $query_length = $result->query_length();
                my $query_descr = $result->query_description();
                my $stats = $result->available_statistics();
                my $hits = $result->num_hits();
                my $hit_count = 0;
                my $score_cutoff = 100;
                my $sig_cutoff = 0.00001;
              HITLOOP: while (my $hits = $result->next_hit()) {
                  my $hit_name = $hits->name();
                  my $hit_length = $hits->length();
                  my $hit_acc = $hits->accession();
                  my $hit_descr = $hits->description();
                  my $hit_score = $hits->score();
                  my $hit_sig = $hits->significance();
                  my $hit_bits = $hits->bits();
                  next HITLOOP unless ($hit_sig < $sig_cutoff and $hit_score > $score_cutoff);
                  $hit_count++;
              }
                my $entry = "$query_name $query_descr\n";
                if ($hit_count == 1) {
                    print SINGLES $entry;
                } elsif ($hit_count == 2) {
                    print DOUBLES $entry;
                } elsif ($hit_count >= 3 and $hit_count <= 10) {
                    print FEW $entry;
                } elsif ($hit_count > 10) {
                    print MANY $entry;
                } else {
                    print ZEROS $entry;
                }
                print COUNTS "$query_name\t$hit_count\n";
            } ## End of the individual blast search  (maybe not needed?)
        } ## End of this search library  (the 7 states to search against)
        close(COUNTS);
        close(ZEROS);
        close(SINGLES);
        close(DOUBLES);
        close(FEW);
        close(MANY);
    }
}

=item C<Parse_Fasta>

    hmmm

=cut
sub Parse_Fasta {
    my $me = shift;
    my %args = @_;
    my $best = 1;
    $best = $args{best} if (defined($args{best}));
    my $sig = 0.0001;
    $sig = $args{sig} if (defiend($args{sig}));
    my $output = qq"$me->{basedir}/outputs/fasta/parsed.txt";
    if (!-d "$me->{basedir}/outputs/fasta") {
        make_path("$me->{basedir}/outputs/fasta");
    }
    my $out = new FileHandle;
    $out->open(">$output");
    ## This works, other attempts I've used for zipped input fail.
    my $f = new FileHandle;
    $f->open("lesspipe $args{input} |");
    ## This works, other attempts at using zipped input fail.
    my $searchio = new Bio::SearchIO(-format => 'fasta', -fh => $f, -best => ${best}, -signif => ${sig});
    print $out "Query Name\tQuery length\tHit ID\tHit Length\tScore\tE\tIdentity\tHit length\tHit Matches\n";
    my $results = 0;
    while (my $result = $searchio->next_result()) {
        $results++;
        while(my $hit = $result->next_hit) {
            my $query_name = $result->query_name();
            my $query_length = $result->query_length();
            my $accession = $hit->accession();
            my $acc2 = $hit->name();
            my $length = $hit->length();
            my $score = $hit->raw_score();
            my $sig = $hit->significance();
            my $ident = $hit->frac_identical();

            my $hit_len;
            my $hit_matches;
            while (my $hsp = $hit->next_hsp) {
                $hit_len = $hsp->length('total');
                my @matches = $hsp->matches(-seq => 'hit');
                $hit_matches = $matches[1];
            }

            print $out "$query_name\t$query_length\t$acc2\t$length\t$score\t$sig\t$ident\t$hit_len\t$hit_matches\n";
            # process the Bio::Search::Hit::HitI object
            #       while(my $hsp = $hit->next_hsp) {
            # process the Bio::Search::HSP::HSPI object
            #       }
        } ## End of each hit of a result.
    } ## End of each result.
    $f->close();
    $out->close();
    return($results);
}

=item C<Split_Align>

    $hpgl->Split_Align();
    Split apart a set of query sequences into $args{number} pieces and align them all separately.

=cut
sub Split_Align_Blast {
    my $me = shift;
    my %args = @_;
    my $query = $args{query};
    my $library = $args{library};
    $args{param} = ' -e 10 ' unless (defined($args{param}));
    $args{tool} = 'blastp' unless (defined($args{tool}));
    $args{number} = 200 unless (defined($args{number}));
    $args{parse} = 0 unless (defined($args{parse}));
    $args{output_type} = 0 unless (defined($args{output_type}));
    $args{num_dirs} = 0 unless (defined($args{num_dirs}));
    $args{best_only} = 0 unless (defined($args{best_only}));
    $args{queue} = 'workstation' unless (defined($args{queue}));
    $args{num_sequences} = 0 unless (defined($args{num_sequences}));
    $args{library} = basename($args{library}, ('.fasta'));
    my $q = basename($args{query}, ('.fasta'));

    ## Number entries 23119
    ## Thus 116 entries in each fasta

    ##  To arrive at this number, I just rounded up (/ 23119 200)
    ## However, we can count the number of fasta entries in the input file
    ##    my $num_per_split = $me->Get_Split(input => $query, number => $number);
    my $num_per_split = $me->Get_Split(%args);
    $args{num_per_split} = $num_per_split;
    print "Going to make $args{number} directories with $num_per_split sequences each.\n";
    my $actual = $me->Make_Directories(%args);
    print "Actually used ${actual} directories to write files.\n";
    $me->Check_Blastdb(%args);
    my $alignment = $me->Make_Align(number => $actual);
    my $out_dir = qq"$me->{basedir}/outputs/blast_out";
    make_path($out_dir) unless (-d "${out_dir}");
    my $out = qq"${out_dir}/blast_out/${q}_vs_$args{library}.txt.gz";
    $me->Concatenate_Searches(depends => $alignment->{pbs_id}, output => $out);
    return($num_per_split);
}


sub Split_Align_Fasta {
    my $me = shift;
    my %args = @_;
    $me->Check_Options(["query","library"]);
    $args{query} = $me->{query} if ($me->{query});
    $args{library} = $me->{library} if ($me->{library});

    ## $args{param} = ' -e 10 ' unless (defined($args{param}));
    ## $args{tool} = 'ggsearch36' unless (defined($args{tool}));
    $args{number} = 200 unless (defined($args{number}));
    $args{parse} = 0 unless (defined($args{parse}));
    $args{output_type} = 0 unless (defined($args{output_type}));
    $args{num_dirs} = 0 unless (defined($args{num_dirs}));
    $args{best_only} = 0 unless (defined($args{best_only}));
    $args{queue} = 'workstation' unless (defined($args{queue}));
    $args{num_sequences} = 0 unless (defined($args{num_sequences}));
    ## $args{library} = basename($args{library}, ('.fasta'));
    my $lib = basename($args{library}, ('.fasta'));
    my $que = basename($args{query}, ('.fasta'));
    my $outdir = qq"$me->{baseidr}/outputs/fasta";
    make_path("${outdir}") unless(-d ${outdir});
    my $output = qq"${outdir}/${que}_vs_${lib}.txt.gz";
    my $concat_job;
    if ($me->{pbs}) {
        my $num_per_split = $me->Get_Split(%args);
        $args{num_per_split} = $num_per_split;
        print "Going to make $args{number} directories with $num_per_split sequences each.\n";
        my $actual = $me->Make_Directories(%args);
        print "Actually used ${actual} directories to write files.\n";
        my $alignment = $me->Make_Fasta_Align(number => $actual);
        $concat_job = $me->Concatenate_Searches(depends => $alignment->{pbs_id}, output => ${output},);
    } else {
        ## If we don't have pbs, force the number of jobs to 1.
        $args{number} = 1;
        my $num_per_split = $me->Get_Split(%args);
        my $actual = $me->Make_Directories(%args);
        my $alignment = $me->Make_Fasta_Align(number => $actual);
        $concat_job = $me->Concatenate_Searches(output=> ${output},);
    }

    my $parse_input = cwd() . qq"/$concat_job->{output}";
    my $comment_string = qq!## I don't know if this will work.!;
    my $job_string = qq?
use HPGL;
my \$h = new HPGL;
\$h->Parse_Search(input => '$parse_input', search_type => 'global_fasta',);
?;
    my $parse_job = $me->Qsub_Perl(job_name => "parse_search",
                              depends => $concat_job->{pbs_id},
                              job_string => $job_string,
                              qsub_shell => '/usr/bin/env perl',
                              comment => $comment_string,
        );
    return($concat_job);
}

sub Make_Fasta_Align {
    my $me = shift;
    my %args = @_;
    my $dep = '';
    $dep = $args{depends} if (defined($args{depends}));
    my $library;
    my $job;
    my $array_end = 1000 + $args{number};
    my $array_string = qq"1000-${array_end}";
    if ($args{parse}) {
        if ($args{output_type} ne '7') {
            print "If one wants to parse the output, it probably should be in the xml format, which may be done by adding --output_type 7 to this script.\n";
            print "This script will still try to parse a blast plain text output, but if it fails, don't say I didn't warn you.\n";
            print "Sleeping for 10 seconds so you can hit Control-C if you want to rerun.\n";
            sleep(10);
        }
    }
    my $job_string = '';
    if ($me->{pbs}) {
        $job_string = qq!
cd $me->{basedir}
mkdir -p $me->{basedir}/status $me->{basedir}/outputs && \\
 $me->{fasta_tool} $me->{fasta_args} -T $me->{cpus} $me->{query} $me->{library} \\
 1>$me->{basedir}/outputs/\${PBS_ARRAYID}.out \\
 2>>$me->{basedir}/split_align_errors.txt
!;
    } else {
        $job_string = qq!
cd $me->{basedir}
mkdir -p $me->{basedir}/status $me->{basedir}/outputs && \\
  $me->{fasta_tool} $me->{fasta_args} -T $me->{cpus} $me->{query} $me->{library} \\
  1>$me->{basedir}/outputs/$me->{fasta_tool}.out \\
  2>>$me->{basedir}/split_align_errors.txt
!;
    }
    my $comment = qq!## Running multiple fasta jobs.!;
    my $fasta_jobs = $me->Qsub(job_name => 'fasta_multi',
                               depends => $dep,
                               job_string => $job_string,
                               comment => $comment,
                               qsub_queue => 'long',
                               qsub_wall => '96:00:00',
                               qsub_args => " $me->{qsub_args} -t ${array_string} ",
        );
    return($fasta_jobs);
}

=item C<Concatenate_Searches>

The function Concatenate_Searches waits until the cluster finishes processing all the blast jobs, then
concatenates all the output files into one gzipped file.
If --parse is on, it will call Parse_Blast()

=cut
sub Concatenate_Searches {
    my $me = shift;
    my %args = @_;
    my $finished = 0;
    my $output = qq"$me->{basedir}/outputs/split_search.txt";
    $output = $me->{output} if (defined($me->{output}));
    $output = $args{output} if (defined($args{output}));
    my $comment_string = qq"## Concatenating the output files into ${output}\n";
    my $job_string = qq!
rm -f ${output} && for i in \$(/bin/ls outputs/*.out); do gzip -c \$i >> ${output}; done
!;
    my $concatenate_job = $me->Qsub(job_name => "concat",
                                    depends => $args{depends},
                                    job_string => $job_string,
                                    comment => $comment_string,
                                    output => qq"${output}",
        );
    return($concatenate_job);
}

=item C<Parse_Search>

    Use Parse_Search() to pass to different parsers.

=cut
sub Parse_Search {
    my $me = shift;
    my %args = @_;
    my @suffixes = ('.txt','.gz');
    $args{output} = basename($args{input}, @suffixes);
    $args{output} = basename($args{output}, @suffixes);
    if ($args{search_type} eq 'blast') {
        $me->Parse_Blast(%args);
    } elsif ($args{search_type} eq 'local_fasta') {
        $me->Parse_Fasta(%args);
    } elsif ($args{search_type} eq 'global_fasta') {
        $me->Parse_Fasta_Global(%args);
    } else {
        $me->Parse_Fasta(%args);
    }
}

=item C<Cleanup>

  Cleanup:  cleans up the mess of temporary files/directories left behind by this.

=cut
sub Cleanup {
    my $me = shift;
    my %args = @_;
    system("rm -rf outputs split status blastdb formatdb.log split_align_submission.sh");
    exit(0);
}

=item C<Parse_Blast>

    Parse_Blast is responsible for parsing blast output, it makes some
    attempt to be flexible vis a vis different formatting.  It prints a
    table of each gene and some relevant statistics about the hits that
    were found.

=cut
sub Parse_Blast {
    my $input = shift;
    my %args = @_;
    my $output = qq"${input}";
    $output =~ s/\.txt\.gz/_parsed\.txt/g;
    print "Writing parsed output to $output\n";
    my $count_table = qq"${input}";
    $count_table =~ s/\.txt\.gz/_counts\.txt/g;
    print "Writing count table to $count_table\n";
    my $counts = new FileHandle;
    $counts->open(">$count_table");
    $counts->autoflush(1);
    print $counts "Query\thits\n";
    my $fh = new FileHandle;
    $fh->open(">$output");
    $fh->autoflush(1);
    $| = 0;  ## I want to watch the output as it happens.
    print "Using /usr/bin/gunzip -c $input to read parsed blast output.\n";
    my $in = new FileHandle;
    $in->open("lesspipe $input |");
    $XML::SAX::ParserPackage = 'XML::SAX::PurePerl';
    my $searchio;
    if ($args{output_type} eq '0') { ## use old blast format
        $searchio = new Bio::SearchIO(-format => 'blast', -fh => $in, -best => 1,);
    } elsif ($args{output_type} eq '7') {  ## use blastxml
        $searchio = new Bio::SearchIO(-format => 'blastxml', -fh => $in, -best => 1,);
    } else { ## I don't know what it is, make searchio guess it
        $searchio = new Bio::SearchIO(-fh => $in, -best => 1,);
        my $test_format = $searchio->format();
        print "I am guessing the format is: $test_format\n";
        $searchio = undef;
        $searchio = new Bio::SearchIO(-format => $test_format, -fh => $in, -best => 1,);
    }
    print $fh "QUERYNAME\tChromosome\tStart\tEnd\t%ID\tScore\tSig\tCompLength\tHit_Ident\thits\n";
    my $result_count = 0;
  RESULT: while(my $result = $searchio->next_result()) {
      my $hit_count = 0;
      $result_count++;
      print "Parsed $result_count of $args{num_sequences} results.\n";
      my $query_name = "";
      while(my $hit = $result->next_hit) {
          $hit_count++;
          $query_name = $result->query_name();
          my $query_length = $result->query_length();
          my $accession = $hit->accession();
          my $acc2 = $hit->name();
          my $acc3 = $hit->description();
          my $length = $hit->length();
          my $score = $hit->raw_score();
          my $sig = $hit->significance();
          my $ident = ($hit->frac_identical() * 100);
          my ($start, $end, $hsp_identical, $hsp_cons);

        HSP: while (my $hsp = $hit->next_hsp) {
            $start = $hsp->start('subject');
            ## Maybe want $hsp->start('subject');
            $end = $hsp->end('subject');
            $hsp_identical = $hsp->frac_identical('total');
            $hsp_identical = sprintf("%.4f", $hsp_identical);
            ## $fun = sprintf("%2d %", $something, $somethingelse);
            $hsp_identical = $hsp_identical * 100;

#           $hsp_cons = $hsp->frac_cons('total');
            last HSP;  ## The 'HSP: ' is a label given to a logical loop
            ## If you then say next HSP; it will immediately stop processing the loop and move to the next iteration of the loop
            ## If instead you say last HSP; it will immediately break out of the loop.
        } ## next hsp

          ## This prints %identity with respect to the overlap of the component
          ## Then the score of the hit, its Evalue, and the component length.
          my $print_string = qq"${query_name}\t${acc2}\t${acc3}\t${start}\t${end}\t${ident}\t${score}\t${sig}\t${query_length}\t${hsp_identical}\t$hit_count\n";
          print $fh $print_string;
          if ($args{best_only}) {
              next RESULT;
          }
      } ## End each hit for a single result
      print $counts "${query_name}\t${hit_count}\n";
#    print "$result had $hit_count hits\n";
  } ## Finish looking at each sequence
    $fh->close();
    $in->close();
    $counts->close();
}

=item C<Parse_Global>

    Parse_Global makes a hard assumption about the structure of the hits in the output of the alignment program.  They should be
    in a global/global search, thus we assume 1 hit / 1 result.

=cut
sub Parse_Fasta_Global {
    my $me = shift;
    my %args = @_;
    $args{best_hit_only} = 0 unless (defined($args{best_hit_only}));
    $args{format} = 'fasta' unless (defined($args{format}));
    $args{max_significance} = 0.001 unless (defined($args{max_significance}));
    $args{many_cutoff} = 10 unless (defined($args{many_cutoff}));
    $args{min_score} = undef unless (defined($args{min_score}));
    $args{check_all_hits} = 0 unless (defined($args{check_all_hits}));
    $args{min_percent} = undef unless (defined($args{min_percent}));
    $args{output} = $args{input} unless($args{output});

    my $in = new FileHandle;
    $in->open("lesspipe $args{input} |");
    my $searchio = new Bio::SearchIO(-format => $args{format}, -fh => $in,
                                     -best_hit_only => $args{best_hit_only} , -max_significance => $args{max_significance},
                                     -check_all_hits => $args{check_all_hits}, -min_score => $args{min_score},
        );
    ##my $outdir = qq"outputs/ggsearch";
    ##print STDERR "WTF: $outdir\n";
    ##make_path("${outdir}") unless (-d ${outdir});
    my $counts = new FileHandle;
    my $parsed = new FileHandle;
    my $singles = new FileHandle;
    my $doubles = new FileHandle;
    my $few = new FileHandle;
    my $many = new FileHandle;
    my $zero = new FileHandle;
    my $all = new FileHandle;
    ##my $count_file = qq"${outdir}/$args{output}.count";
    ##my $parsed_file = qq"${outdir}/$args{output}.parsed";
    ##my $single_file = qq"${outdir}/$args{output}_singles.txt";
    ##my $double_file = qq"${outdir}/$args{output}_doubles.txt";
    ##my $few_file = qq"${outdir}/$args{output}_few.txt";
    ##my $many_file = qq"${outdir}/$args{output}_many.txt";
    ##my $zeroes_file = qq"${outdir}/$args{output}_zero.txt";
    ##my $all_file = qq"${outdir}/$args{output}_all.txt";
    my $count_file = qq"$args{output}.count";
    my $parsed_file = qq"$args{output}.parsed";
    my $single_file = qq"$args{output}_singles.txt";
    my $double_file = qq"$args{output}_doubles.txt";
    my $few_file = qq"$args{output}_few.txt";
    my $many_file = qq"$args{output}_many.txt";
    my $zeroes_file = qq"$args{output}_zero.txt";
    my $all_file = qq"$args{output}_all.txt";

    $counts->open(">$count_file");
    $parsed->open(">$parsed_file");
    $singles->open(">$single_file");
    $doubles->open(">$double_file");
    $few->open(">$few_file");
    $many->open(">$many_file");
    $zero->open(">$zeroes_file");
    $all->open(">$all_file");

    my $seq_count = 0;
    print $parsed "Query Name\tQuery length\tHit ID\tHit Length\tScore\tE\tIdentity\tHit length\tHit Matches\n";
    while (my $result = $searchio->next_result()) {
        $seq_count++;
        my $query_name = $result->query_name();
        my $entry = qq"${query_name}\t";
        my $hit_count = 0;
        HITLOOP: while(my $hit = $result->next_hit) {

            ##my $query_name = $result->query_name();
            my $query_length = $result->query_length();
            my $accession = $hit->accession();
            my $acc2 = $hit->name();
            my $length = $hit->length();
            my $score = $hit->raw_score();
            my $sig = $hit->significance();
            my $ident = $hit->frac_identical();
            my $hit_len;
            my $hit_matches;
            if (defined($args{min_percent})) {
                next HITLOOP if ($ident < $args{min_percent});
            }
            ## If we pass the min_percet threshold, increment the number of hits.
            $hit_count++;
            while (my $hsp = $hit->next_hsp) {
                $hit_len = $hsp->length('total');
                my @matches = $hsp->matches(-seq => 'hit');
                $hit_matches = $matches[1];
            }
            $entry .= "${acc2}:${ident}:${sig} ";
            print $parsed "$query_name\t$query_length\t$acc2\t$length\t$score\t$sig\t$ident\t$hit_len\t$hit_matches\n";
        } ## End iterating over every hit for a sequence.
        $entry .= "\n";
        print $all $entry;
        if ($hit_count == 1) {
            print $singles $entry;
        } elsif ($hit_count == 2) {
            print $doubles $entry;
        } elsif ($hit_count >= 3 and $hit_count <= $args{many_cutoff}) {
            print $few $entry;
        } elsif ($hit_count > $args{many_cutoff}) {
            print $many $entry;
        } else {
            print $zero $entry;
        }
        print $counts "${query_name}\t${hit_count}\n";
    } ## End iterating over every result in searchio
    $counts->close();
    $parsed->close();
    $singles->close();
    $doubles->close();
    $few->close();
    $many->close();
    $zero->close();
    $all->close();
    $in->close();
    return($seq_count);
}


=item C<Check_Blastdb>

    Check_Blastdb makes sure that there is an appropriately formatted
    blastdb for the library and type of blast performed.

=cut
sub Check_Blastdb {
    my $me = shift;
    my %args = @_;
    my $formatdb_ret = "not run.";
    if (!-d qq"$me->{basedir}/blastdb" and !-f qq"$me->{basedir}/blastdb") {
        make_path("$me->{basedir}/blastdb");
    }
    my $found_lib = 0;
    $args{peptide} = ' -p T ';
    if ($args{tool} eq 'blastn' or $args{tool} eq 'tblastx') {
        $args{peptide} eq ' -p F ';
    }

    if ($args{tool} =~ m/blast/ and $args{peptide} eq 'T') {
        print "Checking for blastdb/$args{library}.psq\n";
        if (-f "blastdb/$args{library}.psq") {
            print "Found a protein library for $args{library}.\n";
            $found_lib++;
        }
    } else {
        print "Checking for blastdb/$args{library}.nsq\n";
        if (-f "blastdb/$args{library}.nsq") {
            print "Found a nucleotide library for $args{library}.\n";
            $found_lib++;
        }
    }
    if (!$found_lib) {
        my $formatdb_command = qq"formatdb $args{peptide} -o T -n blastdb/$args{library} -s -i $args{library}.fasta";
        print "The formatdb command run is: $formatdb_command\n";
        $formatdb_ret = qx"$formatdb_command";
    }
    return($formatdb_ret);
}

=item C<Get_Split>

    Get_Split takes the number of sequences in the query library and divides by the number of jobs to run,
    takes the ceiling of that, and prints that many sequences to each file in split/ starting at 1000.
    (So, as long as you have <= 8,999 jobs PBS won't get argsused by the difference between 099 and 100
    because they will be 1099 and 1100 instead.)

=cut
sub Get_Split {
    my $me = shift;
    my %args = @_;
    my $in = new Bio::SeqIO(-file => $args{query},);
    my $seqs = 0;
    while (my $in_seq = $in->next_seq()) {
        $seqs++;
    }
    my $ret = ceil($seqs / $args{number});
    if ($seqs < $args{number}) {
        print "There are fewer sequences than the chosen number of split files, resetting that to $seqs.\n";
        ##$splits = $seqs;
        $args{number} = $seqs;
    }
    return($ret);
}

=item C<Make_Directories>

    Create sequential directories for each job and actually write the sequences using Bio::SeqIO

=cut
sub Make_Directories {
    my $me = shift;
    my %args = @_;
    my $num_per_split = $args{num_per_split};
    my $splits = $args{number};
    ## I am choosing to make directories starting at 1000
    ## This way I don't have to think about the difference from
    ## 99 to 100 (2 characters to 3) as long as no one splits more than 9000 ways...
    my $dir = 1000;

    for my $c ($dir .. ($dir + $splits)) {
        ## print "Making directory: split/$c\n";
        make_path("split/$c") if (!-d "split/$c" and !-f "split/$c");
    }

    my $in = new Bio::SeqIO(-file => $args{query},);
    my $count = 0;
    my $outfile = new FileHandle;
    $outfile->open(">>split/${dir}/in.fasta");
    while (my $in_seq = $in->next_seq()) {
        my $id = $in_seq->id();
        my $seq = $in_seq->seq();
        print $outfile qq!>$id
$seq
!;
        $count++;
        $args{num_sequences}++;
        if ($count >= $num_per_split) {
            $count = 0;
            $dir++;
            my $last = $dir;
            $last--;
            if ($count == 1) {
                print "Writing $num_per_split entries to files 1000 to $last\n";
            }
            ## You might be wondering why this num_dirs is here.
            ## Imagine if you have a query library of 10,004 sequences and you try to write them to 200 files.
            ## This script will write 51 sequences / file, but will therefore not quite finish writing all 200 files.
            ## As a result, when we go back at the end to clean up,
            ## the loop which tries to detect how many jobs are finished will never exit
            ## because it won't ever reach 200.
            ## Therefore this keeps track of the last file written and uses that instead.
            $args{num_dirs} = ($last - 999); ## This - 1000 is there because we start at job 1000, but we start counting at 0
        } ## End for each iteration of $num_per_split files
    } ## End while reading the fasta
    $outfile->close();
    my $actual_number_dirs_used = $args{num_dirs};
    return($actual_number_dirs_used);
}

=item C<Make_Align>

    Perform the alignment with qsub, this writes a file
    'split_align_submission.sh' and passes it to qsub with the options
    specified above.

=cut
sub Make_Blast_Align {
    my $me = shift;
    my %args = @_;
    my $library;
    my $job;
    my $array_end = 1000 + $args{number};
    my $array_string = qq"1000-${array_end}";
    my $dir = getcwd;
    if ($args{parse}) {
        if ($args{output_type} ne '7') {
            print "If one wants to parse the output, it probably should be in the xml format, which may be done by adding --output_type 7 to this script.\n";
            print "This script will still try to parse a blast plain text output, but if it fails, don't say I didn't warn you.\n";
            print "Sleeping for 10 seconds so you can hit Control-C if you want to rerun.\n";
            sleep(10);
        }
    }
    my $string = qq?
cd ${dir}
mkdir -p ${dir}/status ${dir}/outputs && \\
 blastall -m $args{output_type} $args{blast_params} -p $args{blast} -d ${dir}/blastdb/${library} \\
 -i ${dir}/split/\${PBS_ARRAYID}/in.fasta \\
 1>${dir}/outputs/\${PBS_ARRAYID}.out 2>>${dir}/split_align_errors.txt
?;
    print "The blast jobs look like this:
$string
";
    ## The next time I do a blast search and look at this code, refactor it to use Qsub.pm
    ## Doing it this way is dumb, but I need to leave in a moment so I don't have time to fix it.
    my $script = new FileHandle;
    $script->open(">split_align_submission.sh");
    print $script "#PBS -t ${array_string} -V -d $dir -q $args{queue} -l walltime=18:00:00,mem=8gb -jeo -e $dir/status/split_align.out -m n\n";
    print $script $string;
    $script->close();
    my $submit = qq"qsub split_align_submission.sh";
    my $cmd = new FileHandle;
    $cmd->open("$submit |");
    while(my $line = <$cmd>) {
        chomp $line;
        print "The submitted job is: $line\n";
        $job = $line;
    }
    $cmd->close();
}

=back

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

    L<Bio::Tools::Run::StandAloneBlast>

=cut

1;
