package CYOA;
use common::sense;
use autodie qw":all";

use Cwd;
use Bio::Seq;
use Bio::SearchIO::blast;
use Bio::SearchIO::fasta;
use Bio::Tools::Run::StandAloneBlast;
use POSIX qw"ceil";

=head1 NAME
    CYOA::Alignment - Perform Sequence alignments with tools like Blast/Fasta

=head1 SYNOPSIS

    use CYOA;
    my $hpgl = new CYOA;
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
sub Run_Parse_Blast {
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

=item C<Split_Align>

    $hpgl->Split_Align();
    Split apart a set of query sequences into $args{number} pieces and align them all separately.

=cut
sub Split_Align_Blast {
    my $me = shift;
    my %args = @_;
    $me->Check_Options(["query","library", "blast_tool"]);
    $args{query} = $me->{query} if ($me->{query});
    $args{library} = $me->{library} if ($me->{library});
    $args{blast_tool} = $me->{blast_tool} if ($me->{blast_tool});

    $args{param} = ' -e 10 ' unless (defined($args{param}));
    $args{blast_tool} = 'blastn' unless (defined($args{blast_tool}));
    $args{number} = 40 unless (defined($args{number}));
    $args{parse} = 0 unless (defined($args{parse}));
    $args{num_dirs} = 0 unless (defined($args{num_dirs}));
    $args{best_only} = 0 unless (defined($args{best_only}));
    $args{queue} = 'workstation' unless (defined($args{queue}));
    $args{num_sequences} = 0 unless (defined($args{num_sequences}));

    ## This might be wrong (tblastx)
    print STDERR qq"A quick reminder because I (atb) get confused easily:
tblastn is a protein fasta query against a nucleotide blast database.
tblastx is a nucleotide fasta query (which is translated on the fly) against a protein blast db.
blastx is a nucleotide fasta query against a protein blast db.
blastn is normal nucleotide/nucleotide
blastp is normal protein/protein.
";
    ## Also set the default blastdb if it didn't get set.
    if ($args{blast_tool} eq 'blastn' or $args{blast_tool} eq 'tblastx') {
        $args{peptide} = 'F';
        $args{library} = 'nt' if (!defined($args{library}));
    } else {
        $args{peptide} = 'T';
        $args{library} = 'nr' if (!defined($args{library}));
    }
    my $lib = basename($args{library}, ('.fasta'));
    my $que = basename($args{query}, ('.fasta'));
    my $outdir = qq"$me->{basedir}/outputs";
    make_path("${outdir}") unless(-d ${outdir});
    my $output = qq"${outdir}/${que}_vs_${lib}.txt";
    my $concat_job;
    $lib = $me->Check_Blastdb(%args);
    if ($me->{pbs}) {
        my $num_per_split = $me->Get_Split(%args);
        $args{num_per_split} = $num_per_split;
        print "Going to make $args{number} directories with $num_per_split sequences each.\n";
        my $actual = $me->Make_Directories(%args);
        print "Actually used ${actual} directories to write files.\n";
        my $alignment = $me->Make_Blast_Job(number => $actual, library => $lib, output_type => 5);  # Forcing the blast output type to be 5: blastxml
        $concat_job = $me->Concatenate_Searches(depends => $alignment->{pbs_id}, output => ${output},);
    } else {
        ## If we don't have pbs, force the number of jobs to 1.
        $args{number} = 1;
        my $num_per_split = $me->Get_Split(%args);
        my $actual = $me->Make_Directories(%args);
        my $alignment = $me->Make_Blast_Job(number => $actual, library => $lib, output_type => 5);
        $concat_job = $me->Concatenate_Searches(output=> ${output},);
    }

    my $parse_input = cwd() . qq"/$concat_job->{output}";
    my $comment_string = qq!## I don't know if this will work.!;
    my $job_string = qq?
use CYOA;
my \$h = new CYOA;
\$h->Parse_Search(input => '$parse_input', search_type => 'blastxml', best => $args{best_only});
?;
    my $parse_job = $me->Qsub_Perl(job_name => "parse_search",
                              depends => $concat_job->{pbs_id},
                              job_string => $job_string,
                              comment => $comment_string,
        );
    return($concat_job);
}

sub Make_Blast_Job {
    my $me = shift;
    my %args = @_;
    my $dep = '';
    $dep = $args{depends} if (defined($args{depends}));
    my $library = $args{library};
    my $job;
    my $array_end = 1000 + $args{number};
    my $array_string = qq"1000-${array_end}";
    ## The fine folks at NCBI changed the numbers of the blast outputs!
    ## 0=pairwise ; 1=query_anchored identies ; 2=query_anchored no identities ;
    ## 3=flat identities ; 4=flat no identities ; 5=blastxml ; 6=tabular
    ## 7=tabular commented ; 8=text asn1 ; 9=binary asn1 ; 10=csv ; 11=blast archive asn1
    if ($args{output_type} ne '5') {
        print STDERR "If one wants to parse the output, it probably should be in the xml format, which may be done by adding --output_type 7 to this script.\n";
        print STDERR "This script will still try to parse a blast plain text output, but if it fails, don't say I didn't warn you.\n";
        print STDERR  "Sleeping for 3 seconds so you can hit Control-C if you want to rerun.\n";
        sleep(3);
    }
    my $job_string = '';
    ## I have been getting null pointers, there is a bug report suggesting this is because of xml output
    $args{output_type} = 1;
    if ($me->{pbs}) {
        $job_string = qq!
cd $me->{basedir}
$me->{blast_tool} -outfmt $args{output_type} \\
 -query $me->{basedir}/split/\${PBS_ARRAYID}/in.fasta \\
 -db ${library} \\
 -out $me->{basedir}/outputs/\${PBS_ARRAYID}.out \\
 1>$me->{basedir}/outputs/\${PBS_ARRAYID}.stdout \\
 2>>$me->{basedir}/split_align_errors.txt
!;
    } else {
        $job_string = qq!
cd $me->{basedir}
$me->{blast_tool} -outfmt $args{output_type} \\
  -query $me->{query} \\
  -db ${library} \\
  -out $me->{basedir}/outputs/$me->{blast_tool}.out \\
  1>$me->{basedir}/outputs/$me->{blast_tool}.stdout \\
  2>>$me->{basedir}/split_align_errors.txt
!;
    }
    my $comment = qq!## Running multiple fasta jobs.!;
    my $blast_jobs = $me->Qsub(job_name => 'blast_multi',
                               depends => $dep,
                               job_string => $job_string,
                               comment => $comment,
                               qsub_queue => 'large',
                               qsub_wall => '96:00:00',
                               qsub_mem => 32,
                               qsub_args => " $me->{qsub_args} -t ${array_string} ",
        );
    return($blast_jobs);
}

sub Merge_Parse_Blast {
    my $me = shift;
    my %args = @_;
    $me->Check_Options(['output']);
    my $concat = $me->Concatenate_Searches();
    $me->{input} = $concat->{output};
    my $job_string = qq?
use CYOA;
my \$h = new CYOA;
\$h->{input} = '$me->{input}';
\$h->Parse_Search(search_type => 'blastxml',);
?;
    my $parse = $me->Qsub_Perl(job_name => 'parse_search',
                               depends => $concat->{pbs_id},
                               job_string => $job_string,
        );
    return($parse);
}

=item C<Parse_Blast>

    Parse_Blast is responsible for parsing blast output, it makes some
    attempt to be flexible vis a vis different formatting.  It prints a
    table of each gene and some relevant statistics about the hits that
    were found.

=cut
sub Parse_Blast {
    my $me = shift;
    my %args = @_;
    $me->Check_Options(["input",]);
    my $input = $me->{input};
    my $output = qq"${input}";
    my $best = $args{best_only};
    $best = 0 if (!defined($best));
    my $search_type = $args{search_type};
    $search_type = 'blastxml' if (!defined($search_type));
    $output =~ s/\.txt\.gz/_parsed\.txt/g;
    print "Writing parsed output to $output\n";
    my $count_table = qq"${input}";
    $count_table =~ s/\.txt\.gz/_counts\.txt/g;
    print "Writing count table to $count_table\n";
    $args{num_sequences} = 'unknown' if (!defined($args{num_sequences}));
    my $counts = new FileHandle;
    $counts->open(">$count_table");
    $counts->autoflush(1);
    print $counts "Query\tquery_accession\thits\n";
    my $fh = new FileHandle;
    $fh->open(">$output");
    $fh->autoflush(1);
    $| = 0;  ## I want to watch the output as it happens.
    my $in = new FileHandle;
    $in->open("less ${input} |");
    $XML::SAX::ParserPackage = 'XML::SAX::PurePerl';
    ##my $searchio = new Bio::SearchIO(-fh => $in);
    ## Perform a test for the output format, this will close the filehandle, so reopen it after
    ##my $guessed_format = $searchio->format();
    ##print STDERR "Guessed the format is: $guessed_format\n";
    ##$guessed_format = 'blastxml';
    ##print STDERR "It turns out that searchio->format() is crap, setting format to 'blastxml'\n";
    ##$searchio = undef;
    ##$in->close();
    ##$in = new FileHandle;
    ##$in->open("less ${input} |");
    my $searchio = new Bio::SearchIO(-fh => $in, -format => $search_type, -best => ${best},);
##    if ($args{output_type} eq '0') { ## use old blast format
##        $searchio = new Bio::SearchIO(-format => 'blast', -fh => $in, -best => 1,);
##    } elsif ($args{output_type} eq '7') {  ## use blastxml
##        $searchio = new Bio::SearchIO(-format => 'blastxml', -fh => $in, -best => 1,);
##    } else { ## I don't know what it is, make searchio guess it
##        $searchio = new Bio::SearchIO(-fh => $in, -best => 1,);
##        my $test_format = $searchio->format();
##        print "I am guessing the format is: $test_format\n";
##        $searchio = undef;
##        $searchio = new Bio::SearchIO(-format => $test_format, -fh => $in, -best => 1,);
    ##    }
    print $fh "QUERYNAME\treal_name\tChromosome\tStart\tEnd\t%ID\tScore\tSig\tCompLength\tHit_Ident\thits\n";
    my $result_count = 0;
  RESULT: while(my $result = $searchio->next_result()) {
      my $hit_count = 0;
      $result_count++;
      print STDERR "Parsed $result_count of $args{num_sequences} results.\n";
      my $query_name = "";
      my $real_name = "";
      while(my $hit = $result->next_hit) {
          $hit_count++;
          ## $query_name = $result->query_name();
          $query_name = $result->query_name();
          $real_name = $result->query_description();
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
          my $print_string = qq"${query_name}\t${real_name}\t${acc2}\t${acc3}\t${start}\t${end}\t${ident}\t${score}\t${sig}\t${query_length}\t${hsp_identical}\t$hit_count\n";
          print $fh $print_string;
          if ($args{best_only}) {
              next RESULT;
          }
      } ## End each hit for a single result
      print $counts "${query_name}\t${real_name}\t${hit_count}\n";
#    print "$result had $hit_count hits\n";
  } ## Finish looking at each sequence
    $fh->close();
    $in->close();
    $counts->close();
}

=item C<Check_Blastdb>

    Check_Blastdb makes sure that there is an appropriately formatted
    blastdb for the library and type of blast performed.

=cut
sub Check_Blastdb {
    my $me = shift;
    my %args = @_;
    my $formatdb_ret = "not run.";
    my $libname = $args{library};
    ## First check for the relevant library in $ENV{BLASTDB}
    ## If it isn't there, make one in basedir/blastdb/
    my $foundlib = 0;
    my $checklib = "";
    my $checklib_zero = "";
    if ($args{blast_tool} =~ m/blast/ and $args{peptide} eq 'T') {
        $checklib_zero = qq"$args{library}.00.psq";
        $checklib = qq"$args{library}.psq";
    } else {
        $checklib_zero = qq"$args{library}.00.nsq";
        $checklib = qq"$args{library}.nsq";
    }
    my $lib = "";
    my $libname = "";
    print "Looking for ${checklib} / ${checklib_zero} in either $ENV{BLASTDB} or $me->{basedir}/blastdb\n";
    ## Start with BLASTDB
    if (-f "$ENV{BLASTDB}/${checklib}" or -f "$ENV{BLASTDB}/${checklib_zero}") {
        $foundlib++;
        $libname = qq"$ENV{BLASTDB}/$args{library}";
        $lib = qq"$ENV{BLASTDB}/$args{library}";
    } elsif (-f "$me->{basedir}/blastdb/${checklib}" or -f "$me->{basedir}/blastdb/${checklib_zero}") { ## Then try basedir
        $foundlib++;
        $libname = qq"$me->{basedir}/blastdb/$args{library}";
        $lib = qq"$me->{basedir}/blastdb/$args{library}";
    }

    if (!$foundlib) {
        if (!-d qq"$me->{basedir}/blastdb") {
            make_path(qq"$me->{basedir}/blastdb");
        }
        $libname = qq"$me->{basedir}/blastdb/$args{library}";
        my $formatdb_command = qq"formatdb -p $args{peptide} -o T -n blastdb/$args{library} -s -i $args{library}";
        print "The formatdb command run is: $formatdb_command\n";
        $formatdb_ret = qx"$formatdb_command";
    }
    return($libname);
}


=back

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

    L<Bio::Tools::Run::StandAloneBlast>

=cut

1;
