package HPGL;
use common::sense;
use autodie qw":all";

use Bio::Seq;
use Bio::SearchIO::blast;
use Bio::Tools::Run::StandAloneBlast;


sub Blast_Parse {
    my $me = shift;
    my %args = @_;
    my $query = $me->{input};
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
    for my $library (@search_libraries) {
        my $query_library = new Bio::SeqIO(-file => $query, -format => 'Fasta');
        my $output_directory = qq"${query}_vs_${library}";
        if (!-d $output_directory) {
            system("mkdir $output_directory");
        }
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

=head1 NAME

split_align_blast

=head1 SYNOPSIS

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

=head1 DESCRIPTION

This script seeks to make running large blast jobs easier.  It will split the input query file
into many smaller files, format a blastdb appropriate for the alignment, submit them to the cluster,
collect the alignments, parse the outputs, and cleanup the mess.

=head1 AUTHOR - atb

Email abelew@gmail.com

=cut

sub Split_Align_Blast {
    my %conf = (
        number => 200,
        input => undef,
        lib => undef,
        blast_params => ' -e 10 ',
        blast => 'blastn',
        parse => 1,
        peptide => 'F',
        output_type => '0',
        output => 'split_align_output.txt',
        formatdb => 'formatdb -p $peptide -o T -n blastdb/$lib -s -i',
        help => undef,
        num_dirs => 0,
        clean => undef,
        best_only => undef,
        queue => 'workstation',
        skip => undef,
        num_sequences => 0,
        );
    my %options = (
        'best_only' => \$conf{best_only}, ## Whether to only print 1 result/gene
        'blast|b:s' => \$conf{blast},  ## Use --blast or -b followed by 'blastp' or 'tblastx' or whatever to set the blast program.
        'blast_params|bp:s' => \$conf{blast_params},
        'clean' => \$conf{clean}, ## Clean up in case of an incomplete run
        'formatdb|f:s' => \$conf{formatdb}, ## the formatdb command used to index a fasta library of sequences
        'help|h' => \$conf{help},
        'input|i:s' => \$conf{input},
        'lib|l:s' => \$conf{lib},
        'number|n:s' => \$conf{number},
        'output|o:s' => \$conf{output},
        'output_type|m:s' => \$conf{output_type}, ## Blast's output type can only be parsed if it is 7, 0 is default
        'parse:s' => \$conf{parse},
        'peptide|p:s' => \$conf{peptide}, ## T for a peptide database, F for nucleotide
        'queue|q:s' => \$conf{queue},  ## Set the pbs queue, right now throughput isn't moving -- for example
        'skip' => \$conf{skip},  ## If you just want to parse the output
        );
    my $argv_result = GetOptions(%options);
    my $job = "";

    if ($conf{clean}) {
        Cleanup();
    }

    if ($conf{skip}) {
        print "Skipping blast, just parsing $ARGV[0]\n";
        Parse_Blast($ARGV[0]);
        exit(0);
    }

    if ($conf{help}) {
        pod2usage();
        exit(0);
    }

    my $library = basename($conf{lib}, ('.fasta'));

    ## Number entries 23119
    ## Thus 116 entries in each fasta

    ##  To arrive at this number, I just rounded up (/ 23119 200)
    ## However, we can count the number of fasta entries in the input file
    my $num_per_split = Get_Split();
    print "Going to make $conf{number} directories with $num_per_split sequences each.\n";
    Make_Directories($num_per_split);
    Check_Blastdb();
    Make_Align();
    Wait_Concat();
}

=head2 Wait_Concat

The function Wait_Concat waits until the cluster finishes processing all the blast jobs, then
concatenates all the output files into one gzipped file.
If --parse is on, it will call Parse_Blast()

=cut

sub Wait_Concat {
    my $finished = 0;
    my %conf;
    while (!$finished) {
        my $num_finished = 0;
        open(FIN, '/bin/ls status/*.finished 2>/dev/null | wc -w |');
        while(my $line = <FIN>) {
            chomp $line;
            $num_finished = $line;
        }
        close(FIN);
        print "$num_finished jobs are completed out of $conf{num_dirs}.\n";
        if ($num_finished >= $conf{num_dirs}) {
            $finished++;
        } else {
            print "Waiting 20 seconds to see if the blast jobs are finished.\n";
            sleep(20);
        }
    }
    print "Concatenating the output files into $conf{output}.gz\n";
    ## waiting a couple seconds to make sure the files get written to disk
    sleep(20);
    open(CONCAT, "rm -f $conf{output}.gz && for i in \$(/bin/ls outputs/*.out); do echo \"concatenating \$i\"; gzip -c \$i >> $conf{output}.gz; done && rm -f finished/* |");
    while(my $line = <CONCAT>) {
        print $line;
    }
    close(CONCAT);
    if ($conf{parse}) {
        Parse_Blast("$conf{output}.gz");
    }
    print "Unless you hit control-C in the next minute, this script will delete everything except the compressed output file.\n";
    sleep(60);
    Cleanup();
}

=head2 Cleanup

Cleanup:  cleans up the mess of temporary files/directories left behind by this.

=cut

sub Cleanup {
    system("rm -rf outputs split status blastdb formatdb.log split_align_submission.sh");
    exit(0);
}

=head2 Parse_Blast

Parse_Blast is responsible for parsing blast output, it makes some attempt to be flexible vis a vis different formatting
It prints a table of each gene and some relevant statistics about the hits that were found.

=cut

sub Parse_Blast {
    my $input = shift;
    my %conf;
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
    open(F,"gunzip -c $input |");
    $XML::SAX::ParserPackage = 'XML::SAX::PurePerl';
    my $searchio;
    if ($conf{output_type} eq '0') { ## use old blast format
        $searchio = new Bio::SearchIO(-format => 'blast', -fh => \*F, -best => 1,);
    } elsif ($conf{output_type} eq '7') {  ## use blastxml
        $searchio = new Bio::SearchIO(-format => 'blastxml', -fh => \*F, -best => 1,);
    } else { ## I don't know what it is, make searchio guess it
        $searchio = new Bio::SearchIO(-fh => \*F, -best => 1,);
        my $test_format = $searchio->format();
        print "I am guessing the format is: $test_format\n";
        $searchio = undef;
        $searchio = new Bio::SearchIO(-format => $test_format, -fh => \*F, -best => 1,);
    }
    print $fh "QUERYNAME\tChromosome\tStart\tEnd\t%ID\tScore\tSig\tCompLength\tHit_Ident\thits\n";
    my $result_count = 0;
  RESULT: while(my $result = $searchio->next_result()) {
      my $hit_count = 0;
      $result_count++;
      print "Parsed $result_count of $conf{num_sequences} results.\n";
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
          if ($conf{best_only}) {
              next RESULT;
          }
      } ## End each hit for a single result
      print $counts "${query_name}\t${hit_count}\n";
#    print "$result had $hit_count hits\n";
  } ## Finish looking at each sequence
    close(F);
    $fh->close();
    $counts->close();
}

=head2 Check_Blastdb

Check_Blastdb makes sure that there is an appropriately formatted blastdb for the library and type of blast performed.

=cut

sub Check_Blastdb {
    my %conf;
    my $lib = $conf{lib};
    my $library = $lib;
    if (!-d 'blastdb') {
        mkdir('blastdb');
    }
    my $found_lib = 0;
    if ($conf{blast} eq 'blastn') {
        $conf{peptide} eq 'F';
    } elsif ($conf{blast} eq 'blastp') {
        $conf{peptide} eq 'T';
    } elsif ($conf{blast} eq 'blastx') {
        $conf{peptide} eq 'T';
    } elsif ($conf{blast} eq 'tblastx') {
        $conf{peptide} eq 'F';
    }

    if ($conf{peptide} eq 'T') {
        print "Checking for blastdb/${library}.psq\n";
        if (-f "blastdb/${library}.psq") {
            print "Found a protein library for $lib.\n";
            $found_lib++;
        }
    } else {
        print "Checking for blastdb/${library}.nsq\n";
        if (-f "blastdb/${library}.nsq") {
            print "Found a nucleotide library for $lib.\n";
            $found_lib++;
        }
    }
    if (!$found_lib) {
        my $formatdb = qq"$conf{formatdb} $lib";
        $formatdb =~ s/\$peptide/$conf{peptide}/g;
        $formatdb =~ s/\$lib/$library/g;
        print "The formatdb command run is: $formatdb\n";
        system($formatdb);
    }
}

=head2 Get_Split

Get_Split takes the number of sequences in the query library and divides by the number of jobs to run,
takes the ceiling of that, and prints that many sequences to each file in split/ starting at 1000.
(So, as long as you have <= 8,999 jobs PBS won't get confused by the difference between 099 and 100
because they will be 1099 and 1100 instead.)

=cut

sub Get_Split {
    my %conf;
    my $input = $conf{input};
    my $splits = $conf{number};
    my $in = new Bio::SeqIO(-file => $input,);
    my $seqs = 0;
    while (my $in_seq = $in->next_seq()) {
        $seqs++;
    }
    my $ret = ceil($seqs / $splits);
    if ($seqs < $splits) {
        print "There are fewer sequences than the chosen number of split files, resetting that to $seqs.\n";
        $splits = $seqs;
        $conf{number} = $seqs;
    }
    return($ret);
}

=head2 Make_Directories

Create sequential directories for each job and actually write the sequences using Bio::SeqIO

=cut

sub Make_Directories {
    my %conf;
    my $num_per_split = shift;
    my $splits = $conf{number};
    ## I am choosing to make directories starting at 1000
    ## This way I don't have to think about the difference from
    ## 99 to 100 (2 characters to 3) as long as no one splits more than 9000 ways...
    my $dir = 1000;

    for my $c ($dir .. ($dir + $splits)) {
##      print "Making directory: split/$c\n";
        system("mkdir -p split/$c");
    }

    my $in = new Bio::SeqIO(-file => $conf{input},);
    my $count = 0;
    while (my $in_seq = $in->next_seq()) {
        my $id = $in_seq->id();
        my $seq = $in_seq->seq();
        open(OUT, ">>split/$dir/in.fasta");
        print OUT ">$id
$seq
";
        close(OUT);
        $count++;
        $conf{num_sequences}++;
        if ($count >= $num_per_split) {
            $count = 0;
            $dir++;
            my $last = $dir;
            $last--;
            print "Wrote $num_per_split entries to $last\n";
            ## You might be wondering why this num_dirs is here.
            ## Imagine if you have a query library of 10,004 sequences and you try to write them to 200 files.
            ## This script will write 51 sequences / file, but will therefore not quite finish writing all 200 files.
            ## As a result, when we go back at the end to clean up,
            ## the loop which tries to detect how many jobs are finished will never exit
            ## because it won't ever reach 200.
            ## Therefore this keeps track of the last file written and uses that instead.
            $conf{num_dirs} = ($last - 999); ## This - 1000 is there because we start at job 1000, but we start counting at 0
        } ## End for each iteration of $num_per_split files
    } ## End while reading the fasta
}

=head2 Make_Align

Perform the alignment with qsub, this writes a file 'split_align_submission.sh' and passes it to qsub
with the options specified above.

=cut

sub Make_Align {
    my %conf;
    my $library;
    my $job;
    my $array_end = 1000 + $conf{number};
    my $array_string = qq"1000-${array_end}";
    use Cwd;
    my $dir = getcwd;
##    my $string = qq?
##cat <<"EOF" | qsub -t ${array_string} -V -d $dir -q throughput -l walltime=18:00:00,mem=8Gb -j eo -e $dir/split_align.out -m n -
##cd ${dir}
##CMD="mkdir -p ${dir}/stats ${dir}/outputs && \
## blastall -e 1000 -p $conf{blast} -d ${dir}/blast/${library} \
## -i ${dir}/split/\${PBS_ARRAYID}/in.fasta \
## -o ${dir}/outputs/\${PBS_ARRAYID}.out \
## && touch ${dir}status/\${PBS_ARRAYID}.finished"
##echo \$CMD
##eval \$CMD
##EOF
##?;
    if ($conf{parse}) {
        if ($conf{output_type} ne '7') {
            print "If one wants to parse the output, it probably should be in the xml format, which may be done by adding --output_type 7 to this script.\n";
            print "This script will still try to parse a blast plain text output, but if it fails, don't say I didn't warn you.\n";
            print "Sleeping for 10 seconds so you can hit Control-C if you want to rerun.\n";
            sleep(10);
        }
    }
    my $string = qq?
cd ${dir}
CMD="mkdir -p ${dir}/status ${dir}/outputs && \
 blastall -m $conf{output_type} $conf{blast_params} -p $conf{blast} -d ${dir}/blastdb/${library} \
 -i ${dir}/split/\${PBS_ARRAYID}/in.fasta \
 1> ${dir}/outputs/\${PBS_ARRAYID}.out 2>>${dir}/split_align_errors.txt \
 && touch ${dir}/status/\${PBS_ARRAYID}.finished"
echo \$CMD
eval \$CMD
?;
    print "The blast jobs look like this:
$string
";
    open(SCRIPT, ">split_align_submission.sh");
    print SCRIPT "#PBS -t ${array_string} -V -d $dir -q $conf{queue} -l walltime=18:00:00,mem=8gb -jeo -e $dir/status/split_align.out -m n\n";
    print SCRIPT $string;
    close(SCRIPT);
    my $submit = qq"qsub split_align_submission.sh";
    open(CMD, "$submit |");
    while(my $line = <CMD>) {
        chomp $line;
        print "The submitted job is: $line\n";
        $job = $line;
    }
    close(CMD);
}

1;
