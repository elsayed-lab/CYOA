package Bio::Adventure::Align_Blast;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';
use Bio::Adventure::Torque;

use Bio::SearchIO::blast;
use Bio::SearchIO::fasta;
use Bio::Seq;
use Bio::Tools::Run::StandAloneBlast;
use Cwd;
use File::Basename;
use File::Path qw"make_path";
use File::Which qw"which";
use POSIX qw"ceil";

my $out_fmt = 5;

=head1 NAME

    Bio::Adventure::Alignment - Perform Sequence alignments with tools like Blast/Fasta

=head1 SYNOPSIS

    use Bio::Adventure;
    my $hpgl = new Bio::Adventure;
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
    my ($class, %args) = @_;
    my $check = which('blastp');
    die("Could not find blast in your PATH.") unless($check);
    my $options = $class->Get_Vars(args => \%args);
    print STDERR qq"Please note that this function calls blastall
    separately for every sequence in the database.  As a result it is
    not recommended fo use with large sequence libraries.  Instead use
    the separate functions 'Run_Blast()' or 'Split_Align_Blast()'
    followed by 'Parse_Blast()' which does these steps separately.";
    my $blast_program = 'blastp';
    ## Search libraries are the various combinations of esmer, nonesmer, and unassigned which are of potential interest
    ## # TODO: hus all 7 of the possible states
    ## These libraries are generated in make_blast.sh
    my @search_libraries = ('nr',);
    #my @search_libraries = ('esmer',
    ##                       'non',
    ##                       'unas',
    ##                       'esmer-non',
    ##                       'esmer-unas',
    ##                       'non-unas',
    ##                       'all',);
    ## States to save: singles, doubles, few (3-10), many (10+)
    my $blast_output = Bio::SearchIO->new(-format => 'blast',);
    my $query = $options->{query};
    my $number_hits = 0;
    for my $library (@search_libraries) {
        my $query_library = Bio::SeqIO->new(-file => ${query}, -format => 'Fasta');
        my $output_directory = qq"$options->{basedir}/outputs/${query}_vs_${library}";
        make_path("${output_directory}") unless (-d $output_directory);
        my $counts = FileHandle->new(">${output_directory}/counts.txt");
        my $zeros = FileHandle->new(">${output_directory}/zeros.fasta");
        my $singles = FileHandle->new(">${output_directory}/singles.fasta");
        my $doubles = FileHandle->new(">${output_directory}/doubles.fasta");
        my $few = FileHandle->new(">${output_directory}/few.fasta");
        my $many = FileHandle->new(">${output_directory}/many.fasta");
        my $seq_count = 0;

        my @params = (
            -program => $blast_program,
            ## Show GI's in definition lines?
            -F => 'T', ## Filter query sequence
            -G => '-1', ## Cost to open a gap
            -E => '-1', ## Cost to extend a gap
            -X => '0',  ## X dropoff value for gapped alignment
            -I => 'T',  ## Default in blast is F
            ## -q => '-3', ## Penalty for nucleotide mismatch (blastn only)
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
            my $search = Bio::Tools::Run::StandAloneBlast->new(@params);
            my $blast_output;
            eval {
                ## $blast_output = $search->blastall($query_library);
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
                    $number_hits = $number_hits++;
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
                    print $singles $entry;
                } elsif ($hit_count == 2) {
                    print $doubles $entry;
                } elsif ($hit_count >= 3 and $hit_count <= 10) {
                    print $few $entry;
                } elsif ($hit_count > 10) {
                    print $many $entry;
                } else {
                    print $zeros $entry;
                }
                print $counts "$query_name\t$hit_count\n";
            }    ## End of the individual blast search  (maybe not needed?)
        }        ## End of this search library  (the 7 states to search against)
        $counts->close();
        $zeros->close();
        $singles->close();
        $doubles->close();
        $few->close();
        $many->close();
    }
    return($number_hits);
}

=item C<Split_Align>

    $hpgl->Split_Align();
    Split apart a set of query sequences into $args{number} pieces and align them all separately.

=cut
sub Split_Align_Blast {
    my ($class, %args) = @_;
    my $check = which('blastp');
    die("Could not find blast in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ["query","library", "blast_tool"],
        param => ' -e 10 ',
        blast_tool => 'blastn',
        number => 40,
        parse => 0,
        num_dirs => 0,
        best_only => 0,
        interactive => 0,
        num_sequences => 0,
    );
    if ($options->{interactive}) {
        print "Run this with: cyoa --task align --method blastsplit --query something.fasta --library blastdb --blast_tool tblastx\n";
    }

    ## This might be wrong (tblastx)
    print STDERR qq"A quick reminder because I (atb) get confused easily:
tblastn is a protein fasta query against a nucleotide blast database.
tblastx is a nucleotide fasta query (which is translated on the fly) against a nucleotide blast db.
blastx is a nucleotide fasta query against a protein blast db.
blastn is normal nucleotide/nucleotide
blastp is normal protein/protein.
";
    sleep(1);
    ## Also set the default blastdb if it didn't get set.
    if ($options->{blast_tool} eq 'blastn' or
            $options->{blast_tool} eq 'tblastn' or
            $options->{blast_tool} eq 'tblastx') {
        $options->{peptide} = 'F';
        $options->{library} = 'nt' if (!defined($options->{library}));
    } else {
        $options->{peptide} = 'T';
        $options->{library} = 'nr' if (!defined($options->{library}));
    }
    my $lib = basename($options->{library}, ('.fasta'));
    my $que = basename($options->{query}, ('.fasta'));
    my $outdir = qq"$options->{basedir}/outputs";
    make_path("${outdir}") unless(-d ${outdir});
    my $output = qq"${outdir}/${que}_vs_${lib}.txt";
    my $concat_job;
    $lib = Bio::Adventure::Align_Blast::Check_Blastdb($class, %args);
    if ($options->{pbs}) {
        my $num_per_split = Bio::Adventure::Align::Get_Split($class, %args);
        $options = $class->Set_Vars(num_per_split => $num_per_split);
        print "Going to make $options->{number} directories with $num_per_split sequences each.\n";
        my $actual = Bio::Adventure::Align::Make_Directories($class, %args);
        print "Actually used ${actual} directories to write files.\n";
        my $alignment = Bio::Adventure::Align_Blast::Make_Blast_Job(
            $class,
            library => $lib,
            number => $actual,
            output_type => ${out_fmt},
        );
        $concat_job = Bio::Adventure::Align::Concatenate_Searches(
            $class,
            depends => $alignment->{pbs_id},
            output => ${output},
        );
    } else {
        ## If we don't have pbs, force the number of jobs to 1.
        print "Not using the cluster.\n";
        $options = $class->Set_Vars(number => 1);
        my $num_per_split = Bio::Adventure::Align::Get_Split($class);
        $options = $class->Set_Vars(num_per_split => $num_per_split);
        print "TESTME: Going to write $num_per_split sequences in 1 file.\n";
        my $actual = Bio::Adventure::Align::Make_Directories($class);
        my $alignment = Bio::Adventure::Align_Blast::Make_Blast_Job(
            $class,
            library => $lib,
            number => $actual,
            output_type => ${out_fmt},
        );
        $concat_job = Bio::Adventure::Align::Concatenate_Searches(
            $class,
            output => ${output},
        );
    }

    my $parse_input = cwd() . qq"/$concat_job->{output}";
    my $comment_string = qq!## I don't know if this will work.!;
    my $jstring = qq?
use Bio::Adventure;
use Bio::Adventure::Align;
use Bio::Adventure::Align_Blast;
my \$h = new Bio::Adventure;
Bio::Adventure::Align::Parse_Search(\$h, input => '$parse_input', search_type => 'blastxml', best => $args{best_only});
?;
    my $parse_job = $class->Submit(
        comment => $comment_string,
        depends => $concat_job->{pbs_id},
        jname => "parse_search",
        jstring => $jstring,
        language => 'perl',

    );
    return($concat_job);
}

sub Make_Blast_Job {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        depends => '',
    );
    my $dep = $options->{depends};
    my $library = $options->{library};
    my $array_end = 1000 + $options->{number};
    my $array_string = qq"1000-${array_end}";

    ## Handle array job types for slurm/torque.
    my $queue_array_string = 'SLURM_ARRAY_TASK_ID';
    if ($options->{pbs} eq 'torque') {
        $queue_array_string = 'PBS_ARRAYID';
    }

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
    my $jstring = '';
    ## I have been getting null pointers, there is a bug report suggesting this is because of xml output
    ## $args{output_type} = 1;
    $args{output_type} = ${out_fmt};
    if ($options->{pbs}) {
        $jstring = qq!
cd $options->{basedir}
$options->{blast_tool} -outfmt $options->{output_type} \\
 -query $options->{basedir}/split/\${${queue_array_string}}/in.fasta \\
 -db ${library} \\
 -out $options->{basedir}/outputs/\${${queue_array_string}}.out \\
 1>$options->{basedir}/outputs/\${${queue_array_string}}.stdout \\
 2>>$options->{basedir}/split_align_errors.txt
!;
    } else {
        $jstring = qq!
cd $options->{basedir}
$options->{blast_tool} -outfmt $options->{output_type} \\
  -query $options->{query} \\
  -db ${library} \\
  -out $options->{basedir}/outputs/$options->{blast_tool}.out \\
  1>$options->{basedir}/outputs/$options->{blast_tool}.stdout \\
  2>>$options->{basedir}/split_align_errors.txt
!;
    }
    my $comment = qq!## Running multiple blast jobs.!;
    my $blast_jobs = $class->Submit(
        comment => $comment,
        jname => 'blast_multi',
        depends => $dep,
        jstring => $jstring,
        mem => 32,
        array_string => $array_string,
    );
    return($blast_jobs);
}

sub Merge_Parse_Blast {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['output'],
    );
    if ($args{interactive}) {
        print "Run with: cyoa --task align --method mergeparse --output filename\n";
    }
    my $concat = Bio::Adventure::Align::Concatenate_Searches($class);
    my $input = $options->{output};
    my $jstring = qq?
use Bio::Adventure;
use Bio::Adventure::Align;
use Bio::Adventure::Align_Blast;
my \$h = Bio::Adventure->new(input => \$input);
my \$final = Bio::Adventure::Align_Blast->Parse_Search(\$h, search_type => 'blastxml',);
?;
    my $parse = $class->Submit(
        depends => $concat->{pbs_id},
        jname => 'parse_search',
        jstring => $jstring,
        language => 'perl',
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
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        best_only => 0,
        search_type => 'blastxml',
        num_sequences => 'unknown',
    );
    if ($args{interactive}) {
        print "Run with: cyoa --task align --method parseblast --intput filename.\n";
    }
    my $input = $options->{input};
    my $output = qq"${input}";
    my $best = $options->{best_only};
    my $search_type = $options->{search_type};
    $output =~ s/\.txt\.gz//g;
    $output .= "_parsed.txt";
    print "Writing parsed output to ${output}\n";
    my $count_table = qq"${output}";
    $count_table =~ s/_parsed\.txt/_counts\.txt/g;
    print "Writing count table to ${count_table}\n";
    my $counts = FileHandle->new(">$count_table");
    $counts->autoflush(1);
    print $counts "Query\tquery_accession\thits\n";
    my $fh = FileHandle->new(">$output");
    $fh->autoflush(1);
    local $| = 0;               ## I want to watch the output as it happens.
    my $in = FileHandle->new("less ${input} |");
    $XML::SAX::ParserPackage = 'XML::SAX::PurePerl';
    my $searchio = Bio::SearchIO->new(-fh => $in);
    ## Perform a test for the output format, this will close the filehandle, so reopen it after
    my $guessed_format = $searchio->format();
    print STDERR "Guessed the format is: $guessed_format\n";
    ##$guessed_format = 'blastxml';
    ##print STDERR "It turns out that searchio->format() is crap, setting format to 'blastxml'\n";
    $in->close();
    $in = FileHandle->new("less ${input} |");

    $searchio = Bio::SearchIO->new(-fh => $in, -format => $search_type, -best => ${best},);
    ##   if ($args{output_type} eq '0') { ## use old blast format
    ##     $searchio = new Bio::SearchIO(-format => 'blast', -fh => $in, -best => 1,);
    ##   } elsif ($args{output_type} eq '7') {  ## use blastxml
    ##     $searchio = new Bio::SearchIO(-format => 'blastxml', -fh => $in, -best => 1,);
    ##   } else { ## I don't know what it is, make searchio guess it
    ##   $searchio = new Bio::SearchIO(-fh => $in, -best => 1,);
    ##   my $test_format = $searchio->format();
    ##   print "I am guessing the format is: $test_format\n";
    ##   $searchio = undef;
    ##   $searchio = new Bio::SearchIO(-format => $test_format, -fh => $in, -best => 1,);
    ## }
    print $fh "QUERYNAME\treal_name\tChromosome\tStart\tEnd\t%ID\tScore\tSig\tCompLength\tHit_Ident\thits\n";
    my $result_count = 0;
  RESULT: while(my $result = $searchio->next_result()) {
        my $hit_count = 0;
        $result_count++;
        print STDERR "Parsed $result_count of $options->{num_sequences} results.\n";
        my $query_name = "";
        my $real_name = "";
        while (my $hit = $result->next_hit) {
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

                ## $hsp_cons = $hsp->frac_cons('total');
                last HSP;      ## The 'HSP: ' is a label given to a logical loop
                ## If you then say next HSP; it will immediately stop processing the loop and move to the next iteration of the loop
                ## If instead you say last HSP; it will immediately break out of the loop.
            }   ## next hsp

            ## This prints %identity with respect to the overlap of the component
            ## Then the score of the hit, its Evalue, and the component length.
            my $print_string = qq"${query_name}\t${real_name}\t${acc2}\t${acc3}\t${start}\t${end}\t${ident}\t${score}\t${sig}\t${query_length}\t${hsp_identical}\t$hit_count\n";
            print $fh $print_string;
            if ($args{best_only}) {
                next RESULT;
            }
        }                       ## End each hit for a single result
        print $counts "${query_name}\t${real_name}\t${hit_count}\n";
        ## print "$result had $hit_count hits\n";
    }   ## Finish looking at each sequence
    $fh->close();
    $in->close();
    $counts->close();
    return($result_count);
}

=item C<Check_Blastdb>

    Check_Blastdb makes sure that there is an appropriately formatted
    blastdb for the library and type of blast performed.

=cut
sub Check_Blastdb {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['library'],
    );
    my $formatdb_ret = "not run.";
    my $libname = $options->{library};
    ## First check for the relevant library in $ENV{BLASTDB}
    ## If it isn't there, make one in basedir/blastdb/
    my $foundlib = 0;
    my $checklib = "";
    my $checklib_zero = "";
    if ($options->{blast_tool} =~ m/blast/ and $options->{peptide} eq 'T') {
        $checklib_zero = qq"$options->{library}.00.psq";
        $checklib = qq"$options->{library}.psq";
    } else {
        $checklib_zero = qq"$options->{library}.00.nsq";
        $checklib = qq"$options->{library}.nsq";
    }
    my $lib = "";
    if (!defined($ENV{BLASTDB})) {
        $ENV{BLASTDB} = "$options->{basedir}/blastdb";
    }
    print "Looking for ${checklib} / ${checklib_zero} in either $ENV{BLASTDB}.\n";
    ## Start with BLASTDB
    if (-f "$ENV{BLASTDB}/${checklib}" or -f "$ENV{BLASTDB}/${checklib_zero}") {
        $foundlib++;
        $libname = qq"$ENV{BLASTDB}/$options->{library}";
        $lib = qq"$ENV{BLASTDB}/$options->{library}";
        print "Found an existing blast database at $libname.\n";
    } elsif (-f "$options->{basedir}/blastdb/${checklib}" or -f "$options->{basedir}/blastdb/${checklib_zero}") { ## Then try basedir
        $foundlib++;
        $libname = qq"$options->{basedir}/blastdb/$options->{library}";
        $lib = qq"$options->{basedir}/blastdb/$options->{library}";
        print "Found an existing blast database at $libname.\n";
    } else {
        print "Did not find an existing blast database.\n";
    }

    if (!$foundlib) {
        if (!-d qq"$options->{basedir}/blastdb") {
            make_path(qq"$options->{basedir}/blastdb");
        }
        $libname = qq"$options->{basedir}/blastdb/$options->{library}";
        my $formatdb_command = qq"formatdb -p $options->{peptide} -o T -n blastdb/$options->{library} -s -i $options->{library}";
        print "The formatdb command run is: ${formatdb_command}\n";
        $formatdb_ret = qx"${formatdb_command}";
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
