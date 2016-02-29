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

=item C<Parse_Fasta>

    hmmm

=cut
sub Parse_Fasta {
    my $me = shift;
    my %args = @_;
    my $best = 1;
    $best = $args{best} if (defined($args{best}));
    my $sig = 0.0001;
    $sig = $args{sig} if (defined($args{sig}));
    my $input = $me->{input};
    $input = $args{input} unless ($input);
    my $output = $me->{output};
    $output = $args{output} unless ($output);
    unless ($output) {
        $output = $input;
        $output =~ s/\..+/\.parsed/g;
        print "TESTME: $output\n";
    }
    my $out = new FileHandle;
    $out->open(">$output");
    ## This works, other attempts I've used for zipped input fail.
    my $f = new FileHandle;
    $f->open("less ${input} |");
    ## This works, other attempts at using zipped input fail.
    my $searchio = new Bio::SearchIO(-format => 'fasta', -fh => $f, -best => ${best}, -signif => ${sig});
    print $out "QueryName\tQueryLength\tLibID\tLibHitLength\tScore\tE\tIdentity\tHitMatches\tQueryStart\tQueryEnd\tLibStart\tLibEnd\tLStrand\n";
    my $results = 0;
    while (my $result = $searchio->next_result()) {
        $results++;
        while(my $hit = $result->next_hit) {
            my $query_name = $result->query_name();
            my $query_length = $result->query_length();
            my $accession = $hit->accession();
            my $library_id = $hit->name();
            my $length = $hit->length();
            my $score = $hit->raw_score();
            my $sig = $hit->significance();
            my $ident = $hit->frac_identical();
            my $lstrand = $hit->strand('hit');

            my $hit_len;
            my $hit_matches;
            my $first_qstart = -1;
            my $first_qend = -1;
            my $first_lstart = -1;
            my $first_lend = -1;
            my $hsp_count = 0;
            while (my $hsp = $hit->next_hsp) {
                $hsp_count++;
                if ($hsp_count == 1) {
                    $first_qstart = $hsp->start('query');
                    $first_qend = $hsp->end('query');
                    $first_lstart = $hsp->start('subject');
                    $first_lend = $hsp->end('subject');
                }
                $hit_len = $hsp->length('total');
                my @matches = $hsp->matches(-seq => 'hit');
                $hit_matches = $matches[1];
            }
            ##print $out "QueryName\tQueryLength\tLibID\tLibHitLength\tScore\tE\tIdentity\tHitMatches\tQueryStart\tQueryEnd\tLibStart\tLibEnd\tQStrand\tLStrand\n";
            print $out "$query_name\t$query_length\t$library_id\t$length\t$score\t$sig\t$ident\t$hit_matches\t$first_qstart\t$first_qend\t$first_lstart\t$first_lend\t$lstrand\n";
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


sub Split_Align_Fasta {
    my $me = shift;
    my %args = @_;
    $me->Check_Options(["query","library"]);
    $args{query} = $me->{query} if ($me->{query});
    $args{library} = $me->{library} if ($me->{library});
    ## $args{param} = ' -e 10 ' unless (defined($args{param}));
    $args{number} = 40 unless (defined($args{number}));
    $args{parse} = 0 unless (defined($args{parse}));
    $args{num_dirs} = 0 unless (defined($args{num_dirs}));
    $args{best_only} = 0 unless (defined($args{best_only}));
    $args{queue} = 'workstation' unless (defined($args{queue}));
    $args{num_sequences} = 0 unless (defined($args{num_sequences}));
    ## $args{library} = basename($args{library}, ('.fasta'));
    my $lib = basename($args{library}, ('.fasta'));
    my $que = basename($args{query}, ('.fasta'));
    my $outdir = qq"$me->{basedir}/outputs/fasta";
    make_path("${outdir}") unless(-d ${outdir});
    my $output = qq"${outdir}/${que}_vs_${lib}.txt.gz";
    my $concat_job;
    if ($me->{pbs}) {
        my $num_per_split = $me->Get_Split(%args);
        $args{num_per_split} = $num_per_split;
        print "Going to make $args{number} directories with $num_per_split sequences each.\n";
        my $actual = $me->Make_Directories(%args);
        print "Actually used ${actual} directories to write files.\n";
        my $alignment = $me->Make_Fasta_Job(number => $actual, split => 1);
        $concat_job = $me->Concatenate_Searches(depends => $alignment->{pbs_id}, output => ${output},);
    } else {
        ## If we don't have pbs, force the number of jobs to 1.
        $args{number} = 1;
        my $num_per_split = $me->Get_Split(%args);
        my $actual = $me->Make_Directories(%args);
        my $alignment = $me->Make_Fasta_Job(number => $actual);
        $concat_job = $me->Concatenate_Searches(output=> ${output},);
    }

    my $parse_input = $concat_job->{output};
    my $comment_string = qq!## I don't know if this will work.!;
    my $job_string = qq?
use CYOA;
my \$h = new CYOA;
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

sub Make_Fasta_Job {
    my $me = shift;
    my %args = @_;
    my $dep = '';
    $dep = $args{depends} if (defined($args{depends}));
    my $split = 0;
    $split = $args{split} if ($args{split});
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
    if ($split) {
        $job_string = qq!
cd $me->{basedir}
$me->{fasta_tool} $me->{fasta_args} -T $me->{cpus} \\
 $me->{basedir}/split/\${PBS_ARRAYID}/in.fasta \\
 $me->{library} \\
 1>$me->{basedir}/outputs/\${PBS_ARRAYID}.out \\
 2>>$me->{basedir}/split_align_errors.txt
!;
    } else {
        $job_string = qq!
cd $me->{basedir}
  $me->{fasta_tool} $me->{fasta_args} -T $me->{cpus} \\
  $me->{query} \\
  $me->{library} \\
  1>$me->{basedir}/outputs/$me->{fasta_tool}.out \\
  2>>$me->{basedir}/split_align_errors.txt
!;
    }
    my $comment = qq!## Running multiple fasta jobs.!;
    my $fasta_jobs = $me->Qsub(job_name => 'fasta_multi',
                               depends => $dep,
                               depends_type => 'array',
                               job_string => $job_string,
                               comment => $comment,
                               qsub_queue => 'long',
                               qsub_wall => '96:00:00',
                               qsub_args => " $me->{qsub_args} -t ${array_string} ",
        );
    return($fasta_jobs);
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
    $in->open("less $args{input} |");
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

1;
