package Bio::Adventure::Align_Fasta;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use Bio::SearchIO::blast;
use Bio::SearchIO::fasta;
use Bio::Seq;
use Bio::Tools::Run::StandAloneBlast;
use Cwd;
use File::Basename qw"dirname basename";
use File::Path qw"make_path remove_tree";
use File::Which qw"which";
use POSIX qw"ceil";

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

=item C<Parse_Fasta>

    hmmm

=cut
sub Parse_Fasta {
    my ($class, %args) = @_;
    my $check = which('fasta36');
    die("Could not find fasta in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        best => 1,
        sig => 0.0001,
    );
        if ($args{interactive}) {
            print "Run this with: cyoa --task align --method fastaparse --input filename.\n";
        }
    my $best = $options->{best};
    my $sig = $options->{sig};
    my $input = $options->{input};
    my $output = $options->{output};
    unless ($output) {
        $output = $input;
        $output =~ s/\..+//g;
        $output .= "_parsed.txt";
    }
    my $out = FileHandle->new(">${output}");
    ## This works, other attempts I've used for zipped input fail.
    my $f = FileHandle->new("less ${input} |");
    ## This works, other attempts at using zipped input fail.
    my $searchio = Bio::SearchIO->new(-format => 'fasta', -fh => $f,
                                      -best => ${best}, -signif => ${sig});
    print $out "QueryName\tQueryLength\tLibID\tLibHitLength\tScore\tE\tIdentity\tHitMatches\tQueryStart\tQueryEnd\tLibStart\tLibEnd\tLStrand\n";
    my $results = 0;
    while (my $result = $searchio->next_result()) {
        $results++;
        while (my $hit = $result->next_hit) {
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
            ## process the Bio::Search::Hit::HitI object
            ## while(my $hsp = $hit->next_hsp) {
            ## process the Bio::Search::HSP::HSPI object
            ## }
        }   ## End of each hit of a result.
    }       ## End of each result.
    $f->close();
    $out->close();
    return($results);
}

sub Split_Align_Fasta {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['query', 'library'],
        interactive => 0,
        number => 40,
        parse => 0,
        num_dirs => 0,
        best_only => 0,
        num_sequences => 0,
    );
    if ($options->{interactive}) {
        print "Run this with: cyoa --task align --method fastasplit --query filename.fasta --library another.fasta.\n";
    }

    ## $args{library} = basename($options->{library}, ('.fasta'));
    my $lib = basename($options->{library}, ('.fasta'));
    my $que = basename($options->{query}, ('.fasta'));
    my $outdir = qq"$options->{basedir}/outputs/fasta";
    make_path("${outdir}") unless(-d ${outdir});
    my $output = qq"${outdir}/${que}_vs_${lib}.txt.gz";
    my $concat_job;
    if ($options->{pbs}) {
        my $num_per_split = Bio::Adventure::Align::Get_Split($class, %args);
        $options = $class->Set_Vars(num_per_split => $num_per_split);
        print "Going to make $options->{number} directories with ${num_per_split} sequences each.\n";
        my $actual = Bio::Adventure::Align::Make_Directories($class, %args);
        print "Actually used ${actual} directories to write files.\n";
        my $alignment = Bio::Adventure::Align_Fasta::Make_Fasta_Job($class, number => $actual, split => 1);
        $concat_job = Bio::Adventure::Align::Concatenate_Searches(
            $class,
            depends => $alignment->{pbs_id},
            output => $output,
        );
    } else {
        ## If we don't have pbs, force the number of jobs to 1.
        $options->{number} = 1;
        my $num_per_split = Bio::Adventure::Align::Get_Split($class, %args);
        $options = $class->Set_Vars(num_per_split => $num_per_split);
        my $actual = Bio::Adventure::Align::Make_Directories($class, %args);
        my $alignment = Bio::Adventure::Align_Fasta::Make_Fasta_Job($class, number => $actual);
        $concat_job = Bio::Adventure::Align::Concatenate_Searches($class, output=> $output,);
    }

    my $parse_input = $output;
    my $comment_string = qq!## I don't know if this will work.!;
    my $jstring = qq?
use Bio::Adventure;
use Bio::Adventure::Align;
use Bio::Adventure::Align_Fasta;
Bio::Adventure::Align::Parse_Search(\$h, input => '$parse_input', search_type => 'global_fasta',);
?;
    my $parse_job = $class->Submit(
        comment => $comment_string,
        depends => $concat_job->{pbs_id},
        jname => "parse_search",
        jstring => $jstring,
        language => 'perl',
        shell => '/usr/bin/env perl',
    );
    return($concat_job);
}

sub Make_Fasta_Job {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        depends => '',
        split => 0,
        output_type => 7,
    );
    my $dep = $options->{depends};
    my $split = $options->{split};
    my $library;
    my $array_end = 1000 + $args{number};
    my $array_string = qq"1000-${array_end}";
    if ($options->{parse}) {
        if ($options->{output_type} ne '7') {
            print "If one wants to parse the output, it probably should be in the xml format, which may be done by adding --output_type 7 to this script.\n";
            print "This script will still try to parse a blast plain text output, but if it fails, don't say I didn't warn you.\n";
            print "Sleeping for 10 seconds so you can hit Control-C if you want to rerun.\n";
            sleep(10);
        }
    }
    my $jstring = '';
    if ($split) {
        ## $options->{cpus} should be replaced with a slot from Torque->new().
        $jstring = qq!
cd $options->{basedir}
$options->{fasta_tool} $options->{fasta_args} -T $options->{cpus} \\
 $options->{basedir}/split/\${PBS_ARRAYID}/in.fasta \\
 $options->{library} \\
 1>$options->{basedir}/outputs/\${PBS_ARRAYID}.out \\
 2>>$options->{basedir}/split_align_errors.txt
!;
    } else {
        $jstring = qq!
cd $options->{basedir}
  $options->{fasta_tool} $options->{fasta_args} -T $options->{cpus} \\
  $options->{query} \\
  $options->{library} \\
  1>$options->{basedir}/outputs/$options->{fasta_tool}.out \\
  2>>$options->{basedir}/split_align_errors.txt
!;
    }
    my $comment = qq!## Running multiple fasta jobs.!;
    my $fasta_jobs = $class->Submit(
        comment => $comment,
        depends_type => 'array',
        depends => $dep,
        jname => 'fasta_multi',
        jstring => $jstring,
        qsub_args => " $options->{qsub_args} -t ${array_string} ",
        queue => 'long',
        walltime => '96:00:00',
    );
    return($fasta_jobs);
}

=item C<Parse_Global>

    Parse_Global makes a hard assumption about the structure of the hits in the output of the alignment program.  They should be
    in a global/global search, thus we assume 1 hit / 1 result.

=cut
sub Parse_Fasta_Global {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        best_hit_only => 0,
        format => 'fasta',
        max_significance => 0.001,
        many_cutoff => 10,
        min_score => undef,
        check_all_hits => 0,
        min_percent => undef,
    );
    my $output = $options->{input};
    my $outdir = dirname($output);
    $output = basename($output, ('.gz', '.xz'));
    $output = basename($output, ('.txt'));

    my $in = FileHandle->new("less $options->{input} |");
    my $searchio = Bio::SearchIO->new(-format => $options->{format},
                                      -fh => $in,
                                      -best_hit_only => $options->{best_hit_only},
                                      -max_significance => $options->{max_significance},
                                      -check_all_hits => $options->{check_all_hits},
                                      -min_score => $options->{min_score},);

    my $counts = FileHandle->new(qq">${outdir}/${output}.count");
    my $parsed = FileHandle->new(qq">${outdir}/${output}.parsed");
    my $singles = FileHandle->new(qq">${outdir}/${output}_singles.txt");
    my $doubles = FileHandle->new(qq">${outdir}/${output}_doubles.txt");
    my $few = FileHandle->new(qq">${outdir}/${output}_few.txt");
    my $many = FileHandle->new(qq">${outdir}/${output}_many.txt");
    my $zero = FileHandle->new(qq">${outdir}/${output}_zero.txt");
    my $all = FileHandle->new(qq">${outdir}/${output}_all.txt");

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
        }                       ## End iterating over every hit for a sequence.
        $entry .= "\n";
        print $all $entry;
        if ($hit_count == 1) {
            print $singles $entry;
        } elsif ($hit_count == 2) {
            print $doubles $entry;
        } elsif ($hit_count >= 3 and $hit_count <= $options->{many_cutoff}) {
            print $few $entry;
        } elsif ($hit_count > $options->{many_cutoff}) {
            print $many $entry;
        } else {
            print $zero $entry;
        }
        print $counts "${query_name}\t${hit_count}\n";
    }                           ## End iterating over every result in searchio
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
