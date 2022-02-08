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
use Cwd 'abs_path';
use File::Basename qw"dirname basename";
use File::Path qw"make_path remove_tree";
use File::Which qw"which";
use POSIX qw"ceil";

=head1 NAME

 Bio::Adventure::Align_Fasta - Perform Sequence alignments with the fasta suite of tools.

=head1 SYNOPSIS

 All the functions which live here work with the fasta36 suite of
 programs.  They write the scripts, invoke them on the cluster,
 collect the results, and (optionally) parse them into simplified
 tables.

=head1 METHODS

=head2 C<Make_Fasta_Job>

 Write a job for the cluster to run a fasta36 job.

 This handles the creation of the files and directories required when
 splitting up a sequence into a bunch of pieces for a split-fasta search.

=item C<Arguments>

 split(0): Split up the fasta jobs into multiple pieces?
 fasta_tool('ggsearch36'): Use this fasta36 tool.
 output_type(undef): Specify a particular output type, ideally one of
  the parseable formats.
 jdepends(): Put this in a dependency chain?
 jmem(8): Expected memory usage.
 modules('fasta'): Load this environment module.

=cut
sub Make_Fasta_Job {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        align_jobs => 1,
        fasta_tool => 'ggsearch36',
        split => 0,
        output_type => undef,
        jdepends => '',
        jmem => 8,
        modules => ['fasta'],);
    my $dep = $options->{jdepends};
    my $split = $options->{split};
    my $output_type = $options->{output_type};
    my $library;
    my $array_end = 1000 + $args{align_jobs};
    my $array_string = qq"1000-${array_end}";
    my $jstring = '';
    my $type_string = '';
    if (defined($output_type)) {
        $type_string = "-m ${output_type}";
    }
    my $output = '';
    if ($split) {
        $output = qq"$options->{basedir}/outputs/\${PBS_ARRAYID}.stdout";
        $jstring = qq!
cd $options->{basedir}
$options->{fasta_tool} $options->{fasta_args} ${type_string} -T $options->{cpus} \\
 $options->{workdir}/split/\${PBS_ARRAYID}/in.fasta \\
 $options->{library} \\
 1>${output} \\
 2>>$options->{basedir}/split_align.stderr
!;
    } else {
        $output = qq"$options->{workdir}/$options->{fasta_tool}.stdout";
        $jstring = qq!
cd $options->{basedir}
  $options->{fasta_tool} $options->{fasta_args} ${type_string} -T $options->{cpus} \\
  $options->{input} \\
  $options->{library} \\
  1>${output} \\
  2>>$options->{basedir}/split_align.stderr
!;
    }
    my $comment = '## Running $options->{align_jobs} fasta job(s).';

    my $fasta_jobs = $class->Submit(
        comment => $comment,
        depends_type => 'array',
        jdepends => $dep,
        jmem => $options->{jmem},
        jname => 'fasta_multi',
        jstring => $jstring,
        jprefix => "91",
        modules => $options->{modules},
        qsub_args => " $options->{qsub_args} -t ${array_string} ",
        output => $output,
        jqueue => 'long',
        jwalltime => '96:00:00',);
    return($fasta_jobs);
}

=head2 C<Parse_Fasta>

 Parse a fasta36 tool result file and print a summary table.

 Given the output from one of the fasta36 programs: ggsearch36,
 fasta36, etc.  This parses it and prints a simplified table of the
 results.

=item C<Arguments>

 input(required): File containing a completed fasta search output.
 best_only(0): Print only the best hit for each search sequence.
 evalue(0.0001): Filter for only hits with a better evalue than this.

=item C<Invocation>

> cyoa --task align --method parsefasta --input fasta_output.txt.xz

=cut
sub Parse_Fasta {
    my ($class, %args) = @_;
    my $check = which('fasta36');
    die("Could not find fasta in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        best_only => 0,
        evalue => 0.0001,);
    my $best_only = $options->{best_only};
    my $evalue = $options->{evalue};
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
    my $searchio = Bio::SearchIO->new(-format => 'fasta',
                                      -fh => $f,
                                      -best => ${best_only},
                                      -signif => ${evalue});
    print $out "QueryName\tQueryLength\tLibID\tLibHitLength\tScore\tE\tIdentity\tHitMatches\tQueryStart\tQueryEnd\tLibStart\tLibEnd\tLStrand\n";
    my $results = 0;
    while (my $result = $searchio->next_result()) {
      $results++;
      my $hit_count = 0;
    HIT: while (my $hit = $result->next_hit) {
        $hit_count++;
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
      HSP: while (my $hsp = $hit->next_hsp) {
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
        print $out "${query_name}\t${query_length}\t${library_id}\t${length}\t${score}\t${sig}\t${ident}\t${hit_matches}\t${first_qstart}\t${first_qend}\t${first_lstart}\t${first_lend}\t${lstrand}\n";
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

=head2 C<Parse_Fasta_Global>

 A separate parser for global:global searches.

 This was originally written in order to write summary tables of
 alignments of two closely related species (like L.major vs. T.cruzi)
 in order to categorize the 1:1 orthologs vs. gene families.

 Parse_Global makes a hard assumption about the structure of the hits
 in the output of the alignment program.  They should be in a
 global/global search, thus we assume 1 hit / 1 result.

 In addition, this will print summaries of hits in a few different
 files depending on how many hits were observed for each query
 sequence: singletons, doubles, triples, few (3-10), and many (10+).

=item C<Arguments>

 best_only(0): Parse out only the best hit?
 fasta_format(0): Format of the search input
 evalue(0.001): Quality of hit filter
 many_cutoff(10): How many hits define 'many' hits?
 min_score(0): Provide a minimum score required for a hit?
 check_all_hits(0): I am not sure what this does.
 min_percent(0): Add a percentage identity filter

=cut
sub Parse_Fasta_Global {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        best_only => 0,
        fasta_format => 'fasta',
        evalue => 0.001,
        many_cutoff => 10,
        min_score => undef,
        check_all_hits => 0,
        min_percent => undef,
        output => '',
        output_counts => '',
        output_singles => '',
        output_doubles => '',
        output_few => '',
        output_many => '',
        output_zero => '',
        output_all => '',);
    my $output = $options->{input};
    my $outdir = dirname($output);
    $output = basename($output, ('.gz', '.xz'));
    $output = basename($output, ('.txt'));

    my $input = $options->{input};
    my $in = FileHandle->new("less ${input} |");
    my $searchio = Bio::SearchIO->new(
        -format => $options->{fasta_format},
        -fh => $in,
        -best_hit_only => $options->{best_only},
        -max_significance => $options->{evalue},
        -check_all_hits => $options->{check_all_hits},
        -min_score => $options->{min_score},);

    my $counts = FileHandle->new(qq">$options->{output_counts}");
    my $parsed = FileHandle->new(qq">$options->{output}");
    my $singles = FileHandle->new(qq">$options->{output_singles}");
    my $doubles = FileHandle->new(qq">$options->{output_doubles}");
    my $few = FileHandle->new(qq">$options->{output_few}");
    my $many = FileHandle->new(qq">$options->{output_many}");
    my $zero = FileHandle->new(qq">$options->{output_zero}");
    my $all = FileHandle->new(qq">$options->{output_all}");

    my $seq_count = 0;
    print $parsed "Query Name\tQuery length\tHit ID\tHit Length\tScore\tE\tIdentity\tHit length\tHit Matches\n";
    while (my $result = $searchio->next_result()) {
        $seq_count++;
        my $query_name = $result->query_name();
        my $entry = qq"${query_name}\t";
        my $hit_count = 0;
      HITLOOP: while(my $hit = $result->next_hit) {
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
          print $parsed "${query_name}\t${query_length}\t${acc2}\t${length}\t${score}\t${sig}\t${ident}\t${hit_len}\t${hit_matches}\n";
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

=head2 C<Split_Align_Fasta>

 Split apart a set of query sequences into $args{align_jobs} pieces and align them all separately.

 This is the main callable function when one wants to perform a bunch
 of fasta36 searches in parallel.

=item C<Arguments>

 input(required): File containing sequences to search _for_.
 library(required): File containing sequences to search _from_.
 align_jobs(40): How many jobs should be spawned?
 align_parse(0): Start up a parsing job upon completion?
 best_only(0): Passed to the parser, only print the best hit?

=item C<Invocation>

> cyoa --task fastasplit --input query.fasta --library library.fasta --best_only 1

=cut
sub Split_Align_Fasta {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        fasta_tool => 'ggsearch36',
        required => ['input', 'library'],
        interactive => 0,
        align_jobs => 40,
        align_parse => 0,
        num_dirs => 0,
        best_only => 0,
        jmem => 8,
        modules => ['fasta']);
    my $lib = basename($options->{library}, ('.fasta'));
    my $que = basename($options->{input}, ('.fasta'));
    my $outdir = qq"$options->{basedir}/outputs/fasta_${que}_${lib}";
    make_path("${outdir}") unless(-d ${outdir});
    my $output = $options->{input};
    $output = basename($output, ('.gz', '.xz'));
    $output = basename($output, ('.txt', '.fasta'));
    my $output_file = qq"${outdir}/${output}.parsed.txt";
    my $counts = qq"${outdir}/${output}.count";
    my $singles = qq"${outdir}/${output}_singles.txt";
    my $doubles = qq"${outdir}/${output}_doubles.txt";
    my $few = qq"${outdir}/${output}_few.txt";
    my $many = qq"${outdir}/${output}_many.txt";
    my $zero = qq"${outdir}/${output}_zero.txt";
    my $all = qq"${outdir}/${output}_all.txt";
    my ($concat_job, $alignment);
    if ($options->{pbs}) {
        my $num_per_split = $class->Bio::Adventure::Align::Get_Split(%args);
        $options = $class->Set_Vars(num_per_split => $num_per_split);
        $options = $class->Set_Vars(workdir => $outdir);
        print "Going to make $options->{align_jobs} directories with ${num_per_split} sequences each.\n";
        my $actual = $class->Bio::Adventure::Align::Make_Directories(
            %args,
            workdir => $outdir);
        print "Actually used ${actual} directories to write files.\n";
        $alignment = $class->Bio::Adventure::Align_Fasta::Make_Fasta_Job(
            %args,
            align_jobs => $actual,
            workdir => $outdir,
            split => 1);
        $concat_job = $class->Bio::Adventure::Align::Concatenate_Searches(
            jdepends => $alignment->{job_id},
            output => $output,
            workdir => $outdir,
            jprefix => '92',);
        ## Make sure that the alignment job gets the final output
        $alignment->{output} = $concat_job->{output};
    } else {
        ## If we don't have a queue, force the number of jobs to 1.
        $options->{align_jobs} = 1;
        my $num_per_split = $class->Bio::Adventure::Align::Get_Split(%args);
        $options = $class->Set_Vars(num_per_split => $num_per_split);
        $options = $class->Set_Vars(workdir => $outdir);
        $alignment = $class->Bio::Adventure::Align_Fasta::Make_Fasta_Job(
            %args,
            workdir => $outdir,
            align_jobs => 1);
        print "Finished alignment with output: $alignment->{output}\n";
    }

    my $parse_input = $alignment->{output};
    my $comment_string = qq!## I don't know if this will work.!;
    my $jstring = qq?
use Bio::Adventure;
use Bio::Adventure::Align;
use Bio::Adventure::Align_Fasta;
my \$result = \$h->Bio::Adventure::Align_Fasta::Parse_Fasta_Global(
  input => '$parse_input',
  output => '$output_file',
  output_counts => '$counts',
  output_singles => '$singles',
  output_doubles => '$doubles',
  output_few => '$few',
  output_many => '$many',
  output_zero => '$zero',
  output_all => '$all',
  search_type => 'global_fasta',);
?;
    my $parse_job = $class->Submit(
        comment => $comment_string,
        input => $parse_input,
        output => $output_file,
        output_counts => $counts,
        output_singles => $singles,
        output_doubles => $doubles,
        output_few => $few,
        output_many => $many,
        output_zero => $zero,
        output_all => $all,
        jdepends => $concat_job->{job_id},
        jmem => $options->{jmem},
        jname => 'parse_search',
        jstring => $jstring,
        jprefix => '93',
        language => 'perl',);
    $parse_job->{align} = $alignment;
    $parse_job->{concat} = $concat_job;
    return($parse_job);
}

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

L<fasta36>

=cut

1;
