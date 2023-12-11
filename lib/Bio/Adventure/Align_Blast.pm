package Bio::Adventure::Align_Blast;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';
use Bio::Adventure::Config;
use Bio::SearchIO::blast;
use Bio::SearchIO::fasta;
use Bio::Seq;
use Bio::Tools::Run::StandAloneBlast;
use Cwd;
use File::Basename;
use File::Path qw"make_path";
use File::Which qw"which";
use POSIX qw"ceil";

=head1 NAME

 Bio::Adventure::Align_Blast - Perform Sequence alignments with the blast suite of tools.

=head1 SYNOPSIS

 The functions in this file work with the blast family.  They write the scripts, invoke them on
 the cluster, collect the results, and (optionally) parse them into simplified tables.

=head1 METHODS

=head2 C<Make_Blast_Job>

 Write a single blast job for the cluster.

 This function is responsible for actually writing out an appropriate blast job.
 Hopefully, by the time this is called, all the appropriate options have been
 found.

 This is one of the older functions in CYOA and could really use some love.

=over

=item C<Arguments>

 blast_format(5): This needs to be one of the parseable formats, which seems
  to change over time.
 blast_tool('blastn'): Which blast tool to use
 library(nr): Choose a blast database.
 query(required): This ought to be changed to input.
 align_jobs(40): How many jobs to create.
 jdepends(''): What job does this depend upon?
 jmem(24): Expected memory usage.
 modules('blastdb', 'blast'): Load these environment modules.

=cut
sub Make_Blast_Job {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        blast_format => 5,
        jdepends => '',
        jmem => 24,);
    my $dep = $options->{jdepends};
    my $library = $options->{library};
    my $array_end = 100 + $options->{align_jobs};
    my $array_string = qq"100-${array_end}";
    my $blast_args = $class->Passthrough_Args(arbitrary => $options->{blast_args});

    ## Handle array job types for slurm/torque.
    my $queue_array_string = 'SLURM_ARRAY_TASK_ID';
    if ($options->{cluster} eq 'torque') {
        $queue_array_string = 'PBS_ARRAYID';
    }

    ## The fine folks at NCBI changed the numbers of the blast outputs!
    ## 0=pairwise ; 1=query_anchored identies ; 2=query_anchored no identities ;
    ## 3=flat identities ; 4=flat no identities ; 5=blastxml ; 6=tabular
    ## 7=tabular commented ; 8=text asn1 ; 9=binary asn1 ; 10=csv ; 11=blast archive asn1
    if ($options->{blast_format} ne '5') {
        print STDERR "If one wants to parse the output, it probably should be in the xml format, which may be done by adding --output_type 7 to this script.\n";
        print STDERR "This script will still try to parse a blast plain text output, but if it fails, don't say I didn't warn you.\n";
        print STDERR  "Sleeping for 3 seconds so you can hit Control-C if you want to rerun.\n";
        sleep(3);
    }
    my $jstring = '';
    if ($options->{cluster}) {
        $jstring = qq!
cd $options->{basedir}
export BLASTDB=$ENV{BLASTDB}
if [[ -f "outputs/split/\${$queue_array_string}/in.fasta" ]]; then
  $options->{blast_tool} -outfmt $options->{blast_format} \\
    -query $options->{basedir}/outputs/split/\${${queue_array_string}}/in.fasta \\
    -db ${library} ${blast_args} \\
    -out $options->{basedir}/outputs/split/\${${queue_array_string}}.out \\
    1>$options->{basedir}/outputs/\${${queue_array_string}}.stdout \\
    2>>$options->{basedir}/split_align.stderr
else
  echo "The input does not exist."
  exit 0
fi
!;
    } else {
        $jstring = qq!
cd $options->{basedir}
$options->{blast_tool} -outfmt $options->{blast_format} \\
  -query $options->{input} \\
  -db ${library} \\
  -out $options->{basedir}/outputs/$options->{blast_tool}.out \\
  1>$options->{basedir}/outputs/$options->{blast_tool}.stdout \\
  2>>$options->{basedir}/split_align.stderr
!;
    }
    my $comment = '## Running multiple blast jobs.';
    my $blast_jobs = $class->Submit(
        comment => $comment,
        jname => 'blast_multi',
        jdepends => $dep,
        jstring => $jstring,
        jmem => $options->{jmem},
        array_string => $array_string,);
    return($blast_jobs);
}

=back

=head2 C<Merge_Parse_Blast>

 Merge multiple blast outputs and parse them into a table.

 This invokes Concatenate_Searches() and Parse_Search() in order to get the final
 parsed table from a series of blast searches.

=over

=item C<Arguments>

 output(required): Output file into which to write the results.
 jmem(8): Expected memory.

=cut
sub Merge_Parse_Blast {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['output'],
        jmem => 8,);
    my $concat = $class->Bio::Adventure::Align::Concatenate_Searches();
    my $input = $options->{output};
    my $jstring = qq!
use Bio::Adventure;
use Bio::Adventure::Align;
use Bio::Adventure::Align_Blast;
\$h = Bio::Adventure->new(input => \$input);
my \$final = \$h->Bio::Adventure::Align_Blast->Parse_Search(search_type => 'blastxml',);
!;
    my $parse = $class->Submit(
        jdepends => $concat->{job_id},
        jmem => $options->{jmem},
        jname => 'parse_search',
        jstring => $jstring,
        language => 'perl',);
    return($parse);
}

=back

=head2 C<Parse_Blast>

 Do the actual parsing of a blast result and print an easy-to-read table.

 Parse_Blast is responsible for parsing blast output, it makes some attempt to be
 flexible vis a vis different formatting.  It prints a table of each gene and
 some relevant statistics about the hits that were found.

=over

=item C<Arguments>

 input(required): Merged blast output.
 best_only(0): Keep only the best hit for each query?
 evalue(undef): Set an optional evalue cutoff.
 search_type('blastxml'): This needs to be a parseable blast format.

=cut
sub Parse_Blast {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        best_only => 0,
        evalue => undef,
        search_type => 'blastxml',
    );
    my $input = $options->{input};
    my $output = $input;
    my $best = $options->{best_only};
    my $search_type = $options->{search_type};
    $output =~ s/\.txt\.xz/_parsed.txt/g;
    print "Writing parsed output to ${output}\n";
    my $count_table = $output;
    $count_table =~ s/_parsed\.txt/_counts\.txt/g;
    print "Writing count table to ${count_table}\n";
    my $counts = FileHandle->new(">$count_table");
    $counts->autoflush(1);
    print $counts "Query\tquery_accession\thits\n";
    my $fh = FileHandle->new(">$output");
    $fh->autoflush(1);
    local $| = 0; ## I want to watch the output as it happens.
    my $in = FileHandle->new("less ${input} |");
    $XML::SAX::ParserPackage = 'XML::SAX::PurePerl';
    my $searchio = Bio::SearchIO->new(-fh => $in);
    ## Perform a test for the output format, this will close the filehandle, so reopen it after
    my $guessed_format = $searchio->format();
    print STDERR "Guessed the format is: ${guessed_format}\n";
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
        print STDERR "Parsed ${result_count} of $options->{state}->{num_sequences} results.\n";
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
            my $print_string = qq"${query_name}\t${real_name}\t${acc2}\t${acc3}\t${start}\t${end}\t${ident}\t${score}\t${sig}\t${query_length}\t${hsp_identical}\t${hit_count}\n";
            print $fh $print_string;
            if ($args{best_only}) {
                next RESULT;
            }
        } ## End each hit for a single result
        print $counts "${query_name}\t${real_name}\t${hit_count}\n";
        ## print "$result had $hit_count hits\n";
    } ## Finish looking at each sequence
    $fh->close();
    $in->close();
    $counts->close();
    return($result_count);
}

=back

=head2 C<Run_Blast_Parse>

 Create a blast database and run blast on it.

 This creates a blast library from 'library' and performs individual
 blast searches for every sequence in nt_query.fasta against it.
 Finally, it writes a summary of the results in a set of tables.

 Run_Parse_Blast() is rather slow because it uses Bio::SearchIO and
 Bio::Tools::Run::StandAloneBlast which in turn calls blastall
 separately for every sequence.  In many cases, the user will likely
 prefer to use the the 'blastsplit' method.

 In addition, this should have the blast parameters as an input
 variable.

=over

=item C<Arguments>

 input(required): Input sequences to search.
 library(required): Library of sequences to search against.
 blast_tool('blastp'): Use this blast tool.
 evalue(0.01): Set an evalue cutoff.
 output('blast_output.txt'): Output file to write.
 modules('blast'): Which environment module to load.

=cut
sub Run_Parse_Blast {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'library'],
        blast_tool => 'blastp',
        evalue => 0.01,
        output => 'blast_output.txt',);
    my %modules = Get_Modules(caller => 1);
    my $loaded = $class->Module_Loader(%modules);
    my $query = $options->{input};
    my $library_path = $class->Bio::Adventure::Index::Check_Blastdb(%args);
    my $library = $options->{library};
    my $number_hits = 0;
    my $libname = basename($library, ('.fasta'));
    ## my $query_library = Bio::SeqIO->new(-file => ${query}, -format => 'Fasta');
    my $output_directory = qq"$options->{basedir}/outputs/${query}_vs_${libname}";
    make_path("${output_directory}") unless (-d $output_directory);
    my $counts = FileHandle->new(">${output_directory}/counts.txt");
    my $zeros = FileHandle->new(">${output_directory}/zeros.fasta");
    my $singles = FileHandle->new(">${output_directory}/singles.fasta");
    my $doubles = FileHandle->new(">${output_directory}/doubles.fasta");
    my $few = FileHandle->new(">${output_directory}/few.fasta");
    my $many = FileHandle->new(">${output_directory}/many.fasta");
    my $seq_count = 0;

    my @params = (
        -program => $options->{blast_tool},
        ## Show GI's in definition lines?
        ##-F => 'T', ## Filter query sequence
        ##-G => '-1', ## Cost to open a gap
        ##-E => '-1', ## Cost to extend a gap
        ##-X => '0',  ## X dropoff value for gapped alignment
        ##-I => 'T',  ## Default in blast is F
        ## -q => '-3', ## Penalty for nucleotide mismatch (blastn only)
        ## And many more
        -db_name => $library,);

    ##while (my $query_seq = $query_library->next_seq()) {
    ##    ## I think this while loop may not be needed.
    ##    $seq_count++;
    ##    my $id = $query_seq->id;
    ##    my $desc = $query_seq->desc;
    ##    my $seq = $query_seq->seq;
    my $search = Bio::Tools::Run::StandAloneBlastPlus->new(@params);
    my $blast_output;
    eval {
        ## $blast_output = $search->blastall($query_library);
        if ($options->{blast_tool} eq 'blastp') {
            $blast_output = $search->blastp(-query => $options->{input},
                                            -outfile => $options->{output});
        } elsif ($options->{blast_tool} eq 'blastn') {
            $blast_output = $search->blastn(-query => $options->{input},
                                            -outfile => $options->{output});
        } elsif ($options->{blast_tool} eq 'tblastx') {
            $blast_output = $search->tblastx(-query => $options->{input},
                                             -outfile => $options->{output});
        } else {
            $blast_output = $search->blastx(-query => $options->{input},
                                            -outfile => $options->{output});
        }
    };
    if ($@) {
        print "Error? $@\n";
    }
    my $blast_results = Bio::SearchIO->new(-format => 'blast',
        -file => $options->{output});

    my $result_count = 0;
    while (my $result = $blast_results->next_result()) {
        $result_count++;
        my $query_name = $result->query_name();
        my $query_length = $result->query_length();
        my $query_descr = $result->query_description();
        my $stats = $result->available_statistics();
        my $hits = $result->num_hits();
        my $hit_count = 0;
        my $score_cutoff = 100;
        my $sig_cutoff = $options->{evalue};
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
      } ## End iterating over each hit of the result.
        my $entry = "${query_name} ${query_descr}\n";
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
        print $counts "${query_name}\t${hit_count}\n";
    ## } ## End of the individual blast search  (maybe not needed?)
    } ## End of this search library  (the 7 states to search against)
    $counts->close();
    $zeros->close();
    $singles->close();
    $doubles->close();
    $few->close();
    $many->close();
    my $unloaded = $class->Module_Reset(env => $loaded);
    return($number_hits);
}

=back

=head2 C<Split_Align_Blast>

 Split up a pile of sequences and search them in parallel.

 Split apart a set of query sequences into $args{align_jobs} pieces and align
 them all separately.  This is the primary function to call when one
 wants to run a whole lot of blast.

=over

=item C<Arguments>

 input(required): Fasta file containing the sequences for searching.
 library(required): Fasta file containing the library of sequences to search.
 e(10): Default parameters for blast (this should be blast_param)
 blast_tool(blastn): Which blast method to invoke?
 align_jobs(40): How many jobs to invoke?
 align_parse(0): Start a parsing job upon completion?
 blast_format(5: blastxml): Which blast format for the output?
 best_only(0): Report only the best hits per search?

=item C<Invocation>

> cyoa --task align --method blastsplit --input query.fasta --library library.fasta

=cut
sub Split_Align_Blast {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'library',],
        param => ' -evalue 10 ',
        blast_tool => 'blastn',
        align_jobs => 40,
        align_parse => 0,
        blast_format => 5,
        num_dirs => 0,
        best_only => 0,
        interactive => 0,
        jmem => 8,);
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
    my $que = basename($options->{input}, ('.fasta'));
    my $outdir = qq"$options->{basedir}/outputs";
    make_path("${outdir}") unless(-d ${outdir});
    my $output = qq"${outdir}/${que}_vs_${lib}.txt";
    my $concat_job;
    $lib = $class->Bio::Adventure::Index::Check_Blastdb(%args);
    if ($options->{cluster}) {
        my $split_info = $class->Bio::Adventure::Align::Get_Split(%args);
        my $num_per_split = $split_info->{num_per_split};
        print "Going to make $options->{align_jobs} directories with ${num_per_split} sequences each.\n";
        my $actual = $class->Bio::Adventure::Align::Make_Directories(%args,
            num_per_split => $split_info);
        print "Actually used ${actual} directories to write files.\n";
        my $alignment = $class->Bio::Adventure::Align_Blast::Make_Blast_Job(
            library => $lib,
            align_jobs => $actual,
            output_type => $options->{blast_format},);
        $concat_job = $class->Bio::Adventure::Align::Concatenate_Searches(
            jdepends => $alignment->{job_id},
            output => ${output},);
    } else {
        ## If we don't have pbs, force the number of jobs to 1.
        print "Not using the cluster.\n";
        ## $options = $class->Set_Vars(align_jobs => 1);
        ## $class->{align_jobs} = 1;
        my $num_per_split = $class->Bio::Adventure::Align::Get_Split();
        $options = $class->Set_Vars(num_per_split => $num_per_split);
        my $actual = $class->Bio::Adventure::Align::Make_Directories(
            num_per_split => $num_per_split);
        my $alignment = $class->Bio::Adventure::Align_Blast::Make_Blast_Job(
            library => $lib,
            align_jobs => $actual,
            output_type => $options->{blast_format},);
        $concat_job = $class->Bio::Adventure::Align::Concatenate_Searches(
            output => ${output},);
    }

    my $parse_input = $concat_job->{output};
    my $parse_output = $parse_input;
    $parse_output =~ s/\.txt\.xz/_parsed.txt/g;
    my $count_output = $parse_input;
    $count_output =~ s/_parsed\.txt/_counts\.txt/g;
    my $comment_string = qq!## I don't know if this will work.!;
    my $jstring = qq?
use Bio::Adventure;
use Bio::Adventure::Align;
use Bio::Adventure::Align_Blast;
my \$h = new Bio::Adventure;
my \$result = Bio::Adventure::Align::Parse_Search(
  \$h, input => '${parse_input}',
  search_type => 'blastxml',
  parsed_output => '$parse_output',
  count_output => '$count_output',
  best => $options->{best_only});
?;
    my $parse_job = $class->Submit(
        comment => $comment_string,
        count_output => $count_output,
        jdepends => $concat_job->{job_id},
        jmem => $options->{jmem},
        jname => 'parse_search',
        jstring => $jstring,
        language => 'perl',
        parsed_output => $parse_output,);
    $concat_job->{parser} = $parse_job;
    return($concat_job);
}

## For the moment, just submit a job which assumes a 'input' directory in the cwd
## I want to think about how I might handle different ways of running orthomcl.
sub OrthoMCL_Pipeline {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],);
    my $job_name = 'orthomcl';
    ## Note that I am cheating and using a pre-defined pipeline for orthomcl
    ## https://github.com/apetkau/orthomcl-pipeline
    ## This also assumes that everything for orthomcl and its pipeline has already been configured.
    my $inputs = $class->Get_Paths($options->{input});
    my $cwd_name = basename(cwd());

    my $comment = '## This is a orthomcl submission script.';
    my $jstring = qq!
orthomcl-pipeline.pl -i input -o output \\
  --nocompliant --yes
!;
    my $job = $class->Submit(
        comment => $comment,
        jcpu => 12,
        jdepends => $options->{jdepends},
        jname => 'orthomcl',
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => 16,
        output => 'output',
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        walltime => '144:00:00',);
    return($job);
}

=back

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

L<Bio::Tools::Run::StandAloneBlast> L<blast>

=cut

1;
