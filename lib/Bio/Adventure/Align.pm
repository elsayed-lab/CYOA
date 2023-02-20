package Bio::Adventure::Align;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use Bio::SearchIO::fasta;
use Bio::Seq;
use Cwd;
use File::Basename;
use File::Find;
use File::Path qw"make_path remove_tree";
use File::Which qw"which";
use POSIX qw"ceil";

=head1 NAME

 Bio::Adventure::Align - Functions shared by multiple sequence alignment methods.

=head1 SYNOPSIS

 These functions set up and clean up blast/fasta/etc searches. (I need to set up hmmer and friends).

=head1 METHODS

=head2 C<Cleanup>

 Cleans up the mess of temporary files/directories left behind by these tools.

=cut
sub Cleanup {
    my ($class, %args) = @_;
    my $ret = qx("rm -rf outputs split blastdb formatdb.log split_align_submission.sh");
    return($ret);
}

=head2 C<Concatenate_Searches>

 Bring together multiple parallel search results.

 This function waits until the cluster finishes processing all the blast jobs, then
 concatenates all the output files into one compressed file .

=cut
sub Concatenate_Searches {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        workdir => '',
        jmem => 8,);
    my $workdir = $options->{workdir};
    my $finished = 0;
    my $output = qq"${workdir}/split_search.txt";
    $output = $options->{output} if (defined($options->{output}));
    $output .= ".xz" unless ($output =~ /\.xz$/);
    my $comment_string = qq"## Concatenating the output files into ${output}
";
    my $jstring = qq!
rm -f ${output} && for i in \$(/bin/ls outputs/split/*.stdout); do xz -9e -c \$i >> ${output}; done
!;
    my $concatenate = $class->Submit(
        comment => $comment_string,
        depends_type => 'array',
        jdepends => $options->{jdepends},
        jname => "concat",
        jmem => $options->{jmem},
        jstring => $jstring,
        jprefix => $options->{jprefix},
        stdout => $output,
        output => $output,);
    return($concatenate);
}

=head2 C<Duplicate_Remove>

 Get rid of duplicate sequences when writing out fasta files.

=cut
sub Duplicate_Remove {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $input = $options->{input};
    my $output = qq"${input}_removed.txt";
    my $in = FileHandle->new("<$input");
    my @entries = ();
    while (my $line = <$in>) {
        chomp $line;
        my ($self, $others) = split(/\t/, $line);
        my @other_list = split(/ /, $others);
        my $entry = {self => $self, others => \@other_list};
        print "Pushing self: ${self} others: @other_list\n";
        push(@entries, $entry);
    }
    $in->close();
    my $out = FileHandle->new(">$output");
    my $num_found = 0;
    while (scalar(@entries) > 0) {
        my $keeper = shift(@entries);
        my $self = $keeper->{self};
        my $others = $keeper->{others};
        ## Write down this first entry
        print "Writing ${self} to out.\n";
        print $out qq"${self}\n";
        ## Now remove all instances of the 'others' from the list and future consideration.
        my $other_length = scalar(@{$others});
      OTHERS: foreach my $other (@{$others}) {
          my ($other_name, $other_ident, $other_e) = split(/:/, $other);
          print "Checking for ${other_name} in entries.\n";
          my @tmp = @entries;
          my $length = scalar(@tmp);
          my @new = ();
          print "The length of entries is now: ${length}\n";
        CHECK: foreach my $tmp_entry (@tmp) {
            my $self_check = $tmp_entry->{self};
            if ($self_check eq $other_name) {
                $num_found = $num_found + 1;
                ## print "FOUND! $self_check $other_name\n";
                next CHECK;
            } else {
                ## print "NOT FOUND! $self_check $other_name\n";
                push(@new, $tmp_entry);
            }
        }
          @entries = @new;
      }
    }
    $out->close();
    return($num_found);
}

=head2 C<Get_Split>

 Split up sequences into even numbers for submitting in parallel.

 Take the number of sequences in the query library and divides by the
 number of jobs to run, get the ceiling of that, and prints that many
 sequences to each file in split/ starting at 1000.  (So, as long as
 you have <= 8,999 jobs PBS won't get confused by the difference
 between 099 and 100 because they will be 1099 and 1100 instead.)

=cut
sub Get_Split {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $in = Bio::SeqIO->new(-file => $options->{input},);
    my $seqs = 0;
    while (my $in_seq = $in->next_seq()) {
        $seqs++;
    }
    my $seqs_per_job = ceil($seqs / $options->{align_jobs});
    if ($seqs < $options->{align_jobs}) {
        print "There are fewer sequences than the chosen number of split files, resetting that to ${seqs}.\n";
        $seqs_per_job = $seqs;
    }
    print "Get_Split: Making $options->{align_jobs} directories with ${seqs_per_job} sequences.\n";
    my $ret = {
        num_per_split => $seqs_per_job,
        seqs => $seqs,
    };
    return($ret);
}

=head2 Ids_to_Sequence

 Writes relatively nicely formatted fasta files given a genome and gff file.

=cut
sub Ids_to_Sequences {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $id_input_file = $options->{id_list};
    my $id_in = FileHandle->new("<${id_input_file}");
    my $chromosomes = $class->Read_Genome(genome => $class->{genome});
    my $gff_info = $class->Read_GFF(gff => $class->{gff});
    ## Now I have the set of start/stops for all IDs
    my %sequences = ();
    ## Create a set of sequences by ID
    foreach my $chr (keys %{$gff_info}) {
        my %inner = %{$gff_info};
        foreach my $id (keys %inner) {
            my $start = $inner{$id}->{start};
            my $end = $inner{$id}->{end};
            $sequences{$id} = $chromosomes->{$chr}->subseq($start, $end);
        }
    }
    my $id_output_file = qq"${id_input_file}.fasta";
    my $id_out = FileHandle->new(">$id_output_file");
    while (my $line = <$id_in>) {
        chomp $line;
        my $seq = join("\n",($sequences{$line} =~ m/.{1,80}/g));
        print $id_out ">${line}
$sequences{$seq}
";
    }
    $id_out->close();
    $id_in->close();
    ## This should generate files which end in _parsed.txt for the library
    my $concat = Bio::Tools::Adventure::Align_Fasta::Split_Align_Fasta(
        $class,
        library => $id_output_file,
        max_significance => 0.0001,
        min_percent => 0.80,
        query => $id_output_file,);
    my $duplicates = Bio::Tools::Adventure::Align::Duplicate_Remove(
        $class,
        input => $id_output_file,);
    return($concat);
}

=head2 C<Make_Directories>

 Create sequential directories for each job and actually write the
 sequences using Bio::SeqIO

=cut
sub Make_Directories {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        num_per_split => 100,
        align_jobs => 40,);
    my $num_per_split = $options->{num_per_split};
    my $splits = $options->{align_jobs};
    my $workdir = $options->{workdir};
    ## I am choosing to make directories starting at 1000
    ## This way I don't have to think about the difference from
    ## 99 to 100 (2 characters to 3) as long as no one splits more than 9000 ways...
    print "Make_Directories: Making $options->{align_jobs} directories with $options->{num_per_split} sequences.\n";
    my $dir = $options->{array_start};

    remove_tree("outputs/split", {verbose => 0 });
    for my $c ($dir .. ($dir + $splits)) {
        ## print "Making directory: split/$c\n";
        make_path("outputs/split/${c}") if (!-d "outputs/split/${c}" and !-f "outputs/split/${c}");
    }

    my $in = Bio::SeqIO->new(-file => $options->{input},);
    my $count = 0;
    while (my $in_seq = $in->next_seq()) {
        my $id = $in_seq->id();
        my $seq = $in_seq->seq();
        $seq = join("\n", ($seq =~ m/.{1,80}/g));
        my $output_file = qq"outputs/split/${dir}/in.fasta";
        my $outfile = FileHandle->new(">>${output_file}");
        my $out_string = qq!>${id}
${seq}
!;
        print $outfile $out_string;
        $outfile->close();
        $count++;
        $options->{state}->{num_sequences}++;
        if ($count >= $num_per_split) {
            $count = 0;
            $dir++;
            my $last = $dir;
            $last--;
            if ($count == 1) {
                print "Writing ${num_per_split} entries to files $options->{array_start} to ${last}\n";
            }
            ## You might be wondering why this num_dirs is here.
            ## Imagine if you have a query library of 10,004 sequences and you try to write them to 200 files.
            ## This script will write 51 sequences / file, but will therefore not quite finish writing all 200 files.
            ## As a result, when we go back at the end to clean up,
            ## the loop which tries to detect how many jobs are finished will never exit
            ## because it won't ever reach 200.
            ## Therefore this keeps track of the last file written and uses that instead.
            my $max_dirs = $options->{array_start} - 1;
            $options->{num_dirs} = ($last - $max_dirs); ## This - 1000 is there because we start at job 1000, but we start counting at 0
        }                      ## End for each iteration of $num_per_split files
    }                          ## End while reading the fasta
    my $actual_number_dirs_used = $options->{num_dirs};
    print "Finished writing files.\n";
    return($actual_number_dirs_used);
}

=head2 C<Parse_Search>

 This passes off to other parsers depending on the input format.

=cut
sub Parse_Search {
    my ($class, %args) = @_;
    my $ret;
    if ($args{search_type} eq 'blastxml') {
        $ret = Bio::Adventure::Align_Blast::Parse_Blast($class, %args);
    } elsif ($args{search_type} eq 'local_fasta') {
        $ret = Bio::Adventure::Align_Fasta::Parse_Fasta($class, %args);
    } elsif ($args{search_type} eq 'global_fasta') {
        $ret = Bio::Adventure::Align_Fasta::Parse_Fasta_Global($class, %args);
    } else {
        $ret = Bio::Adventure::Align_Fasta::Parse_Fasta($class, %args);
    }
    return($ret);
}

sub Pairwise_Similarity_Matrix {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        cutoff => 500,
        jcpus => 8,
        modules => ['blast'],
        required => ['input',],);
    my $input = $options->{input};
    if (-d $input) {
        ## Then this should be a directory containing fasta files, concatenate them into a single
        ## sequence database for use.
        my @directory = ($input);
        my @file_list = ();
        find(
            sub {
                push(@file_list, $File::Find::name) if ($File::Find::name =~ /\.faa$/);
            }, @directory);
        my $write_input = FileHandle->new(">${input}.faa");
        for my $f (@file_list) {
            my $read = FileHandle->new("<${f}");
            while (my $line = <$read>) {
                chomp $line; ## Just in case some files are \r\n, chomp should help.
                print $write_input "${line}\n";
            }
            $f->close();
        }
        $write_input->close();
    }

    my $pre = $input;
    $pre =~ s/\.faa//g;
    my $pre_dir = dirname($input);

    my $result_data = {};
    my $e_values = {};
    my $bit_values = {};
    my $scores = {};
    my $evalue_fh = FileHandle->new(qq">${pre_dir}/${pre}_pairwise_evalues.tsv");
    my $bit_fh = FileHandle->new(">${pre_dir}/${pre}_pairwise_bits.tsv");
    my $score_fh = FileHandle->new(">${pre_dir}/${pre}_pairwise_scores.tsv");
    my $output = FileHandle->new(">${pre_dir}/${pre}_out.txt");
    my $prefix = basename($options->{input});

    my @params = (
        -create => 1,
        -db_data => $input,
        -db_name => $prefix,
        -num_alignments => $options->{cutoff},
        -num_descriptions => $options->{cutoff},
        -num_threads => $options->{jcpus},
        -program => 'blastp',
        );

    my $search = Bio::Tools::Run::StandAloneBlastPlus->new(@params);
    my @parameters = $search->get_parameters;

    my $outfile = qq"outputs/${prefix}.txt";
    my $blast_report = $search->blastp(-query => $input, -outfile => $outfile);
    my $search_output = Bio::SearchIO->new(-file => $outfile, -format => 'blast');
    my @hit_lst = ();
    my $number_hits = 0;
    my $number_result = 0;
  RESULTS: while (my $result = $search_output->next_result) {
      $number_result++;
      print "Result: ${number_result}\n";
      my $query_name = $result->query_name();
      my $query_length = $result->query_length();
      my $query_descr = $result->query_description();
      my $stats = $result->available_statistics();
      my $hits = $result->num_hits();
      my $datum_ref = {
          description => $query_descr,
          name => $query_name,
          stats => $stats,
          length => $query_length,
          hit_data => {},
      };
      $result_data->{$query_name} = $datum_ref;
      my $hit_count = 0;
    HITLOOP: while (my $hits = $result->next_hit()) {
        $number_hits = $number_hits++;
        my $hit_name = $hits->name();
        ## The hit_name should cross reference to one of the two accession columns in xref_aoh
        ## from the beginning of this function.
        my $longname = '';
        my $hit_length = 0;
        if (defined($hits->length())) {
            $hit_length = $hits->length();
        }
        my $hit_acc = $hits->accession();
        my $hit_descr = $hits->description();
        my $hit_score = $hits->score();
        my $hit_sig = $hits->significance();
        my $hit_bits = $hits->bits();

        ## Create our values matrices
        $e_values->{$query_name}->{$hit_name} = $hit_sig;
        $bit_values->{$query_name}->{$hit_name} = $hit_bits;
        $scores->{$query_name}->{$hit_name} = $hit_score;

        #my %hit_datum = (
        #    query_name => $query_name,
        #    query_length => $query_length,
        #    query_descr => $query_descr,
        #    stats => $stats,
        #    length => $hit_length,
        #    acc => $hit_acc,
        #    description => $hit_descr,
        #    score => $hit_score,
        #    sig => $hit_sig,
        #    bit => $hit_bits,
        #    longname => $longname,);
        #push (@hit_lst, \%hit_datum);
        #$result_data->{$query_name}->{hit_data}->{$hit_name} = \%hit_datum;
    } ## End looking at this for this result.
  } ## End iterating over the results for this sequence.

    my @rows = keys %{$e_values};
    my $header_string = "ID\t";
    for my $r (@rows) {
        $header_string .= qq"${r}\t";
    }
    $header_string =~ s/\t$/\n/;
    print $evalue_fh $header_string;
    print $bit_fh $header_string;
    print $score_fh $header_string;

    my ($evalue_row_string, $bit_row_string, $score_row_string) = '';
    for my $r (@rows) {
        $evalue_row_string = qq"${r}\t";
        $bit_row_string = qq"${r}\t";
        $score_row_string = qq"${r}\t";
        for my $c (@rows) {
            my $pair_evalue = $e_values->{$c}->{$r};
            $pair_evalue = 0 if ($c eq $r);
            $pair_evalue = 1 if (!defined($pair_evalue));
            $evalue_row_string .= qq"${pair_evalue}\t";
            my $pair_bit = $bit_values->{$c}->{$r};
            $pair_bit = 0 if (!defined($pair_bit));
            $bit_row_string .= qq"${pair_bit}\t";
            my $pair_score = $scores->{$c}->{$r};
            $pair_score = 0 if (!defined($pair_score));
            $score_row_string .= qq"${pair_score}\t";
        }
        $evalue_row_string =~ s/\t$/\n/;
        print $evalue_fh $evalue_row_string;
        $bit_row_string =~ s/\t$/\n/;
        print $bit_fh $bit_row_string;
        $score_row_string =~ s/\t$/\n/;
        print $score_fh $score_row_string;
    }
    $evalue_fh->close();
    $bit_fh->close();
    $score_fh->close();
}

sub Map_Accession {
    my %args = @_;
    my $in = Bio::SeqIO->new(-file => $args{fasta}, -format => 'fasta');
    my $map = {};
  SEQS: while (my $seq = $in->next_seq) {
      my $id = $seq->id;
      my $short_id = $id;
      $short_id =~ s/(^.*)\..*$/$1/g;
      my $test = $seq->desc;
      my $taxonomy = $test;
      ## In theory I should be able to make a relatively simple capturing regex to
      ## pull the taxonomy information.  I am insufficiently intelligent to do this.
      my @stuff_taxonomy = split(/\[/, $taxonomy);
      my $taxonomy_info = $stuff_taxonomy[$#stuff_taxonomy];
      $taxonomy_info =~ s/\]$//g;
      my ($genus, $species, $etc) = split(/\s+/, $taxonomy_info);
      $species = '' if (!defined($species));
      my $new = qq"${genus}_${species}_${short_id}";
      $map->{$id} = $new;
  }
    return($map);
}

sub OrthoFinder {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jmem => 24,
        jprefix => '50',
        modules => ['orthofinder'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});

    my $jname = qq'orthofinder';

    my $outdir = qq"outputs/$options->{jprefix}orthofinder";
    my $error_file = qq"${outdir}/orthofinder.stderr";
    my $stdout_file = qq"${outdir}/orthofinder.stdout";
    my $comment = qq!## Attempting to run orthofinder using a faa directory in $options->{input}.
!;
    my $jstring = qq!
orthofinder -f $options->{input} \\
  -o ${outdir} \\
  2>>${error_file} 1>>${stdout_file}
!;
    my $ortho = $class->Submit(
        comment => $comment,
        input => $options->{input},
        jdepends => $options->{jdepends},
        jname => ${jname},
        jprefix => $options->{jprefix},
        jstring => $jstring,
        stderr => $error_file,
        stdout => $stdout_file,
        jmem => $options->{jmem},
        modules => $options->{modules},);
    return($ortho);
}

=head1 AUTHOR - atb

    Email  <abelew@gmail.com>

=head1 SEE ALSO

    L<Bio::Adventure::Align_Blast> L<Bio::Adventure::Align_Fasta>

=cut

1;
