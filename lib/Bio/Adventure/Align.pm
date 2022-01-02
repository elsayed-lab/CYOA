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
        jmem => 8,);
    my $workdir = $options->{workdir};
    my $finished = 0;
    my $output = qq"${workdir}/split_search.txt";
    $output = $options->{output} if (defined($options->{output}));
    $output .= ".xz" unless ($output =~ /\.xz$/);
    my $comment_string = qq"## Concatenating the output files into ${output}
";
    my $jstring = qq!
rm -f ${output} && for i in \$(/bin/ls ${workdir}/*.out); do xz -9e -c \$i >> ${output}; done
!;
    my $concatenate = $class->Submit(
        comment => $comment_string,
        depends_type => 'array',
        jdepends => $options->{jdepends},
        jname => "concat",
        jmem => $options->{jmem},
        jstring => $jstring,
        jprefix => $options->{jprefix},
        output => qq"${output}",);
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
    my $ret = ceil($seqs / $options->{align_jobs});
    if ($seqs < $options->{align_jobs}) {
        print "There are fewer sequences than the chosen number of split files, resetting that to ${seqs}.\n";
        $options = $class->Set_Vars(align_jobs => $seqs);
    }
    print "Get_Split: Making $options->{align_jobs} directories with $ret sequences.\n";
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
    my $dir = 1000;

    remove_tree("${workdir}/split", {verbose => 0 });
    for my $c ($dir .. ($dir + $splits)) {
        ## print "Making directory: split/$c\n";
        make_path("${workdir}/split/${c}") if (!-d "${workdir}/split/${c}" and !-f "${workdir}/split/${c}");
    }

    my $in = Bio::SeqIO->new(-file => $options->{input},);
    my $count = 0;
    while (my $in_seq = $in->next_seq()) {
        my $id = $in_seq->id();
        my $seq = $in_seq->seq();
        $seq = join("\n", ($seq =~ m/.{1,80}/g));
        my $output_file = qq"${workdir}/split/${dir}/in.fasta";
        my $outfile = FileHandle->new(">>${workdir}/split/${dir}/in.fasta");
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
                print "Writing ${num_per_split} entries to files 1000 to ${last}\n";
            }
            ## You might be wondering why this num_dirs is here.
            ## Imagine if you have a query library of 10,004 sequences and you try to write them to 200 files.
            ## This script will write 51 sequences / file, but will therefore not quite finish writing all 200 files.
            ## As a result, when we go back at the end to clean up,
            ## the loop which tries to detect how many jobs are finished will never exit
            ## because it won't ever reach 200.
            ## Therefore this keeps track of the last file written and uses that instead.
            $options->{num_dirs} = ($last - 999); ## This - 1000 is there because we start at job 1000, but we start counting at 0
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

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

L<Bio::Adventure::Align_Blast> L<Bio::Adventure::Align_Fasta>

=cut

1;
