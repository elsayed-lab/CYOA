package CYOA;
use common::sense;
use autodie qw":all";

use Cwd;
use Bio::Seq;
use Bio::SearchIO::fasta;
use POSIX qw"ceil";

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
rm -f ${output}.gz && for i in \$(/bin/ls outputs/*.out); do gzip -c \$i >> ${output}; done
!;
    my $concatenate_job = $me->Qsub(job_name => "concat",
                                    depends => $args{depends},
                                    depends_type => 'array',
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
    if ($args{search_type} eq 'blastxml') {
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
    system("rm -rf outputs split blastdb formatdb.log split_align_submission.sh");
    exit(0);
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

    remove_tree("split", {verbose => 1 });
    for my $c ($dir .. ($dir + $splits)) {
        ## print "Making directory: split/$c\n";
        make_path("split/$c") if (!-d "split/$c" and !-f "split/$c");
    }

    my $in = new Bio::SeqIO(-file => $args{query},);
    my $count = 0;
    while (my $in_seq = $in->next_seq()) {
        my $id = $in_seq->id();
        my $seq = $in_seq->seq();
        my $outfile = new FileHandle;
        my $output_file = qq"split/${dir}/in.fasta";
        $outfile->open(">>split/${dir}/in.fasta");
        my $out_string = qq!>$id
$seq
!;
        print $outfile $out_string;
        $outfile->close();
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
    my $actual_number_dirs_used = $args{num_dirs};
    return($actual_number_dirs_used);
}

=back

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

    L<Bio::Tools::Run::StandAloneBlast>

=cut

## Culling duplicates requires a few steps:
## 1.  Extract the sequences of the genes of interest in any available method
## 2.  Perform some search of the set of genes against itself
## 3.  Extract a single family member from each resulting set
sub Ids_to_Sequences {
    my $me = shift;
    my %args = @_;
    $me->Get_Options(['id_list', 'genome', 'gff']);
    my $id_input_file = $me->{id_list};
    my $id_in = new FileHandle;
    $id_in->open("<$id_input_file");
    my $chromosomes = $me->Read_Genome(genome => $me->{genome});
    my $gff_info = $me->Read_GFF(gff => $me->{gff});
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
    my $id_out = new FileHandle;
    my $id_output_file = qq"${id_input_file}.fasta";
    $id_out->open(">$id_output_file");
    while (my $line = <$id_in>) {
        chomp $line;
        print $id_out ">$line
$sequences{$line}
";
    }
    $id_out->close();
    $id_in->close();
    ## This should generate files which end in _parsed.txt for the library
    my $concat = $me->Split_Align_Fasta(query => $id_output_file, library => $id_output_file,  max_significance => 0.0001, min_percent => 0.80);
    $me->Duplicate_Remove(input => $id_output_file);
}

sub Duplicate_Remove {
    my $me = shift;
    my %args = @_;
    my $input = $args{input};
    my $output = qq"${input}_removed.txt";
    open(IN, "<$input");
    my @entries = ();
    while(my $line = <IN>) {
        chomp $line;
        my ($self, $others) = split(/\t/, $line);
        my @other_list = split(/ /, $others);
        my $entry = {self => $self, others => \@other_list};
        print "Pushing self: $self others: @other_list\n";
        push(@entries, $entry);
    }
    close(IN);
    open(OUT, ">$output");
    while (scalar(@entries) > 0) {
        my $keeper = shift(@entries);
        my $self = $keeper->{self};
        my $others = $keeper->{others};
        ## Write down this first entry
        print "Writing $self to out.\n";
        print OUT qq"${self}\n";
        ## Now remove all instances of the 'others' from the list and future consideration.
        my $other_length = scalar(@{$others});
        print "TESTME: $other_length @{$others}\n";
      OTHERS: foreach my $other (@{$others}) {
          my ($other_name, $other_ident, $other_e) = split(/:/, $other);
          print "Checking for $other_name in entries.\n";
          my @tmp = @entries;
          my $length = scalar(@tmp);
          my @new = ();
          print "The length of entries is now: $length\n";
        CHECK: foreach my $tmp_entry (@tmp) {
            my $self_check = $tmp_entry->{self};
            ##print "TEST: $self_check $other_name\n";
            if ($self_check eq $other_name) {
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
    close(OUT);
}

1;
