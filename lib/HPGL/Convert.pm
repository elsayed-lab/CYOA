=head2
    Gff2Gtf()
=cut
sub Gff2Gtf {
    my $in = shift;
    use File::Basename;
    use Bio::FeatureIO;

    my ($name, $path, $suffix) = fileparse($in, qr/\.gff/);
    my $outFile = $path . $name . ".gtf";

    my $inGFF = new Bio::FeatureIO('-file' => "$in",
                                   '-format' => 'GFF',
                                   '-version' => 3);
    my $outGTF = new Bio::FeatureIO('-file' => ">$outFile",
                                    '-format' => 'GFF',
                                    '-version' => 2.5);

    while (my $feature = $inGFF->next_feature() ) {
        $outGTF->write_feature($feature);
    }
}

=head2
    Sam2Bam()
=cut
sub Sam2Bam {
    my $me = shift;
    $me->Check_Options(["species"]);
    my %args = @_;
    ## A little logic to see if the invocation of this is for an extent .sam file
    ## Or calling it on an existing .fastq(.gz), in which case one must assume
    ## It is being called on one or more files in the bowtie_out/ directory
    ## And which start with $basename, include something like -trimmed-v0M1.sam
    my $basename = $me->{basename};
    my @input_list = ();
    if ($args{input}) {
        push(@input_list, $args{input});
    } elsif (-r $me->{input} and $me->{input} =~ /\.sam$/) {
        push(@input_list, $me->{input});
        my $sam = $me->Samtools(sam => \@input_list);
    } elsif (-r $me->{input} and $me->{input} =~ /\.fastq$/) {
        if (-d "bowtie_out") {
            find({ wanted => sub { push(@input_list, "bowtie_out/$_") if ($_ =~ /\.sam$/); }, follow => 1 }, 'bowtie_out/');
            my $sam = $me->Samtools(sam => \@input_list);
        } else {
            foreach my $k (%{$me->{bt_args}}) {
                my $output_string = "bowtie_out/${basename}-${k}.sam";
                push(@input_list, $output_string);
            }
            my $bt = $me->Bowtie();
            my $sam = $me->Samtools(depends => $bt->{pbs_id}, sam => \@input_list);
        }
    } else {
        die("I don't know what to do without a .fastq file nor a .sam file.\n");
    }
}


=head2
    Samtools()
=cut
sub Samtools {
    my $me = shift;
    my %args = @_;
    my $basename = $me->{basename};
    my $depends = "";
    if ($args{depends}) {
        $depends = $args{depends};
    }
    my $input = $args{input};
    my $output = $input;
    $output =~ s/\.sam$/\.bam/g;
    my $sorted = $input;
    $sorted =~ s/\.sam$//g;
    $sorted = qq"${sorted}-sorted";
    print "Converting to a compressed/sorted bam file.\n";
    my $job_string = qq!samtools view -u -t $me->{libdir}/genome/$me->{species}.fasta -S $input 1>$output
samtools sort -l 9 $output $sorted
rm $output && rm $input && mv ${sorted}.bam $output && samtools index $output
bamtools stats -in $output 2>${output}.stats 1>&2
!;
    my $comment = qq!## Converting the text sam to a compressed, sorted, indexed bamfile.
## Also printing alignment statistics to ${output}.stats!;
    my $samtools = $me->Qsub(job_name => "sam",
                             depends => $depends,
                             job_string => $job_string,
                             input => $input,
                             output => $output,
                             comment => $comment,
			     prescript => $args{prescript},
			     postscript => $args{postscript},
        );
    return($samtools);
}

1;
