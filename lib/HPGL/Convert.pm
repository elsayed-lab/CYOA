package HPGL;

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

sub Gb2Gff {
    my %options = (
        input => undef,
        base => undef,
        );
    my $opt = GetOptions(
        "input|i:s" => \$options{input},
        );
    if (!defined($options{input})) {
        my $term = new Term::ReadLine('>');
        my $attribs = $term->Attribs;
        $attribs->{completion_suppress_append} = 1;
        my $OUT = $term->OUT || \*STDOUT;
        $options{input} = $term->readline("Please provide an input genbank file: ");
        $options{input} =~ s/\s+$//g;
    }
    my @suffix = (".gb", ".genbank");
    $options{base} = basename($options{input}, @suffix);

    my $seqio = new Bio::SeqIO(-format => 'genbank', -file => $options{input});
    my $count = 0;
    my $fasta = new Bio::SeqIO(-file => qq">${options{base}}.fasta", -format => 'fasta', -flush => 0);
    my $gffout = new Bio::Tools::GFF(-file => ">${options{base}}_all.gff", -gff_version => 3);
    my $gene_gff = new Bio::Tools::GFF(-file => ">${options{base}}.gff", -gff_version => 3);
    my $cds_gff = new Bio::Tools::GFF(-file => ">${options{base}}_cds.gff", -gff_version => 3);
    my $inter_gffout = new Bio::Tools::GFF(-file => ">${options{base}}_interCDS.gff", -gff_version => 3);
    my $cds_fasta = new Bio::SeqIO(-file => qq">${options{base}}_cds.fasta", -format => 'Fasta');
    my $pep_fasta = new Bio::SeqIO(-file => qq">${options{base}}_pep.fasta", -format => 'Fasta');
    while(my $seq = $seqio->next_seq) {
        $count++;
        $fasta->write_seq($seq);
        my @feature_list = ();
        # defined a default name

        foreach my $feature ($seq->top_SeqFeatures()) {
            $gffout->write_feature($feature);
            my $feat_count = 0;
        }
        for my $feat_object ($seq->get_SeqFeatures) {
            my $feat_count = 0;
            if ($feat_object->primary_tag eq "CDS") {
                $cds_gff->write_feature($feat_object);
                $feat_count++;
                my $start = $feat_object->start;
                my $end = $feat_object->end;
                my $len = $feat_object->length;
                my $id = $feat_object->{_gsf_tag_hash}->{locus_tag}->[0];
                my $seq = $feat_object->spliced_seq->seq;
                my $pep = $feat_object->spliced_seq->translate->seq;
                my $size = {
                    id => $id,
                    start => $start,
                    end => $end,
                    len => $len,
                    feat => $feat_object,
                };
                push(@feature_list, $size);
                my $seq_object = new Bio::PrimarySeq(-id => $id, -seq => $seq);
                my $pep_object = new Bio::PrimarySeq(-id => $id, -seq => $pep);
                $cds_fasta->write_seq($seq_object);
                ##print "TESTME: $seq_object\n";
                my $ttseq = $seq_object->seq;
                ##print "TESTME: $id $ttseq\n";
                $pep_fasta->write_seq($pep_object);
            } elsif ($feat_object->primary_tag eq "gene") {
                $gene_gff->write_feature($feat_object);
            }
        } ## End looking at every feature and putting them into the @feature_list
        ## Now make a hash from it and fill in the inter-cds data

        my %inter_features = ();
        my ($p_st, $p_en, $c_st, $c_en, $inter_start, $inter_end) = 0;
      INTER: for my $c (0 .. $#feature_list) {
          my $p_size = {};
          if ($c == 0) {
              $p_size = $feature_list[$#feature_list];
          } else {
              $p_size = $feature_list[$c - 1];
          }

          my $c_size = $feature_list[$c];
          next INTER if (!defined($p_size->{start}));
          next INTER if (!defined($p_size->{end}));
          next INTER if (!defined($c_size->{start}));
          next INTER if (!defined($c_size->{end}));

          my $p_start = $p_size->{start};
          my $p_end = $p_size->{end};
          my $c_start = $c_size->{start};
          my $c_end = $c_size->{end};
          my $interstart = 0;
          my $interend = 0;

          if ($p_start > $p_end) {
              $interstart = $p_start;
          } else {
              $interstart = $p_end;
          }
          if ($c_start < $c_end) {
              $interend = $c_start;
          } else {
              $interend = $c_end;
          }
          my $tmp_start = $interstart;
          my $tmp_end = $interend;
          if ($tmp_start > $tmp_end) {
              $interstart = $tmp_end;
              $interend = $tmp_start;
          }
          ## I need a little logic to catch the first feature
          ## As it is currently written, it grabs the position of the last feature
          ## and uses that number rather than the beginning of the chromosome.
          ## Therefore we end up with a first feature which contains the entire chromosome.
          if ($c <= 2) {  ## Arbitrarily set it to check the first 2 features
              if ($interend >= 100000) {  ## If the first feature take >= 1/2 the chromosome
                  $interend = $interstart;
                  $interstart = 1;
              }
          }

          my $new_feature = $c_size->{feat};
          my $inter_location = new Bio::Location::Atomic(-start => $interstart, -end => $interend, -strand => 1);
          $new_feature->{_location} = $inter_location;
          $inter_gffout->write_feature($new_feature);
      }
    }
}

1;
