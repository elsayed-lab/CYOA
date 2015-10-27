
=head2
    HT_Multi()
=cut
sub HT_Multi {
    my $me = shift;
    my %args = @_;
    my @jobs = ();
    $me->Check_Options(["species"]);
    my $basename = $me->{basename};
    my $species = $me->{species};
    my %gff_types = (
        misc => "transcript",
        rintron => "exon",
        nmd => "exon",
        linc => "gene",
        pseudo => "exon",
        antisense => "exon",
        exon => "exon",
        fiveputr => "five_prime_UTR",
        threeputr => "three_prime_UTR",
        sn => "exon",
        mi => "exon",
        sno => "exon",
        rrna => "gene",
        interCDS => "CDS",
        operons => "gene",
    );
    my $htseq_runs = 0;
    my $input = $args{input};
##    my $aligntype = $args{aligntype};
    my $suffix = $args{suffix};
    my $out_prefix = $input;
    $out_prefix =~ s/\.bam$//g;
    foreach my $type (keys %gff_types) {
        my $gff = qq"$me->{libdir}/genome/${species}_${type}.gff";
        my $output = qq"${out_prefix}_${type}.count";
        if (-r "$gff") {
            print "Found $gff, performing htseq with it.\n";
            my $ht = $me->HTSeq(htseq_type => $gff_types{$type},
                                htseq_gff => $gff,
                                jobname => "ht${type}${suffix}",
                                depends => $args{depends},
                                input => $input,
                                output => $output,
				prescript => $args{prescript},
				postscript => $args{postscript},
                                suffix => $args{suffix},
                );
            push(@jobs, $ht);
            $htseq_runs++;
        } else {
            print "Did not find $gff, skipping.\n";
        }
    } ## End foreach type

    ## Also perform a whole genome count
    my $gff = qq"$me->{libdir}/genome/${species}.gff";
    my $output = $input;
    $output =~ s/\.bam$/\.count/g;
    if (-r "$gff") {
        my $ht = $me->HTSeq(htseq_type => "all",
                            htseq_gff => $gff,
                            jobname => "htall${suffix}",
                            depends => $args{depends},
                            input => $input,
                            output => $output,
                            prescript => $args{prescript},
                            postscript => $args{postscript},
            );
        push(@jobs, $ht);
    }
    return(\@jobs);
}

=head2
    HTSeq()
=cut
sub HTSeq {
    my $me = shift;
    my %args = @_;
    $me->Check_Options(["species", "htseq_stranded", "htseq_identifier"]);
    my $basename = $me->{basename};
    my $species = $me->{species};
    my $stranded = $me->{htseq_stranded};
    my $identifier = $me->{htseq_identifier};
    my $htseq_jobname = qq"htseq";
    $htseq_jobname = $args{jobname} if ($args{jobname});
    my $depends = "";
    $depends = $args{depends} if ($args{depends});
    my $libtype = "genome";
    $libtype = $args{libtype} if ($args{libtype});
    my $type = "exon";
    $type = $args{htseq_type} if ($args{htseq_type});
    my $gff = qq"$me->{libdir}/${libtype}/${species}.gff";
    $gff = $args{htseq_gff} if ($args{htseq_gff});
    my $input = "${basename}.bam";
    $input = $args{input} if ($args{input});
    my $output = $input;
    if ($args{suffix}) {
        my $suffix = $args{suffix};
        $output =~ s/\.bam/_${suffix}\.count/g;
    } else {
        $output =~ s/\.bam/\.count/g;
    }
    $output = $args{output} if ($args{output});
    my $error = $input;
    $error =~ s/\.bam/\.error/g;
    if (!-r "$gff") {
        print "Unable to read $gff, please fix this and try again.\n";
    }
    $type = 'gene' if ($type eq 'all');
    my $job_string = qq!htseq-count -q -f bam -s $stranded -t $type -i $identifier $input $gff 1>$output 2>$error && pxz $output
!;
    my $comment = qq!## Counting the number of hits in $input for each feature found in $gff
## Is this stranded? $stranded.  The attribute type (3rd column of the gff file) is $type
## and the comment field (10th column) used to name the feature is $identifier!;
    my $htseq = $me->Qsub(job_name => $htseq_jobname,
                          qsub_mem => 6,
                          depends => $depends,
                          job_string => $job_string,
                          comment => $comment,
                          input => $input,
                          output => $output,
			  prescript => $args{prescript},
			  postscript => $args{postscript},
        );
    return($htseq);
}

1;
