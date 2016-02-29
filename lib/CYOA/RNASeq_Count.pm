package CYOA;

=head1 NAME
    CYOA::RNASeq_Count - Perform Sequence alignments counting with HTSeq

=head1 SYNOPSIS

    use CYOA;
    my $hpgl = new CYOA;
    $hpgl->HT_Multi();

=head2 Methods

=over 2

=item C<HT_Multi>

    some stuff here

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
                                jobname => "ht${type}${suffix}_${species}",
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
                            jobname => "htall_${species}",
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

sub HTSeq {
    my $me = shift;
    my %args = @_;
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
    my $multi_test = 0;
    foreach my $type (keys %gff_types) {
        my $gff = qq"$me->{libdir}/genome/$me->{species}_${type}.gff";
        if (-r "$gff") {
            $multi_test++;
        }
    }
    my $jobs = undef;
    if ($multi_test) {
        ## $jobs = $me->HT_Multi(%args);
        print "HT_Multi needs to be rethought.  Using single for now.\n";
        $jobs = $me->HT_Single(%args);
    } else {
        $jobs = $me->HT_Single(%args);
    }
    return($jobs);
}

=head2
    HTSeq_Single()
=cut
sub HT_Single {
    my $me = shift;
    my %args = @_;
    $me->Check_Options(["species", "htseq_stranded", "htseq_args"]);
    my $basename = $me->{basename};
    my $species = $me->{species};
    my $stranded = $me->{htseq_stranded};
    my $htseq_args = "";
    if (defined($me->{htseq_args}->{$species})) {
        $htseq_args = $me->{htseq_args}->{$species};
    } else {
        $htseq_args = $me->{htseq_args}->{all};
    }
    my $htseq_jobname = qq"hts_${species}";
    $htseq_jobname = $args{jobname} if ($args{jobname});
    my $depends = "";
    $depends = $args{depends} if ($args{depends});
    my $libtype = "genome";
    $libtype = $args{libtype} if ($args{libtype});
    my $gff = qq"$me->{libdir}/${libtype}/${species}.gff";
    $gff = $args{htseq_gff} if ($args{htseq_gff});
    my $gtf = $gff;
    $gtf =~ s/\.gff/\.gtf/g;
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
    
    if (!-r "${gff}" and !-r "${gtf}") {
        die("Unable to read ${gff} nor ${gtf}, please fix this and try again.\n");
    }
    my $annotation = $gtf;
    if (!-r "${gtf}") {
        $annotation = $gff;
    }
    my $job_string = qq!htseq-count -q -f bam -s ${stranded} ${htseq_args} \\
  ${input} ${annotation} \\
  1>${output} 2>${error} && \\
    xz -9e ${output}
!;
    my $comment = qq!## Counting the number of hits in ${input} for each feature found in ${annotation}
## Is this stranded? ${stranded}.  The defaults of htseq are:
## $me->{htseq_args}->{default}
!;
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

=back

=head1 AUTHOR - atb

Email <abelew@gmail.com>

=head1 SEE ALSO

    L<htseq-count>

=cut

1;
