
=head2
    Finishedp_Fastqc()
=cut
sub Finishedp_Fastqc {
    my $me = shift;
    my %args = @_;
    if (-d qq"$me->{hpgl}-unfiltered_fastqc" or -d qq"$me->{hpgl}-filtered_fastqc") {
	return(1);
    } else {
	return(0);
    }
}

=head2
    Finishedp_Trimomatic()
=cut
sub Finishedp_Trimomatic {
    my $me = shift;
    my %args = @_;
    if (-r qq"$me->{basename}-trimmed.fastq" or (-r qq"$me->{basename}_forward-trimmed.fastq" and -r qq"$me->{basename}_reverse-trimmed.fastq")) {
	return(1);
    } else {
	return(0);
    }
}

=head2
    Finishedp_Graphs()
=cut
sub Finishedp_Graphs {
    my $me = shift;
    my %args = @_;
    if (-r qq"$me->{basename}_ntdist.svg") {
	return(1);
    } else {
	return(0);
    }
}

=head2
    Finishedp_Tophat()
=cut
sub Finishedp_Tophat {
    my $me = shift;
    my %args = @_;
    my $subdir = ".";
    $subdir = $args{species} if ($args{species});
    if (-r qq"$subdir/tophat_out/accepted_hits.bam") {
	return(1);
    } else {
	return(0);
    }
}

=head2
    Finishedp_HTseq()
=cut
sub Finishedp_Htseq {
    my $me = shift;
    my %args = @_;
    my $count_table = qq"$me->{hpgl}.count";
    if ($args{species}) {
	$count_table = qq"$me->{hpgl}_$args{species}.count";
    }
    if (-r $count_table) {
	return(1);
    } else {
	return(0);
    }
}

=head2
    Finishedp_Rawseq()
=cut
sub Finishedp_Rawseq {
    my $me = shift;
    my %args = @_;
    if (-r qq"$me->{hpgl}_forward.fastq" && -r qq"$me->{hpgl}_reverse.fastq") {
	return(1);
    } else {
	return(0);
    }
}

1;
