package Bio::Adventure::Status;
## LICENSE: gplv2
## ABSTRACT:  Kitty!
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";

=head2
    Finishedp_Fastqc()
=cut
sub Finishedp_Fastqc {
    my ($class, %args) = @_;
    if (-d qq"$class->{hpgl}-unfiltered_fastqc" or -d qq"$class->{hpgl}-filtered_fastqc") {
	return(1);
    } else {
	return(0);
    }
}

=head2
    Finishedp_Trimomatic()
=cut
sub Finishedp_Trimomatic {
    my ($class, %args) = @_;
    if (-r qq"$class->{basename}-trimmed.fastq" or (-r qq"$class->{basename}_forward-trimmed.fastq" and -r qq"$class->{basename}_reverse-trimmed.fastq")) {
	return(1);
    } else {
	return(0);
    }
}

=head2
    Finishedp_Graphs()
=cut
sub Finishedp_Graphs {
    my ($class, %args) = @_;
    if (-r qq"$class->{basename}_ntdist.svg") {
	return(1);
    } else {
	return(0);
    }
}

=head2
    Finishedp_Tophat()
=cut
sub Finishedp_Tophat {
    my ($class, %args) = @_;
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
    my ($class, %args) = @_;
    my $count_table = qq"$class->{hpgl}.count";
    if ($args{species}) {
	$count_table = qq"$class->{hpgl}_$args{species}.count";
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
    my ($class, %args) = @_;
    if (-r qq"$class->{hpgl}_forward.fastq" && -r qq"$class->{hpgl}_reverse.fastq") {
	return(1);
    } else {
	return(0);
    }
}

1;
