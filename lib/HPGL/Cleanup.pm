package HPGL;
use common::sense;
use autodie qw":all";

=head1 NAME
    HPGL::Cleanup - Delete the various directories/files created by
    the RNASeq tools invoked by HPGL.pm

=head1 SYNOPSIS

    use HPGL;
    my $hpgl = new HPGL;
    $hpgl->Cleanup();

=over4

=item C<Cleanup>

    delete the various directories/files from the RNASeq tools invoked by HPGL.pm

=cut
sub Cleanup {
    my $me = shift;
    my %args = @_;
    my $depends = "";
    if ($args{depends}) {
        $depends = $args{depends};
    }
    my $input_files = $me->{input};
    $input_files =~ s/:/ /g;
    my $unzipped = $input_files;
    $unzipped =~ s/\.gz//g;
    my $trimmed_input = $unzipped;
    $trimmed_input =~ s/\.fastq/-trimmed\.fastq/g;
    my $trimmed_paired = $unzipped;
    $trimmed_paired =~ s/\.fastq/-trimmed_paired\.fastq/g;
    my $trimmed_unpaired = $unzipped;
    $trimmed_unpaired =~ s/\.fastq/-trimmed_unpaired\.fastq/g;
    my $job_string = qq!rm -rf outputs scripts sequences $input_files $unzipped $trimmed_input $trimmed_paired $trimmed_unpaired!;
    my $clean = $me->Qsub(job_name => "clean",
                          qsub_queue => 'throughput',
                          depends => $depends,
                          job_string => $job_string,
			  prescript => $args{prescript},
			  postscript => $args{postscript},
        );
    return($clean);
}

=back

=head1 AUTHOR - atb

Email <abelew@gmail.com>

=head1 SEE ALSO

    L<HPGL>

=cut

1;
