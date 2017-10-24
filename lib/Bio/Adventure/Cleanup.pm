package Bio::Adventure::Cleanup;
## LICENSE: gplv2
## ABSTRACT:  Kitty!
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

=head1 NAME
    Bio::Adventure::Cleanup - Delete the various directories/files created by
    the RNASeq tools invoked by Bio::Adventure.pm

=head1 SYNOPSIS

    use Bio::Adventure;
    my $hpgl = new Bio::Adventure;
    $hpgl->Cleanup();

=over4

=item C<Cleanup>

    delete the various directories/files from the RNASeq tools invoked by Bio::Adventure.pm

=cut
sub Cleanup {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars();
    my $depends = $options->{depends};
    my $input_files = $class->{input};
    $input_files =~ s/:/ /g;
    my $unzipped = $input_files;
    $unzipped =~ s/\.gz//g;
    my $trimmed_input = $unzipped;
    $trimmed_input =~ s/\.fastq/-trimmed\.fastq/g;
    my $trimmed_paired = $unzipped;
    $trimmed_paired =~ s/\.fastq/-trimmed_paired\.fastq/g;
    my $trimmed_unpaired = $unzipped;
    $trimmed_unpaired =~ s/\.fastq/-trimmed_unpaired\.fastq/g;
    my $jstring = qq!rm -rf outputs scripts sequences $input_files $unzipped $trimmed_input $trimmed_paired $trimmed_unpaired!;
    print "Execute: $jstring\n";
    return(1);
}

=back

=head1 AUTHOR - atb

Email <abelew@gmail.com>

=head1 SEE ALSO

    L<Bio::Adventure>

=cut

1;
