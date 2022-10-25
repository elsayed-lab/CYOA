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

Bio::Adventure::Cleanup - Delete the various directories/files created by the
tools invoked by Bio::Adventure.

=head1 SYNOPSIS

use Bio::Adventure;
my $hpgl = new Bio::Adventure;
$hpgl->Cleanup();

=head1 METHODS

=head2 C<Cleanup>

Delete the various directories/files from the RNASeq tools invoked by Bio::Adventure.pm.

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

sub Cleanup_Phage_Assembly {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        jprefix => '99',
        jname => 'cleanup_phage');

    my $job_name = $class->Get_Job_Name();
    my $input_paths = $class->Get_Paths($options->{input});

    my $jstring = qq!
## Rando fastq files
rm -f \$(find . -type f -name '*.fastq')
## Core dumps
rm -f \$(find . -type f -name core)
## Hisat indexes
rm -f \$(find . -type f -name '*.ht2')
## tmp files from the various prediction tools
rm -f \$(find . -type f -name '*.tmp.*')
## Trinotate junk
rm -r \$(find . -type d -name '*_dir*')
## Recompress random fastq files
xz -9e -f \$(find . -name '*.fastq')
!;

    my $comment = '## Cleanup some of the mess.';
    my $clean = $class->Submit(
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => $options->{jname},
        jprefix => $options->{jprefix},
        jstring => $jstring,);
    return($clean);
}

=head1 AUTHOR - atb

Email <abelew@gmail.com>

=head1 SEE ALSO

    L<Bio::Adventure>

=cut

1;
