package Bio::Adventure::Phylogeny;
## LICENSE: gplv2
## ABSTRACT:  Kitty!
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use File::Basename;

=head2 C<Run_Gubbins>

Invoke gubbins to search for optimal trees given an outgroup, starting tree, and MSA.

=over

=item I<input> * Input fasta file.

=item I<outgroup> * Name of an outgroup taxon.

=item I<starting_tree> * Name of start tree.

=back

=cut
sub Run_Gubbins {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'outgroup', 'starting_tree'],
        modules => ['gubbins'],
        cpus => 32,);
    my $cpus = $options->{cpus};
    my $bname = basename($options->{input});
    my $gub_dir = qq"outputs/gubbins_${bname}";
    $gub_dir = $options->{gub_dir} if ($options->{gub_dir});
    my $jstring = qq!mkdir -p ${gub_dir} && run_gubbins.py \\
  $options->{input} \\
  --threads ${cpus} !;
    if ($options->{outgroup}) {
        $jstring .= qq! --outgoup $options->{outgroup} !;
    }
    if ($options->{starting_tree}) {
        $jstring .= qq! --starting_tree $options->{starting_tree} !;
    }

    my $gubbins = $class->Submit(
        basename => 'gubbins',
        comment => 'This should remove recombination events from bacterial genomes.',
        jdepends => $options->{jdepends},
        job_output => 'filtered_polymorphic_sites.fasta',
        jqueue => 'large',
        jprefix => '30',
        jmem => 50,
        modules => $options->{modules},
        cpus => $options->{cpus},
        jstring => $jstring,
        jname => 'gub',);
    return($gubbins);
}

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

L<gubbins>

=cut


1;
