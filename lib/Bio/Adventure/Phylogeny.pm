package Bio::Adventure::Phylogeny;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use File::Basename;

=head2 C<Run_Gubbins>

 Invoke gubbins to search for optimal trees given an outgroup, starting tree, and MSA.
 10.1093/nar/gku1196

=over

=item C<Arguments>

 input(required): Input fasta file.
 outgroup(required): Name of an outgroup taxon.
 starting_tree(required) Name of start tree.

=cut
sub Run_Gubbins {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        jcpu => 8,
        modules => ['gubbins'],
        required => ['input', 'outgroup', 'starting_tree'],);
    my $bname = basename($options->{input});
    my $gub_dir = qq"outputs/gubbins_${bname}";
    $gub_dir = $options->{gub_dir} if ($options->{gub_dir});
    my $jstring = qq!mkdir -p ${gub_dir} && run_gubbins.py \\
  $options->{input} \\
  --threads $options->{jcpu} \\
!;
    if ($options->{outgroup}) {
        $jstring .= qq! --outgoup $options->{outgroup} \\!;
    }
    if ($options->{starting_tree}) {
        $jstring .= qq! --starting_tree $options->{starting_tree} \\!;
    }
    $jstring .= qq"  2>${gub_dir}/gubbins.stderr 1>${gub_dir}/gubbins.stdout";
    my $gubbins = $class->Submit(
        basename => 'gubbins',
        comment => 'This should remove recombination events from bacterial genomes.',
        jdepends => $options->{jdepends},
        job_output => 'filtered_polymorphic_sites.fasta',
        jqueue => 'large',
        jprefix => '30',
        jmem => 50,
        modules => $options->{modules},
        jcpu => $options->{jcpu},
        jstring => $jstring,
        jname => 'gub',);
    return($gubbins);
}

=back

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

L<gubbins>

=cut


1;
