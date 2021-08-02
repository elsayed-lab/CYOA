package Bio::Adventure::PopGen;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use File::Basename;
use Math::Round qw":all";

=head1 NAME

Bio::Adventure::PopGen - Invoke some tools used for high-throughput population genetics.
Because I do not fully understand the interplay between the various SNP tools and
population genetics, this may contain some inappropriate ideas.

For the moment, my goal is to get a set of inputs ready for plink, run plink, and
do some downstream analyses.

=head1 SYNOPSIS

=head1 METHODS

=head2 C<Angsd

=cut
sub Angsd_Filter {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        modules => ['angsd'],
        required => ['input'],
        pval => 0.01,);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $input_paths = $class->Get_Paths($options->{input});
    my $output_dir = qq"outputs/$options->{jprefix}angsd_$input_paths->{dirname}";
    my $jstring = qq?
angsd -b $options->{input} -doHWE 1 -GL 1 -doMajorMinor 1 -doMaf 2 -snp_pval 1e-2 -P 5 -dosnpstat 1 -out ${output_dir}
?;
    my $angsd_filter = $class->Submit(
        jname => qq"angsdfilter_$input_paths->{filename}",
        jprefix => $options->{jprefix},
        jstring => $jstring,);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($angsd_filter);
}

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

=cut

1;
