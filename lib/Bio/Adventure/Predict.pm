package Bio::Adventure::Predict;
## LICENSE: gplv2
## ABSTRACT:  Kitty!
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';
use feature 'try';
no warnings 'experimental::try';

use Bio::SeqFeature::Generic;
use Capture::Tiny qw":all";
use Cwd qw"abs_path getcwd cwd";
use File::Basename;
use File::Spec;
use File::Path qw"make_path";
use File::Which qw"which";
use File::ShareDir qw":ALL";

sub RhoTermPredict {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input',],
        modules => 'rhotermpredict',
        jprefix => '51',
        );
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $input_paths = $class->Get_Paths($options->{input});
    my $input_full = $input_paths->{fullpath};
    my $output_dir = qq"outputs/$options->{jprefix}rhotermpredict_$input_paths->{dirname}";
    my $output_file = qq"${output_dir}/predictions_coordinates_seqname.csv";
    my $info_file = qq"${output_dir}/info_about_predictions_seqname.csv";
    my $jstring = qq?mkdir -p ${output_dir}
start=\$(pwd)
cd ${output_dir}
cp $options->{$input} .
echo $input_paths->{filename} | RhoTermPredict_algorithm.py

?;

    my $rhoterm = $class->Submit(
        jdepends => $options->{jdepends},
        jname => 'rhotermpredict',
        jprefix => $options->{jprefix},
        jstring => $jstring,
        output_file => $output_file,
        output_info => $info_file,);
    return($rhoterm);
}

1;
