package Bio::Adventure::Visualization;
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

use Capture::Tiny qw":all";
use Cwd qw"abs_path getcwd cwd";
use File::Basename;
use File::Spec;
use File::Path qw"make_path";
use File::Which qw"which";
use File::ShareDir qw":ALL";

=head2 C<CGView>

Invoke cgview on an input genbank file in order to create a pretty
picture of the genome.

cgview has a few options which are of interest:

linear: distinguish between a linear and circular genome. (linear)
reading_frame: Show the reading frames (true)
orfs: Print the ORFs (true)
feature_labels: Print the features (true)
gene_labels: Print the genes. (false)

This invocation of cgview has it print both a svg and png copy of the sequence.

=cut
sub CGView {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input',],
        modules => ['cgview'],
        jprefix => '16',
        linear => 1,
        reading_frame => 1,
        gene_labels => 0,
        feature_labels => 1,
        orfs => 1,
        imagemap => 0,);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('cgview');
    die("Could not find cgview in your PATH.") unless($check);

    my $job_name = $class->Get_Job_Name();
    ## Reminder, this gives me: filename, directory, dirname, fullpath
    my $job_paths = $class->Get_Paths($options->{input});
    my $output_directory = qq"outputs/$options->{jprefix}cgview";
    my $output_base = basename($job_paths->{filename}, ('.gbk'));
    my $xml_input = qq"$options->{input}";
    my $xml_output = qq"${output_directory}/${output_base}.xml";
    my $output_file = qq"${output_directory}/${output_base}";
    my $out_paths = $class->Get_Paths($output_file);

    ## Some potential options when making the input xml file:
    ## -linear add a glyph for linear genomes
    ## -reading_frames Separate ORFs by frame
    ## A bunch more
    my $option_string = '-size large-v2 ';
    if ($options->{linear}) {
        $option_string .= '-linear T ';
    }
    if ($options->{reading_frame}) {
        $option_string .= '-reading_frame T ';
    }
    if ($options->{orfs}) {
        $option_string .= ' -orfs T ';
    }
    if ($options->{feature_labels}) {
        $option_string .= '-feature_labels T ';
    }
    if ($options->{gene_labels}) {
        $option_string .= '-gene_labels T ';
    }

    my $jstring = qq!
mkdir -p ${output_directory}
cgview_xml_builder.pl -sequence ${xml_input} \\
  -output ${xml_output} \\
  ${option_string} \\
  2>${output_directory}/cgview.stderr 1>${output_directory}/cgview.stdout

cgview -i ${xml_output} \\
  -o ${output_file}.png -f png \\
  2>>${output_directory}/cgview.stderr 1>>${output_directory}/cgview.stdout
cgview -i ${xml_output} \\
  -o ${output_file}.svg -f svg \\
  2>>${output_directory}/cgview.stderr 1>>${output_directory}/cgview.stdout
!;
    if ($options->{imagemap}) {
        $jstring = qq!
${jstring}
cgview -i ${xml_output} \\
  -s ${output_directory}/imagemap \\
  2>>${output_directory}/cgview.stderr 1>>${output_directory}/cgview.stdout
!;
    }

    my $cgview = $class->Submit(
        jdepends => $options->{jdepends},
        jname => qq"cgview_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        modules => $options->{modules},
        output => $output_file,
        xml_output => $xml_output,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($cgview);
}

1;

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

=cut
