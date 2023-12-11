package Bio::Adventure::Visualization;
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

 Invoke cgview on an input genbank file and make a pretty picture.
 10.1093/nar/gkn179

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
        jprefix => '16',
        linear => 1,
        reading_frame => 1,
        gene_labels => 0,
        feature_labels => 1,
        orfs => 1,
        imagemap => 0,);
    my $job_name = $class->Get_Job_Name();
    ## Reminder, this gives me: filename, directory, dirname, fullpath
    my $job_paths = $class->Get_Paths($options->{input});
    my $job_filename = $job_paths->[0]->{filename};
    my $output_directory = qq"outputs/$options->{jprefix}cgview";
    my $output_base = basename($job_filename, ('.gbk'));
    my $xml_input = qq"$options->{input}";
    my $xml_output = qq"${output_directory}/${output_base}.xml";
    my $output_file = qq"${output_directory}/${output_base}";

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
    my $comment = '## Run cgview given an assembly genbank file.';
    my $jstring = qq!
mkdir -p ${output_directory}
cgview_xml_builder.pl -sequence ${xml_input} \\
  -output ${xml_output} \\
  ${option_string} \\
  2>${output_directory}/cgview.stderr 1>${output_directory}/cgview.stdout
if [[ \$? \!= "0" ]]; then
  echo "Unable to create initial xml file."
  exit 0
fi

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
        comment => $comment,
        jcpu => 1,
        jdepends => $options->{jdepends},
        jname => qq"cgview_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        output => $output_file,
        output_xml => $xml_output,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($cgview);
}

1;

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

=cut
