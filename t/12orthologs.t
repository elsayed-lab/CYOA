# -*-Perl-*-
use Test::More qw"no_plan";
use Bio::Adventure;
use Cwd;
use File::Copy qw"cp";
use File::Path qw"remove_tree make_path rmtree";
use File::ShareDir qw"dist_file";
use String::Diff qw"diff";

my $start = getcwd();
my $new = 'test_output';
mkdir($new);
chdir($new);

my $groupa_strep = dist_file('Bio-Adventure', 'genome/mgas_5005.gb.xz');
my $groupb_strep = dist_file('Bio-Adventure', 'genome/sagalactiae_cjb111.gb.xz');
my $cyoa = Bio::Adventure->new(
    cluster => 0,
    basedir => cwd(),
    libdir => cwd());

## First extract the amino acid sequences of our two streptococci genbank files.
my $groupa_convert = $cyoa->Bio::Adventure::Convert::Gb2Gff(
    input => $groupa_strep, output_dir => '.',);
ok($groupa_convert, 'Converted the group A strep genbank to fasta/gff/etc.');
my $groupb_convert = $cyoa->Bio::Adventure::Convert::Gb2Gff(
    input => $groupb_strep, output_dir => '.',);
ok($groupb_convert, 'Converted the group B strep genbank to fasta/gff/etc.');

## Now invoke orthofinder
my $orthos = $cyoa->Bio::Adventure::Align::OrthoFinder(
    input => qq"$groupa_convert->{output_pep_fasta}:$groupb_convert->{output_pep_fasta}");
