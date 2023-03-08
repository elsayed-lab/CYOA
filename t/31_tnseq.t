# -*-Perl-*-
use Test::More qw"no_plan";
use Bio::Adventure;
use Cwd;
use File::Copy qw"copy move";
use File::Path qw"remove_tree make_path rmtree";
use File::ShareDir qw"dist_file";
use String::Diff qw"diff";

my $start = getcwd();
my $new = 'test_output';
mkdir($new);
chdir($new);
my $input = 'consolidated.fastq.xz';
my $input_file = dist_file('Bio-Adventure', $input);
if (!-r $input) {
    ok(copy($input_file, $input), 'Copying data.');
}
make_path('genome/indexes'); ## Make a directory for the phix indexes.

my $genome = 'sagalactiae_cjb111';
my $cyoa = Bio::Adventure->new(
    cluster => 0,
    basedir => cwd(),
    libdir => cwd(),
    species => $genome,
    gff_tag => 'locus_tag',
    gff_type => 'gene',
    stranded => 'no',);

my $groupb_strep = dist_file('Bio-Adventure', qq"genome/${genome}.gb.xz");
if (!-r qq"genome/${genome}.fasta") {
    ok(copy($groupb_strep, qq"genome/${genome}.gb.xz"), 'Copying group b genome.');
}
my $groupb_convert = $cyoa->Bio::Adventure::Convert::Gb2Gff(
    input => qq"genome/${genome}.gb.xz", output_dir => 'genome',);
ok($groupb_convert, 'Converted the group B strep genbank to fasta/gff/etc.');
my $moved = move(qq"genome/${genome}.fsa", qq"genome/${genome}.fasta");

my $tpp = $cyoa->Bio::Adventure::TNSeq::Transit_TPP(
    input => $input,
    species => $genome,
    protocol => 'Sassetti',
    jprefix => '31',);
