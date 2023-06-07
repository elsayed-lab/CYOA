# -*-Perl-*-
use strict;
use Test::More qw"no_plan";
use Bio::Adventure;
use Cwd;
use File::Copy qw"cp mv copy move";
use File::Path qw"remove_tree make_path rmtree";
use File::ShareDir qw"dist_file module_dir dist_dir";
use String::Diff qw"diff";
use Test::File::ShareDir::Dist { 'Bio-Adventure' => 'share/' };
my $start_dir = dist_dir('Bio-Adventure');
my $input = 'consolidated.fastq.xz';
my $input_file = qq"${start_dir}/${input}";
my $phix_fasta = qq"${start_dir}/genome/phix.fastq";
my $phix_gff = qq"${start_dir}/genome/phix.gff";

my $start = getcwd();
my $new = 'test_output';
mkdir($new);
chdir($new);
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
