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

my $input_r1 = dist_file('Bio-Adventure', 'r1.fastq.xz');
ok(cp($input_r1, 'r1.fastq.xz'), 'Copying r1.');
my $input_r2 = dist_file('Bio-Adventure', 'r2.fastq.xz');
ok(cp($input_r2, 'r2.fastq.xz'), 'Copying r2.');

my $cyoa = Bio::Adventure->new(cluster => 0, basedir => cwd());
my $assemble = $cyoa->Bio::Adventure::Pipeline::Phage_Assemble(
    input => 'r1.fastq.xz:r2.fastq.xz',);
