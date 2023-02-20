# -*-Perl-*-
use Test::More qw"no_plan";
use Bio::Adventure;
use Cwd;
use File::Copy qw"cp";
use File::Path qw"remove_tree";
use File::ShareDir qw"dist_file";
use String::Diff qw"diff";

my $start = getcwd();
##my $new = tempdir(CLEANUP => 0, TEMPLATE => 'test_XXXX');
my $new = 'test_output';
mkdir($new);
chdir($new);
my $cyoa = Bio::Adventure->new(cluster => 0, basedir => cwd());

my $input_file = dist_file('Bio-Adventure', 'test_forward.fastq.gz');
my $trimmer = $cyoa->Bio::Adventure::Trim::Cutadapt(input => $input_file);
ok($trimmer, 'Run Trimomatic');

chdir($start);
