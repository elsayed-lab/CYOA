# -*-Perl-*-
use Test::More qw"no_plan";
use Bio::Adventure;
use Cwd;
use File::Path qw"remove_tree";
use File::Copy qw"cp";
use String::Diff qw" diff ";
use File::ShareDir qw"dist_file";

my $start = getcwd();
##my $new = tempdir(CLEANUP => 0, TEMPLATE => 'test_XXXX');
my $new = 'test_output';
mkdir($new);
chdir($new);
my $cyoa = Bio::Adventure->new(cluster => 0, basedir => cwd());

my $input_file = dist_file('Bio-Adventure', 'test_forward.fastq.gz');
my $trimmer = $cyoa->Bio::Adventure::Trim::Trimomatic_Single(input => $input_file);
ok($trimmer, 'Run Trimomatic');

my $csv_file = $trimmer->{stats}->{output};

ok(my $actual = $cyoa->Last_Stat(input => $csv_file),
   'Collect Trimomatic Statistics');

my $expected = 'test_forward,10000,9316,684';
unless(ok($expected eq $actual,
          'Are the trimomatic results the expected value?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

chdir($start);
