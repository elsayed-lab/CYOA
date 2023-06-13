# -*-Perl-*-
use strict;
use Test::More qw"no_plan";
use Bio::Adventure;
use Cwd;
use File::Copy qw"cp mv";
use File::Path qw"remove_tree make_path rmtree";
use File::ShareDir qw"dist_file module_dir dist_dir";
use String::Diff qw"diff";
use Test::File::ShareDir::Dist { 'Bio-Adventure' => 'share/' };
my $start_dir = dist_dir('Bio-Adventure');
my $input_file = qq"${start_dir}/test_forward.fastq.gz";

my $start = getcwd();
##my $new = tempdir(CLEANUP => 0, TEMPLATE => 'test_XXXX');
my $new = 'test_output';
mkdir($new);
chdir($new);
my $cyoa = Bio::Adventure->new(cluster => 0, basedir => cwd());

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
