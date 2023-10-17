# -*-Perl-*-
use strict;
use Test::More qw"no_plan";
use Bio::Adventure;
use Cwd;
use File::Basename qw"basename";
use File::Copy qw"cp mv";
use File::Path qw"remove_tree make_path rmtree";
use File::ShareDir qw"dist_file module_dir dist_dir";
use String::Diff qw"diff";
use Test::File::ShareDir::Dist { 'Bio-Adventure' => 'share/' };
my $start_dir = dist_dir('Bio-Adventure');
my $input_forward = qq"${start_dir}/test_forward.fastq.gz";
my $input_forward_local = basename($input_forward);

my $start = getcwd();
my $new = 'test_output';
mkdir($new);
chdir($new);
if (!-r $input_forward_local) {
    ok(cp($input_forward, $input_forward_local), 'Copied r1.');
}

my $cyoa = Bio::Adventure->new(basedir => cwd());
my $trimmer = $cyoa->Bio::Adventure::Trim::Trimomatic_Single(input => $input_forward_local);
ok($trimmer, 'Submitted trimomatic SE job.');
my $status = $cyoa->Wait(job => $trimmer);
ok($status->{State} eq 'COMPLETED', 'The trimomatic jobs completed.');
my $csv_file = $trimmer->{stats}->{output};
ok(my $actual = $cyoa->Last_Stat(input => $csv_file),
   'Collect Trimomatic Statistics');
my $expected = 'test_forward,10000,9390,610';
unless(ok($expected eq $actual,
          'Are the trimomatic results the expected value?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

chdir($start);
