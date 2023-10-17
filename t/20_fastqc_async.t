# -*-Perl-*-
use strict;
use Test::More qw"no_plan";
use Bio::Adventure;
use Cwd;
use File::Basename qw"basename";
use File::Copy qw"cp mv";
use File::Path qw"remove_tree make_path rmtree";
use File::ShareDir qw"dist_file module_dir dist_dir";
use String::Diff qw" diff_fully diff diff_merge diff_regexp ";
use Test::File::ShareDir::Dist { 'Bio-Adventure' => 'share/' };

my $start_dir = dist_dir('Bio-Adventure');
my $input_file = qq"${start_dir}/test_forward.fastq.gz";
my $input_local = basename($input_file);

## Create a temporary directory for working with this data
my $start = getcwd();
my $new = 'test_output';
mkdir($new);
chdir($new);
if (!-r $input_local) {
    ok(cp($input_file, $input_local), 'Copying data.');
}

my $cyoa = Bio::Adventure->new(basedir => cwd());

my $fastqc = $cyoa->Bio::Adventure::QA::Fastqc(input => $input_local);
ok($fastqc, 'Run Fastqc');
my $status = $cyoa->Wait(job => $fastqc);
ok($status->{State} eq 'COMPLETED', 'Fastqc completed.);
ok(-r 'scripts/01fqc_test_forward.sh',
   'Fastqc script exists?');
ok(my $actual = $cyoa->Last_Stat(input => 'outputs/fastqc_stats.csv'),
   'Collect Fastqc Statistics');
my $expected = 'fqcst,10000,0,pass,warn,pass,pass,pass,warn,fail,0';
unless(ok($expected eq $actual,
          'Are the fastqc results the expected value?')) {
    my($old, $new) = diff($expected, $actual);
    diag("$old\n$new\n");
}

## Go back to the top level.
chdir($start);
