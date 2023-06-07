# -*-Perl-*-
use strict;
use Cwd;
use File::Copy qw"cp mv";
use File::Path qw"remove_tree make_path rmtree";
use File::Path qw"remove_tree make_path rmtree";
use File::ShareDir qw"dist_file module_dir dist_dir";
use File::Temp qw"tempdir";
use String::Diff qw" diff_fully diff diff_merge diff_regexp ";
use Test::More qw"no_plan";
use Test::File::ShareDir::Dist { 'Bio-Adventure' => 'share/' };

use Bio::Adventure;
my $start_dir = dist_dir('Bio-Adventure');
my $input_file = qq"${start_dir}/test_forward.fastq.gz";

## Create a temporary directory for working with this data
my $start = getcwd();
##my $new = tempdir(CLEANUP => 0, TEMPLATE => 'test_XXXX');
my $new = 'test_output';
mkdir($new);
chdir($new);

my $cyoa = Bio::Adventure->new(cluster => 0, basedir => cwd());

ok($cyoa->Bio::Adventure::QA::Fastqc(input => $input_file),
   'Run Fastqc');

ok(-r 'scripts/01fqc_test_forward_Bio-Adventure.sh',
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
