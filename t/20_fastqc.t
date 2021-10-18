# -*-Perl-*-
use Test::More qw"no_plan";
use Bio::Adventure;
use File::Path qw"remove_tree";
use String::Diff qw" diff_fully diff diff_merge diff_regexp ";
use File::Temp qw"tempdir";
use File::ShareDir qw":ALL";
use Cwd;

## Create a temporary directory for working with this data
my $start = getcwd();
##my $new = tempdir(CLEANUP => 0, TEMPLATE => 'test_XXXX');
my $new = 'test_output';
mkdir($new);
chdir($new);

my $input_file = dist_file('Bio-Adventure', 'test_forward.fastq.gz');
my $cyoa = Bio::Adventure->new(cluster => 0, basedir => cwd());

ok($cyoa->Bio::Adventure::QA::Fastqc_Single(input => $input_file),
   'Run Fastqc');

ok(-r 'scripts/01fqc_test_forward_Bio-Adventure.sh',
   'Fastqc script exists?');

ok(qx"fastqc --help", 'Can run fastqc?');

ok(my $actual = $cyoa->Last_Stat(input => 'outputs/fastqc_stats.csv'),
   'Collect Fastqc Statistics');

my $expected = 'fqc_test_forward_Bio-Adventure,10000,0,pass,warn,pass,pass,pass,warn,fail,0';
unless(ok($expected eq $actual,
          'Are the fastqc results the expected value?')) {
    my($old, $new) = diff($expected, $actual);
    diag("$old\n$new\n");
}

## Go back to the top level.
chdir($start);
