use Test::More qw"no_plan";
use Bio::Adventure;
use File::Path qw"remove_tree";
use File::Copy qw"cp";
use String::Diff qw" diff ";

my $cyoa = Bio::Adventure->new();

ok(Bio::Adventure::RNASeq_Trim::Trimomatic_Single($cyoa, input => qq"share/test_forward.fastq.gz"),
   'Run Trimomatic');

##my $fqc_dir = qx"find outputs/";
##diag($fqc_dir);
##my $log = qx'cat outputs/test_forward-trimomatic.out';
##diag($log);

ok(my $actual = $cyoa->Last_Stat(input => 'outputs/trimomatic_stats.csv'),
   'Collect Trimomatic Statistics');

my $expected = 'test_forward,10000,9762,238';
unless(ok($expected eq $actual,
          'Are the trimomatic results the expected value?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}
