# -*-Perl-*-
use Test::More qw"no_plan";
use Bio::Adventure;
use File::Path qw"remove_tree";
use File::Copy qw"cp";

my $cyoa = Bio::Adventure->new();

ok(cp("t/data/test_forward.fastq.gz", "test_forward.fastq.gz"),
    'Copying data.');

ok(Bio::Adventure::RNASeq_Trim::Trimomatic_Single($cyoa, input => qq"test_forward.fastq.gz"),
   'Run Trimomatic');

ok(my $actual = $cyoa->Last_Stat(input => 'outputs/trimomatic_stats.csv'),
   'Collect Trimomatic Statistics');

my $expected = 'test_forward,10000,9762,238';
unless(ok($expected eq $actual,
          'Are the trimomatic results the expected value?')) {
    my($old, $new) = String::Diff::diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}
