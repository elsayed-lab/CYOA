# -*-Perl-*-
use Test::More qw"no_plan";
use Bio::Adventure;
use File::Path qw"remove_tree";
use String::Diff qw( diff_fully diff diff_merge diff_regexp );

my $cyoa = Bio::Adventure->new(pbs => 0);

ok(Bio::Adventure::RNASeq_QA::Fastqc_Single($cyoa, input => qq"t/data/test_forward.fastq.gz"),
   'Run Fastqc');

ok(my $actual = $cyoa->Last_Stat(input => 'outputs/fastqc_stats.csv'),
   'Collect Fastqc Statistics');

my $expected = 'fqc_test_forward,10000,0,pass,warn,pass,pass,pass,warn,fail,pass';
unless(ok($expected eq $actual,
          'Are the fastqc results the expected value?')) {
    my($old, $new) = String::Diff::diff($expected, $actual);
    diag("$old\n$new\n");
}
