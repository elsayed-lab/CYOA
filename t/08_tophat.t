# -*-Perl-*-
use Test::More qw"no_plan";
use Bio::Adventure;
use File::Path qw"remove_tree";
use File::Copy qw"cp";
use String::Diff qw( diff_fully diff diff_merge diff_regexp );

my ($expected, $actual, $test_num) = "01";
my $cyoa = new Bio::Adventure();

ok(cp('t/data/test_forward.fastq.gz', 'test_forward.fastq.gz'),
   $test_num);

mkdir('t/data/genome/indexes'); ## Make a directory for the phix indexes.
ok(Bio::Adventure::RNASeq_Map::Tophat($cyoa, input => qq"test_forward.fastq.gz", pbs => 0,
                                      htseq_id => 'ID', htseq_type => 'CDS',
                                      species => 'phix', libdir => 't/data'),
   'Run Tophat');

ok($actual = $cyoa->Last_Stat(input => 'outputs/tophat_stats.csv'),
   'Collect Tophat Statistics');

$expected = qq'CYOA2,phix,10000,10000,40,9960,25000,accepted_hits.count.xz';
unless(ok($expected eq $actual,
          'Are the tophat results the expected value?')) {
    my($old, $new) = String::Diff::diff($expected, $actual);
    diag("$old\n$new\n");
}

$expected = qq"phiX174p01\t1
phiX174p02\t0
phiX174p03\t0
phiX174p04\t0
phiX174p05\t0
phiX174p06\t12
phiX174p07\t0
phiX174p08\t0
phiX174p09\t0
phiX174p10\t0
phiX174p11\t0
__no_feature\t0
__ambiguous\t27
__too_low_aQual\t0
__not_aligned\t0
__alignment_not_unique\t0
";
$actual = qx"xzcat outputs/tophat_phix/accepted_hits.count.xz";
unless(ok($expected eq $actual,
          'Is the resulting count table as expected?')) {
    my($old, $new) = String::Diff::diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}
