# -*-Perl-*-
use Test::More qw"no_plan";
use Bio::Adventure;
use File::Path qw"remove_tree";
use File::Copy qw"cp";
use String::Diff qw( diff_fully diff diff_merge diff_regexp );
use Archive::Extract;

my ($expected, $actual);
my $test_num = "01";
my $cyoa = new Bio::Adventure();
mkdir('t/data/genome/indexes'); ## Make a directory for the phix indexes.

ok(cp("t/data/test_forward.fastq.gz", "test_forward.fastq.gz"),
   "${test_num}: Copy test files.");
$test_num++;


ok(Bio::Adventure::RNASeq_Map::BWA($cyoa, input => qq"test_forward.fastq.gz", pbs => 0,
                                   htseq_id => 'ID', htseq_type => 'CDS',
                                   species => 'phix', libdir => 't/data',),
   "${test_num}: Run Bwa");
$test_num++;

ok($actual = $cyoa->Last_Stat(input => 'outputs/bwa_stats.csv'),
   "${test_num}: Collect bwa Statistics");
$test_num++;

$expected = qq"CYOA2,10000,10000,44,49,0,";
unless(ok($actual eq $expected,
          "${test_num}: Are the bwa results the expected value?")) {
    my ($old, $new) = String::Diff::diff($expected, $actual);
    diag("expected:\n$old\nactual:\n$new\n");
}
$test_num++;

$expected = qq"phiX174p01\t5
phiX174p02\t0
phiX174p03\t0
phiX174p04\t0
phiX174p05\t0
phiX174p06\t13
phiX174p07\t0
phiX174p08\t0
phiX174p09\t0
phiX174p10\t0
phiX174p11\t0
__no_feature\t0
__ambiguous\t31
__too_low_aQual\t0
__not_aligned\t9951
__alignment_not_unique\t0
";
$actual = qx"xzcat outputs/bwa_phix/CYOA2_mem.count.xz";
unless(ok($expected eq $actual,
          "${test_num}: Check bwa count tables mem version.")) {
    my ($old, $new) = String::Diff::diff($expected, $actual);
    diag("expected:\n$old\nactual:\n$new\n");
}
$test_num++;

$expected = qq"phiX174p01\t3
phiX174p02\t0
phiX174p03\t0
phiX174p04\t0
phiX174p05\t0
phiX174p06\t13
phiX174p07\t0
phiX174p08\t0
phiX174p09\t0
phiX174p10\t0
phiX174p11\t0
__no_feature\t0
__ambiguous\t28
__too_low_aQual\t0
__not_aligned\t9956
__alignment_not_unique\t0
";
$actual = qx"xzcat outputs/bwa_phix/CYOA2_aln.count.xz";
unless(ok($expected eq $actual,
          "${test_num}: Check bwa count tables aln version.")) {
    my ($old, $new) = String::Diff::diff($expected, $actual);
    diag("expected:\n$old\nactual:\n$new\n");
}
