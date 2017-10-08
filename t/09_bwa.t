# -*-Perl-*-
use Test::More qw"no_plan";
use Bio::Adventure;
use File::Path qw"remove_tree";
use File::Copy qw"cp";
use String::Diff qw"diff";
use Archive::Extract;

my ($expected, $actual);
my $test_num = "01";
my $cyoa = new Bio::Adventure();
mkdir('share/genome/indexes'); ## Make a directory for the phix indexes.

ok(cp('share/test_forward.fastq.gz', 'test_forward.fastq.gz'),
   "${test_num}: Copy test files.");
$test_num++;


ok(Bio::Adventure::RNASeq_Map::BWA(
       $cyoa,
       htseq_id => 'ID',
       htseq_type => 'CDS',
       input => qq"test_forward.fastq.gz",
       libdir => 't/data',
       pbs => 0,
       species => 'phix',
   ),
   "${test_num}: Run Bwa");
$test_num++;

##ok($actual = $cyoa->Last_Stat(input => 'outputs/bwa_stats.csv'),
##   "${test_num}: Collect bwa Statistics");
##$test_num++;
##$expected = qq'CYOA,10000,10000,44,49,0,';
##unless(ok($actual eq $expected,
##          "${test_num}: Are the bwa results the expected value?")) {
##    my ($old, $new) = diff($expected, $actual);
##    diag("expected:\n$old\nactual:\n$new\n");
##}
##$test_num++;

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
";
$actual = qx"xzcat outputs/bwa_phix/CYOA_mem.count.xz | head";
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
";
$actual = qx"xzcat outputs/bwa_phix/CYOA_aln.count.xz | head";
unless(ok($expected eq $actual,
          "${test_num}: Check bwa count tables aln version.")) {
    my ($old, $new) = String::Diff::diff($expected, $actual);
    diag("expected:\n$old\nactual:\n$new\n");
}
