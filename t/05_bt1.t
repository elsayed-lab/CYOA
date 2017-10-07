# -*-Perl-*-
use Test::More qw"no_plan";
use Bio::Adventure;
use File::Path qw"remove_tree";
use File::Copy qw"cp";
use String::Diff qw( diff_fully diff diff_merge diff_regexp );

my $cyoa = Bio::Adventure->new();
ok(cp("share/test_forward.fastq.gz", "test_forward.fastq.gz"),
    'Copying data.');

mkdir('share/genome/indexes'); ## Make a directory for the phix indexes.
ok(Bio::Adventure::RNASeq_Map::Bowtie($cyoa,
                                      input => qq"test_forward.fastq.gz",
                                      species => 'phix',
                                      htseq_id => 'ID',
                                      htseq_type => 'CDS',
                                      libdir => 'share'),
   'Run Bowtie1');

my $fqc_dir = qx"find outputs/";
diag($fqc_dir);
my $script = qx'cat scripts/10test_forward.sh';
diag($script);
my $log = qx'cat outputs/bowtie_phix/CYOA2-v0M1.err';
diag($log);

ok(my $actual = $cyoa->Last_Stat(input => 'outputs/bowtie_stats.csv'),
   'Collect Bowtie1 Statistics');

my $expected = qq"CYOA2,v0M1,10000,10000,30,9970,0,33333.3333333333,CYOA2-v0M1.count.xz";
ok($actual eq $expected,
   'Are the bowtie results the expected value?');
##my($old, $new) = String::Diff::diff($expected, $actual);
##print STDERR "$old\n";
##print STDERR "---\n";
##print STDERR "$new\n";
##print STDERR "---\n";

$expected = qq"phiX174p01\t1
phiX174p02\t0
phiX174p03\t0
phiX174p04\t0
phiX174p05\t0
phiX174p06\t8
phiX174p07\t0
phiX174p08\t0
phiX174p09\t0
phiX174p10\t0
phiX174p11\t0
__no_feature\t0
__ambiguous\t21
__too_low_aQual\t0
__not_aligned\t9970
__alignment_not_unique\t0
";
$actual = qx"xzcat outputs/bowtie_phix/CYOA2-v0M1.count.xz";
unless(ok($expected eq $actual,
          'Is the resulting count table as expected?')) {
    my($old, $new) = String::Diff::diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}
