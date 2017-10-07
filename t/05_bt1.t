use Test::More qw"no_plan";
use Bio::Adventure;
use File::Path qw"remove_tree";
use File::Copy qw"cp";
use String::Diff qw"diff";

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

ok(my $actual = $cyoa->Last_Stat(input => 'outputs/bowtie_stats.csv'),
   'Collect Bowtie1 Statistics');

my $expected = qq"CYOA2,v0M1,10000,10000,30,9970,0,33333.3333333333,CYOA2-v0M1.count.xz";
ok($actual eq $expected,
   'Are the bowtie results the expected value?');

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
";

## WHAT IN THE SHIT!?!
my $ls = qx"/bin/ls outputs/bowtie_phix/*.xz";
chomp($ls);
diag($ls);
$actual = qx"xzcat ${ls} | head";
diag($actual);

unless(ok($expected eq $actual,
          'Is the resulting count table as expected?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}
