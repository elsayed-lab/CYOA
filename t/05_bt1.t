# -*-Perl-*-
use Test::More qw"no_plan";
use Bio::Adventure;
use File::Path qw"remove_tree";
use File::Copy qw"cp";
use String::Diff qw"diff";
use File::ShareDir qw":ALL";
use Cwd;

my $start = getcwd();
my $new = 'test_output';
mkdir($new);
chdir($new);

my $input_file = dist_file('Bio-Adventure', 'test_forward.fastq.gz');
my $cyoa = Bio::Adventure->new(cluster => 0, basedir => cwd());
ok(cp($input_file, "test_forward.fastq.gz"),
    'Copying data.');

mkdir('share/genome/indexes'); ## Make a directory for the phix indexes.
ok($cyoa->Bio::Adventure::Map::Bowtie(
       input => qq"test_forward.fastq.gz",
       htseq_id => 'ID',
       htseq_type => 'CDS',
       libdir => 'share',
       species => 'phix',),
   'Run Bowtie1');

ok(my $actual = $cyoa->Last_Stat(input => 'outputs/bowtie_stats.csv'),
   'Collect Bowtie1 Statistics');

my $expected = qq"CYOA,v0M1,10000,10000,30,9970,0,33333.3333333333,CYOA-v0M1.count.xz";
my $old_expected = qq"CYOA,v0M1,0,10000,30,9970,0,33333.3333333333,CYOA-v0M1.count.xz";
## old bowtie provides different numbers and I am not chasing them down.
unless(ok(($expected eq $actual || $old_expected eq $actual),
          'Are the bowtie stats as expected?')) {
    my ($old, $new) = diff($expected, $actual);
    my $stupid;
    ($old, $stupid) = diff($old_expected, $actual);
    diag("--\n${old}\n--\n${new}\n--\n${stupid}\n");
}

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

$actual = qx"xzcat outputs/bowtie_phix/CYOA-v0M1.count.xz | head";
unless(ok($expected eq $actual,
          'Is the resulting count table as expected?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

chdir($start);
