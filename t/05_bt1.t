# -*-Perl-*-
use Test::More qw"no_plan";
use Bio::Adventure;
use Cwd;
use File::Copy qw"cp";
use File::Path qw"remove_tree make_path rmtree";
use File::ShareDir qw"dist_file";
use String::Diff qw"diff";

my $start = getcwd();
my $new = 'test_output';
mkdir($new);
chdir($new);
make_path('genome/indexes'); ## Make a directory for the phix indexes.
my $input_file = dist_file('Bio-Adventure', 'test_forward.fastq.gz');
my $phix_fasta = dist_file('Bio-Adventure', 'genome/phix.fasta');
my $phix_gff = dist_file('Bio-Adventure', 'genome/phix.gff');
if (!-r 'test_forward.fastq.gz') {
    ok(cp($input_file, 'test_forward.fastq.gz'), 'Copying data.');
}

if (!-r 'genome/phix.fasta') {
    ok(cp($phix_fasta, 'genome/phix.fasta'), 'Copying phix fasta file.');
    ## my $uncompressed = qx"gunzip genome/phix.fastq.gz && mv genome/phix.fasta.gz genome/phix.fasta";
}

if (!-r 'genome/phix.gff') {
    ok(cp($phix_gff, 'genome/phix.gff'), 'Copying phix gff file.');
    ## my $uncompressed = qx"gunzip genome/phix.gff.gz && mv genome/phix.gff.gz genome/phix.gff";
}

my $cyoa = Bio::Adventure->new(cluster => 0, basedir => cwd());
my $bt1 = $cyoa->Bio::Adventure::Map::Bowtie(
    input => qq"test_forward.fastq.gz",
    htseq_id => 'gene_id',
    htseq_type => 'gene',
    libdir => '.',
    species => 'phix',);
ok($bt1, 'Run Bowtie1');
my $sam_file = $bt1->{samtools}->{output};
my $htseq_file = $bt1->{htseq}->[0]->{output};

ok(-f $sam_file, 'The sorted bamfile was created.');
ok(-f $htseq_file, 'The count table was created.');

my $actual = $cyoa->Last_Stat(input => 'outputs/bowtie_stats.csv');
ok($actual, 'Collect Bowtie1 Statistics');

my $expected = qq"test_output,v0M1,10000,10000,30,9970,0,33333.3333333333,test_output-v0M1.count.xz";
my $old_expected = qq"test_output,v0M1,0,10000,30,9970,0,33333.3333333333,test_output-v0M1.count.xz";
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
phiX174p11\t0
__no_feature\t0
__ambiguous\t21
__too_low_aQual\t0
__not_aligned\t9970
__alignment_not_unique\t0
";

$actual = qx"less ${htseq_file}";
unless(ok($expected eq $actual, 'Is the resulting count table as expected?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

chdir($start);
