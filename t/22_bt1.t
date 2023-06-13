# -*-Perl-*-
use strict;
use Test::More qw"no_plan";
use Bio::Adventure;
use Cwd;
use File::Copy qw"cp mv";
use File::Path qw"remove_tree make_path rmtree";
use File::ShareDir qw"dist_file module_dir dist_dir";
use String::Diff qw"diff";
use Test::File::ShareDir::Dist { 'Bio-Adventure' => 'share/' };
my $start_dir = dist_dir('Bio-Adventure');
my $input_file = qq"${start_dir}/test_forward.fastq.gz";
my $phix_fasta = qq"${start_dir}/genome/phix.fastq";
my $phix_gff = qq"${start_dir}/genome/phix.gff";

my $start = getcwd();
my $new = 'test_output';
mkdir($new);
chdir($new);
make_path('genome/indexes'); ## Make a directory for the phix indexes.

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

my $cyoa = Bio::Adventure->new(
    cluster => 0,
    basedir => cwd(),
    libdir => cwd(),
    stranded => 'no');
my $bt1 = $cyoa->Bio::Adventure::Map::Bowtie(
    input => qq"test_forward.fastq.gz",
    gff_tag => 'gene_id',
    gff_type => 'gene',
    jprefix => '22',
    species => 'phix',);
ok($bt1, 'Run Bowtie1');
my $sam_file = $bt1->{samtools}->{output};
my $htseq_file = $bt1->{htseq}->[0]->{output};

ok(-f $sam_file, 'The sorted bamfile was created.');
ok(-f $htseq_file, 'The count table was created.');

my $actual = $cyoa->Last_Stat(input => 'outputs/bowtie_stats.csv');
ok($actual, 'Collect Bowtie1 Statistics');
my $expected = qq"test_output-v0M1.stderr,v0M1,0,10000,0,9970,0,outputs/22bowtie_phix/test_output-v0M1_sno_gene_gene_id.count.xz";
unless(ok($expected eq $actual, 'Are the bowtie stats as expected?')) {
    print "IN UNLESS?\n";
    my ($old, $new) = diff($expected, $actual);
    print "--Expected--\n${old}\n--Actual--\n${new}\n";
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
    diag("--Expected--\n${old}\n--Actual--\n${new}\n");
}

chdir($start);
