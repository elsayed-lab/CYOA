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

my $start = getcwd();
my $new = 'test_output';
mkdir($new);
chdir($new);

make_path('genome/indexes'); ## Make a directory for the phix indexes.
my $input_file = qq"${start_dir}/test_forward.fastq.gz";
my $phix_fasta = qq"${start_dir}/genome/phix.fasta";
my $phix_gff = qq"${start_dir}/genome/phix.gff";
if (!-r 'test_forward.fastq.gz') {
    ok(cp($input_file, 'test_forward.fastq.gz'), 'Copying data.');
}
if (!-r 'genome/phix.fasta') {
    ok(cp($phix_fasta, 'genome/phix.fasta'), 'Copying phix fasta file.');
}
if (!-r 'genome/phix.gff') {
    ok(cp($phix_gff, 'genome/phix.gff'), 'Copying phix gff file.');
}

my $cyoa = Bio::Adventure->new(
    basedir => cwd(),
    gff_tag => 'gene_id',
    gff_type => 'gene',
    libdir => cwd(),
    species => 'phix',
    stranded => 'no',);
my $hisat = $cyoa->Bio::Adventure::Map::Hisat2(
    input => 'test_forward.fastq.gz',);
ok($hisat, 'Submit hisat2 jobs to the cluster.');
my $status = $cyoa->Wait(job => $hisat);
ok($status->{State} eq 'COMPLETED', 'The hisat2 jobs completed.');
my $sam_file = $hisat->{samtools}->{output};
my $htseq_file = $hisat->{htseq}->[0]->{output};
my $stats_file = $hisat->{stats}->{output};

ok(-f $sam_file, 'The sorted bamfile was created.');
ok(-f $htseq_file, 'The count table was created.');
ok(-f $stats_file, 'The hisat stats were recorded.');

my $actual = $cyoa->Last_Stat(input => $stats_file);
ok($actual, 'Collect Hisat Statistics');

my $expected = qq"hisat2_phix_genome_test_output.stderr,10000,49,9951,0";
unless(ok($expected eq $actual, 'Are the hisat stats as expected?')) {
    my ($old, $new) = diff($expected, $actual);
    diag("--Expected--\n${old}\n--Actual--\n${new}\n");
}

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

$actual = qx"less ${htseq_file}";
unless(ok($expected eq $actual, 'Is the resulting count table as expected?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--Expected--\n${old}\n--Actual--\n${new}\n");
}

chdir($start);
