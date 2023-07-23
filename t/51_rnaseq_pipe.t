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

my ($actual, $expected) = '';
## Invoke the pipeline, keep it within our test directory with basedir.
my $cyoa = Bio::Adventure->new(
    cluster => 0,
    basedir => cwd(),
    gff_tag => 'gene_id',
    gff_type => 'gene',
    libdir => cwd(),
    species => 'phix',);
my $rnaseq = $cyoa->Bio::Adventure::Pipeline::Process_RNAseq(
    input => 'test_forward.fastq.gz',);

## Check the trimomatic output.
my $test_file = $rnaseq->{'01trim'}->{stderr};
my $comparison = ok(-f $test_file, qq"Checking trimomatic output: ${test_file}");
print "Passed.\n" if ($comparison);
$expected = qq"Using PrefixPair: 'TACACTCTTTCCCTACACGACGCTCTTCCGATCT' and 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
Using Long Clipping Sequence: 'GGAATTCTCGGGTGCCAAGGAACTCCAA'
Using Long Clipping Sequence: 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA'
Using Long Clipping Sequence: 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
Using Long Clipping Sequence: 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
Using Long Clipping Sequence: 'TACACTCTTTCCCTACACGACGCTCTTCCGATCT'
Using Long Clipping Sequence: 'AGATCGGAAGAGCACACGTCTGAAC'
ILLUMINACLIP: Using 1 prefix pairs, 6 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences
Input Read Pairs: 108212 Both Surviving: 95613 (88.36%) Forward Only Surviving: 10173 (9.40%) Reverse Only Surviving: 865 (0.80%) Dropped: 1561 (1.44%)
TrimmomaticPE: Completed successfully
";
$actual = qx"tail ${test_file}";
#$comparison = ok($expected eq $actual, 'Checking trimomatic result:');
#if ($comparison) {
#    print "Passed.\n";
#} else {
#    my ($e, $a) = diff($expected, $actual);
#    diag("-- expected\n${e}\n-- actual\n${a}\n");
#}

## Look at the fastqc outputs
$test_file = $rnaseq->{'02fastqc'}->{txtfile};
$comparison = ok(-f $test_file, qq"Checking fastqc output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"less ${test_file}";
$expected = qq"PASS\tBasic Statistics\tr1-trimmed.fastq
PASS\tPer base sequence quality\tr1-trimmed.fastq
WARN\tPer tile sequence quality\tr1-trimmed.fastq
PASS\tPer sequence quality scores\tr1-trimmed.fastq
WARN\tPer base sequence content\tr1-trimmed.fastq
WARN\tPer sequence GC content\tr1-trimmed.fastq
PASS\tPer base N content\tr1-trimmed.fastq
WARN\tSequence Length Distribution\tr1-trimmed.fastq
WARN\tSequence Duplication Levels\tr1-trimmed.fastq
WARN\tOverrepresented sequences\tr1-trimmed.fastq
PASS\tAdapter Content\tr1-trimmed.fastq
";


my $sam_file = $rnaseq->{'03hisat'}->{samtools}->{output};
my $htseq_file = $rnaseq->{'03hisat'}->{htseq}->[0]->{output};
my $stats_file = $rnaseq->{'03hisat'}->{stats}->{output};

ok(-f $sam_file, 'The sorted bamfile was created.');
ok(-f $htseq_file, 'The count table was created.');
ok(-f $stats_file, 'The rnaseq stats were recorded.');

$actual = $cyoa->Last_Stat(input => $stats_file);
ok($actual, 'Collect Rnaseq Statistics');

$expected = qq"hisat2_phix_genome_test_output.stderr,9316,35,9281,0";
unless(ok($expected eq $actual, 'Are the rnaseq stats as expected?')) {
    my ($old, $new) = diff($expected, $actual);
    diag("--Expected--\n${old}\n--Actual--\n${new}\n");
}

$expected = qq"phiX174p01\t5
phiX174p02\t3
phiX174p03\t0
phiX174p04\t0
phiX174p05\t0
phiX174p06\t10
phiX174p07\t0
phiX174p08\t0
phiX174p09\t0
phiX174p10\t0
phiX174p11\t0
__no_feature\t0
__ambiguous\t22
__too_low_aQual\t0
__not_aligned\t9281
__alignment_not_unique\t0
";

$actual = qx"less ${htseq_file}";
unless(ok($expected eq $actual, 'Is the resulting count table as expected?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--Expected--\n${old}\n--Actual--\n${new}\n");
}

chdir($start);
