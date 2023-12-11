# -*-Perl-*-
use strict;
use Test::More qw"no_plan";
use Bio::Adventure;
use Cwd;
use File::Basename qw"basename";
use File::Copy qw"cp mv";
use File::Path qw"remove_tree make_path rmtree";
use File::ShareDir qw"dist_file module_dir dist_dir";
use String::Diff qw"diff";
use Test::File::ShareDir::Dist { 'Bio-Adventure' => 'share/' };
my $start_dir = dist_dir('Bio-Adventure');
my $input_file = qq"${start_dir}/test_forward.fastq.gz";
my $phix_fasta = qq"${start_dir}/genome/phix.fasta";
my $phix_gff = qq"${start_dir}/genome/phix.gff";
my $input_local = basename($input_file);
my $phix_fasta_local = qw'genome/' . basename($phix_fasta);
my $phix_gff_local = qw'genome/' . basename($phix_gff);

my $start = getcwd();
my $new = 'test_output';
mkdir($new);
chdir($new);
make_path('genome/indexes'); ## Make a directory for the phix indexes.

if (!-r $input_local) {
    ok(cp($input_file, $input_local), 'Copying data.');
}
if (!-r $phix_fasta_local) {
    ok(cp($phix_fasta, $phix_fasta_local), 'Copying phix fasta file.');
}
if (!-r $phix_gff_local) {
    ok(cp($phix_gff, $phix_gff_local), 'Copying phix gff file.');
}

my $cyoa = Bio::Adventure->new(
    basedir => cwd(),
    libdir => cwd(),
    stranded => 'no');
my $bt1 = $cyoa->Bio::Adventure::Map::Bowtie(
    input => $input_local,
    gff_tag => 'gene_id',
    gff_type => 'gene',
    jprefix => '22',
    species => 'phix',);
ok($bt1, 'Submitted bowtie1 jobs.');
my $status = $cyoa->Wait(job => $bt1);
ok($status->{State} eq 'COMPLETED', 'The bowtie jobs completed.');
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
