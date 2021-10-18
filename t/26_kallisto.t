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

## In order to make a decoy aware transcriptome file,
## I would like to have the genome and transcript file sitting in the
## same place.  Thus I will copy the genome and gff files and invoke
## gff2fasta on them.
my $phix_genome = 'genome/phix.fasta';
my $phix_annot = 'genome/phix.gff';
my $phix_transcripts = 'genome/phix_cds_nt.fasta';
if (!-r $phix_genome) {
    ok(cp($phix_fasta, $phix_genome), 'Copying phix fasta file.');
    ## my $uncompressed = qx"gunzip genome/phix.fastq.gz && mv genome/phix.fasta.gz genome/phix.fasta";
}
if (!-r $phix_annot) {
    ok(cp($phix_gff, $phix_annot), 'Copying phix gff file.');
    ## my $uncompressed = qx"gunzip genome/phix.gff.gz && mv genome/phix.gff.gz genome/phix.gff";
}

my $cyoa = Bio::Adventure->new(cluster => 0, basedir => cwd());
my $gff2fasta = $cyoa->Bio::Adventure::Convert::Gff2Fasta(
    input => $phix_genome, gff => $phix_annot, libdir => '.');
ok($gff2fasta, 'Run gff2fasta.');

my $index = $cyoa->Bio::Adventure::Map::Kallisto_Index(
    input => $phix_transcripts,
    species => 'phix');
my $kallisto_file = 'outputs/46kallisto_phix/abundance.tsv';
my $kallisto = $cyoa->Bio::Adventure::Map::Kallisto(
    input => 'test_forward.fastq.gz',
    libdir => '.',
    species => 'phix',);
ok($kallisto, 'Run Kallisto');
ok(-f $kallisto_file, 'The kallisto abundance.tsv file was created.');

$expected = qq"target_id\tlength\teff_length\test_counts\ttpm
NC_001422_phiX174p01_1\t1406\t1125.2\t8.50998\t189149
NC_001422_phiX174p01_2\t136\t96.9989\t0\t0
NC_001422_phiX174p02_3\t890\t735.894\t0.281695\t9573.51
NC_001422_phiX174p02_4\t136\t96.9989\t0\t0
NC_001422_phiX174p03_5\t312\t99.6032\t1.20833\t303402
NC_001422_phiX174p03_6\t51\t14.352\t0\t0
NC_001422_phiX174p04_7\t171\t131.999\t0\t0
NC_001422_phiX174p05_8\t261\t179.863\t1\t139048
NC_001422_phiX174p06_9\t3618\t3648.51\t29.5203\t202354
NC_001422_phiX174p07_10\t634\t637.649\t2.84075\t111419
NC_001422_phiX174p08_11\t276\t236.999\t0\t0
NC_001422_phiX174p09_12\t117\t77.9989\t0\t0
NC_001422_phiX174p10_13\t528\t488.999\t0\t0
NC_001422_phiX174p11_14\t987\t909.775\t1.63894\t45054.3
";

$actual = qx"less ${kallisto_file}";
unless(ok($expected eq $actual, 'Is the resulting count table as expected?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

chdir($start);
