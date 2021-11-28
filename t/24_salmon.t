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

my $index = $cyoa->Bio::Adventure::Index::Salmon_Index(
    input => $phix_transcripts,
    species => 'phix');
my $salmon = $cyoa->Bio::Adventure::Map::Salmon(
    input => 'test_forward.fastq.gz',
    libdir => '.',
    jprefix => 24,
    species => 'phix',);
ok($salmon, 'Run salmon');
my $salmon_file = $salmon->{output};
ok(-f $salmon_file, 'The salmon quant.sf file was created.');

my $stats_file = $salmon->{stats}->{output};
ok(-f $stats_file, "The salmon stats file was created: ${stats_file}");
my $actual = $cyoa->Last_Stat(input => $stats_file);
ok($actual, 'Collect salmon Statistics');

my $expected = qq"test_output,phix,44,44,44,0,0.4318181818181818";
unless(ok($expected eq $actual, 'Are the mapping stats from salmon expected?')) {
    my ($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

$expected = qq"Name\tLength\tEffectiveLength\tTPM\tNumReads
NC_001422_phiX174p01_1\t1406\t1165.522\t79003.686307\t7.540
NC_001422_phiX174p01_2\t136\t4.000\t0.000000\t0.000
NC_001422_phiX174p02_3\t890\t639.000\t0.000000\t0.000
NC_001422_phiX174p03_5\t312\t62.000\t287513.563515\t1.460
NC_001422_phiX174p03_6\t51\t2.000\t0.000000\t0.000
NC_001422_phiX174p04_7\t171\t6.000\t0.000000\t0.000
NC_001422_phiX174p05_8\t261\t24.000\t508825.648290\t1.000
NC_001422_phiX174p06_9\t3618\t3330.751\t124657.101888\t34.000
NC_001422_phiX174p07_10\t634\t359.172\t0.000000\t0.000
NC_001422_phiX174p08_11\t276\t32.000\t0.000000\t0.000
NC_001422_phiX174p09_12\t117\t3.000\t0.000000\t0.000
NC_001422_phiX174p10_13\t528\t277.000\t0.000000\t0.000
NC_001422_phiX174p11_14\t987\t743.819\t0.000000\t0.000
";

$actual = qx"less ${salmon_file}";
unless(ok($expected eq $actual, 'Is the resulting count table as expected?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

chdir($start);
