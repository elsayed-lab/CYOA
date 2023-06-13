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
my $phix_fasta = qq"${start_dir}/genome/phix.fasta";
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
    ok(cp($phix_fasta, 'genome/phix.fasta'), 'Copying genome.');
}
if (!-r 'genome/phix.gff') {
    ok(cp($phix_gff, 'genome/phix.gff'), 'Copying gff.');
}

## In order to make a decoy aware transcriptome file,
## I would like to have the genome and transcript file sitting in the
## same place.  Thus I will copy the genome and gff files and invoke
## gff2fasta on them.

#if (!-r $phix_fasta) {
#    ok(cp($phix_fasta, $phix_genome), 'Copying phix fasta file.');
#    ## my $uncompressed = qx"gunzip genome/phix.fastq.gz && mv genome/phix.fasta.gz genome/phix.fasta";
#}
#if (!-r $phix_annot) {
#    ok(cp($phix_gff, $phix_annot), 'Copying phix gff file.');
    ## my $uncompressed = qx"gunzip genome/phix.gff.gz && mv genome/phix.gff.gz genome/phix.gff";
#}

my $cyoa = Bio::Adventure->new(
    cluster => 0,
    basedir => cwd(),
    libdir => cwd(),);
my $gff2fasta = $cyoa->Bio::Adventure::Convert::Gff2Fasta(
    input => 'genome/phix.fasta', gff => 'genome/phix.gff', libdir => '.',
    gff_type => 'gene', gff_tag => 'gene_id');
ok($gff2fasta, 'Run gff2fasta.');
my $phix_transcripts = 'genome/phix_cds_nt.fasta';
ok(cp('phix_gene_gene_id_nt.fasta', $phix_transcripts));


my $index = $cyoa->Bio::Adventure::Index::Salmon_Index(
    input => $phix_transcripts,
    decoy => 0,
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

my $expected = qq"lib_format_counts.json,phix,44,44,44,0,0.4318181818181818";
unless(ok($expected eq $actual, 'Are the mapping stats from salmon expected?')) {
    my ($old, $new) = diff($expected, $actual);
    diag("--Expected--\n${old}\n--Actual--\n${new}\n");
}

$expected = qq"Name	Length	EffectiveLength	TPM	NumReads
chr_NC_001422_id_phiX174p01_start_3981_end_5386	1406	1165.522	79003.686307	7.540
chr_NC_001422_id_phiX174p01_start_1_end_136	136	4.000	0.000000	0.000
chr_NC_001422_id_phiX174p02_start_4497_end_5386	890	639.000	0.000000	0.000
chr_NC_001422_id_phiX174p03_start_5075_end_5386	312	62.000	287513.563515	1.460
chr_NC_001422_id_phiX174p03_start_1_end_51	51	2.000	0.000000	0.000
chr_NC_001422_id_phiX174p04_start_51_end_221	171	6.000	0.000000	0.000
chr_NC_001422_id_phiX174p05_start_133_end_393	261	24.000	508825.648290	1.000
chr_NC_001422_id_phiX174p06_start_358_end_3975	3618	3330.751	124657.101888	34.000
chr_NC_001422_id_phiX174p07_start_358_end_991	634	359.172	0.000000	0.000
chr_NC_001422_id_phiX174p08_start_568_end_843	276	32.000	0.000000	0.000
chr_NC_001422_id_phiX174p09_start_848_end_964	117	3.000	0.000000	0.000
chr_NC_001422_id_phiX174p10_start_2395_end_2922	528	277.000	0.000000	0.000
chr_NC_001422_id_phiX174p11_start_2931_end_3917	987	743.819	0.000000	0.000
";

$actual = qx"less ${salmon_file}";
unless(ok($expected eq $actual, 'Is the resulting count table as expected?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

chdir($start);
