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

my $cyoa = Bio::Adventure->new(
    cluster => 0,
    basedir => cwd(),
    libdir => cwd(),
    species => 'phix',);
my $gff2fasta = $cyoa->Bio::Adventure::Convert::Gff2Fasta(
    input => $phix_genome, gff => $phix_annot,
    gff_type => 'gene', gff_tag => 'gene_id',);
ok($gff2fasta, 'Run gff2fasta.');
ok(cp('phix_gene_gene_id_nt.fasta', $phix_transcripts));
my $index = $cyoa->Bio::Adventure::Index::Kallisto_Index(
    input => $phix_transcripts,);
## There is a problem here, in all the other tests so far I have not needed
## to redefine the species as I do here, but for some crazy ass reason
## it gets lost in this invocation and I do not know why.
my $kallisto = $cyoa->Bio::Adventure::Map::Kallisto(
    input => 'test_forward.fastq.gz',
    species => 'phix',
    jprefix => '26',);
ok($kallisto, 'Run Kallisto');

##my $kallisto_file = 'outputs/46kallisto_phix/abundance.tsv';
my $kallisto_file = $kallisto->{output};
ok(-f $kallisto_file, qq"The kallisto otuput file was created: ${kallisto_file}.");

my $expected = qq"target_id	length	eff_length	est_counts	tpm
chr_NC_001422_id_phiX174p01_start_3981_end_5386	1406	1125.2	8.50998	189149
chr_NC_001422_id_phiX174p01_start_1_end_136	136	96.9989	0	0
chr_NC_001422_id_phiX174p02_start_4497_end_5386	890	735.894	0.281695	9573.51
chr_NC_001422_id_phiX174p02_start_1_end_136	136	96.9989	0	0
chr_NC_001422_id_phiX174p03_start_5075_end_5386	312	99.6032	1.20833	303402
chr_NC_001422_id_phiX174p03_start_1_end_51	51	14.352	0	0
chr_NC_001422_id_phiX174p04_start_51_end_221	171	131.999	0	0
chr_NC_001422_id_phiX174p05_start_133_end_393	261	179.863	1	139048
chr_NC_001422_id_phiX174p06_start_358_end_3975	3618	3648.51	29.5203	202354
chr_NC_001422_id_phiX174p07_start_358_end_991	634	637.649	2.84075	111419
chr_NC_001422_id_phiX174p08_start_568_end_843	276	236.999	0	0
chr_NC_001422_id_phiX174p09_start_848_end_964	117	77.9989	0	0
chr_NC_001422_id_phiX174p10_start_2395_end_2922	528	488.999	0	0
chr_NC_001422_id_phiX174p11_start_2931_end_3917	987	909.775	1.63894	45054.3
";

my $actual = qx"less ${kallisto_file}";
unless(ok($expected eq $actual, 'Is the resulting count table as expected?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

chdir($start);
