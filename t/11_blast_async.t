# -*-Perl-*-
use strict;
use Test::More qw"no_plan";
use Cwd;
use File::Basename qw"basename";
use File::Copy qw"cp mv";
use File::Path qw"remove_tree make_path rmtree";
use File::ShareDir qw"dist_file module_dir dist_dir";
use String::Diff qw"diff";
use Test::File::ShareDir::Dist { 'Bio-Adventure' => 'share/' };
use Bio::SeqIO;
use Bio::Adventure;

my $start = getcwd();
my $new = 'test_output';
mkdir($new);
chdir($new);

my $start_dir = dist_dir('Bio-Adventure');

my $phix_fasta = qq"${start_dir}/genome/phix.fasta";
my $phix_gff = qq"${start_dir}/genome/phix.gff";
my $phix_local = 'genome/phix.fasta';
my $gff_local = 'genome/phix.gff';
my $cds_local = 'phix_gene_gene_id_nt.fasta';

make_path('genome'); ## Make a directory for the phix indexes.
if (!-r $phix_local) {
    ok(cp($phix_fasta, $phix_local), 'Copying phix fasta file.');
}

if (!-r $gff_local) {
    ok(cp($phix_gff, $gff_local), 'Copying phix gff file.');
}
my $cyoa = Bio::Adventure->new(basedir => cwd());
if (!-r $cds_local) {
    my $gff2fasta = $cyoa->Bio::Adventure::Convert::Gff2Fasta(
        input => $phix_local, gff => $gff_local,
        gff_tag => 'gene_id', gff_type => 'gene',
        libdir => 'share');
    ok($gff2fasta, 'Run gff2fasta.');
}

## Perform a check to see if the desired blast library exists,
## create it if it does not.
my $phix_library = $cyoa->Bio::Adventure::Index::Check_Blastdb(
    input => $phix_local, type => 'nucl');
my $actual = basename($phix_library);
ok($actual eq 'phix', 'The database name is "phix".');

$ENV{BLASTDB} = 'blastdb';
## Run a standalone blast search of the phix CDSes vs. the phix genome.
my $standalone = $cyoa->Bio::Adventure::Align_Blast::Run_Parse_Blast(
    input => $cds_local,
    blast_tool => 'blastn',
    output => 'blast_result.txt',
    library => $phix_library);

my $run_blast = $cyoa->Bio::Adventure::Align_Blast::Split_Align_Blast(
    input => $cds_local,
    library => $phix_local,
    number => 1,
    parse => 0,);
my $status = $cyoa->Wait(job => $run_blast);
ok($status->{Status} eq 'COMPLETED', 'The blast job has completed.');
##ok($run_blast, 'Run Split_Align_Blast.');
my $expected = qq"QUERYNAME	real_name	Chromosome	Start	End	%ID	Score	Sig	CompLength	Hit_Ident	hits
Query_1	chr_NC_001422_id_phiX174p01_start_3981_end_5386	NC_001422		3981	5386	100	1406	0	1406	100	1
Query_2	chr_NC_001422_id_phiX174p01_start_1_end_136	NC_001422		1	136	100	136	7.66468e-71	136	100	1
Query_3	chr_NC_001422_id_phiX174p02_start_4497_end_5386	NC_001422		4497	5386	100	890	0	890	100	1
Query_4	chr_NC_001422_id_phiX174p02_start_1_end_136	NC_001422		1	136	100	136	7.66468e-71	136	100	1
Query_5	chr_NC_001422_id_phiX174p03_start_5075_end_5386	NC_001422		5075	5386	100	312	2.68412e-168	312	100	1
Query_6	chr_NC_001422_id_phiX174p03_start_1_end_51	NC_001422		1	51	100	51	4.41014e-24	51	100	1
Query_7	chr_NC_001422_id_phiX174p04_start_51_end_221	NC_001422		51	221	100	171	3.41394e-90	171	100	1
Query_8	chr_NC_001422_id_phiX174p05_start_133_end_393	NC_001422		133	393	100	261	4.99259e-140	261	100	1
Query_9	chr_NC_001422_id_phiX174p06_start_358_end_3975	NC_001422		358	3975	100	3618	0	3618	100	1
";
my $parsed_txt = $run_blast->{parser}->{parsed_output};
my $count_txt = $run_blast->{parser}->{count_output};
ok(-r $parsed_txt, 'The parsed blast output file was created.');
ok(-r $count_txt, 'The count-table blast output was created.');
$actual = qx"head ${parsed_txt}";
unless(ok($expected eq $actual,
          'Did we parse the expected blast output?')) {
    my($old, $new) = diff($expected, $actual);
    diag("$old\n$new\n");
}
