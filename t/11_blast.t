# -*-Perl-*-
use strict;
use Test::More qw"no_plan";
use Cwd;
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
my $cyoa = Bio::Adventure->new(cluster => 0, basedir => cwd());
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
ok($phix_library, 'Run Check_Blastdb()');
ok($phix_library eq 'blastdb/phix', 'The database name is "phix".');
$ENV{BLASTDB} = 'blastdb';
## Run a standalone blast search of the phix CDSes vs. the phix genome.
my $standalone = $cyoa->Bio::Adventure::Align_Blast::Run_Parse_Blast(
    input => $cds_local,
    blast_tool => 'blastn',
    output => 'blast_result.txt',
    library => $phix_library);

##my $run_blast = $cyoa->Bio::Adventure::Align_Blast::Split_Align_Blast(
##    input => $cds_local,
##    library => $phix_local,
##    number => 1,
##    parse => 0,);
##ok($run_blast, 'Run Split_Align_Blast.');
