# -*-Perl-*-
use Test::More qw"no_plan";
use Cwd;
use File::Copy qw"cp";
use File::Path qw"remove_tree make_path rmtree";
use File::ShareDir qw"dist_file";
use String::Diff qw"diff";
use Bio::SeqIO;
use Bio::Adventure;

my $start = getcwd();
my $new = 'test_output';
mkdir($new);
chdir($new);

my $phix_fasta = dist_file('Bio-Adventure', 'genome/phix.fasta');
my $phix_gff = dist_file('Bio-Adventure', 'genome/phix.gff');
my $phix_local = 'genome/phix.fasta';
my $gff_local = 'genome/phix.gff';
my $cds_local = 'genome/phix_cds_nt.fasta';

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
        libdir => 'share');
    ok($gff2fasta, 'Run gff2fasta.');
}

## Perform a check to see if the desired blast library exists,
## create it if it does not.
my $phix_library = $cyoa->Bio::Adventure::Index::Check_Blastdb(
    input => $phix_local, type => 'nucl');
ok($phix_library, 'Run Check_Blastdb()');
ok($phix_library eq 'blastdb/phix', 'The database name is "phix".');

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
