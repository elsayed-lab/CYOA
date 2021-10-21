# -*-Perl-*-
use Test::More qw"no_plan";
use Cwd;
use File::Copy qw"cp";
use File::Path qw"remove_tree make_path rmtree";
use File::ShareDir qw"dist_file";
use String::Diff qw"diff";
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

my $gff2fasta = $cyoa->Bio::Adventure::Convert::Gff2Fasta(
    input => $phix_local, gff => $gff_local,
    libdir => 'share');
ok($gff2fasta, 'Run gff2fasta.');

my $run_fasta = $cyoa->Bio::Adventure::Align_Fasta::Split_Align_Fasta(
    input => $cds_local,
    library => $phix_local,
    fasta_tool => 'fasta36',
    number => 1,
    parse => 0,);
ok($run_fasta, 'Run Split_Align_Fasta.');

my $parsed_file = 'outputs/fasta_phix_cds_nt_phix/phix_cds_nt_vs_phix.parsed.txt';
## Caveat: every fasta36 run will give slightly different E-values due to the usage of rand().
my $expected = qq"Name
1406
1406
136
136
890
890
136
136
312
";

$actual = qx"less ${parsed_file} | awk '{print \$2}' | head";
unless(ok($expected eq $actual, 'Is the resulting table of hits expected?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

chdir($start);
