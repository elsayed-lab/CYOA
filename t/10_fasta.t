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

my $cyoa = Bio::Adventure->new(cluster => 0, basedir => cwd());

my $gff2fasta = $cyoa->Bio::Adventure::Convert::Gff2Fasta(
    input => $phix_fasta, gff => $phix_gff,
    genome => $phix_fasta, libdir => 'share');
ok($gff2fasta, 'Run gff2fasta.');

my $run_fasta = $cyoa->Bio::Adventure::Align_Fasta::Split_Align_Fasta(
    query => 'phix_cds_nt.fasta',
    library => 'phix_cds_nt.fasta',
    number => 1,
    parse => 0,);
ok($run_fasta, 'Run Split_Align_Fasta.');
