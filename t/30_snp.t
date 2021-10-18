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

if (!-r 'genome/phix.fasta') {
    ok(cp($phix_fasta, 'genome/phix.fasta'), 'Copying phix fasta file.');
    ## my $uncompressed = qx"gunzip genome/phix.fastq.gz && mv genome/phix.fasta.gz genome/phix.fasta";
}

if (!-r 'genome/phix.gff') {
    ok(cp($phix_gff, 'genome/phix.gff'), 'Copying phix gff file.');
    ## my $uncompressed = qx"gunzip genome/phix.gff.gz && mv genome/phix.gff.gz genome/phix.gff";
}

my $cyoa = Bio::Adventure->new(cluster => 0, basedir => cwd());

my $variant = $cyoa->Bio::Adventure::SNP::Align_SNP_Search(
    input => qq"test_forward.fastq.gz",
    htseq_id => 'ID',
    htseq_type => 'CDS',
    libdir => 'share',
    species => 'phix_3N',
    vcf_cutoff => 1,);

ok($variant, 'Ran variant search.');
