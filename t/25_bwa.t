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

my $cyoa = Bio::Adventure->new(
    cluster => 0,
    basedir => cwd(),
    libdir => cwd(),
    stranded => 'no',
    gff_type => 'gene',
    gff_id => 'gene_id',
    species => 'phix',);

my $index = $cyoa->Bio::Adventure::Index::BWA_Index(
    input => $phix_fasta,);
## Check that the indexes were created:
ok(-f $index->{output_sa}, "The .sa index file was created: $index->{output_sa}");
ok(-f $index->{output_pac}, "The .sa index file was created: $index->{output_pac}");
ok(-f $index->{output_bwt}, "The .sa index file was created: $index->{output_bwt}");
ok(-f $index->{output_ann}, "The .sa index file was created: $index->{output_ann}");
ok(-f $index->{output_amb}, "The .sa index file was created: $index->{output_amb}");
ok(-f $index->{output_fa}, "The .sa index file was created: $index->{output_fa}");

my $bwa = $cyoa->Bio::Adventure::Map::BWA(
    input => qq'test_forward.fastq.gz',
    jprefix => 25,);
ok($bwa, 'Run Bwa.');

## Some files of interest:
my $htseq_mem = $bwa->{htseq_mem}->[0]->{output};
my $htseq_aln = $bwa->{htseq_aln}->[0]->{output};
my $reporter_sam = $bwa->{reporter}->{output};
my $bwa_output_files = $bwa->{output};
my $aln_bam = $bwa->{samtools_aln}->{output};
ok(-f $htseq_mem, "htseq output from the mem alignment was created: ${htseq_mem}");
ok(-f $aln_bam, "samtools converted the sam to a bam file: ${aln_bam}");

$expected = qq"phiX174p01\t5
phiX174p02\t0
phiX174p03\t0
phiX174p04\t0
phiX174p05\t0
phiX174p06\t13
phiX174p07\t0
phiX174p08\t0
phiX174p09\t0
phiX174p10\t0
phiX174p11\t0
__no_feature\t0
__ambiguous\t31
__too_low_aQual\t0
__not_aligned\t9951
__alignment_not_unique\t0
";
$actual = qx"less ${htseq_mem}";
unless(ok($expected eq $actual, "Check bwa count tables mem version.")) {
    my ($old, $new) = diff($expected, $actual);
    diag("expected:\n$old\nactual:\n$new\n");
}

$expected = qq"phiX174p01\t3
phiX174p02\t0
phiX174p03\t0
phiX174p04\t0
phiX174p05\t0
phiX174p06\t13
phiX174p07\t0
phiX174p08\t0
phiX174p09\t0
phiX174p10\t0
phiX174p11\t0
__no_feature\t0
__ambiguous\t28
__too_low_aQual\t0
__not_aligned\t9956
__alignment_not_unique\t0
";

$actual = qx"less ${htseq_aln}";
unless(ok($expected eq $actual, "Check bwa count tables aln version.")) {
    my ($old, $new) = diff($expected, $actual);
    diag("expected:\n$old\nactual:\n$new\n");
}
