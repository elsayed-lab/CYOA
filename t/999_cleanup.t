# -*-Perl-*-
use Test::More qw"no_plan";
use Bio::Adventure;
use File::Path qw"remove_tree";

ok(remove_tree('scripts'),
    'Removed scripts/');
ok(remove_tree('t/data/genome/indexes'),
    'Removed indexes directory.');

my @files = ('phix.fasta.pdata','phix.fasta.pdata','phix_3N.gff.pdata',
             'phix_3N.fasta.pdata', 'phix.gff', 'phix.fasta', 'phix_cds_nt.fasta',
             'phix_cds_aa.fasta', 'split_align_errors.txt', 'test_forward.fastq',
             'test_forward.fastq.gz', 'test_forward.fastq.xz', 'test_forward-trimmed.fastq',
             'test_forward.log', 'xzal.log', 'xzun.log');

for my $f (@files) {
    ok(unlink($f),
       "Removed ${f}.");
}

