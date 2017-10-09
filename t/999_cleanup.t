# -*-Perl-*-
use Test::More qw"no_plan";
use Bio::Adventure;
use File::Path qw"remove_tree";

ok(remove_tree('scripts'),
    'Removed scripts/');
ok(remove_tree('share/genome/indexes'),
    'Removed indexes directory.');

my @all_files = (
    'phix.fasta.pdata','phix_3N.gff.pdata', 'phix_3N.fasta.pdata', 'phix.gff',
    'phix.fasta', 'phix_cds_nt.fasta', 'phix_cds_aa.fasta', 'test_forward.fastq',
    'test_forward.fastq.xz', 'test_forward-trimmed.fastq', 'test_forward.log',
    'xzal.log', 'xzun.log',
    );
for my $f (@all_files) {
    ok(unlink($f),
       "Removed ${f}.");
}

my @not_travis_files = ('split_align_errors.txt');
unless ($ENV{TRAVIS}) {
    for my $g (@not_travis_files) {
        ok(unlink($g),
           "Removed ${g}.");
    }
}
