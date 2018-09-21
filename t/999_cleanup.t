# -*-Perl-*-
use Test::More qw"no_plan";
use Bio::Adventure;
use File::Path qw"remove_tree";

ok(remove_tree('scripts'),
    'Removed scripts/');
ok(remove_tree('share/genome/indexes'),
    'Removed indexes directory.');

my @all_files = (
    'test_forward-trimmed.fastq',
    );
for my $f (@all_files) {
    ok(unlink($f),
       "Removed ${f}.");
}

my @not_travis_files = (
    'split_align_errors.txt', 'phix.fasta', 
    'phix.gff', 'phix_cds_nt.fasta', 'phix_cds_aa.fasta',
    );
unless ($ENV{TRAVIS}) {
    for my $g (@not_travis_files) {
        ok(unlink($g),
           "Removed ${g}.");
    }
}
