# -*-Perl-*-
use Test::More qw"no_plan";
use Bio::Adventure;
use Bio::SeqIO;
use File::Path qw"remove_tree";
use File::Copy qw"cp";
use String::Diff qw"diff";

## Travis does not have bcftools in its repository.
if ($ENV{TRAVIS}) {
    print STDERR "Skipping these tests as bcftools is not on travis.\n";
    ok(1);
} else {
    my $cyoa = Bio::Adventure->new();
    ok(cp('share/test_forward.fastq.gz', 'test_forward.fastq.gz'),
       'Copy the data files.');
    
    mkdir('share/genome/indexes'); ## Make a directory for the phix indexes.
    ok(Bio::Adventure::SNP::Align_SNP_Search(
           $cyoa,
           input => qq"test_forward.fastq.gz",
           htseq_id => 'ID',
           htseq_type => 'CDS',
           libdir => 'share',
           species => 'phix_3N',
           vcf_cutoff => 1,
       ),
       'Run SNP_Search_Align');


    ## There are a couple other places which do not agree with the canonical genome it appears.
    my $original = Bio::SeqIO->new(-file => 'share/genome/phix_modified.fasta');
    my $new = Bio::SeqIO->new(-file => 'outputs/vcfutils_phix_3N/CYOA_phix_3N_modified.fasta');
    my $expected = $original->next_seq();
    my $actual = $new->next_seq();
    $expected = $expected->seq;
    $actual = $actual->seq;
    unless(ok($expected eq $actual,
              'Did we get back the old genome?')) {
        my($old, $new) = diff($expected, $actual);
        diag("$old\n$new\n");
    }
}
