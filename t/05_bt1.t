# -*-Perl-*-
BEGIN {
    use Test::More qw"no_plan";
    use Bio::Adventure;
    use File::Path qw"remove_tree";
    use File::Copy qw"cp";
    use String::Diff qw( diff_fully diff diff_merge diff_regexp );
}

my $cyoa = Bio::Adventure->new();
diag("Copying data file to cwd().");
ok(cp("t/data/test_forward.fastq.gz", "test_forward.fastq.gz"));

mkdir('t/data/genome/indexes'); ## Make a directory for the phix indexes.
diag("Does bowtie execute?");
ok(Bio::Adventure::RNASeq_Map::Bowtie($cyoa,
                                      input => qq"test_forward.fastq.gz",
                                      species => 'phix',
                                      htseq_id => 'ID',
                                      htseq_type => 'CDS',
                                      libdir => 't/data'),
   'Run Bowtie1');

diag("Can I collect bowtie statistics into bowtie_stats.csv?");
ok(my $actual = $cyoa->Last_Stat(input => 'outputs/bowtie_stats.csv'),
   'Collect Bowtie1 Statistics');

diag("Does the last entry of bowtie_stats.csv match the expected output?");
my $expected = qq"CYOA2,v0M1,10000,10000,30,9970,0,33333.3333333333,CYOA2-v0M1.count.xz";
ok($actual eq $expected,
   'Are the bowtie results the expected value?');
##my($old, $new) = String::Diff::diff($expected, $actual);
##print STDERR "$old\n";
##print STDERR "---\n";
##print STDERR "$new\n";
##print STDERR "---\n";

$expected = qq"phiX174p01\t1
phiX174p02\t0
phiX174p03\t0
phiX174p04\t0
phiX174p05\t0
phiX174p06\t8
phiX174p07\t0
phiX174p08\t0
phiX174p09\t0
phiX174p10\t0
phiX174p11\t0
__no_feature\t0
__ambiguous\t21
__too_low_aQual\t0
__not_aligned\t9970
__alignment_not_unique\t0
";
$actual = qx"xzcat outputs/bowtie_phix/CYOA2-v0M1.count.xz";
diag('Is the htseq output what was expected for phix?');
ok($expected eq $actual);
