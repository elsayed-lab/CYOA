# -*-Perl-*-
BEGIN {
    use Test::More qw"no_plan";
    use CYOA;
    use File::Path qw"remove_tree";
    use File::Copy qw"cp";
    use String::Diff qw( diff_fully diff diff_merge diff_regexp );
}
diag("Copying data file to cwd().");
ok(cp("t/data/test_forward.fastq.gz", "test_forward.fastq.gz"));
my $cyoa = new CYOA(input => qq"test_forward.fastq.gz", pbs => 0, species => 'phix', libdir => 't/data');
mkdir('t/data/genome/indexes'); ## Make a directory for the phix indexes.
diag("Does bowtie execute?");
ok($cyoa->Bowtie(),     'Run Bowtie1');
diag("Can I collect bowtie statistics into bowtie_stats.csv?");
ok(my $actual_csv = $cyoa->Last_Stat(input => 'outputs/bowtie_stats.csv'),      'Collect Bowtie1 Statistics');
diag("Does the last entry of bowtie_stats.csv match the expected output?");
my $expected_csv = qq"test_forward,v0M1,10000,10000,30,9970,0,33333.3333333333,test_forward-v0M1.count.xz";
ok($actual_csv eq $expected_csv,      'Are the bowtie results the expected value?');
my $expected_htseq = qq"phiX174p01\t1
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
my $actual_htseq = qx"xzcat outputs/bowtie/test_forward-v0M1.count.xz";
diag("Is the htseq output what was expected for phix?");
ok($expected_htseq eq $actual_htseq);

##my($old, $new) = String::Diff::diff($expected_htseq, $actual_htseq);
##print STDERR "$old\n";
##print STDERR "---\n";
##print STDERR "$new\n";
##print STDERR "---\n";


