# -*-Perl-*-
BEGIN {
    use Test::More qw"no_plan";
    use CYOA;
    use File::Path qw"remove_tree";
    use File::Copy qw"cp";
    use String::Diff qw( diff_fully diff diff_merge diff_regexp );
}
diag("Copy data file to cwd().");
ok(cp("t/data/test_forward.fastq.gz", "test_forward.fastq.gz"));
my $cyoa = new CYOA(input => qq"test_forward.fastq.gz", pbs => 0, species => 'phix', libdir => 't/data');
mkdir('t/data/genome/indexes'); ## Make a directory for the phix indexes.
diag("Does topthat execute?");
ok($cyoa->Tophat(),     'Run Tophat');
diag("Can I collect tophat statistics into tophat_stats.csv?");
ok(my $actual_csv = $cyoa->Last_Stat(input => 'outputs/tophat_stats.csv'),      'Collect Tophat Statistics');
diag("Does the last entry of tophat_stats.csv match the expected output?");
my $expected_csv = qq"test_forward,0,10000,40,9960,25000,accepted_hits.count.xz";
ok($actual_csv eq $expected_csv,      'Are the tophat results the expected value?');
my $expected_htseq = qq"phiX174p01\t1
phiX174p02\t0
phiX174p03\t0
phiX174p04\t0
phiX174p05\t0
phiX174p06\t12
phiX174p07\t0
phiX174p08\t0
phiX174p09\t0
phiX174p10\t0
phiX174p11\t0
__no_feature\t0
__ambiguous\t27
__too_low_aQual\t0
__not_aligned\t0
__alignment_not_unique\t0
";
my $actual_htseq = qx"xzcat outputs/tophat/accepted_hits.count.xz";
diag("Is the htseq output what was expected for phix?");
ok($expected_htseq eq $actual_htseq);

##my($old, $new) = String::Diff::diff($expected_htseq, $actual_htseq);
##print STDERR "$old\n";
##print STDERR "---\n";
##print STDERR "$new\n";
##print STDERR "---\n";


