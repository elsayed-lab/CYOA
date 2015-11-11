BEGIN {
    use Test::More qw"no_plan";
    use HPGL;
    use File::Path qw"remove_tree";
    use File::Copy qw"cp";
    use String::Diff qw( diff_fully diff diff_merge diff_regexp );
    use Archive::Extract;
}
diag("Copy data file to cwd().");
ok(cp("t/data/test_forward.fastq.gz", "test_forward.fastq.gz"));
my $hpgl = new HPGL(input => qq"test_forward.fastq.gz", pbs => 0, species => 'phix', libdir => 't/data');
mkdir('t/data/genome/indexes'); ## Make a directory for the phix indexes.
diag("Does bwa execute?");
ok($hpgl->BWA(),     'Run Bwa');
diag("Can I collect bwa statistics into bwa_stats.csv?");
ok(my $actual_csv = $hpgl->Last_Stat(input => 'outputs/bwa_stats.csv'),      'Collect Tophat Statistics');
diag("Does the last entry of bwa_stats.csv match the expected output?");
my $expected_csv = qq"test_forward,0,10000,44,49,0,";
ok($actual_csv eq $expected_csv,      'Are the tophat results the expected value?');

my $expected_htseq_mem = qq"phiX174p01\t5
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
my $expected_htseq_aln = qq"phiX174p01\t3
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
my $actual_htseq_mem = qx"xzcat outputs/bwa/test_forward_mem.count.xz";
diag("Is the htseq output for the mem alignment expected?");
ok($expected_htseq_mem eq $actual_htseq_mem);
my $actual_htseq_aln = qx"xzcat outputs/bwa/test_forward_aln.count.xz";
diag("Is the htseq output for the aln alignment expected?");
ok($expected_htseq_aln eq $actual_htseq_aln);

##my($old, $new) = String::Diff::diff($expected_htseq, $actual_htseq);
##print STDERR "$old\n";
##print STDERR "---\n";
##print STDERR "$new\n";
##print STDERR "---\n";


