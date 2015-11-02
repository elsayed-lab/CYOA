BEGIN {
    use Test::More qw"no_plan";
    use HPGL;
    use File::Path qw"remove_tree";
    use File::Copy qw"cp";
}
diag("Copying data file to cwd().");
ok(cp("t/data/test_forward.fastq.gz", "test_forward.fastq.gz"));
my $hpgl = new HPGL(input => qq"test_forward.fastq.gz", pbs => 0);
diag("Does trimomatic execute?");
ok($hpgl->Trimomatic_Single(),     'Run Trimomatic');
diag("Can I collect trimomatic statistics into trimomatic_stats.csv?");
ok(my $trim_result = $hpgl->Last_Stat(input => 'outputs/trimomatic_stats.csv'),      'Collect Trimomatic Statistics');
diag("Does the last entry of trimomatic_stats.csv match the expected output?");
ok($trim_result eq 'test_forward,10000,9762,238',      'Are the trimomatic results the expected value?');
