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
ok(my $fastq_result = $hpgl->Last_Stat(input => 'stats/trimomatic_stats.csv'),      'Collect Trimomatic Statistics');
diag("Does the last entry of trimomatic_stats.csv match the expected output?");
ok($fastq_result eq 'test_forward,10000,9762,238',      'Are the trimomatic results the expected value?');
diag("Can I remove the trimomatic outputs?");
ok(unlink('test_forward-trimmed.fastq'),      'Remove the trimmed sequences');
diag("Can I clean up the mess?");
ok(remove_tree('outputs') and remove_tree('stats') and remove_tree('status') and remove_tree('sequences') and remove_tree('scripts'))
