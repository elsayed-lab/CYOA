BEGIN {
    use Test::More qw"no_plan";
    use HPGL;
    use File::Path qw"remove_tree";
}
my $hpgl = new HPGL(input => qq"t/data/test_forward.fastq.gz", pbs => 0);
diag("Does fastqc execute?");
ok($hpgl->Fastqc_Single(),     'Run Fastqc');
diag("Can I collect fastqc statistics into fastqc_stats.csv?");
ok(my $fastq_result = $hpgl->Fastqc_Last_Stat(),      'Collect Fastqc Statistics');
diag("Does the last entry of fastqc_stats.csv match the expected output?");
ok($fastq_result eq 'test_forward,10000,0,pass,warn,pass,pass,pass,warn,fail,pass',      'Are the Fastqc results the expected value?');
diag("Can I remove the fastqc output tree?");
ok(remove_tree('test_forward-unfiltered_fastqc'),      'Remove the Fastqc output directory');
diag("Can I remove the fastqc_stats.csv file?");
ok(unlink('fastqc_stats.csv'),      'Remove the fastqc_stats.csv output file.');
