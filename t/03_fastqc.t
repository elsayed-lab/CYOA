# -*-Perl-*-
BEGIN {
    use Test::More qw"no_plan";
    use HPGL;
    use File::Path qw"remove_tree";
}
my $hpgl = new HPGL(input => qq"t/data/test_forward.fastq.gz", pbs => 0);
diag("Does fastqc execute?");
ok($hpgl->Fastqc_Single(),     'Run Fastqc');
diag("Can I collect fastqc statistics into fastqc_stats.csv?");
ok(my $fastq_result = $hpgl->Last_Stat(input => 'outputs/fastqc_stats.csv'),      'Collect Fastqc Statistics');
diag("Does the last entry of fastqc_stats.csv match the expected output?");
ok($fastq_result eq 'test_forward,10000,0,pass,warn,pass,pass,pass,warn,fail,pass',      'Are the Fastqc results the expected value?');
