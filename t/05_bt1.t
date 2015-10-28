BEGIN {
    use Test::More qw"no_plan";
    use HPGL;
    use File::Path qw"remove_tree";
    use File::Copy qw"cp";
}
diag("Copying data file to cwd().");
ok(cp("t/data/test_forward.fastq.gz", "test_forward.fastq.gz"));
my $hpgl = new HPGL(input => qq"test_forward.fastq.gz", pbs => 0, species => 'phix', libdir => 't/data');
diag("Does bowtie execute?");
ok($hpgl->Bowtie(),     'Run Bowtie1');
#diag("Can I collect bowtie statistics into bowtie_stats.csv?");
#ok(my $bt1_result = $hpgl->Last_Stat(input => 'stats/bowtie_stats.csv'),      'Collect Bowtie1 Statistics');
#diag("Does the last entry of bowtie_stats.csv match the expected output?");
#ok($bt1_result eq 'asdsadadasdasd',      'Are the bowtie results the expected value?');
#diag("Can I remove the bowtie outputs?");
#ok(unlink('test_forward.fastq'),      'Remove the sequences');
#diag("Can I clean up the mess?");
#ok(remove_tree('outputs') and remove_tree('stats') and remove_tree('status') and remove_tree('sequences') and remove_tree('scripts'))
