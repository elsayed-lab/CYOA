# -*-Perl-*-
BEGIN {
    use Test::More qw"no_plan";
    use Bio::Adventure;
    use File::Path qw"remove_tree";
    use File::Copy qw"cp";
}
my $cyoa = Bio::Adventure->new();

diag("Copying data file to cwd().");
ok(cp("t/data/test_forward.fastq.gz", "test_forward.fastq.gz"));

diag("Does trimomatic execute?");
ok(Bio::Adventure::RNASeq_Trim::Trimomatic_Single($cyoa, input => qq"test_forward.fastq.gz"),
   'Run Trimomatic');

diag("Can I collect trimomatic statistics into trimomatic_stats.csv?");
ok(my $actual = $cyoa->Last_Stat(input => 'outputs/trimomatic_stats.csv'),
   'Collect Trimomatic Statistics');

my $expected = 'test_forward,10000,9762,238',
diag("Does the last entry of trimomatic_stats.csv match the expected output?");
ok($expected eq $actual,
   'Are the trimomatic results the expected value?');
