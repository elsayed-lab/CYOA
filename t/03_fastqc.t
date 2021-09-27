# -*-Perl-*-
use Test::More qw"no_plan";
use Bio::Adventure;
use File::Path qw"remove_tree";
use String::Diff qw" diff_fully diff diff_merge diff_regexp ";

my $cyoa = Bio::Adventure->new(cluster => 0);

ok($cyoa->Bio::Adventure::QA::Fastqc_Single(input => qq"share/test_forward.fastq.gz"),
   'Run Fastqc');

ok(-r 'scripts/01fqc_test_forward_share.sh',
   'Fastqc script exists?');

ok(qx"fastqc --help",
   'Can run fastqc?');

ok(my $actual = $cyoa->Last_Stat(input => 'outputs/fastqc_stats.csv'),
   'Collect Fastqc Statistics');

my $expected = 'fqc_test_forward_share,10000,0,pass,warn,pass,pass,pass,warn,fail,0';
unless(ok($expected eq $actual,
          'Are the fastqc results the expected value?')) {
    my($old, $new) = diff($expected, $actual);
    diag("$old\n$new\n");
}
