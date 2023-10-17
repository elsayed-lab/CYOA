# -*-Perl-*-
use strict;
use Test::More qw"no_plan";
use Bio::Adventure;
use Cwd;
use File::Copy qw"cp mv";
use File::Path qw"remove_tree";
use File::Basename qw"basename";
use File::ShareDir qw"dist_file module_dir dist_dir";
use String::Diff qw"diff";
use Test::File::ShareDir::Dist { 'Bio-Adventure' => 'share/' };

my $start_dir = dist_dir('Bio-Adventure');
my $input_file = qq"${start_dir}/test_forward.fastq.gz";
my $input_local = basename($input_file);
my $start = getcwd();
my $new = 'test_output';
mkdir($new);
chdir($new);
my $cyoa = Bio::Adventure->new(basedir => cwd());

if (!-r $input_local) {
    ok(cp($input_file, $input_local), 'Copying data.');
}

my $trimmer = $cyoa->Bio::Adventure::Trim::Cutadapt(input => $input_local);
ok($trimmer, 'Submit cutadapt job.');
my $status = $cyoa->Wait(job => $trimmer);
ok($status->{State} eq 'COMPLETED', 'Cutadapt completed on cluster.');
ok(-r 'scripts/12cutadapt_test_forward.sh',
   'Cutadapt script exists?');
ok(my $actual = $cyoa->Last_Stat(input => 'outputs/cutadapt_stats.csv'),
   'Collect cutadapt Statistics');
my $expected = 'cutst,10000,2240,11,7760,2229';
unless(ok($expected eq $actual,
          'Are the cutadapt results the expected value?')) {
    my($old, $new) = diff($expected, $actual);
    diag("$old\n$new\n");
}

chdir($start);
