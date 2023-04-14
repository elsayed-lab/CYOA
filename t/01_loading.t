# -*-Perl-*-
use Test::More qw"no_plan";
use File::Path qw"remove_tree";
my $test_directory = 'test_output';
if (-d $test_directory) {
    ok(remove_tree($test_directory), 'Removed previous test directory.');
}

use_ok(Bio::Adventure);
my $cyoa = new Bio::Adventure;
ok($cyoa, 'Able to load the parent Adventure class.');

## I might remove this, on my shell it is sending the output to less and therefore
## making me hit 'q' before it will continue, which is annoying, though I guess it would
## be nice if I actually wanted to look at the help.
my $help = $cyoa->Help();
ok($help == 0, 'help works');
