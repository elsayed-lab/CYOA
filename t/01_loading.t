# -*-Perl-*-
use Test::More qw"no_plan";
use File::Path qw"remove_tree";

if (-d 'scripts') {
    ok(remove_tree('scripts'),
       'Removed scripts/');
}
if (-d 'outputs') {
    ok(remove_tree('outputs'),
       'Removed outputs/');
}

use_ok(Bio::Adventure);
my $cyoa = new Bio::Adventure;
my $help = $cyoa->Help();
ok($help == 0, 'help works');
