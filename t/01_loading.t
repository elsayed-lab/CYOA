# -*-Perl-*-
BEGIN {

}
use Test::More qw"no_plan";
use_ok(CYOA);
my $cyoa = new CYOA;
my $help = $cyoa->Help();
ok($help == 0, 'help works');
