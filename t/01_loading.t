# -*-Perl-*-
BEGIN {

}
use Test::More qw"no_plan";
use_ok(HPGL);
my $hpgl = new HPGL;
my $help = $hpgl->Help();
ok($help == 0, 'help works');
