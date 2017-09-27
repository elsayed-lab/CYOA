# -*-Perl-*-
BEGIN {

}
use Test::More qw"no_plan";
use_ok(Bio::Adventure);
my $cyoa = new Bio::Adventure;
my $help = $cyoa->Help();
ok($help == 0, 'help works');
