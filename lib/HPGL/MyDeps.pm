package MyDeps;
use CPAN;
use strict;
use Test::More;

use vars qw($VERSION);
$VERSION='20151101';

#use vars qw($VERSION);
#our @ISA = qw(Exporter);
#our @EXPORT = qw(deps);    # Symbols to be exported by default

our @use_deps = (
    'Bio::Seq',
    'Bio::SeqIO',
    'Bio::Tools::Run::StandAloneBlast',
    'Error',
    'Getopt::Long',
    'Log::Log4perl',
    );

sub Test {
    foreach my $d (@use_deps) {
	my $response = use_ok($d);
	diag("Testing usability of $d\n");
	if ($response != 1) {
	    diag("$d appears to be missing.  Please run fixdeps.pl\n");
	}
    }
}

sub Res {
    $ENV{FTP_PASSIVE} = 1;
    foreach my $module (@use_deps) {
        print "Loading $module\n";
        my $load_return = eval("use $module; 1");
        if (defined($module)) {
            my $version = $module->VERSION;
            print "Its version is: $version\n";
        }
        if (!$load_return) {
            my $dep_count = 0;
            while ($dep_count <= 50) {
                $dep_count++;
                my $missing = $@;
                if ($missing =~ /Global/) {
                    print "syntax error found in $module.\n";
                    print "$@";
                } elsif ($missing =~ /^Can\'t locate/) {
                    my @response = split(/ /, $missing);
                    $missing = $response[2];
                    $missing =~ s/\//::/g;
                    $missing =~ s/\.pm//g;
                    print "Going to install $missing\n";
                    CPAN::Shell->force("install", $missing);
                }
            }
        }
    }
}

1;
