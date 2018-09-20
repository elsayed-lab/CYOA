# -*-Perl-*-
use Test::More qw(no_plan);
use Bio::Adventure::MyDeps;

our @use_deps = (
    'local::lib',
    'autodie',
    'AppConfig',
    'Archive::Extract',
    'Bio::DB::Universal',
    'Bio::Tools::GFF',
    'Bio::Seq',
    'Bio::SeqIO',
    'Bio::Tools::Run::StandAloneBlast',
    'Cwd',
    'Data::Dumper',
    'Digest::MD5',
    'Error',
    'File::Basename',
    'File::Copy',
    'File::Find',
    'File::Path',
    'File::Spec',
    'File::Which',
    'FileHandle',
    'Getopt::Long',
    'IO::Uncompress::UnXz',
    'List::MoreUtils',
    'Log::Log4perl',
    'LWP',
    'Modern::Perl',
    'PerlIO',
    'Pod::Usage',
    'String::Approx',
    'String::Diff',
    'Term::ReadLine',
    'Test::More',
    );


$ENV{FTP_PASSIVE} = 1;
foreach my $module (@use_deps) {
    print "Loading $module\n";
    my $load_return = eval("use $module; 1");
    if (defined($module)) {
        ##my $version = $module->VERSION;
        ##print "Its version is: $version\n";
        print "Library: $module is installed.\n";
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
