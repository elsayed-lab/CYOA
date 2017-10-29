use strict;
use warnings;
use ExtUtils::MakeMaker;
use File::ShareDir::Install;

install_share module => 'Bio::Adventure' => 'share';

my $eumm_version = version->parse($ExtUtils::MakeMaker::VERSION);
my %prereq = (
    'Term::ProgressBar' => 2.18,
    'autodie' => 2.26,
    'AppConfig' => 1.71,
    'Archive::Extract' => 0.80,
    'Bio::FeatureIO' => 1.6,
    'Bio::DB::GFF' => 0,
    'Bio::DB::Universal' => 0,
    'Bio::DB::Sam' => 0,  ## I am not sure how to use this with travis.
    'Bio::SeqIO' => 0,
    'Bio::Tools::Run::StandAloneBlast' => 0,
    'Cwd' => 0.05,
    'Data::Dumper' => 2.022,
    'Digest::MD5' => 0.08,
    'Error' => 0.17,
    'File::Basename' => 0,
    'File::Copy' => 0.06,
    'File::Find' => 0.34,
    'File::Path' => 1.02,
    'File::Spec' => 0.073,
    'File::Which' => 1.21,
    'FileHandle' => 0,
    'Getopt::Long' => 2.49,
    'IO::String' => 1.08,
    'IO::Uncompress::UnXz' => 2.07,
    'IPC::System::Simple' => 1.25,
    'List::MoreUtils' => 0,
    'Log::Log4perl' => 1.49,
    'Math::Round' => 0.07,
    'Modern::Perl' => 2017,
    'Net::Amazon::S3' => 0.80,
    'PerlIO' => 0,
    'PerlIO::gzip' => 0,
    'Pod::Usage' => 1.69,
    'Storable' => 0.07,
    'String::Diff' => 0.07,
    'Term::ReadLine' => 1.03,
    'Term::UI' => 0.46,
    'Test::Directory' => 0,
    'Test::More' => 0,
    'Text::CSV' => 1.91,
    'TryCatch' => 0,
);

my @args = (
    NAME => 'Bio::Adventure',
    VERSION_FROM => 'lib/Bio/Adventure.pm',
    ABSTRACT => 'Choose your own adventure in sequence processing.',
    AUTHOR => 'abelew@gmail.com',
    LICENSE => 'perl',
    EXE_FILES => [ qw"bin/cyoa bin/cyoa.pl" ],
    test => { TESTS => 't/*.t' },
);

if ($eumm_version >= 6.46) {
    push @args,
    META_MERGE     => {
	requires   => { perl => '5.010001' },
	resources  => {
                license    => 'http://dev.perl.org/licenses/',
                bugtracker => 'https://github.com/abelew/CYOA/issues',
                repository => 'https://github.com/abelew/CYOA.git',
	},
    };
}

my %eumm_req = (
    'POSIX' => 0,
    'Time::Local' => 0,
    %prereq,
    );


if ($eumm_version >= 6.64) {
    push @args, PREREQ_PM => \%prereq, TEST_REQUIRES => \%eumm_req;
}
elsif ($eumm_version >= 6.5503) {
    push @args, PREREQ_PM => \%prereq, BUILD_REQUIRES => \%eumm_req;
}
else {
    push @args, PREREQ_PM => \%eumm_req;
}

WriteMakefile(@args);

package MY;
use File::ShareDir::Install 'postamble';
