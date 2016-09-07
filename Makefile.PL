use strict;
use warnings;
use ExtUtils::MakeMaker;

my $eumm_version = version->parse($ExtUtils::MakeMaker::VERSION);
my $prereq = {
    'autodie' => 0,
    'common::sense' => 0, ## likely going to be removed in favor of modern::perl
    'local::lib' => 0,
    'AppConfig' => 0,
    'Archive::Extract' => 0,
    'Bio::FeatureIO' => 0,
    'Bio::DB::GFF' => 0,
    'Bio::DB::Universal' => 0,
##    'Bio::DB::Sam' => 0,  ## I am not sure how to use this with travis.
    'Bio::SeqIO' => 0,
    'Cwd' => 0,
    'Data::Dumper' => 0,
    'Digest::MD5' => 0,
    'Error' => 0,
    'File::Basename' => 0,
    'File::Copy' => 0,
    'File::Find' => 0,
    'File::Path' => 0,
    'File::Spec' => 0,
    'File::Which' => 0,
    'FileHandle' => 0,
    'Getopt::Long' => 0,
    'IO::String' => 0,
    'IO::Uncompress::UnXz' => 0,
    'IPC::System::Simple' => 0,
    'List::MoreUtils' => 0,
    'Modern::Perl' => 0,
    'Net::Amazon::S3' => 0,
    'PerlIO' => 0,
    'PerlIO::gzip' => 0,
    'Pod::Usage' => 0,
    'Storable' => 0,
    'String::Diff' => 0,
    'Term::ReadLine' => 0,
    'Term::UI' => 0,
    'Test::More' => 0,
    'Text::CSV' => 0,
};

my @args = (
    NAME => 'CYOA',
    VERSION_FROM => 'lib/CYOA.pm',
    ABSTRACT => 'Choose your own adventure in sequence processing.',
    AUTHOR => 'abelew@gmail.com',
    LICENSE => 'gplv2',
    EXE_FILES => [ qw"bin/cyoa" ],
    test => {TESTS => 't/*.t'},);

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
    %$prereq,
    );


if ($eumm_version >= 6.64) {
    push @args, PREREQ_PM => $prereq, TEST_REQUIRES => \%eumm_req;
}
elsif ($eumm_version >= 6.5503) {
    push @args, PREREQ_PM => $prereq, BUILD_REQUIRES => \%eumm_req;
}
else {
    push @args, PREREQ_PM => \%eumm_req;
}

WriteMakefile(@args);
