use strict;
use warnings;
use ExtUtils::MakeMaker;
use File::ShareDir::Install;

install_share 'share';

WriteMakefile(
    NAME => 'Bio-Adventure',
    AUTHOR => 'Ashton Trey Belew <abelew@gmail.com>',
    VERSION_FROM => 'lib/Bio/Adventure.pm',
    ABSTRACT_FROM => 'lib/Bio/Adventure.pm',
    LICENSE => 'perl',
    PREREQ_PM => {
        'autodie' => '0',
        'AppConfig' => '0',
        'Archive::Extract' => '0',
        'Bio::FeatureIO' => '0',
        'Term::ProgressBar' => '0',
        'Bio::DB::GFF' => '0',
        'Bio::DB::Universal' => '0',
        'Bio::DB::Sam' => '0',
        'Bio::FeatureIO' => '0',
        'Bio::SeqIO' => '0',
        'Bio::Tools::Run::StandAloneBlast' => '0',
        'Bio::Tools::Run::WrapperBase' => '0',
        'Cwd' => '0',
        'Data::Dumper' => '0',
        'Data::Table' => '0',
        'Data::Table::Excel' => '0',
        'Digest::MD5' => '0',
        'Env::Modulecmd' => '0',
        'Error' => '0',
        'File::Basename' => '0',
        'File::Copy' => '0',
        'File::Find' => '0',
        'File::Path' => '0',
        'File::Spec' => '0',
        'File::Which' => '0',
        'FileHandle' => '0',
        'Getopt::Long' => '0',
        'IO::String' => '0',
        'IO::Uncompress::UnXz' => '0',
        'IPC::System::Simple' => '0',
        'List::MoreUtils' => '0',
        'Log::Log4perl' => '0',
        'Math::Round' => '0',
        'Modern::Perl' => '0',
        'LWP::Protocol::https' => '0',
        'Parse::CSV' => '0',
        'PerlIO' => '0',
        'PerlIO::gzip' => '0',
        'Pod::Usage' => '0',
        'Storable' => '0',
        'String::Approx' => '0',
        'String::Diff' => '0',
        'Switch' => '0',
        'Template' => '0',
        'Term::ReadLine' => '0',
        'Term::UI' => '0',
        'Test::Directory' => '0',
        'Test::More' => '0',
        'Text::CSV' => '0',
        'Text::CSV_XS::TSV' => '0',
     },
     TEST_REQUIRES => {
         'Test::More' => '0.47'
     },
     EXE_FILES => [
         'script/cyoa.pl'
     ],
    );

package MY;
use File::ShareDir::Install qw"postamble";
