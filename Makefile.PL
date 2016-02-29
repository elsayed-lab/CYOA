use ExtUtils::MakeMaker;

WriteMakefile(
    NAME => 'CYOA',
    VERSION_FROM => 'lib/CYOA.pm',  ## Finds $VERSION
    PREREQ_PM => {
		   Log::Log4perl => 0,
		   Error => 0,
                   Bio::DB::GFF => 0,
		   Bio::DB::Universal => 0,
		   AppConfig => 0,
		   Getopt::Long => 0,
                   Cwd => 0,
                   File::Basename => 0,
                   Term::ReadLine => 0,
                   Archive::Extract => 0,
                   PerlIO => 0,
                   Digest::MD5 => 0,
                   Net::Amazon::S3 => 0,
                   common::sense => 0,
                   autodie => 0,
                   },
##    ABSTRACT_FROM => 'lib/HPGL.pm',
    EXE_FILES    => [ qw(bin/cyoa) ],
    AUTHOR => 'Ashton Trey Belew <abelew@gmail.com>'
    );

