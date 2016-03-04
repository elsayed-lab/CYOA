use ExtUtils::MakeMaker;

WriteMakefile(
    NAME => 'CYOA',
    VERSION_FROM => 'lib/CYOA.pm',  ## Finds $VERSION
    PREREQ_PM => {
		   AppConfig => 0,
                   Archive::Extract => 0,
                   autodie => 0,
                   Bio::FeatureIO => 0,
                   Bio::DB::GFF => 0,
		   Bio::DB::Universal => 0,
                   Bio::DB::Sam => 0,
                   Bio::SeqIO => 0,
                   common::sense => 0,
                   Cwd => 0,
                   Digest::MD5 => 0,
		   Error => 0,
                   File::Basename => 0,
                   File::Find => 0,
                   File::Path => 0,
                   File::Spec => 0,
                   File::Which => 0,
                   Getopt::Long => 0,
                   IO::String => 0,
                   IO::Uncompress::UnXz => 0,
                   IPC::System::Simple => 0,
                   Net::Amazon::S3 => 0,
                   PerlIO => 0,
                   Pod::Usage => 0,
                   Storable => 0,
                   Term::ReadLine => 0,
                   Term::UI => 0,
                   Text::CSV => 0,
                   },
##    ABSTRACT_FROM => 'lib/HPGL.pm',
    EXE_FILES    => [ qw(bin/cyoa) ],
    AUTHOR => 'Ashton Trey Belew <abelew@gmail.com>'
    );

