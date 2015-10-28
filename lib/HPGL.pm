package HPGL;
use AppConfig qw":argcount :expand";
use Archive::Extract;
use autodie qw":all";
use common::sense;
use Bio::DB::Universal;
use Bio::Root::RootI;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Cwd qw"cwd";
use Data::Dumper;
use Digest::MD5 qw"md5 md5_hex md5_base64";
use File::Basename;
use File::Find;
use File::Which qw"which";
use File::Path;
use FileHandle;
use Getopt::Long;
use HPGL::SeqMisc;
use Log::Log4perl;
use Log::Log4perl::Level;
use Net::Amazon::S3;
use PerlIO;
use Pod::Usage;
use Term::ReadLine;
use warnings qw"all";

use HPGL::Compress;
use HPGL::RNASeq_QA;
use HPGL::RNASeq_Trim;
use HPGL::RNASeq_Aligners;
use HPGL::RNASeq_Count;
use HPGL::Convert;
use HPGL::PBS;

our $AUTOLOAD;
require Exporter;
our @ISA = qw"Exporter";
our @EXPORT = qw"";
use vars qw"$VERSION";
$VERSION="20151101";
$ENV{COMPRESSION} = "xz";
$ENV{XZ_OPTS} = "-9e";
$ENV{XZ_DEFAULTS} = "-9e";
$ENV{GZIP} = "--best";
my $term = new Term::ReadLine('>');
my $attribs = $term->Attribs;
$attribs->{completion_suppress_append} = 1;
my $OUT = $term->OUT || \*STDOUT;

=head1 NAME

    HPGL - A perl library to make preprocessing high-throughput data easier.

=head1 SYNOPSIS

    HPGL.pm HPGL::Job.pm Some perl libraries to submit jobs to the
    umiacs cluster.  This can take a variety of options:

  --debug|d   : Get it to print some debugging information.
  --btmulti   : Perform multiple bowtie option sets.
  --tpmulti   : Perform multiple tophat option sets (not implemented).
  --bwmulti   : Perform multiple bwa option sets (not implemented).
  --help      : Print some help information.
  --input|i   : Input file(s).
  --pbs|p     : Use Torque?
  --species:s : Species used for alignments.
  --stranded  : Are the libraries stranded?
  --identifier: The identifier tag in the annotation gff

    use HPGL;
    my $hpgl = new HPGL;
    $hpgl->{species} = 'mmusculus';  ## There had better be a mmusculus.[fasta&gff] in $hpgl->{libdir}
    my $hpgl = new HPGL(species => 'scerevisiae', libdir => '/home/bob/libraries');
    $hpgl->{input} = 'test.fastq';  ## Also via with script.pl --input=test.fastq

    ## Run Fastqc on untrimmed sequence.
    my $bp = $hpgl->Fastqc()
    ## Generates test-trimmed.fastq from test.fastq
    my $trim = $hpgl->Trimomatic();
    ## Graph statistics from test-trimmed.fastq
    my $bp = $hpgl->Biopieces_Graph(depends => $trim->{pbs_id});
    ## Run bowtie1, convert the output to sorted/indexed bam, and run htseq-count against mmusculus
    ## bowtie1 outputs go (by default) into bowtie_out/
    my $bt = $hpgl->Bowtie();
    ## Run htseq-count with a different gff file, use 'mRNA' as the 3rd column of the gff, and 'locus_tag' as the identifier in the last column.
    my $ht = $hpgl->HTSeq(depends => $bt->{pbs_id}, input => $bt->{output}, htseq_gff => 'other.gff', htseq_type => 'mRNA', htseq_identifier => 'locus_tag'
    ## Run tophat
    my $tp = $hpgl->TopHat(depends => $trim->{pbs_id});
    ## or BWA
    my $bw = $hpgl->BWA(depends => $trim->{pbs_id});

=head1 DESCRIPTION

    This library should write out PBS compatible job files for torque and
    submit them to appropriate queues on the cluster.  It should also
    collect the outputs and clean up the mess.

=head1 AUTHOR - atb

Email abelew@gmail.com

=cut

=head2
    new() instantiates a new HPGL object.
    It has a plethora of options which are either pulled from
    GetOpt::Long or via a default hash reference -- I probably should
    move some of those options from the default reference and give
    them flags on the command line.

    Since each tool used herein also can take an annoyingly large set
    of options, they all pull %args from @_ and use that for more
    information.

    There are a few things which are commonly passed, but not necessarily in new()
    The following code snippets describe a couple:
=cut
sub new {
    my ($class, %args) = @_;
    my $me = bless {}, $class;
    foreach my $key (keys %args) {
        $me->{$key} = $args{$key} if (defined($args{$key}));
    }
    $me->{appconfig} = new AppConfig({
        CASE => 1,
        CREATE => 1,
        PEDANTIC => 0,
        ## ERROR => eval(),
        GLOBAL => {
            EXPAND => EXPAND_ALL,
            EXPAND_ENV => 1,
            EXPAND_UID => 1,
            DEFAULT => "unset",
            ARGCOUNT => 1,
        },});
    ## First fill out a set of default configuration values
    ## The following for loop will include all of these in conf_specification_temp
    ## which will in turn be passed to GetOpts
    ## As a result, _all_ of these variables may be overridden on the command line.
    $me->{shell} = "/usr/bin/bash" if (!defined($me->{shell}));
    $me->{align_trimmed} = 1 if (!defined($me->{align_trimmed}));
    $me->{aligner} = "bowtie" if (!defined($me->{aligner}));
    $me->{aws_access_key_id} = 'PQOD56HDPQOXJPYQYHH9' if (!defined($me->{aws_access_key_id}));
    $me->{aws_hostname} = 'gembox.cbcb.umd.edu' if (!defined($me->{aws_hostname}));
    $me->{aws_secret_access_key} = 'kNl27xjeVSnChzEww9ziq1VkgUMNmNonNYjxWkGw' if (!defined($me->{aws_secret_access_key}));
    $me->{bam_small} = "21-26" if (!defined($me->{bam_small}));
    $me->{bam_mid} = "27-32" if (!defined($me->{bam_mid}));
    $me->{bam_big} = "33-40" if (!defined($me->{bam_big}));
    $me->{base} = undef if (!defined($me->{base}));
    $me->{basedir} = cwd() if (!defined($me->{basedir}));
    $me->{bt_args} = {
        v0M1 => "--best -v 0 -M 1",
        v1M1 => "--best -v 1 -M 1",
        v2M1 => "--best -v 2 -M 1",
    } if (!defined($me->{bt_args}));
    $me->{btmulti} = 0 if (!defined($me->{btmulti}));
    $me->{config_file} = qq"$ENV{HOME}/.config/hpgl.conf" if (!defined($me->{config_file}));
    $me->{debug} = 1 if (!defined($me->{debug}));
    $me->{depends} = {} if (!defined($me->{depends}));
    $me->{help} = undef if (!defined($me->{help}));
    $me->{hpgl} = undef if (!defined($me->{hpgl}));
    $me->{hpglid} = undef if (!defined($me->{hpglid}));
    $me->{htseq_identifier} = "ID" if (!defined($me->{htseq_identifier}));
    $me->{htseq_stranded} = "no" if (!defined($me->{htseq_stranded}));
    $me->{input} = undef if (!defined($me->{input}));
    $me->{jobs} = [] if (!defined($me->{jobs}));
    $me->{jobids} = {} if (!defined($me->{jobids}));
    $me->{len_min} = 21 if (!defined($me->{len_min}));
    $me->{len_max} = 40 if (!defined($me->{len_max}));
    $me->{libdir} = "$ENV{HOME}/libraries" if (!defined($me->{libdir}));
    $me->{libtype} = "genome" if (!defined($me->{libtype}));
    $me->{paired} = 0 if (!defined($me->{paired}));
    $me->{pbs} = 1 if (!defined($me->{pbs}));
    $me->{qsub_args} = "-j oe -V -m n" if (!defined($me->{qsub_args}));
    $me->{qsub_queue} = "throughput" if (!defined($me->{qsub_queue}));
    $me->{qsub_queues} = ["throughput","workstation","long","large"] if (!defined($me->{qsub_queues}));
    $me->{qsub_shell} = "bash" if (!defined($me->{qsub_shell}));
    $me->{qsub_mem} = 6 if (!defined($me->{qsub_mem}));
    $me->{qsub_wall} = "10:00:00" if (!defined($me->{qsub_wall}));
    $me->{qsub_cpus} = "4" if (!defined($me->{qsub_cpus}));
    $me->{qsub_depends} = "depend=afterok:" if (!defined($me->{qsub_depends}));
    $me->{qsub_loghost} = "ibissub00.umiacs.umd.edu" if (!defined($me->{qsub_loghost}));
    $me->{qsub_logdir} = qq"$me->{qsub_loghost}:$me->{basedir}/outputs" if (!defined($me->{qsub_logdir}));
    $me->{species} = undef if (!defined($me->{species}));
    $me->{suffixes} = [".fastq",".gz",".xz", ".fasta", ".sam", ".bam", ".count"] if (!defined($me->{suffixes}));
    $me->{trimmer} = "trimmomatic" if (!defined($me->{trimmer}));
    $me->{type} = "rnaseq" if (!defined($me->{type}));

    ## Now that the set of defaults has been created, allow it to be changed via a config file
    ## All the stuff above may therefore be overwritten by entries in the config file.
    ## which defaults to ~/.config/hpgl.conf
    my ($open, %data, $config_option);
    if (-r $me->{config_file}) {
        $open = $me->{appconfig}->file($me->{config_file});
        %data = $me->{appconfig}->varlist("^.*");
        for $config_option (keys %data) {
            $me->{$config_option} = $data{$config_option};
            undef $data{$config_option};
        }
    }

    ## Options set in the config file are trumped by those placed on the comamnd line.
    my (%conf, %conf_specification, %conf_specification_temp);
    foreach my $default (keys %{$me}) {
        $conf_specification{"$default:s"} = \$conf{$default};
    }
    ## For some command line options, one might want shortcuts etc, set those here:
    %conf_specification_temp = (
        "debug|d" => \$conf{debug},
        "btmulti:i" => \$conf{btmulti},
        "help" => \$conf{help},
        "hpgl|h:s" => \$conf{hpgl},
        "input|i:s" => \$conf{input},
        "pbs|p:i" => \$conf{pbs},
        "species|s:s" => \$conf{species},
        "stranded:s" => \$conf{htseq_stranded},
        "identifier:s" => \$conf{htseq_identifier},
        );
    ## This makes both of the above groups command-line changeable
    foreach my $name (keys %conf_specification_temp) {
        $conf_specification{$name} = $conf_specification_temp{$name};
    }
    undef(%conf_specification_temp);
    ## Now pull out the command line options.
    my $argv_result = GetOptions(%conf_specification);
    unless ($argv_result) {
        Help();
    }

    ## Finally merge the options from the config file to those
    ## in the default hash and those from the command line.
    foreach my $opt (keys %conf) {
        if (defined($conf{$opt})) {
            $me->{$opt} = $conf{$opt};
        }
    }
    undef(%conf);

    Help() if (defined($me->{help}));

    ##    $me->Check_Options(["input",]);
    if ($me->{input}) {
        my $base = $me->{input};
        my @suffixes = @{$me->{suffixes}};
        $base = basename($base, @suffixes);
        $base = basename($base, @suffixes);
        $me->{basename} = $base;
        if (!defined($me->{hpglid})) {
            my $tmp = $me->{input};
            $tmp =~ s/^(hpgl\d+).*/$1/g;
            $me->{hpglid} = $tmp;
        }
    }

    if ($me->{pbs}) {
        my $qsub_path = which 'qsub';
        $me->{pbs} = 0 unless ($qsub_path);
    }

    return($me);
}

=head2
    Help()
=cut
sub Help {
    my $me = shift;
    my $fh = \*STDERR;
    ##    my $usage = pod2usage(-output => $fh, -verbose => 99, -sections => "SYNOPSIS");
    use Pod::Find qw(pod_where);
    pod2usage(-input => pod_where({-inc => 1}, __PACKAGE__), -verbose => 2, -output => $fh, -sections => "NAME|SYNOPSIS|DESCRIPTION|VERSION", -exitval => 'NOEXIT');
    ##my $usage = pod2usage(-verbose => 2, -output => $fh, -sections => "NAME|SYNOPSIS|DESCRIPTION|VERSION", -exitval => 'NOEXIT');
    ##print STDERR "Ran pod2usage\n";
    ##print Dumper $usage;
    return(0);
}


=head2
    Check_Options()
=cut
sub Check_Options {
    my $me = shift;
    my $needed_list = shift;
    my @options = @{$needed_list};
    foreach my $option (@options) {
        if (!defined($me->{$option})) {
            my $query = qq"The option:${option} was missing, please fill it in: ";
            $me->{$option} = $term->readline($query);
            $me->{$option} =~ s/\s+$//g;
        }
    }
}

=head2
    Get_Input()
=cut
sub Get_Input {
    my $me = shift;
    my $input = $me->{input};
    my $id = $me->{hpgl};
    my $actual = "";
    ## There are a few problems with how I send input to these scripts
    ## Sometimes I put in --hpgl hpgl0415 when I mean --hpgl HPGL0415
    ## Sometimes I put in --input hpgl0415.fastq  when I mean hpgl0415.fastq.(gz|xz)
    ## Sometimes I put in -i hpgl0415.fastq.(gz|xz) when I mean hpgl0415.fastq
    ## Sometimes I put in -i hpgl0415.fastq when I mean hpgl0415-trimmed.fastq(.gz|.xz)
    ## Sometimes I put in -i hpgl0415_forward.fastq:hpgl0415_reverse.fastq

    ## So, this function should make this unambiguous and consistent no matter what I type
    ## First: lower-case whatever I typed because uppercase characters are obnoxious
    my @input_list = undef;
    if ($input =~ /:/) {
        @input_list = split(/:/, $input);
    } else {
        push(@input_list, $input);
    }

    foreach my $in (@input_list) {
        my $low = $in;
        $low =~ tr/[A-Z]/[a-z]/;
        unless ($in eq $low) {
            system(qq"mv ${in}.fastq ${low}.fastq");
            if (-r "${in}_forward.fastq") {
                system(qq"mv ${in}_forward.fastq ${low}_forward.fastq");
            }
            if (-r "${in}_reverse.fastq") {
                system(qq"mv ${in}_reverse.fastq ${low}_reverse.fastq");
            }
        }
##        if (-r "${in}" and $in =~ /\.gz$/) {
##            print "gunzip of ${in}\n";
##            system("nice gunzip ${in}");
##        }
##        if (-r "${in}" and $in =~ /\.xz$/) {
##            print "xz -d ${in}\n";
##           system("nice xz -d ${in}");
##        }
    }
    $input =~ tr/[A-Z]/[a-z]/;
##    $input =~ s/\.gz//g;
##    $input =~ s/\.xz//g;

    return($input);
}

sub Last_Stat {
    my $me = shift;
    my %args = @_;
    my $input = new FileHandle;
    my $input_filename = $args{input};
    $input->open("<$input_filename");
    my ($line, $last);
    while ($line = <$input>) {
        chomp $line;
        ## I am doing this due to terminal newlines on the file, there is probably a much smrtr way.
        $last = $line unless ($line =~ /^$/);
    }
    $input->close();
    return($last);
}

1;

