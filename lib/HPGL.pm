package HPGL;
use AppConfig qw":argcount :expand";
use Archive::Extract;
use autodie qw":all";
use common::sense;
use Bio::DB::Universal;
use Bio::Root::RootI;
use Cwd qw"cwd";
use Digest::MD5 qw"md5 md5_hex md5_base64";
use File::Basename;
use File::Find;
use Getopt::Long;
use HPGL::SeqMisc;
use HPGL::Ceph;
use HPGL::RNASeq;
use HPGL::Aligners;
use Log::Log4perl;
use Log::Log4perl::Level;
use Net::Amazon::S3;
use PerlIO;
use Pod::Usage;
use Term::ReadLine;

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

=head1 SYNOPSIS

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

=cut
sub new {
    my ($class, %args) = @_;
    my $me = bless {}, $class;
    foreach my $key (keys %args) {
        $me->{$key} = $args{$key} if ($args{$key});
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
    $me->{shell} = "/usr/bin/bash";
    $me->{align_trimmed} = 1;
    $me->{aligner} = "bowtie";
    $me->{aws_access_key_id} = 'PQOD56HDPQOXJPYQYHH9';
    $me->{aws_hostname} = 'gembox.cbcb.umd.edu';
    $me->{aws_secret_access_key} = 'kNl27xjeVSnChzEww9ziq1VkgUMNmNonNYjxWkGw';
    $me->{bam_small} = "21-26";
    $me->{bam_mid} = "27-32";
    $me->{bam_big} = "33-40";
    $me->{base} = undef;
    $me->{basedir} = cwd();
    $me->{bt_args} = {
        v0M1 => "--best -v 0 -M 1",
        v1M1 => "--best -v 1 -M 1",
        v2M1 => "--best -v 2 -M 1",
    };
    $me->{btmulti} = 0;
    $me->{config} = qq"$ENV{HOME}/.config/hpgl.conf",
    $me->{debug} = 1;
    $me->{depends} = {};
    $me->{help} = undef;
    $me->{hpgl} = undef;
    $me->{hpglid} = undef;
    $me->{htseq_identifier} = "ID";
    $me->{htseq_stranded} = "no";
    $me->{input} = undef;
    $me->{jobs} = [];
    $me->{jobids} = {};
    $me->{len_min} = 21;
    $me->{len_max} = 40;
    $me->{libdir} = "$ENV{HOME}/libraries";
    $me->{libtype} = "genome";
    $me->{paired} = 0;
    $me->{pbs} = 1;
    $me->{qsub_args} = "-j oe -V -m n";
    $me->{qsub_queue} = "throughput";
    $me->{qsub_queues} = ["throughput","workstation","long","large"];
    $me->{qsub_shell} = "$bash";
    $me->{qsub_mem} = 6;
    $me->{qsub_wall} = "10:00:00";
    $me->{qsub_cpus} = "4";
    $me->{qsub_depends} = "depend=afterok:";
    $me->{qsub_loghost} = "ibissub00.umiacs.umd.edu";
    $me->{qsub_logdir} = qq"$options{qsub_loghost}:$options{basedir}/outputs";
    $me->{species} = undef;
    $me->{suffixes} = (".fastq",".gz",".xz", ".fasta", ".sam", ".bam", ".count");
    $me->{trimmer} = "trimmomatic";
    $me->{type} = "rnaseq";

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
        "debug|d" => \$options{debug},
        "btmulti:i" => \$options{btmulti},
        "help" => \$options{help},
        "hpgl|h:s" => \$options{hpgl},
        "input|i:s" => \$options{input},
        "pbs|p:i" => \$options{pbs},
        "species|s:s" => \$options{species},
        "stranded:s" => \$options{htseq_stranded},
        "identifier:s" => \$options{htseq_identifier},
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
    my $base = $me->{input};
    my @suffixes = @{$me->{suffixes}};
    $base = basename($base, @suffixes);
    $base = basename($base, @suffixes);
    $me->{base} = $base;
    if (!defined($me->{hpglid})) {
        my $tmp = $me->{input};
        $tmp =~ s/^(hpgl\d+).*/$1/g;
        $me->{hpglid} = $tmp;
    }
    return($me);
}

=head2
    Help()
=cut
sub Help {
    pod2usage();
    exit(0);
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

sub AUTOLOAD {
    my $me = shift;
    my $type = ref($me) or Callstack(die => 1, message => "$me is not an object");
    my $name = $AUTOLOAD;
    $name =~ s/.*://;   # strip fully-qualified portion
    if (@_) {
        return $me->{$name} = shift;
    }
    else {
        return $me->{$name};
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

1;

