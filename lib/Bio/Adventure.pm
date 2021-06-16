package Bio::Adventure;
## LICENSE: gplv2
## ABSTRACT:  Kitty!
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
use vars qw"$VERSION";

use AppConfig qw":argcount :expand";
use Archive::Extract;
##use Bio::DB::Sam;
use Bio::SeqIO;
use Bio::DB::Universal;
use Bio::Root::RootI;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Carp qw"croak carp confess cluck longmess";
use Cwd qw"cwd";
use Digest::MD5 qw"md5 md5_hex md5_base64";
use Env qw"COMPRESSION XZ_OPTS XZ_DEFAULTS HOME";
use File::Basename;
use File::Find;
use File::Spec;
use File::Which qw"which";
use File::Path qw"make_path remove_tree";
use File::Temp;
use FileHandle;
use Getopt::Long qw"GetOptionsFromArray";
use IO::String;
use JSON -convert_blessed_universally;
use Log::Log4perl;
##use Log::Log4perl::Level;
use PerlIO;
use Pod::Usage;
use Storable qw"store retrieve nstore";
use Term::ReadLine;
use Term::UI;

use Bio::Adventure::Align;
use Bio::Adventure::Align_Blast;
use Bio::Adventure::Align_Fasta;
use Bio::Adventure::Annotation;
use Bio::Adventure::Assembly;
use Bio::Adventure::Cleanup;
use Bio::Adventure::Count;
use Bio::Adventure::Compress;
use Bio::Adventure::Convert;
use Bio::Adventure::Local;
use Bio::Adventure::Map;
use Bio::Adventure::Prepare;
use Bio::Adventure::Phylogeny;
use Bio::Adventure::Riboseq;
use Bio::Adventure::QA;
use Bio::Adventure::SeqMisc;
use Bio::Adventure::Slurm;
use Bio::Adventure::SNP;
use Bio::Adventure::TNSeq;
use Bio::Adventure::Torque;
use Bio::Adventure::Trim;

has task => (is => 'rw', default => 'rnaseq');
has method => (is => 'rw');
has thief => (is => 'rw', default => 'Bilbo Baggins');
has options => (is => 'rw');
has option_file => (is => 'rw');
has sbatch_path => (is => 'rw', default => My_Which('sbatch'));
has qsub_path => (is => 'rw', default => My_Which('qsub'));
has bash_path => (is => 'rw', default => My_Which('bash'));
##has menus => (is => 'ro', default => Get_Menus());

our $AUTOLOAD;
##our @EXPORT_OK = qw"";
$VERSION = '20151101';
$COMPRESSION = 'xz';
$XZ_OPTS = '-9e';
$XZ_DEFAULTS = '-9e';
$ENV{LESSOPEN} = '| lesspipe %s';

=head1 NAME

    Bio::Adventure - A perl library to make preprocessing high-throughput data easier.

=head1 SYNOPSIS

    Bio::Adventure.pm  Methods to simplify creating job scripts and submitting them to a
    computing cluster.  Bio::Adventure makes heavy use of command line options and as a result
    accepts many options.  The simplest way to access these methods is via the cyoa script, but
    everything may be called directly in perl.

    use Bio::Adventure;
    my $hpgl = new Bio::Adventure;
    $hpgl->{species} = 'mmusculus';  ## There had better be a mmusculus.[fasta&gff] in $hpgl->{libdir}
    my $hpgl = new Bio::Adventure(species => 'scerevisiae', libdir => '/home/bob/libraries');
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

    This library should write out PBS/slurm compatible job files and
    submit them to the appropriate computing cluster.  It should also
    collect the outputs and clean up the mess.

=head2 Methods

=over 4

=item C<Help>

    Help() always gives 0.
    Before it returns, it will hopefully print some useful information
    regarding ways to invoke Bio::Adventure.pm.

=cut
sub Help {
    my $class = shift @_;
    my $fh = \*STDOUT;
    ##    my $usage = pod2usage(-output => $fh, -verbose => 99, -sections => "SYNOPSIS");
    use Pod::Find qw(pod_where);
    pod2usage(-input => pod_where({-inc => 1}, __PACKAGE__),
              -verbose => 2, -output => $fh,
              -sections => "NAME|SYNOPSIS|DESCRIPTION|VERSION",
              -exitval => 'NOEXIT');
    return(0);
}

sub BUILDARGS {
    my ($class, %args) = @_;
    my $attribs = {};
    ## First populate the class with a series of hopefully useful default values.
    my $defaults = Get_Defaults();
    foreach my $k (keys %{$defaults}) {
        $attribs->{$k} = $defaults->{$k};
    }
    ## Add a little check for slurm/torque
    my $queue_test = My_Which('sbatch');
    if ($queue_test) {
        $attribs->{pbs} = 'slurm';
    } else {
        $queue_test = My_Which('qsub');
        if ($queue_test) {
            $attribs->{pbs} = 'torque';
        }
    }

    if ($#ARGV > 0) {
        ## Make a log of command line arguments passed.
        my $arg_string = '';
        foreach my $a (@ARGV) {
            $arg_string .= "$a ";
        }
        make_path("outputs/", {verbose => 0}) unless (-r qq"outputs/");
        my $out = FileHandle->new(">>outputs/log.txt");
        my $d = qx'date';
        chomp $d;
        print $out "# Started CYOA at ${d} with arguments: $arg_string.\n";
        $out->close();
    }

    ## If one desires, use a configuration file to replace/augment those values.
    my $appconfig = AppConfig->new({CASE => 1,
                                    CREATE => 1,
                                    PEDANTIC => 0,
                                    ## ERROR => eval(),
                                    GLOBAL => {EXPAND => EXPAND_ALL,
                                               EXPAND_ENV => 1,
                                               EXPAND_UID => 1,
                                               DEFAULT => 'unset',
                                               ARGCOUNT => 1,
                                    },});
    if (-r $defaults->{config_file}) {
        my $open = $appconfig->file($defaults->{config_file});
        my %data = $appconfig->varlist("^.*");
        for my $config_option (keys %data) {
            $attribs->{$config_option} = $data{$config_option};
            undef $data{$config_option};
        }
    }

    ## Now pull an arbitrary set of command line arguments
    ## Options set in the config file are trumped by those placed on the comamnd line.
    my (%conf, %conf_specification, %conf_specification_temp);
    foreach my $default (keys %{$defaults}) {
        ## This line tells Getopt::Long to use as a string argument anything in the set of defaults.
        my $def_string = qq"${default}:s";
        $conf_specification{$def_string} = \$conf{$default};
    }
    ## For some command line options, one might want shortcuts etc, set those here:
    %conf_specification_temp = ("de|d" => \$conf{debug},
                                "hp|h:s" => \$conf{hpgl},
                                "in|i:s" => \$conf{input},
                                "pb|p:i" => \$conf{pbs},
                                "sp|s:s" => \$conf{species},
        );
    ## This makes both of the above groups command-line changeable
    foreach my $name (keys %conf_specification_temp) {
        $conf_specification{$name} = $conf_specification_temp{$name};
    }
    undef(%conf_specification_temp);

    ## Use Getopt::Long to acquire all the command line options.
    ## Now pull out the command line options.
    my $argv_result = GetOptions(%conf_specification);
    ##unless ($argv_result) {
    ##    Help();
    ##}

    ## Finally merge the options from the config file to those
    ## in the default hash and those from the command line.
    foreach my $opt (keys %conf) {
        if (defined($conf{$opt})) {
            $attribs->{$opt} = $conf{$opt};
        }
    }
    undef(%conf);

    ## Take a moment to create a simplified job basename
    ## Eg. if the input file is 'hpgl0523_forward-trimmed.fastq.gz'
    ## Take just hpgl0523 as the job basename
    my $job_basename;
    if ($attribs->{input}) {
        $job_basename = $attribs->{input};
        ## Start by pulling apart any colon/comma separated inputs
        if ($job_basename =~ /:|\,/) {
            my @tmp = split(/:|\,/, $job_basename);
            ## Remove likely extraneous information
            $job_basename = $tmp[0];
        }
        $job_basename = basename($job_basename, @{$attribs->{suffixes}});
        $job_basename = basename($job_basename, @{$attribs->{suffixes}});
    } else {
        $job_basename = basename($attribs->{basedir}, @{$attribs->{suffixes}});
    }
    $job_basename =~ s/_forward.*//g;
    $job_basename =~ s/_R1.*//g;
    $job_basename =~ s/-trimmed.*//g;
    $attribs->{job_basename} = $job_basename;

    ## Remove backslashes from arbitrary arguments
    if ($attribs->{arbitrary}) {
        $attribs->{arbitrary} =~ s/\\//g;
    }

    ## These are both problematic when (n)storing the data.
    $attribs->{menus} = Get_Menus();
    ## Extract a function reference for running a job here.
    my $methods_to_run = Get_TODOs(task => $attribs->{task}, method => $attribs->{method});
    $attribs->{methods_to_run} = $methods_to_run;
    $args{options} = $attribs;
    return \%args;
}

sub Get_Term {
    my $term = Term::ReadLine->new('>');
    my $attribs = $term->Attribs;
    $attribs->{completion_suppress_append} = 1;
    my $OUT = $term->OUT || \*STDOUT;
    $Term::UI::VERBOSE = 0;
    $term->ornaments(0);
    return($term);
}

sub Check_Input {
    my ($class, %args) = @_;
    my $file_list;
    if (ref($args{files}) eq 'SCALAR' || ref($args{files}) eq '') {
        if ($args{files} =~ /:/) {
            my @tmp = split(/:/, $args{files});
            $file_list = \@tmp;
        } else {
            $file_list->[0] = $args{files};
        }
    } elsif (ref($args{files}) eq 'ARRAY') {
        $file_list = $args{files};
    } else {
        my $unknown_class = ref($args{files});
        warn("I do not know type: ${unknown_class}.");
    }
    my $found = {};
    foreach my $file (@{$file_list}) {
        $found->{$file} = 0;
        my $first_test = $file;
        my $second_test = basename($file, $class->{options}->{suffixes});
        my $third_test = basename($second_test, $class->{options}->{suffixes});
        $found->{$file} = $found->{$file} + 1 if (-r $first_test);
        $found->{$file} = $found->{$file} + 1 if (-r $second_test);
        $found->{$file} = $found->{$file} + 1 if (-r $third_test);
    }
    for my $f (keys %{$found}) {
        if ($found->{$f} == 0) {
            die("Unable to find a file corresponding to $f.");
        }
    }
    return($found);
}

sub Get_Basename {
    my ($class, $string) = @_;
    my ($in1, $in2) = "";
    if ($string =~ /:/) {
        ($in1, $in2) = split(/:/, $string);
    } else {
        $in1 = $string;
    }
    $in1 = basename($in1, $class->{options}->{suffixes});
    $in1 = basename($in1, $class->{options}->{suffixes});
    return($in1);
}

=item C<Get_Defaults>

    This is the primary source for the configuration data in Bio::Adventure.
    It every option in it is accessible as a command line argument.  Those command line
    arguments are described below:

=over 4

=item C<General Options>

  basedir - Base directory from which all other directories inherit. (cwd)
  config_file - Optional configuration file to hold chosen default options. (~/.config/cyoa.conf)
  input - Generic input file, used by most functions.
  output - Generic output file, used by most functions.
  genome - Location of the chosen genome file, used by many tasks.
  species - Species name for mapping and other tasks.
  gff - Location of the chosen gff file, used by many tasks.
  libdir - Directory containing the various indexes for mapping. (~/libraries)
  libtype - We may have ribosomal or genomic indexes, or indeed other types if required. (genome)
  task - Top-level task to perform (trimming, conversion, mapping, etc).
  method - Which function to actually call?
  raw_dir - Location of the raw data when it might lie elsewhere.
  type - Used in a few places, genomic vs ribosomal, tnseq vs rnaseq, etc.
  suffixes - Set of suffixes to remove from files when looking for basenames. (.fastq, .gz, .fasta, .sam, .bam, .count, .csfast, .qual)
  debug - Print debugging information? (false)
  help - Print help information.

=item C<Cluster Specific Options>

  depends - A string describing the set of dependencies for a job.
  jname - A machine-readable name for each job.
  jprefix - A prefix number for each job, to make reading successions easier.
  jstring - The string comprising the actual job to be run.
  job_id - The job ID returned by Torque/Slurm for a job.
  queue_array_string - When running an array of jobs, this holds the requisite formatted string.
  language - Choose if the cluster options should be appropriate for a shell or perl script. (bash)
  pbs - Use a computing cluster or just invoke scripts one at a time?
  qsub_args - Default arguments when using torque. (-j oe -V -m n)
  sbatch_depends - Dependency string on a slurm cluster. (afterok:)
  sbatch_dependsarray - Dependency string on a slurm for array jobs. (afterok:)
  qsub_depends - Dependency string on a torque cluster. (depend=afterok:)
  qsub_dependsarray - Dependency string on a torque cluster for array jobs. (depend=afterokarray:)
  loghost - Host to send logs for torque. (localhost)
  mem - Memory in Gb to request from the cluster. (6)
  queue - Default qos/queue to send jobs. (workstation)
  queues - Set of available queues. (throughput, workstation, long, large)
  shell - Default shell to request when using torque/slurm. (/usr/bin/bash)
  wall - Amount of walltime to request from the cluster. (10:00:00)

=item C<Alignment Arguments>

  best_only - When parsing blast/fasta output, return only the best hit? (false)
  align_jobs - Split alignments into how many jobs? (40)
  blast_format - Which blast output format should we use? 5 is blastxml (5)
  fasta_format - Which fastq output format should we use? 0 is tabular. (0)
  align_parse - Control whether to parse an alignment's output. (true)
  blast_params - Default options for blast searches. (-e 10)
  peptide - Explicitly tell blast if this is (or not) a peptide search. (F)
  blast_tool - Choose which member of the blast family to use.
  library - Choose the name of the blast/fasta database to search against.
  evalue - E-value cutoff when performing blast/fasta searches.  This is not currently entirely used. (1)
  identity - Percent identical cutoff when doing blast/fasta.  Not used as it should be. (70)
  fasta_args - Default arguments when invoking fasta36 programs. (-b 20 -d 20)
  fasta_tool - Which fasta tool to invoke? (ggsearch36)

=item C<Mapping Arguments>

  mapper - Set the default mapping program. (salmon)
  bt_args - Some default arguments when using bowtie. (--very-sensitive -L 14)
  bt2_args - Similar for bowtie2.
  bt_type - Choose a specific type of bowtie mapping to perform. (v0M1)
  btmulti - Or perform bowtie with all of the chosen options in bt_args? (false)
  stranded - Choose options for stranded libraries (primarily kallisto) (false)

=item C<Counting Arguments>

  htseq_args - Some default arguments for htseq. (--order=name --idattr=gene_id --minaqual=10 --type exon --stranded=yes --mode=union)
  htseq_stranded - Use htseq in stranded mode? (no)
  htseq_type - Which gff type (3rd column) to use in htseq? (exon)
  htseq_id - Which gff ID (tag in the last column) to use in htseq? (gene_id)
  mi_genome - Database of miRNA sequences for counting small RNA sequences.
  mature_fasta - Database of miRNA mature sequences when counting small RNA sequences.
  mirbase_data - Location of the full miRbase data when counting small RNA.

=item C<Conversion Arguments>

  bamfile - Where to put output bam data?  Used in a few places.
  feature_type - Primarily used when looking at gff files, ergo htseq and converting gff to other formats.
  taxid - Used for making gaf files.

=item C<Preparation Arguments>

  csv_file - A csv file containing the metadata for multiple samples. (all_samples.csv)
  sampleid - Chosen ID for a given sample.

=item C<TNSeq Arguments>

  index_file - File containing a key/value pairing of sample IDs and TNSeq indexes. (indexes.txt)
  runs - How many gumbel runs to apply in essentiality? (1000)

=item C<Riboseq Arguments>

  orientation - When looking at reads vs. starts/ends, count start->end or end->start? (start)
  riboasite,psite,esite - Count A/P/E site reads. (1,1,1)
  ribominsize/maxsize - Range for riboseq sizes. (24,36)
  ribominpos/maxpos - Relative positions when counting vs. start/stop codons. (-30, 30)
  riboanchor - Anchor reads with respect to the start/end. (start)
  ribosizes - Set of sizes to count. (25,26,27,28,29,30,31,32,33,34)

=item C<Trimming Arguments>

  maxlength - Maximum length to keep when trimming sequence data (usually TNSeq). (42)
  minlength - Minimum length to keep when trimming sequence data (usually TNSeq). (8)
  phred - Score range to work in, also used for some QC tasks. (33)
  qual - cutadapt specific quality string.

=item C<QC Arguments>

  paired - Primarily Fastqc, boolean check if the data is paired. (false)

=item C<Phylogeny Arguments>

  outgroup - Name of a chosen phylogenetic outgroup.
  starting_tree - Name of a starting tree for improving phylogenetic trees.

=item C<SNP Arguments>

  trim - Trim raw sequence with options more appropriate for SNP calling.
  varfilter - Filter variants by variance. (true)
  vcf_cutoff - Minimum reads for a variant call. (10)
  vcf_minpct - Minimum percent agreement for a variant call. (0.8)

=back

=cut
sub Get_Defaults {
    my $defaults = {
        ## General options
        arbitrary => undef, ## Pass arbitrary options
        basedir => cwd(), ## The base directory when invoking various shell commands
        config_file => qq"${HOME}/.config/cyoa.conf", ## A config file to read to replace these values.
        directories => undef, ## Apply a command to multiple input directories.
        input => undef, ## Generic input argument
        output => undef, ## Generic output argument
        genome => undef, ## Which genome to use?
        species => undef, ## Chosen species
        gff => undef, ## A default gff file!
        libdir => "${HOME}/libraries", ## Directory of libraries for mapping rnaseq reads
        libtype => 'genome', ## Type of rnaseq mapping to perform (genomic/rrna/contaminants)
        task => undef,
        method => undef,
        raw_dir => undef, ## Directory containing raw reads
        type => undef, ## Type!
        suffixes => ['.fastq', '.gz', '.xz', '.fasta', '.sam', '.bam', '.count', '.csfasta', '.qual'], ## Suffixes to remove when invoking basename
        debug => 0, ## Debugging?
        help => undef, ## Ask for help?

        ## Cluster options
        depends => "",      ## A flag for PBS telling what each job depends upon
        jname => "",        ## A name for a job!
        jprefix => "",      ## An optional prefix number for the job
        jstring => "",      ## The string of text printed as a job
        jobids => {},       ## And a hash of jobids
        jobs => [],         ## A list of jobs
        language => 'bash', ## Shell or perl script?
        pbs => undef,       ## Use pbs?
        qsub_args => '-j oe -V -m n', ## What arguments will be passed to qsub by default?
        cpus => '4',                  ## Number of to request in jobs
        qsub_depends => 'depend=afterok:', ## String to pass for dependencies
        qsub_dependsarray => 'depend=afterokarray:', ## String to pass for an array of jobs
        sbatch_depends => 'afterok:',
        sbatch_dependsarray => 'afterok:', ## String to pass for an array of jobs
        loghost => 'localhost',            ## Host to which to send logs
        mem => 24,                          ## Number of gigs of ram to request
        queue => 'workstation',            ## What queue will jobs default to?
        queues => ['throughput','workstation','long','large'], ## Other possible queues
        shell => '/usr/bin/bash', ## Default qsub shell
        wall => '10:00:00',       ## Default time to request

        ## Alignment options
        best_only => 0,
        align_jobs => 40, ## How many blast/fasta alignment jobs should we make?
        align_blast_format => 5, ## Which alignment type should we use?
        ## (5 is blastxml for blast, 0 is tabular for fasta36)
        align_parse => 1,
        blast_params => ' -e 10 ',
        peptide => 'F', ## The flag for blast deciding whether or not to perform a peptide formatdb/alignment
        blast_tool => undef, ## Valid choices are blastn blastp tblastn tblastx and I'm sure some others I can't remember
        library => undef,   ## The library to be used for fasta36/blast searches
        evalue => 0.001,        ## Alignment specific: Filter hits by e-value.
        identity => 70, ## Alignment specific: Filter hits by sequence identity percent.
        fasta_args => ' -b 20 -d 20 ', ## Default arguments for the fasta36 suite
        fasta_tool => 'ggsearch36',    ## Which fasta36 program to run

        ## Mapping options
        mapper => 'salmon',     ## Use this aligner if none is chosen
        bt_args => { def => '', ## A series of default bowtie1 arguments
                     v0M1 => '--best -v 0 -M 1',
                     v1M1 => '--best -v 1 -M 1',
                     v1M1l10 => '--best -v 1 -M 1 -y -l 15',
                     v2M1 => '--best -v 2 -M 1',
        },
        bt2_args => ' --very-sensitive -L 14 ', ## A boolean to decide whether to use multiple bowtie1 argument sets
        bt_type => 'v0M1',
        btmulti => 0, ## Use the following configuration file to overwrite options for these scripts
        stranded => 0,

        ## Counting options
        htseq_args => {         ##hsapiens => " ",
            ## mmusculus => " -i ID ",
            ## lmajor => " -i ID ",
            default => " --order=name --idattr=gene_id --minaqual=10 --type=exon --stranded=yes --mode=union ",
            ## all => " -i ID ",
            ## all => " ", ## Options chosen by specific species, this should be removed.
        },
        htseq_stranded => 'no', ## Use htseq stranded options?
        htseq_type => 'exon', ## The identifier flag passed to htseq (probably should be moved to feature_type)
        htseq_id => 'gene_id',
        mapper => 'hisat2', ## What was used to map the input
        mi_genome => undef, ## miRbase genome to search
        mature_fasta => undef, ## Database of mature miRNA sequences
        mirbase_data => undef, ## Database of miRNA annotations from mirbase

        ## Conversion options
        bamfile => undef, ## Used for mimapping for now, but should also be used for tnseq/riboseq/etc
        taxid => '353153', ## Default taxonomy ID
        tag => 'gene_id',

        ## Preparation options
        csv_file => 'all_samples.csv',
        sampleid => undef, ## An hpgl identifier

        ## TNSeq options
        index_file => 'indexes.txt', ## The default input file when demultiplexing tnseq data
        ta_offset => 0, ## When counting TAs, this is either 0 or 2 depending on if the TA was removed.
        runs => 1000,

        ## Ribosome Profiling options
        orientation => 'start', ## What is the orientation of the read with respect to start/stop codon (riboseq) FIXME: Rename this
        riboasite => 1,         ## Count the A site ribosomes (riboseq)?
        ribopsite => 1,         ## Count the P site ribosomes (riboseq)?
        riboesite => 1,         ## Count the P site ribosomes (riboseq)?
        ribominsize => 24, ## Minimum size to search for ribosome positions (riboseq)
        ribomaxsize => 36, ## Maximum size to search for ribosome positions (riboseq)
        ribominpos => -30, ## Minimum position for counting (riboseq)
        ribomaxpos => 30,  ## Maximum position for counting (riboseq)
        ribocorrect => 1,  ## Correct ribosome positions for biases
        riboanchor => 'start', ## When correcting, use the start or end position as an anchor
        ribosizes => '25,26,27,28,29,30,31,32,33,34', ## Use these sizes for riboseq reads

        ## Trimming options
        maxlength => 42,
        minlength => 8,
        phred => 33,
        qual => undef,          ## Quality string for cutadapt

        ## QC Options
        paired => 0,            ## Paired reads?

        ## Phylogen Options
        outgroup => undef,
        starting_tree => undef,

        ## SNP
        varfilter => 1,       ## Do a varFilter in variant searches
        vcf_cutoff => 10,     ## Minimum depth cutoff for variant searches
        vcf_minpct => 0.8,    ## Minimum percent agreement for variant searches.
        gff_tag => 'gene_id', ## Which ID tag to use when snp searching.
        gff_type => 'gene',   ## Which gff type to use when snp searching.

        ## State
        state => {
            num_sequences => 0,
        },
    };

    $defaults->{logdir} = qq"$defaults->{basedir}/outputs/logs";
    $defaults->{blast_format} = qq"formatdb -p $defaults->{peptide} -o T -s -i -n species";
    return($defaults);
}

sub Get_Menus {
    my $menus = {
        Alignment => {
            name => 'alignment',
            message => 'Hari Seldon once said violence is the last refuge of the incompetent.  Go to page 6626070.',
            choices => {
                    '(blastsplit): Split the input sequence into subsets and align with blast.' => \&Bio::Adventure::Align_Blast::Split_Align_Blast,
                    '(fastasplit): Split the input sequence into subsets and align with fasta36.' => \&Bio::Adventure::Align_Fasta::Split_Align_Fasta,
                    '(concat): Merge split searches into a single set of results.' => \&Bio::Adventure::Align::Concatenate_Searches,
                    '(fastaparse): Parse fasta36 output into a reasonably simple table of hits.' => \&Bio::Adventure::Align_Fasta::Parse_Fasta,
                    '(blastparse): Parse blast output into a reasonably simple table of hits.' => \&Bio::Adventure::Align_Blast::Parse_Blast,
                    '(fastamerge): Merge and Parse fasta36 output into a reasonably simple table of hits.' => \&Bio::Adventure::Align_Fasta::Merge_Parse_Fasta,
                    '(blastmerge): Merge and Parse blast output into a reasonably simple table of hits.' => \&Bio::Adventure::Align_Blast::Merge_Parse_Blast,
            },
        },
        Annotation => {
            name => 'annotation',
            message => 'How come Aquaman can control whales?  They are mammals!  Makes no sense.',
            choices =>  {
                '(aragorn): Search for tRNAs with aragorn.' => \&Bio::Adventure::Annotation::Aragorn,
                '(glimmer): Use glimmer to search for ORFs.' => \&Bio::Adventure::Annotation::Glimmer,
                '(resfinder): Search for antimicrobial resistance genes.' => \&Bio::Adventure::Annotation::Resfinder,
                '(trnascan): Search for tRNA genes with trnascan.' => \&Bio::Adventure::Annotation::tRNAScan,    
            },
        },
        Assembly => {
            name => 'assembly',
            message => 'The wise man fears the wrath of a gentle heart. Go to page 314159.',
            choices => {
                '(extract_trinotate): Extract the most likely hits from Trinotate.' => \&Bio::Adventure::Assembly::Extract_Trinotate,
                '(extend_kraken): Extend a kraken2 database with some new sequences.' => \&Bio::Adventure::Assembly::Extend_Kraken_DB,
                '(kraken2): Taxonomically classify reads.' => \&Bio::Adventure::Assembly::Kraken,
                '(transdecoder):  Run transdecoder on a putative transcriptome.' => \&Bio::Adventure::Assembly::Transdecoder,
                '(trinotate): Perform de novo transcriptome annotation with trinotate.' => \&Bio::Adventure::Assembly::Trinotate,
                '(trinity): Perform de novo transcriptome assembly with trinity.' => \&Bio::Adventure::Assembly::Trinity,
                '(trinitypost): Perform post assembly analyses with trinity.' => \&Bio::Adventure::Assembly::Trinity_Post,
                '(velvet): Perform de novo genome assembly with velvet.' => \&Bio::Adventure::Assembly::Velvet,
                '(shovill): Perform the shovill pre/post processing with spades.' => \&Bio::Adventure::Assembly::Shovill,
            },
        },
        Conversion => {
            name => 'convert',
            message => 'And it rained a fever. And it rained a silence. And it rained a sacrifice. And it rained a miracle. And it rained sorceries and saturnine eyes of the totem.  Go to page 2584981.',
            choices => {
                '(sam2bam): Convert a sam mapping to compressed/sorted/indexed bam.' => \&Bio::Adventure::Convert::Sam2Bam,
                '(gb2gff): Convert a genbank flat file to gff/fasta files.' => \&Bio::Adventure::Convert::Gb2Gff,
                '(gff2fasta): Convert a gff file to a fasta file.' => \&Bio::Adventure::Convert::Gff2Fasta,
            },
        },
        Counting => {
            name => 'count',
            message => 'Once men turned their thinking over to machines in the hope that this would set them free. But that only permitted other men with machines to enslave them.  Go to page 27812',
            choices => {
                '(htseq): Count mappings with htseq-count.' =>  \&Bio::Adventure::Count::HTSeq,
                '(htmulti): Use different option sets for counting with htseq.' => \&Bio::Adventure::Count::HT_Multi,
                '(mimap): Count mature miRNA species.' => \&Bio::Adventure::Count::Mi_Map,
                '(countstates): Count ribosome positions.' => \&Bio::Adventure::Riboseq::Count_States,
            },
        },
        Mapping => {
            name => 'map',
            message => 'The world is dark and full of terrors, take this and go to page 6022140.',
            choices => {
                '(bowtie): Map trimmed reads with bowtie1 and count with htseq.' => \&Bio::Adventure::Map::Bowtie,
                '(bt2): Map trimmed reads with bowtie2 and count with htseq.' => \&Bio::Adventure::Map::Bowtie2,
                '(hisat): Map trimmed reads with hisat2 and count with htseq.' => \&Bio::Adventure::Map::Hisat2,
                '(btmulti): Map trimmed reads and count using multiple bowtie1 option sets.' => \&Bio::Adventure::Map::BT_Multi,
                '(bwa): Map reads with bwa and count with htseq.' => \&Bio::Adventure::Map::BWA,
                '(kallisto): Pseudo-align and count reads using kallisto.' => \&Bio::Adventure::Map::Kallisto,
                '(salmon): Pseudo-align and count reads using salmon.' => \&Bio::Adventure::Map::Salmon,
                '(star): Pseudo-align and count reads using STAR.' => \&Bio::Adventure::Map::STAR,
                '(mimap): Attempt to map reads explicitly to mature miRNA species.' => \&Bio::Adventure::Count::Mi_Map,
                '(rrnabowtie): Map rRNA reads using bowtie1.' => \&Bio::Adventure::Map::Bowtie_RRNA,
                '(rsem): Quantify reads using rsem.' => \&Bio::Adventure::Map::RSEM,
                '(tophat): Map reads using tophat2 and count with htseq.' => \&Bio::Adventure::Map::Tophat,
                '(indexbt1): Create bowtie1 compatible indexes.' => \&Bio::Adventure::Map::BT1_Index,
                '(indexbt2): Create bowtie2 compatible indexes.' => \&Bio::Adventure::Map::BT2_Index,
                '(indexhisat): Create hisat2 compatible indexes.' => \&Bio::Adventure::Map::HT2_Index,
                '(indexbwa): Create bwa compatible indexes.' => \&Bio::Adventure::Map::BWA_Index,
                '(indexkallisto): Create kallisto compatible indexes.' => \&Bio::Adventure::Map::Kallisto_Index,
                '(indexrsem): Create rsem indexes.' => \&Bio::Adventure::Map::RSEM_Index,
                '(indexsalmon): Create salmon indexes.' => \&Bio::Adventure::Map::Salmon_Index,
            },
        },
        Phylogeny => {
            name => 'phylogeny',
            message => '',
            choices => {
                '(gubbins): Run Gubbins with an input msa.' => \&Bio::Adventure::Phylogeny::Run_Gubbins,
            },
        },
        Pipeline => {
            name => 'pipeline',
            message => 'When Mr. Bilbo Baggins announced he would shortly be celebrating his eleventyfirst birthday, there was much talk and excitement in Hobbiton.  Go to page 1618033',
            choices => {
                '(priboseq): Perform a preset pipeline of ribosome profiling tasks.' => \&Bio::Adventure::Pipeline_Riboseq,
                '(ptnseq): Perform a preset pipeline of TNSeq tasks.' => \&Bio::Adventure::Pipeline_TNSeq,
                '(pbt1): Perform a preset group of bowtie1 tasks.' => \&Bio::Adventure::Pipeline_Bowtie,
                '(pbt2): Use preset bowtie2 tasks.' => \&Bio::Adventure::Pipeline_Bowtie2,
                '(ptophat): Use preset tophat tasks.' => \&Bio::Adventure::Pipeline_Tophat,
                '(pbwa): Try preset bwa tasks.' => \&Bio::Adventure::Pipeline_BWA,
                '(pkallisto): Try preset kallisto tasks.' => \&Bio::Adventure::Pipeline_Kallisto,
            },
        },
        Prepare => {
            name => 'preparation',
            message => 'Whan that Aprille withe her shoures sote, the droughte of Marche hath perced to the rote.  go to Cantebury.',
            choices => {
                '(read_samples): Read samples using a csv file to determine the raw data locations.' => \&Bio::Adventure::Prepare::Read_Samples,
                '(copyraw): Copy data from the raw data archive to scratch.' => \&Bio::Adventure::Prepare::Copy_Raw,
                '(fastqdump): Download data from sra.' => \&Bio::Adventure::Prepare::Fastq_Dump,
            },
        },
        QA => {
            name => 'qa',
            message => 'There is a time when the operation of the machine becomes so odious that you cannot take part.',
            choices => {
                '(biopieces): Use biopieces to graph some metrics of the data.' => \&Bio::Adventure::QA::Biopieces_Graph,
                '(cutadapt): Perform adapter trimming with cutadapt.' => \&Bio::Adventure::Trim::Cutadapt,
                '(fastqc): Use fastqc to check the overall quality of the raw data.' => \&Bio::Adventure::QA::Fastqc,
                '(racer): Perform sequence correction with hitec/RACER.' => \&Bio::Adventure::Trim::Racer,
                '(trimomatic): Perform adapter trimming with Trimomatic.' => \&Bio::Adventure::Trim::Trimomatic,
            },
        },
        RiboSeq => {
            name => 'riboseq',
            message => 'Awake Awake Fear Fire Foes!  Go to page 5291772',
            choices => {
                '(biopieces): Make some plots of the demultiplexed/trimmed reads.' => \&Bio::Adventure::QA::Biopieces_Graph,
                '(cutadapt): Use cutadapt to remove the adapters.' => \&Bio::Adventure::Trim::Cutadapt,
                '(rrnabowtie): Map rRNA reads using bowtie1.' => \&Bio::Adventure::Map::Bowtie_RRNA,
                '(btmulti): Use bowtie1 to find putative ribosomal positions.' => \&Bio::Adventure::Map::BT_Multi,
                '(calibrate): Calibrate the positions of the a/p/e sites of the ribosomes.' => \&Bio::Adventure::Riboseq::Calibrate,
                '(countstates): Count the positions of a/p/e/etc sites across the genome/transcriptome.' => \&Bio::Adventure::Riboseq::Count_States,
                '(graphreads): Plot the coverage of ribosomes across the genome by ORF.' => \&Bio::Adventure::Riboseq::Graph_Reads,
            },
        },
        SNP => {
            name => 'snp',
            message => qq"When my god comes back I'll be waiting for him with a shotgun.  And I'm keeping the last shell for myself. (inexact quote)  Go to page 667408",
            choices => {
                '(trim): Trim sequence with an additional rule to remove the first 10 nucleotides.' => \&Bio::Adventure::Trim::Trimomatic,
                '(bwa): Map reads with bwa and count with htseq.' => \&Bio::Adventure::Map::BWA,
                '(bowtie): Map trimmed reads with bowtie1 and count with htseq.' => \&Bio::Adventure::Map::Bowtie,
                '(bt2): Map trimmed reads with bowtie2 and count with htseq.' => \&Bio::Adventure::Map::Bowtie2,
                '(hisat): Map trimmed reads with hisat2 and count with htseq.' => \&Bio::Adventure::Map::Hisat2,
                '(snpsearch): Perform a search for variant positions against a reference genome. (bam input)' => \&Bio::Adventure::SNP::SNP_Search,
                '(snpratio): Count the variant positions by position and create a new genome. (bcf input)' => \&Bio::Adventure::SNP::SNP_Ratio,
                '(snp): Perform alignments and search for variants. (fastq input)' => \&Bio::Adventure::SNP::Align_SNP_Search,
                '(snippy): Invoke snippy. (fastq and genbank inputs)' => \&Bio::Adventure::SNP::Snippy,
            },
        },
        TNSeq => {
            name => 'tnseq',
            message => 'You have enterred a world of jumping DNA, be ware and go to page 42.',
            choices => {
                '(sortindex): Demultiplex raw reads based on the peculiar TNSeq indexes.' => \&Bio::Adventure::TNSeq::Sort_Indexes,
                '(cutadapt): Use cutadapt to remove the odd tnseq adapters.' => \&Bio::Adventure::Trim::Cutadapt,
                '(tacheck): Make certain that all reads have a leading or terminal TA.' => \&Bio::Adventure::TNSeq::TA_Check,
                '(biopieces): Make some plots of the demultiplexed/trimmed reads.' => \&Bio::Adventure::QA::Biopieces_Graph,
                '(essentialityta): Count the hits/TA in preparation for essentiality.' => \&Bio::Adventure::TNSeq::Essentiality_TAs,
                '(runessentiality): Run the essentiality suite of tools.' => \&Bio::Adventure::TNSeq::Run_Essentiality,
                '(gumbel): Run the essentiality suite of tools on the ta counts.' => \&Bio::Adventure::TNSeq::Run_Essentiality,
                '(tpp): Run the transit preprocessing script.' => \&Bio::Adventure::TNSeq::Transit_TPP,
            },
        },
        Test => {
            name => 'test',
            message => 'All happy families are happy in the same way. Go to page 5670367.',
            choices => {
                '(test): Run a test job' => \&Test_Job,
            },
        },
    };
    return($menus);
}

sub Get_Job_Name {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars();
    my $name = 'unknown';
    $name = $options->{input} if ($options->{input});
    if ($name =~ /\:|\s+|\,/) {
        my @namelst = split(/\:|\s+|\,/, $name);
        $name = $namelst[0];
    }
    $name = basename($name, @{$options->{suffixes}});
    $name = basename($name, @{$options->{suffixes}});
    $name =~ s/\-trimmed//g;
    return($name);
}

sub Get_TODOs {
    my %args = @_;
    my $todo_list = ();
    my $possible_todos = {
        "aragorn+" => \$todo_list->{todo}{'Bio::Adventure::Annotation::Aragorn'},
        "biopieces+" => \$todo_list->{todo}{'Bio::Adventure::QA::Biopieces_Graph'},
        "blastmerge+" => \$todo_list->{todo}{'Bio::Adventure::Align_Blast::Merge_Blast_Parse'},
        "blastparse+" => \$todo_list->{todo}{'Bio::Adventure::Align_Blast::Blast_Parse'},
        "blastsplitalign+" => \$todo_list->{todo}{'Bio::Adventure::Align_Blast::Split_Align_Blast'},
        "bowtie+" => \$todo_list->{todo}{'Bio::Adventure::Map::Bowtie'},
        "bowtierrna+" => \$todo_list->{todo}{'Bio::Adventure::Map::Bowtie_RRNA'},
        "rrnabowtie+" => \$todo_list->{todo}{'Bio::Adventure::Map::Bowtie_RRNA'},
        "bt2+" => \$todo_list->{todo}{'Bio::Adventure::Map::Bowtie2'},
        "btmulti+" => \$todo_list->{todo}{'Bio::Adventure::Map::BT_Multi'},
        "bwa+" => \$todo_list->{todo}{'Bio::Adventure::Map::BWA'},
        "calibrate+" => \$todo_list->{todo}{'Bio::Adventure::Riboseq::Calibrate'},
        "copyraw+" => \$todo_list->{todo}{'Bio::Adventure::Prepare::Copy_Raw'},
        "countstates+" => \$todo_list->{todo}{'Bio::Adventure::Riboseq::Count_States'},
        "concat+" => \$todo_list->{todo}{'Bio::Adventure::Align::Concatenate_Searches'},
        "cutadapt+" => \$todo_list->{todo}{'Bio::Adventure::Trim::Cutadapt'},
        "essentialitytas+" => \$todo_list->{todo}{'Bio::Adventure::TNSeq::Essentiality_TAs'},
        "extendkraken+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Extend_Kraken_DB'},
        "extracttrinotate+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Extract_Trinotate'},
        "splitalignfasta+" => \$todo_list->{todo}{'Bio::Adventure::Align_Fasta::Split_Align_Fasta'},
        "fastado+" => \$todo_list->{todo}{'Bio::Adventure::Align_Fasta::Do_Fasta'},
        "fastasplitalign+" => \$todo_list->{todo}{'Bio::Adventure::Align_Fasta::Split_Align_Fasta'},
        "fastamerge+" => \$todo_list->{todo}{'Bio::Adventure::Align_Fasta::Merge_Parse_Fasta'},
        "fastaparse+" => \$todo_list->{todo}{'Bio::Adventure::Align_Fasta::Parse_Fasta'},
        "fastqct+" => \$todo_list->{todo}{'Bio::Adventure::QA::Fastqc'},
        "fastqdump+" => \$todo_list->{todo}{'Bio::Adventure::Prepare::Fastq_Dump'},
        "gb2gff+" => \$todo_list->{todo}{'Bio::Adventure::Convert::Gb2Gff'},
        "gff2fasta+" => \$todo_list->{todo}{'Bio::Adventure::Convert::Gff2Fasta'},
        "glimmer+" => \$todo_list->{todo}{'Bio::Adventure::Annotation::Glimmer'},    
        "graphreads+" => \$todo_list->{todo}{'Bio::Adventure::Riboseq::Graph_Reads'},
        "gubbins+" => \$todo_list->{todo}{'Bio::Adventure::Phylogeny::Run_Gubbins'},
        "gumbel+" => \$todo_list->{todo}{'Bio::Adventure::TNSeq::Run_Essentiality'},
        "hisat+" => \$todo_list->{todo}{'Bio::Adventure::Map::Hisat2'},
        "htmulti+" => \$todo_list->{todo}{'Bio::Adventure::Count::HT_Multi'},
        "hisat+" => \$todo_list->{todo}{'Bio::Adventure::Map::Hisat2'},
        "htseq+" => \$todo_list->{todo}{'Bio::Adventure::Count::HTSeq'},
        "indexhisat+" => \$todo_list->{todo}{'Bio::Adventure::Map::HT2_Index'},
        "indexbt1+" => \$todo_list->{todo}{'Bio::Adventure::Map::BT1_Index'},
        "indexbt2+" => \$todo_list->{todo}{'Bio::Adventure::Map::BT2_Index'},
        "indexbwa+" => \$todo_list->{todo}{'Bio::Adventure::Map::BWA_Index'},
        "indexkallisto+" => \$todo_list->{todo}{'Bio::Adventure::Map::Kallisto_Index'},
        "indexrsem+" => \$todo_list->{todo}{'Bio::Adventure::Map::RSEM_Index'},
        "indexsalmon+" => \$todo_list->{todo}{'Bio::Adventure::Map::Salmon_Index'},
        "kallisto+" => \$todo_list->{todo}{'Bio::Adventure::Map::Kallisto'},
        "kraken+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Kraken'},
        "mergeparse+" => \$todo_list->{todo}{'Bio::Adventure::Align_Blast::Merge_Parse_Blast'},
        "mimap+" => \$todo_list->{todo}{'Bio::Adventure::MiRNA::Mi_Map'},
        "pbt1+" => \$todo_list->{todo}{'Bio::Adventure::RNAseq_Pipeline_Bowtie'},
        "pbt2+" => \$todo_list->{todo}{'Bio::Adventure::RNAseq_Pipeline_Bowtie2'},
        "pbwa+" => \$todo_list->{todo}{'Bio::Adventure::RNAseq_Pipeline_BWA'},
        "pkallisto+" => \$todo_list->{todo}{'Bio::Adventure::RNAseq_Pipeline_Kallisto'},
        "ptophat+" => \$todo_list->{todo}{'Bio::Adventure::RNAseq_Pipeline_Tophat'},
        "ptnseq+" => \$todo_list->{todo}{'Bio::Adventure::TNseq_Pipeline'},
        "priboseq+" => \$todo_list->{todo}{'Bio::Adventure::Riboseq_Pipeline'},
        "parseblast+" => \$todo_list->{todo}{'Bio::Adventure::Align_Blast::Parse_Blast'},
        "posttrinity+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Trinity_Post'},
        "racer+" => \$todo_list->{todo}{'Bio::Adventure::Trim::Racer'},
        "readsample+" => \$todo_list->{todo}{'Bio::Adventure::Prepare::Read_Samples'},
        "resfinder+" => \$todo_list->{todo}{'Bio::Adventure::Annotation::Resfinder'},    
        "rsem+" => \$todo_list->{todo}{'Bio::Adventure::Map::RSEM'},
        "runessentiality+" => \$todo_list->{todo}{'Bio::Adventure::TNSeq::Run_Essentiality'},
        "sam2bam+" => \$todo_list->{todo}{'Bio::Adventure::Convert::Sam2Bam'},
        "salmon+" => \$todo_list->{todo}{'Bio::Adventure::Map::Salmon'},
        "shovill+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Shovill'},
        "snippy+" => \$todo_list->{todo}{'Bio::Adventure::SNP::Snippy'},
        "snp+" => \$todo_list->{todo}{'Bio::Adventure::SNP::Align_SNP_Search'},
        "snpsearch+" => \$todo_list->{todo}{'Bio::Adventure::SNP::SNP_Search'},
        "snpratio+" => \$todo_list->{todo}{'Bio::Adventure::SNP::SNP_Ratio'},
        "snpgenome+" => \$todo_list->{todo}{'Bio::Adventure::SNP::Make_Genome'},
        "sortindexes+" => \$todo_list->{todo}{'Bio::Adventure::TNSeq::Sort_Indexes'},
        "splitalign+" => \$todo_list->{todo}{'Bio::Adventure::Align::Split_Align'},
        "star+" => \$todo_list->{todo}{'Bio::Adventure::Map::STAR'},
        "tacheck+" => \$todo_list->{todo}{'Bio::Adventure::TNSeq::TA_Check'},
        "test+" => \$todo_list->{todo}{'Bio::Adventure::PBS::Test_Job'},
        "tophat+" => \$todo_list->{todo}{'Bio::Adventure::Map::Tophat'},
        "tpp+" => \$todo_list->{todo}{'Bio::Adventure::TNSeq::Transit_TPP'},
        "transdecoder+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Transdecoder'},
        "trimomatic+" => \$todo_list->{todo}{'Bio::Adventure::Trim::Trimomatic'},
        "trinity+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Trinity'},
        "trinitypost+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Trinity_Post'},
        "trinotate+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Trinotate'},
        "tritrypdownload+" => \$todo_list->{todo}{'Bio::Adventure::Convert::TriTryp_Download'},
        "tritryp2text+" => \$todo_list->{todo}{'Bio::Adventure::Convert::TriTryp2Text'},
        "trnascan+" => \$todo_list->{todo}{'Bio::Adventure::Annotation::tRNAScan'},    
        "variantgenome+" => \$todo_list->{todo}{'Bio::Adventure::SNP::Make_Genome'},
        "velvet+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Velvet'},
        "help+" => \$todo_list->{todo}{'Bio::Adventure::Adventure_Help'},
    };

    if (defined($args{method})) {
        my $m = '--' . $args{method};
        my @method = ($m,);
        my $array_result = GetOptionsFromArray(\@method, %{$possible_todos});
    }

    my @methods_to_run = ();
    my $set_of_todos = $todo_list->{todo};
    foreach my $job (keys %{$set_of_todos}) {
        my $thing = $set_of_todos->{$job};
        if ($thing) {
            if (!defined($args{task})) {
                print "Going on an adventure to ${job}.\n";
            } else {
                print "Going on an adventure to ${job} in the $args{task} context.\n";
            }
            my $final = \&${job};
            push(@methods_to_run, $final);
        }
    }
    return(\@methods_to_run);
}

=item C<Get_Input>

    Get_Input() attempts to standardize the inputs passed to Bio::Adventure.
    It returns a stringified and standard representation of the likely
    input(s) to the script.

    There are a few problems with how I send input to these scripts:
    Sometimes I put in --hpgl hpgl0415 when I mean --hpgl HPGL0415,
    Sometimes I put in --input hpgl0415.fastq  when I mean hpgl0415.fastq.(gz|xz),
    Sometimes I put in -i hpgl0415.fastq.(gz|xz) when I mean hpgl0415.fastq,
    Sometimes I put in -i hpgl0415.fastq when I mean hpgl0415-trimmed.fastq(.gz|.xz)
    Sometimes I put in -i hpgl0415_forward.fastq:hpgl0415_reverse.fastq

    So, this function should make this unambiguous and consistent no matter what I type

=cut
sub Get_Input {
    my $class = shift;
    my $input = $class->{input};
    my $id = $class->{hpgl};
    my $actual = "";

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
            qx(qq"mv ${in}.fastq ${low}.fastq");
            if (-r "${in}_forward.fastq") {
                qx(qq"mv ${in}_forward.fastq ${low}_forward.fastq");
            }
            if (-r "${in}_reverse.fastq") {
                qx(qq"mv ${in}_reverse.fastq ${low}_reverse.fastq");
            }
        }
        ## if (-r "${in}" and $in =~ /\.gz$/) {
        ##   print "gunzip of ${in}\n";
        ##   system("nice gunzip ${in}");
        ## }
        ## if (-r "${in}" and $in =~ /\.xz$/) {
        ##   print "xz -d ${in}\n";
        ##   system("nice xz -d ${in}");
        ## }
    }
    $input =~ tr/[A-Z]/[a-z]/;
    ## $input =~ s/\.gz//g;
    ## $input =~ s/\.xz//g;
    return($input);
}

=item C<Get_Vars>

  Handle the peculiar mix of instance options held in $class->{options},
  the set of arguments passed to %args, a list of required and potentially
  missing options, and whatever default values one wishes to set.

  my $options = $class->Get_Vars(args => \%args, required => ['input', 'genome'], bob => 'jane');

=cut
sub Get_Vars {
    my ($class, %args) = @_;
    my $arglist = $args{args};
    my $reqlist = $args{required};
    my %remaining = ();
    my $options = $class->{options};

    ## Use this for loop to fill in default values from the class.
    for my $varname (keys %{$class}) {
        next if ($varname eq 'options');
        if (!defined($options->{$varname})) {
            $options->{$varname} = $class->{$varname};
        }
    }

    ## Use this loop to fill in default values from the argument hash.
    for my $varname (keys %args) {
        next if ($varname eq 'required');
        next if ($varname eq 'args');
        if (!defined($options->{$varname})) {
            $options->{$varname} = $args{$varname};
        }
    }

    ## Options in arglist take precedence.
    for my $k (keys %{$arglist}) {
        $options->{$k} = $arglist->{$k};
    }

    ## Check for the definition of required arguments last.
    my $needed_list = [];
    if (defined($args{required})) {
        $needed_list = $args{required};
    }
    my @options = @{$needed_list};
    foreach my $option (@options) {
        if (!defined($options->{$option})) {
            my $query = qq"The option: ${option} was missing, please fill it in:";
            if (!defined($options->{term})) {
                $options->{term} = Get_Term();
            }
            my $response = $options->{term}->readline($query);
            $response =~ s/\s+$//g;
            $response =~ s/\@|\*|\+//g;
            $options->{$option} = $response;
        }
    }

    ## Final sanity check(s)
    for my $k (keys %{$options}) {
        next unless (defined($options->{$k}));
        $options->{$k} =~ s/^\~/${HOME}/g;

        ## Finally, fill in a few special cases here.
        if ($k eq 'input') {
            my $uncomp = basename($options->{input}, (".gz", ".bz2", ".xz"));
            if (!-r $options->{input}) {
                if (-r $uncomp) {
                    $options->{input} = $uncomp;
                } else {
                    ## Do nothing, we must assume that the file will exist when required.
                }
            }
        }
        if ($k eq 'jname') {
            if ($options->{jname} eq '') {
                my $name = 'unknown';
                $name = $options->{input} if ($options->{input});
                if ($name =~ /\:|\s+|\,/) {
                    my @namelst = split(/\:|\s+|\,/, $name);
                    $name = $namelst[0];
                }
                $name = basename($name, (".gz", ".xz", ".bz2", ".bai", ".fai"));
                $name = basename($name, (".fasta", ".fastq", ".bam", ".sam", ".count"));
                $name =~ s/_forward//g;
                $options->{jname} = $name;
            }

        }
    }
    ## End special cases.

    return($options);
}

=item C<Set_Vars>

  Handle the peculiar mix of instance options held in $class->{options},
  the set of arguments passed to %args, a list of required and potentially
  missing options, and whatever default values one wishes to set.

  $options = $class->Set_Vars(exclude => 'bob');

=cut
sub Set_Vars {
    my ($class, %args) = @_;
    my $options = $class->{options};
    foreach my $k (keys %args) {
        $options->{$k} = $args{$k};
    }
    $class->{options} = $options;
    return($options);
}

sub My_Which {
    my $executable = shift;
    my $result = which($executable);
    if (!defined($result)) {
        return(undef);
    } else {
        return($result);
    }
}

=item C<Last_Stat>

    Last_Stat() reads the final line of an input file and returns it
    as a string.

    This is useful because many of these tools append to .csv files
    with summary information (alignment statistics, sequence sizes,
    etc) and the tests keep a set of expected outputs for the most
    recent run.

=cut
sub Last_Stat {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $input_filename = $options->{input};
    my $input = FileHandle->new("<${input_filename}");
    ##$input->open("<$input_filename");
    my ($line, $last);
    while ($line = <$input>) {
        chomp $line;
        ## I am doing this due to terminal newlines on the file, there is probably a much smrtr way.
        $last = $line unless ($line =~ /^$/);
    }
    $input->close();
    return($last);
}

=item C<Read_Genome_Fasta>

    Read a fasta file and return the chromosomes.

=cut
sub Read_Genome_Fasta {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $chromosomes = {};
    if (!defined($options->{genome})) {
        if (defined($options->{fasta})) {
            $options = $class->Set_Vars(genome => $options->{fasta});
        }
    }
    my $fasta_name = basename($options->{genome}, ['.fasta']);
    my $data_file = qq"$options->{basedir}/${fasta_name}.pdata";
    if (-r $data_file) {
        $chromosomes = retrieve($data_file);
    } else {
        my $input_genome = Bio::SeqIO->new(-file => $options->{genome}, -format => 'Fasta');
        while (my $genome_seq = $input_genome->next_seq()) {
            next unless(defined($genome_seq->id));
            my $id = $genome_seq->id;
            my $sequence = $genome_seq->seq;
            my $length = $genome_seq->length;
            my $empty_forward = [];
            my $empty_reverse = [];
            for my $c (0 .. $length) {
                $empty_forward->[$c] = 0;
                $empty_reverse->[$c] = 0;
            }
            $chromosomes->{$id}->{forward} = $empty_forward;
            $chromosomes->{$id}->{sequence} = $sequence;
            $chromosomes->{$id}->{reverse} = $empty_reverse;
            $chromosomes->{$id}->{obj} = $genome_seq;
        }                       ## End reading each chromosome
        store($chromosomes, $data_file);
    }                           ## End checking for a .pdata file
    return($chromosomes);
}

=item C<Read_Genome_GFF>

    Read a GFF file and extract the annotation information from it.

=cut
sub Read_Genome_GFF {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['gff'],
        feature_type => 'exon',
        id_tag => 'gene_id');
    my $feature_type = $options->{feature_type};
    my $id_tag = $options->{id_tag};
    my $annotation_in = Bio::Tools::GFF->new(-file => "$options->{gff}", -gff_version => 3);
    my $gff_out = { stats => { chromosomes => [],
                               lengths => [],
                               feature_names => [],
                               cds_starts => [],
                               cds_ends => [],
                               inter_starts => [],
                               inter_ends => [],
                    },};
    my $gff_name = basename($args{gff}, ['.gff']);
    my $data_file = qq"$options->{basedir}/${gff_name}.pdata";
    print "Reading $options->{gff}, seeking features tagged (last column ID) $id_tag,
 type (3rd column): $feature_type.\n";
    if (-r $data_file && !$options->{debug}) {
        $gff_out = retrieve($data_file);
    } else {
        print "Starting to read gff: $args{gff}\n" if ($class->{debug});
        my $hits = 0;
        my @chromosome_list = @{$gff_out->{stats}->{chromosomes}};
        my @feature_names = @{$gff_out->{stats}->{feature_names}};
        my @cds_starts = @{$gff_out->{stats}->{cds_starts}};
        my @inter_starts = @{$gff_out->{stats}->{inter_starts}};
        my @cds_ends = @{$gff_out->{stats}->{cds_ends}};
        my @inter_ends = @{$gff_out->{stats}->{inter_ends}};
        my $start = 1;
        my $old_start = 1;
        my $end = 1;
        my $old_end = 1;
      LOOP: while(my $feature = $annotation_in->next_feature()) {
          next LOOP unless ($feature->{_primary_tag} eq $feature_type);
          $hits++;
          my $location = $feature->{_location};
          $old_start = $start;
          $old_end = $end;
          $start = $location->start();
          $end = $location->end();
          my $strand = $location->strand();
          my @ids = $feature->each_tag_value($id_tag);
          my $id = "";
          my $gff_chr = $feature->{_gsf_seq_id};
          my $gff_string = $annotation_in->gff_string($feature);
          foreach my $i (@ids) {
              $i =~ s/^cds_//g;
              $i =~ s/\-\d+$//g;
              $id .= "$i ";
          }
          $id =~ s/\s+$//g;
          my @gff_information = split(/\t+/, $gff_string);
          my $description_string = $gff_information[8];
          my $orf_chromosome = $gff_chr;
          ## Add the chromosome to the list of chromosomes if it is not there already.
          push(@chromosome_list, $orf_chromosome) if ($orf_chromosome !~~ @chromosome_list);
          push(@feature_names, $id);
          push(@cds_starts, $start);
          push(@inter_starts, $old_start+1);
          push(@cds_ends, $end);
          ##push(@inter_ends, $old_start-1);
          push(@inter_ends, $start-1);
          my $annot = {
              id => $id,
              start => $start, ## Genomic coordinate of the start codon
              end => $end, ## And stop codon
              strand => $strand,
              description_string => $description_string,
              chromosome => $gff_chr,
          };
          $gff_out->{$gff_chr}->{$id} = $annot;
      } ## End looking at every gene in the gff file
        $gff_out->{stats}->{chromosomes} = \@chromosome_list;
        $gff_out->{stats}->{feature_names} = \@feature_names;
        $gff_out->{stats}->{cds_starts} = \@cds_starts;
        $gff_out->{stats}->{cds_ends} = \@cds_ends;
        $gff_out->{stats}->{inter_starts} = \@inter_starts;
        $gff_out->{stats}->{inter_ends} = \@inter_ends;
        print STDERR "Not many hits were observed, do you have the right feature type?  It is: ${feature_type}\n" if ($hits < 1000);
        if (-f $data_file) {
            unlink($data_file);
        }
        store($gff_out, $data_file);
    } ## End looking for the gff data file
    return($gff_out);
}

sub Reorder_Fasta {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $genome = $class->Read_Genome_Fasta(%args);

    my @chromosomes = sort { $a cmp $b } keys %$genome;
    my $out = FileHandle->new(">$options->{output}");
    foreach my $c (@chromosomes) {
        my $start = $genome->{c};
        my $end = join("\n", ($start =~ m/.{1,80}/g));
        print $out ">$c
$end
";
    }
    $out->close();
}

sub Submit {
    my ($class, %args) = @_;
    ## Adding this to avoid overwriting data in $class, not sure if this is correct.
    my $runner = $class;
    my $opts = $runner->Get_Vars(args => \%args);
    my %options = %{$opts};

    ## If we are invoking an indirect job, we need a way to serialize the options
    ## in order to get them passed to the eventual interpreter
    my $option_file = "";
    if ($options{language} eq 'perl') {
        ## I think this might be required as per:
        ## https://metacpan.org/pod/release/AMS/Storable-2.21/Storable.pm#CODE_REFERENCES
        $Storable::Deparse = 1;
        $Storable::Eval = 1;
        $option_file = File::Temp->new(
            TEMPLATE => 'optionsXXXX',
            DIR => $options{basedir},
            SUFFIX => '.pdata',
            );
        my $opt = \%options;
        ## Code references are invalid for these things...
        $opt->{menus} = undef;
        $opt->{methods_to_run} = undef;
        ## Why is it that periodically I get this error?
        ## The result of B::Deparse::coderef2text was empty - maybe you're trying to serialize an XS function?
        my $stored = nstore(\%options, $option_file);
        ##my $stored = store(\%options, $option_file);
        ##$utf8_encoded_json_text = encode_json $perl_hash_or_arrayref;
        ##$perl_hash_or_arrayref  = decode_json $utf8_encoded_json_text;
        ##my $json = JSON->new->allow_nonref->convert_blessed;
        ##my $stored_text = $json->encode($opt);
        ##my $stored_file = FileHandle->new(">${option_file}");
        ##print $stored_file $stored_text;
        ##$stored_file->close();
        $args{option_file} = $option_file;
    }
    my $job;
    if ($runner->{sbatch_path}) {
        $job = Bio::Adventure::Slurm->new();
    } elsif ($runner->{qsub_path}) {
        $job = Bio::Adventure::Torque->new();
    } elsif ($runner->{bash_path}) {
        ## I should probably have something to handle gracefully bash jobs.
        $job = Bio::Adventure::Local->new();
    } else {
        die("Could not find sbatch, qsub, nor bash.");
    }
    my $result = $job->Submit(%args);
    return($result);
}

sub Pipeline_Riboseq {
    my ($class, %args) = @_;
    my $fastqc_job = Bio::Adventure::QA::Fastqc($class, %args);
    my $cutadapt_job = Bio::Adventure::Trim::Cutadapt($class, %args);
    $args{depends} = $cutadapt_job->{pbs_id};
    my $biopieces = Bio::Adventure::QA::Biopieces_Graph($class, %args);
    my $rrna_job = Bio::Adventure::Map::Bowtie_RRNA($class, %args);
    $args{depends} = $rrna_job->{pbs_id};
    my $bt_jobs = Bio::Adventure::Map::Bowtie($class, %args);
    my $ret = {
        fastqc => $fastqc_job,
        cutadapt => $cutadapt_job,
        biopieces => $biopieces,
        rrna => $rrna_job,
        bt => $bt_jobs};
    return($ret);
}

sub Pipeline_RNAseq_Bowtie {
    my ($class, %args) = @_;
    $args{aligner} = 'bowtie';
    my $rnaseq_jobs = $class->Pipeline_RNAseq(%args);
    return($rnaseq_jobs);
}

sub Pipeline_RNAseq_Bowtie2 {
    my ($class, %args) = @_;
    $args{aligner} = 'bowtie2';
    my $rnaseq_jobs = $class->Pipeline_RNAseq(%args);
    return($rnaseq_jobs);
}

sub Pipeline_RNAseq_Tophat {
    my ($class, %args) = @_;
    $args{aligner} = 'tophat';
    my $rnaseq_jobs = $class->Pipeline_RNAseq(%args);
    return($rnaseq_jobs);
}

sub Pipeline_RNAseq_BWA {
    my ($class, %args) = @_;
    $args{aligner} = 'bwa';
    my $rnaseq_jobs = $class->Pipeline_RNAseq(%args);
    return($rnaseq_jobs);
}

sub Pipeline_RNAseq_Kallisto {
    my ($class, %args) = @_;
    $args{aligner} = 'kallisto';
    my $rnaseq_jobs = $class->Pipeline_RNAseq(%args);
    return($rnaseq_jobs);
}

sub Pipeline_RNAseq {
    my ($class, %args) = @_;
    my $fastq_job = Bio::Adventure::QA::Fastqc($class, %args);
    my $trim_job = Bio::Adventure::Trim::Trimomatic($class, %args);
    $args{depends} = $trim_job->{pbs_id};
    my $biopieces_job = Bio::Adventure::QA::Biopieces_Graph($class, %args);
    my $rrna_job = Bio::Adventure::Map::Bowtie_RRNA($class, %args);
    $args{depends} = $rrna_job->{pbs_id};
    my $align_jobs;
    if ($args{aligner} eq 'bowtie') {
        $align_jobs = Bio::Adventure::Map::Bowtie($class, %args);
    } elsif ($args{aligner} eq 'bowtie2') {
        $align_jobs = Bio::Adventure::Map::Bowtie2($class, %args);
    } elsif ($args{aligner} eq 'tophat') {
        $align_jobs = Bio::Adventure::Map::Tophat($class, %args);
    } elsif ($args{aligner} eq 'bwa') {
        $align_jobs = Bio::Adventure::Map::BWA($class, %args);
    } elsif ($args{aligner} eq 'kallisto') {
        $align_jobs = Bio::Adventure::Map::Kallisto($class, %args);
    } else {
        $align_jobs = Bio::Adventure::Map::Tophat($class, %args);
    }

    my $ret = {
        fastqc => $fastq_job,
        trim => $trim_job,
        bioieces => $biopieces_job,
        rrna => $rrna_job,
        mapping => $align_jobs,
    };
    return($ret);
}

sub Pipeline_TNseq {
    my ($class, %args) = @_;
    my $fastqc_job = Bio::Adventure::QA::Fastqc($class, %args);
    $args{type} = 'tnseq';
    my $cutadapt_job = Bio::Adventure::Trim::Cutadapt($class, %args);
    $args{depends} = $cutadapt_job->{pbs_id};
    my $biopieces_job = Bio::Adventure::QA::Biopieces_Graph($class, %args);
    my $bt_jobs = Bio::Adventure::Map::Bowtie($class, %args);
    my $ret = {
        fastqc => $fastqc_job,
        cutadapt => $cutadapt_job,
        biopieces => $biopieces_job,
        bt => $bt_jobs,
    };
    return($ret);
}

sub Adventure_Help {
    my ($class, %args) = @_;
    my %methods = %{$class->{methods}};
    print qq"The command line program 'cyoa' has a series of shortcuts intended to make it easy to use and flexible.
The following comprises the set of strings you may feed it as 'methods':\n";
    my $c = 0;
    for my $k (sort keys %methods) {
        my $sep = '\t';
        if (($c % 3) == 0) {
            $sep = '\n';
        } else {
            $sep = '\t';
        }
        $c++;
        print "${k}${sep}";
    }
    print "\n";
    print qq"cyoa uses GetOptions, so you can shortcut all the 'methods', so:
'cyoa --task ri --method cut --input test.fastq' calls the cutadapt with
  options suitable for ribosome profiling data.  Conversely:
'cyoa --task rna --method top --input test.fastq' calls tophat assuming
  rnaseq data.
'cyoa --method blastsplit --query test.fasta --library nr --blast_tool blastp'
  Splits the test.fasta into a bunch of pieces (settable with --number), calls
  blastp on them, merges the outputs, and parses the result into a table of hits.
";
    $class->Help();
    return(0);
}

## Make empty functions for the stuff provided by Bio::Adventure::


=back

=head1 AUTHOR - atb

Email abelew@gmail.com

=head1 SEE ALSO

    L<Bio::Seq> L<autodie> L<local::lib> L<Bio::Adventure::Aligners>
    L<Bio::Adventure::Assembly> L<Bio::Adventure::Count> L<Bio::Adventure::Aligners> L<Bio::Adventure::QA>
    L<Bio::Adventure::Trim> L<Bio::Adventure::Align_Blast> L<Bio::Adventure::Align_Fasta> L<Bio::Adventure::Align>
    L<Bio::Adventure::Cleanup> L<Bio::Adventure::Compress> L<Bio::Adventure::Convert> L<Bio::Adventure::PBS>
    L<Bio::Adventure::Prepare> L<Bio::Adventure::Riboseq> L<Bio::Adventure::SeqMisc> L<Bio::Adventure::SNP>
    L<Bio::Adventure::TNSeq>

=cut

1;
