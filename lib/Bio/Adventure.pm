package Bio::Adventure;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
use vars qw"$VERSION";
use feature qw"try";
no warnings qw"experimental::try";

use AppConfig qw":argcount :expand";
use Bio::SeqIO;
use Bio::DB::Universal;
use Bio::Root::RootI;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Capture::Tiny qw":all";
use Carp qw"croak carp confess cluck longmess";
use Cwd qw"abs_path getcwd cwd";
use Data::Dumper;
use Env qw"COMPRESSION XZ_OPTS XZ_DEFAULTS HOME";
use Env::Modulecmd;
use File::Basename;
use File::Find;
use File::Spec;
use File::Which qw"which";
use File::Path qw"make_path remove_tree";
use File::Temp;
use FileHandle;
use Getopt::Long qw"GetOptionsFromArray";
use IO::String;
use Log::Log4perl;
use PerlIO;
use Pod::Usage;
use Signal::StackTrace qw"TTIN";
use Storable qw"lock_store lock_retrieve";
use Term::ReadLine;
use Term::UI;

use Bio::Adventure::Config;
use Bio::Adventure::Slurm;
use Bio::Adventure::Local;

use Bio::Adventure::Align;
use Bio::Adventure::Align_Blast;
use Bio::Adventure::Align_Fasta;
use Bio::Adventure::Annotation;
use Bio::Adventure::Annotation_Genbank;
use Bio::Adventure::Assembly;
use Bio::Adventure::Cleanup;
use Bio::Adventure::Count;
use Bio::Adventure::Compress;
use Bio::Adventure::Convert;
use Bio::Adventure::Feature_Prediction;
use Bio::Adventure::Index;

use Bio::Adventure::Map;
use Bio::Adventure::Metadata;
use Bio::Adventure::Phage;
use Bio::Adventure::Phylogeny;
use Bio::Adventure::PopGen;
use Bio::Adventure::Prepare;
use Bio::Adventure::Resistance;
use Bio::Adventure::Riboseq;
use Bio::Adventure::QA;
use Bio::Adventure::SeqMisc;
use Bio::Adventure::Structure;

use Bio::Adventure::SNP;
use Bio::Adventure::TNSeq;
use Bio::Adventure::Torque;
use Bio::Adventure::Trim;
use Bio::Adventure::Visualization;
use Bio::Adventure::Pipeline;

=head1 NAME

    Bio::Adventure - A perl library to make preprocessing high-throughput data easier.

=head1 SYNOPSIS

    Bio::Adventure.pm  Methods to simplify creating job scripts and submitting them to a
    computing cluster.  Bio::Adventure makes heavy use of command line options and as a result
    accepts many options.  The simplest way to access these methods is via the cyoa script, but
    everything may be called directly in perl.

    use Bio::Adventure;
    my $hpgl = Bio::Adventure->new(cluster => 0);  ## Force it to run locally.
    $hpgl->{species} = 'mmusculus';  ## There had better be a mmusculus.[fasta&gff] in $hpgl->{libdir}
    my $hpgl = Bio::Adventure->new(species => 'scerevisiae', libdir => '/home/bob/libraries');
    $hpgl->{input} = 'test.fastq';  ## Also via with cyoa --input=test.fastq

    ## Run Fastqc on untrimmed sequence.
    my $bp = $hpgl->Fastqc()
    ## Generates test-trimmed.fastq from test.fastq
    my $trim = $hpgl->Trimomatic();

    Conversely, just run 'cyoa' from the command line and follow the prompts and/or provide
    all the parameters with the appropriate command line options:

    > cyoa --task map --method hisat --species hg38_100 --input r1-trimmed.fastq.xz:r2-trimmed.fastq.xz
    > cyoa --task align --method blast --library nr --input transcripts.fasta

=head1 DESCRIPTION

    This library should write out PBS/slurm compatible job files and
    submit them to the appropriate computing cluster.  It should also
    collect the outputs and clean up the mess.

=head2 Options

  Most options for cyoa are defined in Adventure.pm, there are a small number are used
  only for specific functions defined within them.  The document L<Bio::Adventure/Arguments.pm>
  will seek to explain how to access them and their various purposes.
  the most likely purpose for each parameter follows:

=cut

## Lets move all the default values here.
has align_blast_format => (is => 'rw', default => 5); ## Which alignment type should we use? (5 is blastxml)
has align_jobs => (is => 'rw', default => 40); ## How many blast/fasta alignment jobs should we make when splitting alignments across nodes?
has align_parse => (is => 'rw', default => 1); ## Parse blast searches into a table?
has arbitrary => (is => 'rw', default => ''); ## Extra arbitrary arguments to pass
has array_start => (is => 'rw', default => 100);
has bamfile => (is => 'rw', default => undef); ## Default bam file for converting/reading/etc.
has basedir => (is => 'rw', default => cwd());  ## This was cwd() but I think that may cause problems.
has bash_path => (is => 'rw', default => scalar_which('bash'));
has best_only => (is => 'rw', default => 0); ## keep only the best search result when performing alignments?
has blast_params => (is => 'rw', default => ' -e 10 '); ## Default blast parameters
has blast_tool => (is => 'rw', default => 'blastn'); ## Default blast tool to use
has bt_default => (is => 'rw', default => '--best'); ## Default bt1 arguments.
has bt_varg => (is => 'rw', default => '-v 0');
has bt_marg => (is => 'rw', default => '-M 0');
has bt_larg => (is => 'rw', default => '-y -l 15');
has bt2_args => (is => 'rw', default => ' --very-sensitive -L 14 '); ## My favorite bowtie2 arguments
has btmulti => (is => 'rw', default => 0); ## Perform multiple bowtie searches?
has bwa_method => (is => 'rw', default => 'mem,aln'); ## Default bwa method to use.
has chosen_tag => (is => 'rw', default => 'ODDS');
has clean => (is => 'rw', default => 0); ## Cleanup after yourself?
has cluster => (is => 'rw', default => undef); ## Are we running on a cluster?
has comment => (is => 'rw', default => undef); ## Set a comment in running slurm/bash/etc scripts.
has compress => (is => 'rw', default => 1); ## Compress output files?
has config => (is => 'rw', default => undef); ## Not sure
has count => (is => 'rw', default => 1); ## Quantify reads after mapping?
has coverage => (is => 'rw', default => undef); ## Provide a coverage cutoff
has coverage_tag => (is => 'rw', default => 'DP');
has csv_file => (is => 'rw', default => 'all_samples.csv'); ## Default csv file to read/write.
has cutoff => (is => 'rw', default => 0.05); ## Default cutoff (looking at your vcftools, e.g. I haven't changed those yet).
has decoy => (is => 'rw', default => 1); ## Add decoys
has debug => (is => 'rw', default => 0); ## Print debugging information.
has directories => (is => 'rw', default => undef); ## Apply a command to multiple input directories.
has download => (is => 'rw', default => 1);
has email => (is => 'rw', default => 'abelew@umd.edu');
has evalue => (is => 'rw', default => 0.001); ## Default e-value cutoff
has fasta_args => (is => 'rw', default => ' -b 20 -d 20 '); ## Default arguments for the fasta36 suite
has fasta_tool => (is => 'rw', default => 'ggsearch36'); ## Which fasta36 program to run?
has filtered => (is => 'rw', default => 'unfiltered');  ## Whether or not Fastqc is running on filtered data.
has freebayes => (is => 'rw', default => 0);
has fsa_input => (is => 'rw'); ## fsa genome output file for creating a genbank file
has gcode => (is => 'rw', default => '11'); ## Choose a genetic code
has genome => (is => 'rw', default => undef); ## Choose a genome to work on.
has genus => (is => 'rw', default => undef); ## Choose a genus when using prokka and potentially others like kraken
has gff => (is => 'rw', default => undef); ## Feature file to read/write
has gff_tag => (is => 'rw', default => 'gene_id'); ## Likely redundant with htseq_id
## Ahh I remember, htseq_type was added to facilitate performing multiple htseq counts on multiple gff files.
## Notably the interCDS for bacterial genomes.
has gff_type => (is => 'rw', default => 'gene'); ## When blank, do it on the whole gff file, otherwise use that suffix.
has gff_cds_parent_type => (is => 'rw', default => 'mRNA');
has gff_cds_type => (is => 'rw', default => 'CDS');
has help => (is => 'rw', default => undef); ## Ask for help?
has hisat_args => (is => 'rw', default => ' --sensitive ');
has host_filter => (is => 'rw', default => 1);  ## When performing an assembly, do a host filter?
has htseq_args => (is => 'rw', default => ' --order=name --idattr=gene_id --minaqual=10 --type=exon --stranded=yes --mode=union '); ## Most likely htseq options
has identity => (is => 'rw', default => 70); ## Alignment specific identity cutoff
has index_file => (is => 'rw', default => 'indexes.txt'); ## File containing indexes:sampleIDs when demultiplexing samples - likely tnseq
has index_hash => (is => 'rw', default => undef);
has input => (is => 'rw', default => undef); ## Generic input argument
has input_abricate => (is => 'rw', default => 'outputs/12abricate_10prokka_09termreorder_08phageterm_07rosalind_plus/abricate_combined.tsv'); ## Used when merging annotation files into a xlsx/tbl/gbk file.
has input_aragorn => (is => 'rw', default => 'outputs/21aragorn/aragorn.txt'); ## Used when merging annotation files into a xlsx/tbl/gbk file.
has input_classifier => (is => 'rw', default => 'outputs/18classifier/ictv_filtered.tsv'); ## Similar taxa detected by tblastx
has input_genbank => (is => 'rw', default => undef); ## Existing genbank file for merging annotations.
has input_glimmer => (is => 'rw', default => 'outputs/16glimmer/glimmer.predict');
has input_fastq => (is => 'rw', default => undef);
has input_interpro => (is => 'rw', default => 'outputs/13_interproscan_10prokka_09termreorder_08phageterm_07rosalind_plus/interproscan.tsv'); ## interpro output file when merging annotations.
has input_phageterm => (is => 'rw', default => 'outputs/08phageterm_07rosalind_plus/direct-term-repeats.fasta'); ## phageterm output file when merging annotations.
has input_phanotate => (is => 'rw', default => 'outputs/16phanotate/phanotate.tsv.xz');
has input_prodigal => (is => 'rw', default => 'outputs/17prodigal/predicted_cds.gff');
has input_prokka_tsv => (is => 'rw', default => undef); ## Prokka tsv file for merging annotations.
has input_trinotate => (is => 'rw', default => '11trinotate_10prokka_09termreorder_08phageterm_07rosalind_plus/Trinotate.tsv'); ## trinotate output, used when merging annotations.
has input_umi => (is => 'rw', default => 'umi.txt');
has interactive => (is => 'rw', default => 0); ## Is this an interactive session?
has introns => (is => 'rw', default => 0); ## Is this method intron aware? (variant searching).
has jobs => (is => 'rw', default => undef); ## List of currently active jobs, possibly not used right now.
has jobids => (is => 'rw', default => undef); ## A place to put running jobids, maybe no longer needed.
has jbasename => (is => 'rw', default => basename(cwd())); ## Job basename
has jcpu => (is => 'rw', default => 2); ## Number of processors to request in jobs
has jgpu => (is => 'rw', default => 0);
has jdepends => (is => 'rw', default => '');  ## Flag to start a dependency chain
has jmem => (is => 'rw', default => 20); ## Number of gigs of ram to request
has jname => (is => 'rw', default => undef); ## Job name on the cluster
has jnice => (is => 'rw', default => 0); ## Set the niceness of a job, if it starts positive, we can set a lower nice to preempt
has jpartition => (is => 'rw', default => 'dpart');
has jprefix => (is => 'rw', default => ''); ## Prefix number for the job
has jqueue => (is => 'rw', default => 'workstation'); ## What queue will jobs default to?
has jqueues => (is => 'rw', default => 'throughput,workstation,long,large'); ## Other possible queues
has jsleep => (is => 'rw', default => '0.5'); ## Set a sleep between jobs
has jstring => (is => 'rw', default => undef); ## String of the job
has jtemplate => (is => 'rw', default => undef);
has jwalltime => (is => 'rw', default => '10:00:00'); ## Default time to request
has kingdom => (is => 'rw', default => undef); ## Taxonomic kingdom, prokka/kraken
has language => (is => 'rw', default => 'bash'); ## What kind of script is this?
has length => (is => 'rw', default => 17); ## kmer length, other stuff too.
has libdir => (is => 'rw', default => "\${HOME}/libraries"); ## Directory containing genomes/gff/indexes
has libpath => (is => 'rw', default => "$ENV{HOME}/libraries");
has library => (is => 'rw', default => undef);  ## The library to be used for fasta36/blast searches
has libtype => (is => 'rw', default => 'genome'); ## Type of sequence to map against, genomic/rRNA/contaminants
has locus_tag => (is => 'rw', default => undef); ## Used by prokka to define gene prefixes
has logdir => (is => 'rw', default => 'outputs/logs'); ## place to dump logs
has loghost => (is => 'rw', default => 'localhost'); ## Host to which to send logs
has mapper => (is => 'rw', default => 'hisat'); ## Use this aligner if none was chosen.
has mature_fasta => (is => 'rw', default => undef); ## Database of mature miRNA sequences to search
has maximum => (is => 'rw', default => undef);  ## catchall maximum threshold
has maxlength => (is => 'rw', default => 42); ## Maximum sequence length when trimming
has method => (is => 'rw', default => undef);
has mi_genome => (is => 'rw', default => undef); ## Set a miRbase genome to hunt for mature miRNAs
has min_depth => (is => 'rw', default => 5); ## Default use: variant searching, depth limit
has min_value => (is => 'rw', default => 0.5);  ## Also variant searching.
has minimum => (is => 'rw', default => undef); ## catchall minimum threshold
has minlength => (is => 'rw', default => 8); ## Minimum length when trimming
has mirbase_data => (is => 'rw', default => undef); ## miRbase annotation dataset.
has modulecmd => (is => 'rw', default => '');
has modules => (is => 'rw', default => undef); ## Environment modules to load
has module_string => (is => 'rw', default => '');
has option_file => (is => 'rw', default => undef);
has orientation => (is => 'rw', default => 'start'); ## Default orientation when normalizing riboseq reads
has outgroup => (is => 'rw', default => undef); ## Outgroup for phylogenetic tools
has output => (is => 'rw', default => undef); ## Generic output argument
has outdir => (is => 'rw', default => undef);
has overlap => (is => 'rw', default => 20);
has paired => (is => 'rw', default => 1); ## Is the input paired?
has pdata => (is => 'rw', default => 'options.pdata');
has phred => (is => 'rw', default => 33); ## Minimum quality score when trimming
has postscript => (is => 'rw', default => undef); ## String to put after a cluter job.
has preprocess_dir => (is => 'rw', default => 'preprocessing');
has prescript => (is => 'rw', default => undef); ## String to put before a cluster job.
has primary_key => (is => 'rw', default => 'locus_tag'); ## Choose a keytype for merging data
has product_columns => (is => 'rw', default => 'trinity_sprot_Top_BLASTX_hit,inter_Pfam,inter_TIGRFAM'); ## When merging annotations, choose the favorites when upgrading an annotation to 'product'
has product_transmembrane => (is => 'rw', default => 'inter_TMHMM'); ## Column containing transmembrane domain data when upgrading annotations to 'product'
has product_signal => (is => 'rw', default => 'inter_signalp'); ## Column containing signal peptide domain data when upgrading annotations to 'product'
has protocol => (is => 'rw', default => 'Sassetti'); ## TNSeq protocol for TPP
has pval => (is => 'rw', default => undef); ## Pvalue cutoffs.
has qsub_args => (is => 'rw', default => '-j oe -V -m n'); ## What arguments will be passed to qsub by default?
has qsub_depends => (is => 'rw', default => 'depend=afterok:'); ## String to pass for dependencies
has qsub_dependsarray => (is => 'rw', default => 'depend=afterokarray:'); ## String to pass for an array of jobs
has qsub_path => (is => 'rw', default => scalar_which('qsub'));
has qual => (is => 'rw', default => undef); ## cutadapt quality string
has query => (is => 'rw', default => undef); ## Used for searches when input is already taken, most likely blast/fasta
has restart => (is => 'rw', default => 0); ## Restart job(s) in the middle of a group
has riboanchor => (is => 'rw', default => 'start'); ## When correcting, use the start or end position as an anchor
has riboasite => (is => 'rw', default => 1); ## Count riboseq A site positions?
has ribopsite => (is => 'rw', default => 1); ## Count riboseq P site positions?
has riboesite => (is => 'rw', default => 1); ## Count riboseq E site positions?
has ribominsize => (is => 'rw', default => 24); ## Minimum size to search for ribosome positions (riboseq)
has ribomaxsize => (is => 'rw', default => 36); ## Maximum size to search for ribosome positions (riboseq)
has ribominpos => (is => 'rw', default => -30); ## Minimum position for counting (riboseq)
has ribomaxpos => (is => 'rw', default => 30); ## Maximum position for counting (riboseq)
has ribocorrect => (is => 'rw', default => 1); ## Correct ribosome positions for biases
has ribosizes => (is => 'rw', default => '25,26,27,28,29,30,31,32,33,34'); ## Use these sizes for riboseq reads
has runs => (is => 'rw', default => 1000); ## Number of runs for bayesian methods.
has sampleid => (is => 'rw', default => undef); ## Identifier to use for a sample.
has samtools_mapped => (is => 'rw', default => 0); ## Extract mapped reads with samtools.
has samtools_unmapped => (is => 'rw', default => 0); ## Extract unmapped reads with samtools.
has sbatch_depends => (is => 'rw', default => 'afterok:');
has sbatch_dependsarray => (is => 'rw', default => 'afterok:'); ## String to pass for an array of jobs
has sbatch_path => (is => 'rw', default => scalar_which('sbatch'));
has search_string => (is => 'rw', default => 'tail');
has shell => (is => 'rw', default => '/usr/bin/bash'); ## Default qsub shell
has species => (is => 'rw', default => undef); ## Primarily for getting libraries to search against
has sra => (is => 'rw', default => 0);
has starting_tree => (is => 'rw', default => undef); ## Starting tree for phylogenetic analyses
## Note 202212: Now most of the sequencing kits used by our sequencer are reverse.
has stranded => (is => 'rw', default => 'reverse'); ## Did this data come from a stranded library kit?
has suffixes => (is => 'rw', default => '.fastq,.gz,.xz,.fasta,.sam,.bam,.count,.csfasta,.qual,.fsa,.faa,.fna,.gbf,.gbk,.tsv,.csv,.gff,.tbl,.ffn,.sf'); ## Suffixes to remove when invoking basename
has ta_offset => (is => 'rw', default => 0); ## When counting TAs, this is either 0 or 2 depending on if the TA was removed.
has task => (is => 'rw', default => 'tnseq');
has taxid => (is => 'rw', default => '353153'); ## Default taxonomy ID, unknown for now.
has test_file => (is => 'rw', default => 'direct-term-repeasts.fasta'); ## There are a few places where testing for the existence of a test file is useful.
has threshold => (is => 'rw', default => 0.05); ## A second cutoff value (looking at you, glimmer.)
has trim => (is => 'rw', default => 1); ## Perform trimming (rnaseq pipeline, trinity)
has type => (is => 'rw', default => undef); ## Possibly superceded by htseq_type
has varfilter => (is => 'rw', default => 1); ## use a varfilter when performing variant searches.
has verbose => (is => 'rw', default => 0); ## Print extra information while running?
has vcf_cutoff => (is => 'rw', default => 5); ## Minimum depth cutoff for variant searches
has vcf_method => (is => 'rw', default => 'freebayes');
has vcf_minpct => (is => 'rw', default => 0.8); ## Minimum percent agreement for variant searches.
## A few variables which are by definition hash references and such
has slots_ignored => (is => 'ro', default => 'slots_ignored,methods_to_run,menus,term,todos,variable_getvars_args,variable_function_overrides,variable_getopt_overrides,variable_current_state');  ## Ignore these slots when poking at the class.
has methods_to_run => (is => 'rw', default => undef); ## Set of jobs to run.
has menus => (is => 'rw', default => undef); ## The menus when in an interactive session.
has term => (is => 'rw', default => undef); ## A fun Readline terminal for tab completion.
has todos => (is => 'rw', default => undef); ## Set of TODOs to perform.
## These last are used by Get_Vars()
has variable_iterations => (is => 'rw', default => 0); ## After the first iteration of Get_Vars()
## We should clear out getopt_overrides so they don't mess things up.
has variable_getvars_args => (is => 'rw', default => undef); ## One step higher precedence than defaults.
has variable_function_overrides => (is => 'rw', default => undef); ## One step higher than getvars_args.
has variable_getopt_overrides => (is => 'rw', default => undef); ## Highest precedence
has variable_current_state => (is => 'rw', default => undef); ## Current variable state

our $AUTOLOAD;
##our @EXPORT_OK = qw"";
$VERSION = '202308';
$COMPRESSION = 'xz';
$XZ_OPTS = '-9e';
$XZ_DEFAULTS = '-9e';
$ENV{LESS} = '--buffers 0 -B';
$ENV{PERL_USE_UNSAFE_INC} = 0;
my $lessopen = Get_Lesspipe();
## Added to test Bio:DB::SeqFeature::Store
if (!defined($ENV{PERL_INLINE_DIRECTORY})) {
    my $this_tmpdir = '/tmp';
    $this_tmpdir = $ENV{TMPDIR} if (defined($ENV{TMPDIR}));
    my $filename = File::Temp::tempnam($this_tmpdir, "inline_$ENV{USER}");
    $ENV{PERL_INLINE_DIRECTORY} = $filename;
    my $made = make_path($filename);
}

sub scalar_which {
    my $exe = $_[0];
    my $path = which($exe);
    return($path);
}

=head2 C<Help>

    Help() always gives 0.
    Before it returns, it will hopefully print some useful information
    regarding ways to invoke Bio::Adventure.pm.

=cut
sub Help {
    my ($class, $args) = @_;
    my $fh = \*STDOUT;
    ##    my $usage = pod2usage(-output => $fh, -verbose => 99, -sections => "SYNOPSIS");
    use Pod::Find qw"pod_where";
    pod2usage(-input => pod_where({-inc => 1}, __PACKAGE__),
              -verbose => 2, -output => $fh,
              -sections => "NAME|SYNOPSIS|DESCRIPTION|VERSION",
              -exitval => 'NOEXIT');
    return(0);
}

=head2 C<BUILD>

This is the main CYOA constructor.  It takes up the variables at the
top of Adventure.pm along with any arguments pulled in from @ARGV and
incorporates them into the class.

Most of the logic in it is working toward setting up GetOpt::Long
options.  Thus, if there is confusion about what is happening, that is
the first place to look.

=cut
sub BUILD {
    my ($class, $args) = @_;
    ## There are a few default variables which we cannot fill in with MOO defaults.
    ## Make a hash of the defaults in order to make pulling command line arguments easier
    my %defaults;

    ## The modulecmd comand is kind of a hard-prerequisite for this to work.
    my $check = which('modulecmd');
    die("Could not find environment modules in your PATH:
$ENV{PATH}.") unless($check);

    my @ignored;
    if (defined($class->{slots_ignored})) {
        @ignored = split(/\,/, $class->{slots_ignored});
    }
    foreach my $k (keys %{$class}) {
        for my $ignore (@ignored) {
            next if ($k eq $ignore);
        }
        $defaults{$k} = $class->{$k};
    }

    ## Make a log of command line arguments passed.
    if ($#ARGV > 0) {
        my $arg_string = '';
        foreach my $a (@ARGV) {
            $arg_string .= "$a ";
        }
        make_path("outputs/", {verbose => 0}) unless (-r qq"outputs/");
        my $out = FileHandle->new(">>outputs/log.txt");
        my $d = qx'date';
        chomp $d;
        print $out "# Started CYOA at ${d} with arguments: ${arg_string}.\n";
        $out->close();
    }

    ## Now pull an arbitrary set of command line arguments
    ## Options set in the config file are trumped by those placed on the command line.
    my (%override, %conf_specification, %conf_specification_temp);
    foreach my $spec (keys %defaults) {
        ## This line tells Getopt::Long to use as a string argument anything in the set of defaults.
        ## Every option which gets set here will get put into $override->{} thing
        my $def_string = qq"${spec}:s";
        $conf_specification{$def_string} = \$override{$spec};
    }
    ## For some command line options, one might want shortcuts etc, set those here:
    %conf_specification_temp = (
        "de|d" => \$override{debug},
        "hp|h:s" => \$override{hpgl},
        "in|i:s" => \$override{input},
        "sp|s:s" => \$override{species},);
    ## This makes both of the above groups command-line changeable
    foreach my $name (keys %conf_specification_temp) {
        $conf_specification{$name} = $conf_specification_temp{$name};
    }
    undef(%conf_specification_temp);
    my $argv_result = GetOptions(%conf_specification);
    ## Finally move the override hash to $class->{variable_getopt_overrides};
    $class->{variable_getopt_overrides} = \%override;

    ## Take a moment to create a simplified job basename
    ## Eg. if the input file is 'hpgl0523_forward-trimmed.fastq.gz'
    ## Take just hpgl0523 as the job basename
    my $job_basename = '';
    my @suffixes = ('.gz', '.xz', '.bz2');
    if (defined($class->{suffixes})) {
        @suffixes = split(/,/, $class->{suffixes});
    }
    if (defined($class->{input})) {
        $job_basename = $class->{input};
        ## Start by pulling apart any colon/comma separated inputs
        if ($job_basename =~ /:|\,/) {
            my @tmp = split(/:|\,/, $job_basename);
            ## Remove likely extraneous information
            $job_basename = $tmp[0];
        }
        $job_basename = basename($job_basename, @suffixes);
        $job_basename = basename($job_basename, @suffixes);
    } elsif (defined($class->{basedir})) {
        $job_basename = basename($class->{basedir}, @suffixes);
    }
    $job_basename =~ s/_forward.*//g;
    $job_basename =~ s/_R1.*//g;
    $job_basename =~ s/-trimmed.*//g;
    $class->{jbasename} = $job_basename;

    ## Remove backslashes from arbitrary arguments
    if ($class->{arbitrary}) {
        $class->{arbitrary} =~ s/\\//g;
    }

    ## These are both problematic when (n)storing the data.
    if (defined($class->{interactive}) && $class->{interactive}) {
        $class->{menus} = Get_Menus();
    }
    $class->{methods_to_run} = Get_TODOs(%{$class->{variable_getopt_overrides}});
    my $path_agrees = Check_Libpath(libdir => $class->{libdir}, libpath => $class->{libpath});
    $class->{libpath} = $path_agrees->{libpath};
    $class->{libdir} = $path_agrees->{libdir};

    ## Check that the module command is available as a bash function.
    ## I was initially going to do this in BUILD(), but I think that might mess up jobs
    ## which set the jprefix parameter.
    my $modulecmd_handle = IO::Handle->new;
    my $modulecmd_pid = open($modulecmd_handle, qq"bash -i -c 'type module' |");
    my $modulecmd_text = '';
    while (my $line = <$modulecmd_handle>) {
        $modulecmd_text .= $line;
    }
    $modulecmd_handle->close();
    $class->{modulecmd} = $modulecmd_text;
    return($args);
}

=head2 C<Check_Input>

Given a set of inputs, do a little checking to try to ensure that they
are usable.

This function is not used in many places and may need culling.  Either
that or it should be improved a little and propagated.

=cut
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
        my $second_test = basename($file, $class->{suffixes});
        my $third_test = basename($second_test, $class->{suffixes});
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

sub Check_Libpath {
    my %args = @_;
    my $provided_libdir = 0;
    if (defined($args{libdir})) {
        $provided_libdir = 1;
    } else {
        $args{libdir} = '\$HOME/libraries';
    }
    my $provided_libpath = 0;
    if (defined($args{libpath})) {
        $provided_libpath = 1;
    } else {
        $args{libpath} = "$ENV{HOME}/libraries";
    }
    ## Make sure that the libdir and libpath agree with one another.
    ## I am intending to use libdir often(always?) as the shell variable $HOME.
    ## But this may change from host to host and I want to be able to
    ## perform operations on stuff in that directory while still having some
    ## degree of flexibility.
    ## Thus, this will check if libdir starts with a \$, if so, extract that variable
    ## from the environment, and check if the libpath agrees with it.

    my $interpolated_libdir = $args{libdir};
    if ($args{libdir} =~ /\$/) {
        $interpolated_libdir =~ /^(.*)?\$\{*([A-Z]+)\}*(.*)$/;
        my $prefix = $1;
        my $varname = $2;
        my $suffix = $3;
        $interpolated_libdir = qq"${prefix}$ENV{$varname}${suffix}";
    }

    my $ret = {
        libdir => $args{libdir},
        libpath => $interpolated_libdir,
    };

    if ($provided_libdir && $provided_libpath) {
        ## Both were passed to the constructure, check that they agree, if not
        ## Then warn and go with libdir.
    }

    return($ret);
}

=head2 C<Get_Lesspipe>

 Do the equivalent of $(eval lesspipe)

 On my computer, lesspipe returns:
 export LESSOPEN="| /usr/bin/lesspipe %s";
 export LESSCLOSE="/usr/bin/lesspipe %s %s";
 So, I want to split on = and pull the pieces.

=cut
sub Get_Lesspipe {
    my $lesspipe = FileHandle->new("lesspipe |");
    my $result = {};
    while (my $line = <$lesspipe>) {
        chomp $line;
        my ($var, $string) = split(/=/, $line);
        $var =~ s/export\s+//g;
        $string =~ s/"|\;//g;
        $ENV{$var} = $string;
        $result->{$var} = $string;
    }
    $lesspipe->close();
    return($result);
}

=head2 C<Get_Paths>

Given a file/directory name, provide a set of paths which will
hopefully prove useful for the various functions in CYOA.  I think
this function is a good candidate for replacing Check_Input() below.

=cut
sub Get_Paths {
    my ($class, @inputs) = @_;
    my %ret = ();
    if ($inputs[0] =~ /:|\,|\s+/) {
        @inputs = split(/:|\,|\s+/, $inputs[0]);
    }
    my $num_inputs = scalar(@inputs);
    if ($num_inputs == 0) {
        die("This requires an input filename.");
    }
    my @outputs = ();
    for my $in (@inputs) {
        my $filename = basename($in);
        my $filebase_compress = basename($in, ('.gz', '.xz', '.bz2'));
        my @exts = ('.fastq', '.fasta', '.fsa', '.faa', 'fna', '.fa', '.ffn',
                    '.tsv', 'gff', '.gff3', '.gbk', '.gbf', '.sqn', 'tbl');
        my $filebase_extension = basename($filebase_compress, @exts);
        my $directory = dirname($in);
        my $dirname = basename($directory);
        ## Check that doing basename(dirname(input)) leaves us with something interesting
        if ($dirname eq '.') {
            $dirname = undef;
        }
        ## This test/path build is because abs_path only works on stuff which exists.
        if (! -e $directory) {
            make_path($directory);
        }
        my $full_path = abs_path($directory);
        $full_path .= "/${filename}";

        ## Make a new job_basename which doesn't have some cruft
        my $jbasename = $filebase_extension;
        $jbasename =~ s/_forward.*//g;
        $jbasename =~ s/_R1.*//g;
        $jbasename =~ s/-trimmed.*//g;
        $jbasename =~ s/_trimmed.*//g;

        my %ret = (
            filename => $filename,
            filebase_compress => $filebase_compress,
            filebase_extension => $filebase_extension,
            directory => $directory,
            dirname => $dirname,
            fullpath => $full_path,
            jbasename => $jbasename);
        push(@outputs, \%ret);
    }

    return(\@outputs);
}

=head2 C<Get_Term>

Used to get a useful interactive terminal when more options are required.
Ideally, this should get us a terminal with tab completion.

=cut
sub Get_Term {
    my $term = Term::ReadLine->new('>');
    my $attribs = $term->Attribs;
    $attribs->{completion_suppress_append} = 1;
    my $OUT = $term->OUT || \*STDOUT;
    $Term::UI::VERBOSE = 0;
    $term->ornaments(0);
    return($term);
}

=head2 C<Get_Job_Name>

This attempts to make a reasonable job name given a filename or
cwd().  In practice it just strips off the suffixes from the input
file(s) and uses that.

=cut
sub Get_Job_Name {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args);
    my $name = 'unknown';
    $name = $options->{input} if ($options->{input});
    if ($name =~ /\:|\s+|\,/) {
        my @namelst = split(/\:|\s+|\,/, $name);
        $name = $namelst[0];
    }
    $name = basename($name, split(/,/, $class->{suffixes}));
    $name = basename($name, split(/,/, $class->{suffixes}));
    $name =~ s/\-trimmed//g;
    return($name);
}

=head2 C<Get_Input>

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

    }
    $input =~ tr/[A-Z]/[a-z]/;
    return($input);
}

=head2 C<Get_Vars>

  Handle the peculiar mix of instance options held in $class->{options},
  the set of arguments passed to %args, a list of required and potentially
  missing options, and whatever default values one wishes to set.

  my $options = $class->Get_Vars(args => \%args, required => ['input', 'genome'], bob => 'jane');

=cut
sub Get_Vars {
    my ($class, %args) = @_;
    ## I finally figured out why things are inconsistent.
    ## I need to separate the default values from the argv values.
    ## Then the order of operations should be:
    ##  1. Check for the default value
    ##  2. Overwrite with %args provided to Get_Vars() if they exist.
    ##  3. Overwrite with args from argv if they exist.
    ## Currently, I am putting argv into the defaults and then letting them
    ## get overwritten.

    ## First, the defaults, we will set this to the current_state if it is defined.
    my %default_vars;
    ## We will grab out the default variables or current state, depending on what is available.
    if (defined($class->{variable_current_state})) {
        %default_vars = %{$class->{variable_current_state}};
    } else {
        foreach my $k (keys %{$class}) {
            my @ignored = split(/,/, $class->{slots_ignored});
            for my $ignore (@ignored) {
                next if ($k eq $ignore);
            }
            $default_vars{$k} = $class->{$k};
        }
    }

    ## Then the set of variables passed to Get_Vars();
    ## These are effectively function defaults
    my %getvars_default_vars = %args;  ## This guy has jprefix...
    ## The parent function's variables.
    my %function_override_vars = ();
    if (defined($getvars_default_vars{args})) {
        %function_override_vars = %{$getvars_default_vars{args}};
    }
    ## Keep count of how often we come here.
    $class->{variable_iterations}++;
    if ($class->{variable_iterations} > 1) {
        $class->{variable_getopt_overrides} = undef;
    }
    ## If we have come here more than 1 time, then we should stop using the getopt
    ## overrides, because that will totally fubar chains of functioncalls.

    my %getopt_override_vars = ();
    if (defined($class->{variable_getopt_overrides})) {
        %getopt_override_vars = %{$class->{variable_getopt_overrides}};
    }

    ## The returned args will be a hopefully useful combination of the above.
    my %returned_vars = ();

    ## Get the variablenames which are required, these are defined by the 'required' tag
    ## in the argument list to this function.
    my @required_varnames;
    if (defined($getvars_default_vars{required})) {
        @required_varnames = @{$getvars_default_vars{required}};
    }
    ## Check to see if the required options are defined in any of the above...
    foreach my $needed (@required_varnames) {
        my $found = {
            sum => 0,
            default => undef,
            this => undef,
            override => undef,
        };
        if (defined($default_vars{$needed})) {
            $found->{sum}++;
            $found->{default} = $default_vars{$needed};
        }
        ## This one should not be needed, because it implies
        ## that the function caller set both required and the option.
        if (defined($getvars_default_vars{$needed})) {
            print "It appears you defined the variable in required => [] and as an option.\n";
            $found->{sum}++;
            $found->{this} = $getvars_default_vars{$needed};
        }
        if (defined($function_override_vars{$needed})) {
            $found->{sum}++;
            $found->{override} = $function_override_vars{$needed};
        }
        if (defined($getopt_override_vars{$needed})) {
            $found->{sum}++;
            $found->{override} = $getopt_override_vars{$needed};
        }

        ## Now see if sum got incremented...
        if ($found->{sum} == 0) {
            ## Nope, ask for the value:
            my $query = qq"The option: ${needed} was missing, please fill it in:";
            if (!defined($default_vars{term})) {
                ## This may need to go back to $class
                $default_vars{term} = Get_Term();
            }
            my $response = $default_vars{term}->readline($query);
            $response =~ s/\s+$//g;
            $response =~ s/\@|\*|\+//g;
            $getopt_override_vars{$needed} = $response;
        }
    } ## Done looking for required varnames.

    ## Use this for loop to fill in default values from the class.
    for my $varname (keys %default_vars) {
        $returned_vars{$varname} = $default_vars{$varname};
    }

    ## The set of arguments passed to the function includes:
    ## $this_function_vars{args} which contains args to the parent.
    ## Now let us iterate over this_function_vars
    my %defined_in_this_funcall = ();
    for my $varname (keys %getvars_default_vars) {
        ## required and args are provided in the function call, so skip them.
        next if ($varname eq 'required' || $varname eq 'args');
        $returned_vars{$varname} = $getvars_default_vars{$varname};
        $defined_in_this_funcall{$varname} = 1;
    }
    ## Pick up options passed in the parent function call, e.g.
    ## my $bob = $class->Function(bob => 'jane');
    PARENT_ARG: for my $varname (keys %function_override_vars) {
        ## required and args are provided in the function call, so skip them.
        next PARENT_ARG if ($varname eq 'required' || $varname eq 'args');
        ## Do not redefine variables which were defined in the current functioncall.
        if ($defined_in_this_funcall{$varname}) {

        }
        #next PARENT_ARG if ($defined_in_this_funcall{$varname});
        $returned_vars{$varname} = $function_override_vars{$varname};
        ## Try to ensure that the shell reverts to the default 'bash'
        ## I think the following 4 lines are no longer needed.
    }
    ## Final loop to pick up options from the commandline or a TERM prompt.
    ## These supercede everything else.
    for my $varname (keys %getopt_override_vars) {
        next if ($varname eq 'required' || $varname eq 'args');
        next unless (defined($getopt_override_vars{$varname}));
        $returned_vars{$varname} = $getopt_override_vars{$varname};
    }

    ## I previously had the cluster check in the BUILD() function.
    ## I think this is the wrong place and should instead be here?
    my $torque_test = which('qsub');
    my $slurm_test = which('sbatch');
    if (defined($returned_vars{cluster})) {
        if ($returned_vars{cluster} ne 'bash') {
            $returned_vars{qsub_path} = $torque_test;
            $returned_vars{sbatch_path} = $slurm_test;
        } else {
            $returned_vars{sbatch_path} = '';
            $returned_vars{qsub_path} = '';
            $returned_vars{cluster} = 'bash';
        }
    } else {
        if ($slurm_test) {
            $returned_vars{cluster} = 'slurm';
            $returned_vars{sbatch_path} = $slurm_test;
        } elsif ($torque_test) {
            $returned_vars{cluster} = 'torque';
            $returned_vars{qsub_path} = $torque_test;
        } else {
            $returned_vars{cluster} = 'bash';
        }
    }

    ## Final sanity check(s)
    for my $varname (keys %returned_vars) {
        next unless (defined($returned_vars{$varname}));
        $returned_vars{$varname} =~ s/^\~/${HOME}/g;
        ## Finally, fill in a few special cases here.
        if ($varname eq 'jname') {
            if ($returned_vars{jname} eq '') {
                my $name = 'unknown';
                $name = $returned_vars{input} if ($returned_vars{input});
                if ($name =~ /\:|\s+|\,/) {
                    my @namelst = split(/\:|\s+|\,/, $name);
                    $name = $namelst[0];
                }
                $name = basename($name, (".gz", ".xz", ".bz2", ".bai", ".fai"));
                $name = basename($name, (".fasta", ".fastq", ".bam", ".sam", ".count"));
                $name =~ s/_forward//g;
                $returned_vars{jname} = $name;
            }
        } ## End checking on job name
    } ## End final iteration over the options keeps.
    ## End special cases.

    ## So at this point we should have the full set of defaults + function + overrides.
    ## However, we need to ensure that child functions get this information, but
    ## unfortunately we are currently calling those functions with:
    ## Bio::Adventure::Something::Something($class, %args);

    ## After the global defaults, these are the lowest priority
    $class->{variable_getvars_args} = \%getvars_default_vars;
    ## Then the function overrides
    $class->{variable_function_overrides} = \%function_override_vars;
    ## And last, the getopt overrides
    $class->{variable_getopt_overrides} = \%getopt_override_vars;
    $class->{variable_current_state} = \%returned_vars;
    return(\%returned_vars);
}

sub Load_Vars {
    my ($class, %args) = @_;
    my $num_changed = 0;
    for my $varname (keys %args) {
        if ($varname eq 'input') {
            local $Storable::Eval = 1;
            my $input_options = lock_retrieve($args{$varname});
            for my $opt (keys %{$input_options}) {
                $num_changed ++;
                $class->{$opt} = $input_options->{$opt};
            }
        } else {
            $num_changed++;
            $class->{$varname} = $args{$varname};
        }
    }
    return($num_changed);
}

=head2 C<Reset_Vars>

  Make sure that the environment is returned to a useful default state when perturbed.

=cut
sub Reset_Vars {
    my ($class, %args) = @_;
    my %original_values = (
        array_string => undef,
        jbasename => basename(cwd()),
        jdepends => '',
        jmem => 12,
        jname => 'undefined',
        jprefix => '01',
        jstring => '',
        jwalltime => '10:00:00',
        language => 'bash',
        prescript => '',
        postscript => '',
        modules => [],
        shell => '/usr/bin/env bash',);
    for my $k (keys %original_values) {
        $class->{$k} = $original_values{$k};
        $class->{variable_current_state}->{$k} = $original_values{$k};
        $class->{variable_function_overrides}->{$k} = $original_values{$k};
        $class->{variable_getvars_args}->{$k} = $original_values{$k};
    }
    return($class);
}


=head2 C<Set_Vars>

  Handle the peculiar mix of instance options held in $class->{options},
  the set of arguments passed to %args, a list of required and potentially
  missing options, and whatever default values one wishes to set.

  $options = $class->Set_Vars(exclude => 'bob');

This function is likely no longer needed because of the way
Get_Options now handles the options hash.

=cut
sub Set_Vars {
    my ($class, %args) = @_;
    my $options;
    if (defined($args{options})) {
        $options = $args{options};
    } else {
        $options = $class->{options};
    }
    my $ref = ref($options);
    foreach my $k (keys %args) {
        $options->{$k} = $args{$k};
    }
    $class->{options} = $options;
    return($options);
}

sub Last_Stat {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args, required => ['input']);
    my $input_filename = $options->{input};
    my $input = FileHandle->new("<${input_filename}");
    my ($line, $last);
    while ($line = <$input>) {
        chomp $line;
        ## I am doing this due to terminal newlines on the file, there is probably a much smrtr way.
        $last = $line unless ($line =~ /^$/);
    }
    $input->close();
    return($last);
}

=head2 C<Module_Loader>

This loads environment-modules into the current ENV.

=cut
sub Module_Loader {
    my ($class, %args) = @_;
    my $action = 'load';
    if ($args{action}) {
        $action = $args{action};
    }

    ## Start with the current ENV and an empty set of things to load and test.
    my %old_env = %ENV;
    my @mod_lst = ();
    my $mod_class = ref($args{modules});
    my @test_lst = ();

    my $test_class = ref($args{exe});
    if ($mod_class eq 'SCALAR') {
        push(@mod_lst, $args{modules});
    } elsif ($mod_class eq 'ARRAY') {
        @mod_lst = @{$args{modules}};
    } elsif ($mod_class eq 'HASH') {
        my %mods = %{$args{modules}};
        @mod_lst = keys %mods;
        for my $k (@mod_lst) {
            push(@test_lst, $mods{$k});
        }
    } elsif (!$mod_class) {
        push(@mod_lst, $args{modules});
    } else {
        print "I do not know this class: $mod_class, $args{modules}\n";
    }

    if ($test_class eq 'SCALAR') {
        push(@test_lst, $args{exe});
    } elsif ($test_class eq 'ARRAY') {
        for my $e (@{$args{exe}}) {
            push(@test_lst, $e);
        }
    } elsif (!$test_class) {
        push(@test_lst, $args{exe});
    } else {
        print "I do not know this class: $test_class, $args{exe}\n";
    }

    my $count = 0;
    my $mod;
    if ($action eq 'unload') {
        @mod_lst = reverse(@mod_lst);
    }
    my @messages = ();
    my @loads = ();
    LOADER: for $mod (@mod_lst) {
        ## Skip one's self.
        my ($stdout, $stderr, @result) = capture {
            try {
                my $test;
                push(@loads, $mod);
                if ($action eq 'load') {
                    $test = Env::Modulecmd::load($mod);
                } else {
                    $test = Env::Modulecmd::unload($mod);
                }
            } catch ($e) {
                ## print "There was an error loading ${mod}, ${e}\n";
                my $message = qq"There was an error loading ${mod}: ${e}";
                push(@messages, $message);
            };
        };
        $count++;
    } ## End iterating over every mod in the list

    ## If we got a set of test programs, check that they are in the PATH:
    if (scalar(@test_lst) > 0) {
      EXELOOP: for my $exe (@test_lst) {
          next EXELOOP unless($exe);
          my $check = which($exe);
          die(qq"Could not find ${exe} in the PATH.") unless ($check);
      }
    }
    $class->{modules} = \@mod_lst;
    return(%old_env);
}

sub Module_Reset {
    my ($class, %args) = @_;
    my %old_env = %args{env};
    for my $k (keys %old_env) {
        $ENV{$k} = $old_env{$k};
    }
}

=head2 C<Passthrough_Args>

Add some flexibility in how arbitrary arguments get passed through to
various tools invoked by cyoa.  Using this function you may have any
number of extra args passed through to the final script by separating
them.  If it sees a single letter, it will prefix a single dash.  If
there are multiple letters it will double-dash it.  Conversely, if
the arbitrary string starts with a dash, it will leave it alone.

Thus:

--arbitrary ':--funkytown=bob;L 10,very-sensitive'

Will get the following arguments passed to your downstream tool:

'--funkytown=bob -L 10 --very-sensitive'.

=cut
sub Passthrough_Args {
    my ($class, %args) = @_;
    my $argstring = $args{arbitrary};
    my $new_string = '';
    for my $arg (split /\,|:|\;/, $argstring) {
        if ($arg =~ /^\w{1}$|^\w{1}\W+/) {
            ## print "Single letter arg passthrough\n";
            $arg = qq" -${arg} ";
        } elsif ($arg =~ /^\w{2}/) {
            ## print "Multi letter arg passthrough\n";
            $arg = qq" --${arg} ";
        } elsif ($arg =~ /^\-/) {
            $arg= qq" ${arg} ";
        }
        $new_string .= $arg;
    }
    return($new_string);
}

=head2 C<Read_Genome_Fasta>

Read a fasta file and return the chromosomes.

=cut
sub Read_Genome_Fasta {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        extra => 0);
    my $chromosomes = {};
    ## FIXME: This should be removed and I need to
    ## ensure that all calls use either fasta or genome.
    if (!defined($options->{genome})) {
        if (defined($options->{fasta})) {
            $options->{genome} = $options->{fasta};
        }
    }
    my $fasta_name = basename($options->{genome}, ['.fasta']);
    my $genome_file = qq"$options->{basedir}/${fasta_name}.pdata";
    ##if (-r $genome_file) {
    ##    $chromosomes = lock_retrieve($genome_file);
    ##} else {
    my $input_genome = Bio::SeqIO->new(-file => $options->{genome}, -format => 'Fasta');
    while (my $genome_seq = $input_genome->next_seq()) {
        next unless(defined($genome_seq->id));
        my $id = $genome_seq->id;
        my $sequence = $genome_seq->seq;
        my $length = $genome_seq->length;
        my $empty_forward = [];
        my $empty_reverse = [];
        if ($options->{extra}) {
            for my $c (0 .. $length) {
                $empty_forward->[$c] = 0;
                $empty_reverse->[$c] = 0;
            }
        }
        $chromosomes->{$id}->{forward} = $empty_forward;
        $chromosomes->{$id}->{sequence} = $sequence;
        $chromosomes->{$id}->{reverse} = $empty_reverse;
        $chromosomes->{$id}->{obj} = $genome_seq;
    } ## End reading each chromosome
    ## my $stored = lock_store($chromosomes, $genome_file);
    ##} ## End checking for a .pdata file
    return($chromosomes);
}

=head2 C<Read_Genome_GFF>

Read a GFF file and extract the annotation information from it.

=cut
sub Read_Genome_GFF {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['gff'],
        gff_type => 'gene',
        gff_tag => 'locus_tag');
    my $annotation_in = Bio::Tools::GFF->new(-file => "$options->{gff}", -gff_version => 3);
    my $gff_out = {
        stats => {
            chromosomes => [],
            lengths => [],
            feature_names => [],
            cds_starts => [],
            cds_ends => [],
            inter_starts => [],
            inter_ends => [],
        },
    };
    my $gff_name = basename($options->{gff}, ['.gff']);
    ##my $genome_file = qq"$options->{basedir}/${gff_name}.pdata";
    print "Reading $options->{gff}, seeking features tagged (last column ID) $options->{gff_tag},
 type (3rd column): $options->{gff_type}.\n";
    ##if (-r $genome_file && !$options->{debug}) {
    ##    $gff_out = lock_retrieve($genome_file);
    ##} else {
    print "Starting to read gff: $options->{gff}\n" if ($class->{debug});
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
      next LOOP unless ($feature->primary_tag eq $options->{gff_type});
      $hits++;
      my $location = $feature->location;
      $old_start = $start;
      $old_end = $end;
      $start = $feature->start();
      $end = $feature->end();
      my $strand = $location->strand();
      my @ids = $feature->each_tag_value($options->{gff_tag});
      my $id = "";
      my $gff_chr = $feature->{_gsf_seq_id};
      ## my $gff_chr = $feature->display_name;
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
      push(@inter_starts, $old_start + 1);
      push(@cds_ends, $end);
      ##push(@inter_ends, $old_start-1);
      push(@inter_ends, $start - 1);
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
    print STDERR "Not many hits were observed, do you have the right feature type?  It is: $options->{gff_type}\n" if ($hits < 1000);
    ##if (-f $genome_file) {
    ##        unlink($genome_file);
    ##    }
    ##    my $stored = lock_store($gff_out, $genome_file);
## } ## End looking for the gff data file
    return($gff_out);
}

=head2 C<Reorder_Fasta>

Some tools which pass back a set of chromosomes do not always return
them in the most sensible order.  This attempts to fix that, and make
them prettier.

=cut
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

=head2 C<Submit>

Pass a job to slurm or torque after doing a little sanity checking.  I
am increasinbly certain that keeping a copy of the options in a .pdata
file is not needed.

=cut
sub Submit {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args);
    my %modules = Get_Modules(caller => 2);
    my $loaded = $class->Module_Loader(%modules);
    my $modulecmd_check = $class->{modulecmd};
    make_path("$options->{logdir}", {verbose => 0}) unless (-r "$options->{logdir}");
    my $module_string = '';
    my @module_lst = ();
    @module_lst = @{$modules{modules}} if (defined($modules{modules}));
    if ($options->{language} eq 'perl') {
        push(@module_lst, 'cyoa');
    }
    if (scalar(@module_lst) > 0) {
        $module_string = 'mod=$( { type -t module || true; } )
if [[ -z "${mod}" ]]; then
  module() {
    eval $(/usr/bin/modulecmd bash $*);
  }
  export -f module
fi
module purge
module add ';
        for my $m (@module_lst) {
            $module_string .= qq" ${m}" if (defined($m));
        }
        $module_string .= ' 2>/dev/null 1>&2';
        $options->{module_string} = $module_string;
    }

    ## If we are invoking an indirect job, we need a way to serialize the options
    ## in order to get them passed to the eventual interpreter
    my $option_file = '';
    if ($options->{language} eq 'perl') {
        my $option_directory;
        if (defined($options->{output_dir})) {
            $option_directory = $options->{output_dir}
        } elsif (defined($options->{input})) {
            $option_directory = dirname($options->{input});
        } else {
            $option_directory = $options->{basedir}
        }
        if (!-d $option_directory) {
            my $created = make_path($option_directory);
        }
        ## I think this might be required as per:
        ## https://metacpan.org/pod/release/AMS/Storable-2.21/Storable.pm#CODE_REFERENCES
        $Storable::Deparse = 1;
        $Storable::Eval = 1;
        if ($options->{pdata}) {
            my $jname = 'unknown';
            $jname = $options->{jname} if (defined($options->{jname}));
            $option_file = File::Temp->new(
                TEMPLATE => qq"${jname}XXXX",
                DIR => $option_directory,
                UNLINK => 0,
                SUFFIX => '.pdata',);
            my $option_filename = $option_file->filename;
            $options->{option_file} = $option_filename;
            $class->{option_file} = $option_filename;
            ## Code references are invalid for these things...
            ## Why is it that periodically I get this error?
            ## The result of B::Deparse::coderef2text was empty - maybe you're trying to serialize an XS function?
            my %saved_options = ();
          SAVED: foreach my $k (keys %{$options}) {
              my $r = ref($options->{$k});
              next SAVED if ($r eq 'ARRAY' || $r eq 'HASH' || $r eq 'GLOB');
              $saved_options{$k} = $options->{$k};
          }
            try {
                my $stored = lock_store(\%saved_options, $option_file);
          } catch ($e) {
              warn "An error occurred when storing the options: $e";
              print "HERE ARE THE OPTIONS:\n";
              print Dumper \%saved_options;
          }
        }
    }
    my $runner;
    if ($options->{cluster} eq 'slurm') {
        $runner = Bio::Adventure::Slurm->new();
    } elsif ($options->{cluster} eq 'torque') {
        $runner = Bio::Adventure::Torque->new();
    } elsif ($options->{cluster} eq 'bash') {
        ## I should probably have something to handle gracefully bash jobs.
        $runner = Bio::Adventure::Local->new();
    } elsif ($options->{cluster} eq '0') {
        ## On occasion I set cluster to 0 which is bash.
        $runner = Bio::Adventure::Local->new();
    } else {
        carp("Could not find sbatch, qsub, nor bash.");
        print "Assuming this is running on a local shell.\n";
        $runner = Bio::Adventure::Local->new();
    }

    ## Add the current options to the runner:
    for my $k (keys %{$options}) {
        $runner->{$k} = $options->{$k};
    }
    my $result = $runner->Submit($class, %args);
    my $unloaded = $class->Module_Reset(env => $loaded);
    $class = $class->Reset_Vars();
    return($result);
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
