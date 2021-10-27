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
use Storable qw"lock_store lock_retrieve";
use Term::ReadLine;
use Term::UI;

## Some things which are not currently in use:
## use Archive::Extract; ## When I was opening tarballs
## use Bio::DB::Sam;  ## For reading bamfiles
## use JSON -convert_blessed_universally;  ## Intended for use instead of Storable
## use Data::Dumper;  ## Often used for testing data structures

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
use Bio::Adventure::Phage;
use Bio::Adventure::Phylogeny;
use Bio::Adventure::Pipeline;
use Bio::Adventure::PopGen;
use Bio::Adventure::Prepare;
use Bio::Adventure::Resistance;
use Bio::Adventure::Riboseq;
use Bio::Adventure::QA;
use Bio::Adventure::SeqMisc;
use Bio::Adventure::Slurm;
use Bio::Adventure::SNP;
use Bio::Adventure::TNSeq;
use Bio::Adventure::Torque;
use Bio::Adventure::Trim;
use Bio::Adventure::Visualization;

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

  The options for cyoa are all defined in Adventure.pm, but usually have preferred domains,
  the most likely purpose for each parameter follows:

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

  cluster - Use a computing cluster or just invoke scripts one at a time?
  jdepends - A string describing the set of dependencies for a job.
  jname - A machine-readable name for each job.
  jprefix - A prefix number for each job, to make reading successions easier.
  jstring - The string comprising the actual job to be run.
  job_id - The job ID returned by Torque/Slurm for a job.
  jmem - Memory in Gb to request from the cluster.
  jwalltime - Amount of walltime to request from the cluster. (10:00:00)
  queue_array_string - When running an array of jobs, this holds the requisite formatted string.
  language - Choose if the cluster options should be appropriate for a shell or perl script. (bash)
  qsub_args - Default arguments when using torque. (-j oe -V -m n)
  sbatch_depends - Dependency string on a slurm cluster. (afterok:)
  sbatch_dependsarray - Dependency string on a slurm for array jobs. (afterok:)
  qsub_depends - Dependency string on a torque cluster. (depend=afterok:)
  qsub_dependsarray - Dependency string on a torque cluster for array jobs. (depend=afterokarray:)
  queue - Default qos/queue to send jobs. (workstation)
  queues - Set of available queues. (throughput, workstation, long, large)
  shell - Default shell to request when using torque/slurm. (/usr/bin/bash)

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

  htseq_args - Some default arguments for htseq. (--order=name --idattr=gene_id --minaqual=10
     --type exon --stranded=yes --mode=union)
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
  species - When converting from sam->bam, use this to pull the appropriate genome fasta file.
  paired - Explicitly set whether a sam file is from paired data and therefore filter unpaired reads.

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

## Lets move all the default values here.
has align_blast_format => (is => 'rw', default => 5); ## Which alignment type should we use? (5 is blastxml)
has align_jobs => (is => 'rw', default => 40); ## How many blast/fasta alignment jobs should we make when splitting alignments across nodes?
has align_parse => (is => 'rw', default => 1); ## Parse blast searches into a table?
has arbitrary => (is => 'rw', default => undef); ## Extra arbitrary arguments to pass
has bamfile => (is => 'rw', default => undef); ## Default bam file for converting/reading/etc.
has basedir => (is => 'rw', default => cwd());  ## This was cwd() but I think that may cause problems.
has bash_path => (is => 'rw', default => My_Which('bash'));
has best_only => (is => 'rw', default => 0); ## keep only the best search result when performing alignments?
has blast_params => (is => 'rw', default => ' -e 10 '); ## Default blast parameters
has blast_tool => (is => 'rw', default => 'blastn'); ## Default blast tool to use
has bt_default => (is => 'rw', default => '--best'); ## Default bt1 arguments.
has bt_varg => (is => 'rw', default => '-v 0');
has bt_marg => (is => 'rw', default => '-M 0');
has bt_larg => (is => 'rw', default => '-y -l 15');
has bt2_args => (is => 'rw', default => ' --very-sensitive -L 14 '); ## My favorite bowtie2 arguments
has btmulti => (is => 'rw', default => 0); ## Perform multiple bowtie searches?
has classifier_input => (is => 'rw', default => 'outputs/18classifier/ictv_filtered.tsv'); ## Similar taxa detected by tblastx
has cluster => (is => 'rw', default => undef); ## Are we running on a cluster?
has comment => (is => 'rw', default => undef); ## Set a comment in running slurm/bash/etc scripts.
has config => (is => 'rw', default => undef); ## Not sure
has coverage => (is => 'rw', default => undef); ## Provide a coverage cutoff
has cpus => (is => 'rw', default => 4); ## Number of processors to request in jobs
has csv_file => (is => 'rw', default => 'all_samples.csv'); ## Default csv file to read/write.
has directories => (is => 'rw', default => undef); ## Apply a command to multiple input directories.
has evalue => (is => 'rw', default => 0.001); ## Default e-value cutoff
has fasta_args => (is => 'rw', default => ' -b 20 -d 20 '); ## Default arguments for the fasta36 suite
has fasta_tool => (is => 'rw', default => 'ggsearch36'); ## Which fasta36 program to run?
has filtered => (is => 'rw', default => 'unfiltered');  ## Whether or not Fastqc is running on filtered data.
has fsa_input => (is => 'rw'); ## fsa genome output file for creating a genbank file
has gcode => (is => 'rw', default => '11'); ## Choose a genetic code
has genbank_input => (is => 'rw', default => undef); ## Existing genbank file for merging annotations.
has genome => (is => 'rw', default => undef); ## Choose a genome to work on.
has genus => (is => 'rw', default => undef); ## Choose a genus when using prokka and potentially others like kraken
has gff => (is => 'rw', default => undef); ## Feature file to read/write
has gff_tag => (is => 'rw', default => 'gene_id'); ## Likely redundant with htseq_id
## Ahh I remember, htseq_type was added to facilitate performing multiple htseq counts on multiple gff files.
## Notably the interCDS for bacterial genomes.
has gff_type => (is => 'rw', default => ''); ## When blank, do it on the whole gff file, otherwise use that suffix.
has help => (is => 'rw', default => undef); ## Ask for help?
has htseq_args => (is => 'rw', default => ' --order=name --idattr=gene_id --minaqual=10 --type=exon --stranded=yes --mode=union '); ## Most likely htseq options
has htseq_id => (is => 'rw', default => 'ID'); ## Default htseq ID tag
has htseq_stranded => (is => 'rw', default => 'no'); ## Use htseq stranded options?
has htseq_type => (is => 'rw', default => 'gene'); ## What type of gff file should be read by htseq?
has identity => (is => 'rw', default => 70); ## Alignment specific identity cutoff
has index_file => (is => 'rw', default => 'indexes.txt'); ## File containing indexes:sampleIDs when demultiplexing samples - likely tnseq
has index_hash => (is => 'rw', default => undef);
has input => (is => 'rw', default => undef); ## Generic input argument
has input_abricate => (is => 'rw', default => 'outputs/12abricate_10prokka_09termreorder_08phageterm_07watson_plus/abricate_combined.tsv'); ## Used when merging annotation files into a xlsx/tbl/gbk file.
has interactive => (is => 'rw', default => 0); ## Is this an interactive session?
has interpro_input => (is => 'rw', default => 'outputs/13_interproscan_10prokka_09termreorder_08phageterm_07watson_plus/interproscan.tsv'); ## interpro output file when merging annotations.
has jobs => (is => 'rw', default => undef); ## List of currently active jobs, possibly not used right now.
has jobids => (is => 'rw', default => undef); ## A place to put running jobids, maybe no longer needed.
has jbasename => (is => 'rw', default => undef); ## Job basename
has jdepends => (is => 'rw', default => undef);  ## Flag to start a dependency chain
has jmem => (is => 'rw', default => 12); ## Number of gigs of ram to request
has jname => (is => 'rw', default => undef); ## Job name on the cluster
has jnice => (is => 'rw', default => 10); ## Set the niceness of a job, if it starts positive, we can set a lower nice to preempt
has jpartition => (is => 'rw', default => 'dpart');
has jprefix => (is => 'rw', default => undef); ## Prefix number for the job
has jqueue => (is => 'rw', default => 'workstation'); ## What queue will jobs default to?
has jqueues => (is => 'rw', default => 'throughput,workstation,long,large'); ## Other possible queues
has jsleep => (is => 'rw', default => '0.5'); ## Set a sleep between jobs
has jstring => (is => 'rw', default => undef); ## String of the job
has jwalltime => (is => 'rw', default => '10:00:00'); ## Default time to request
has kingdom => (is => 'rw', default => undef); ## Taxonomic kingdom, prokka/kraken
has language => (is => 'rw', default => 'bash'); ## What kind of script is this?
has length => (is => 'rw', default => 17); ## kmer length, other stuff too.
has libdir => (is => 'rw', default => "$ENV{HOME}/libraries"); ## Directory containing genomes/gff/indexes
has library => (is => 'rw', default => undef);  ## The library to be used for fasta36/blast searches
has libtype => (is => 'rw', default => 'genome'); ## Type of sequence to map against, genomic/rRNA/contaminants
has locus_tag => (is => 'rw', default => undef); ## Used by prokka to define gene prefixes
has logdir => (is => 'rw', default => 'outputs/logs'); ## place to dump logs
has loghost => (is => 'rw', default => 'localhost'); ## Host to which to send logs
has mapper => (is => 'rw', default => 'hisat'); ## Use this aligner if none was chosen.
has mature_fasta => (is => 'rw', default => undef); ## Database of mature miRNA sequences to search
has maxlength => (is => 'rw', default => 42); ## Maximum sequence length when trimming
has method => (is => 'rw', default => undef);
has mi_genome => (is => 'rw', default => undef); ## Set a miRbase genome to hunt for mature miRNAs
has minlength => (is => 'rw', default => 8); ## Minimum length when trimming
has mirbase_data => (is => 'rw', default => undef); ## miRbase annotation dataset.
has modules => (is => 'rw', default => undef); ## Environment modules to load
has option_file => (is => 'rw', default => undef);
has orientation => (is => 'rw', default => 'start'); ## Default orientation when normalizing riboseq reads
has outgroup => (is => 'rw', default => undef); ## Outgroup for phylogenetic tools
has output => (is => 'rw', default => undef); ## Generic output argument
has outdir => (is => 'rw', default => undef);
has paired => (is => 'rw', default => 0); ## Is the input paired?
has phageterm_input => (is => 'rw', default => 'outputs/08phageterm_07watson_plus/direct-term-repeats.fasta'); ## phageterm output file when merging annotations.
has phred => (is => 'rw', default => 33); ## Minimum quality score when trimming
has postscript => (is => 'rw', default => undef); ## String to put after a cluter job.
has prescript => (is => 'rw', default => undef); ## String to put before a cluster job.
has primary_key => (is => 'rw', default => 'locus_tag'); ## Choose a keytype for merging data
has product_columns => (is => 'rw', default => 'trinity_sprot_Top_BLASTX_hit,inter_Pfam,inter_TIGRFAM'); ## When merging annotations, choose the favorites when upgrading an annotation to 'product'
has product_transmembrane => (is => 'rw', default => 'inter_TMHMM'); ## Column containing transmembrane domain data when upgrading annotations to 'product'
has product_signal => (is => 'rw', default => 'inter_signalp'); ## Column containing signal peptide domain data when upgrading annotations to 'product'
has prokka_tsv_input => (is => 'rw', default => undef); ## Prokka tsv file for merging annotations.
has protocol => (is => 'rw', default => 'Sassetti'); ## TNSeq protocol for TPP
has pval => (is => 'rw', default => undef); ## Pvalue cutoffs.
has qsub_args => (is => 'rw', default => '-j oe -V -m n'); ## What arguments will be passed to qsub by default?
has qsub_depends => (is => 'rw', default => 'depend=afterok:'); ## String to pass for dependencies
has qsub_dependsarray => (is => 'rw', default => 'depend=afterokarray:'); ## String to pass for an array of jobs
has qsub_path => (is => 'rw', default => My_Which('qsub'));
has qual => (is => 'rw', default => undef); ## cutadapt quality string
has query => (is => 'rw', default => undef); ## Used for searches when input is already taken, most likely blast/fasta
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
has sbatch_depends => (is => 'rw', default => 'afterok:');
has sbatch_dependsarray => (is => 'rw', default => 'afterok:'); ## String to pass for an array of jobs
has sbatch_path => (is => 'rw', default => My_Which('sbatch'));
has shell => (is => 'rw', default => '/usr/bin/bash'); ## Default qsub shell
has species => (is => 'rw', default => undef); ## Primarily for getting libraries to search against
has starting_tree => (is => 'rw', default => undef); ## Starting tree for phylogenetic analyses
has stranded => (is => 'rw', default => 0); ## Did this data come from a stranded library kit?
has suffixes => (is => 'rw', default => '.fastq,.gz,.xz,.fasta,.sam,.bam,.count,.csfasta,.qual'); ## Suffixes to remove when invoking basename
has ta_offset => (is => 'rw', default => 0); ## When counting TAs, this is either 0 or 2 depending on if the TA was removed.
has task => (is => 'rw', default => undef);
has taxid => (is => 'rw', default => '353153'); ## Default taxonomy ID, unknown for now.
has test_file => (is => 'rw', default => 'direct-term-repeasts.fasta'); ## There are a few places where testing for the existence of a test file is useful.
has trinotate_input => (is => 'rw', default => '11trinotate_10prokka_09termreorder_08phageterm_07watson_plus/Trinotate.tsv'); ## trinotate output, used when merging annotations.
has type => (is => 'rw', default => undef); ## Possibly superceded by htseq_type
has varfilter => (is => 'rw', default => 1); ## use a varfilter when performing variant searches.
has verbose => (is => 'rw', default => 0); ## Print extra information while running?
has vcf_cutoff => (is => 'rw', default => 10); ## Minimum depth cutoff for variant searches
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
$VERSION = '20151101';
$COMPRESSION = 'xz';
$XZ_OPTS = '-9e';
$XZ_DEFAULTS = '-9e';
$ENV{LESSOPEN} = '| lesspipe %s';
$ENV{LESS} = '--buffers 0';

=over

=item C<Help>

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

sub BUILD {
    my ($class, $args) = @_;
    ## There are a few default variables which we cannot fill in with MOO defaults.
    ## Make a hash of the defaults in order to make pulling command line arguments easier
    my %defaults;
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

    if (defined($class->{cluster})) {
        if (!$class->{cluster}) {
            $class->{sbatch_path} = 0;
            $class->{qsub_path} = 0;
            $class->{cluster} = 'bash';
        }
    }
    ## Figure out what kind of cluster we are using, if any.
    my $queue_test = My_Which('sbatch');
    if ($queue_test) {
        $class->{cluster} = 'slurm';
    } else {
        $queue_test = My_Which('qsub');
        if ($queue_test) {
            $class->{cluster} = 'torque';
        }
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
        print $out "# Started CYOA at ${d} with arguments: $arg_string.\n";
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
        "pb|p:i" => \$override{pbs},
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
    if (defined($class->interactive) && $class->interactive) {
        $class->{menus} = Get_Menus();
    }
    $class->{methods_to_run} = Get_TODOs(%{$class->{variable_getopt_overrides}});
    return $args;
}


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

        my %ret = (
            filename => $filename,
            filebase_compress => $filebase_compress,
            filebase_extension => $filebase_extension,
            directory => $directory,
            dirname => $dirname,
            fullpath => $full_path,);
        push(@outputs, \%ret);
    }

    my $final = \@outputs;
    ## When there is just one file, just return its information.
    if (scalar(@outputs) == 1) {
        $final = $outputs[0];
    }
    return($final);
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
                '(abricate): Search for Resistance genes across databases.' => \&Bio::Adventure::Resistance::Abricate,
                '(aragorn): Search for tRNAs with aragorn.' => \&Bio::Adventure::Annotation::Aragorn,
                '(classify_virus): Use ICTV data to classify viral sequences/contigs.' => \&Bio::Adventure::Phage::Classify_Phage,
                '(extend_kraken): Extend a kraken2 database with some new sequences.' => \&Bio::Adventure::Annotation::Extend_Kraken_DB,
                '(glimmer): Use glimmer to search for ORFs.' => \&Bio::Adventure::Annotation::Glimmer,
                '(interproscan): Use interproscan to analyze ORFs.' => \&Bio::Adventure::Annotation::Interproscan,
                '(kraken2): Taxonomically classify reads.' => \&Bio::Adventure::Annotation::Kraken,
                '(phageterm): Invoke phageterm to hunt for likely phage ends.' => \&Bio::Adventure::Phage::Phageterm,
                '(phastaf): Invoke phastaf to attempt classifying phage sequence.' => \&Bio::Adventure::Phage::Phastaf,
                '(terminasereorder): Reorder an assembly based on the results of a blast search.' => \&Bio::Adventure::Phage::Terminase_ORF_Reorder,
                '(prodigal): Run prodigal on an assembly.' => \&Bio::Adventure::Annotation::Prodigal,
                '(prokka): Invoke prokka to annotate a genome.' => \&Bio::Adventure::Annotation::Prokka,
                '(resfinder): Search for antimicrobial resistance genes.' => \&Bio::Adventure::Resistance::Resfinder,
                '(rgi): Search for resistance genes with genecards and peptide fasta input.' => \&Bio::Adventure::Resistance::Rgi,
                '(trnascan): Search for tRNA genes with trnascan.' => \&Bio::Adventure::Annotation::tRNAScan,
                '(trainprodigal): Train prodgial using sequences from a species/strain.' => \&Bio::Adventure::Annotation::Train_Prodigal,
                '(mergeannotations): Merge annotations into a genbank file.' => \&Bio::Adventure::Annotation::Merge_Annotations_Make_Gbk,
            },
        },
        Assembly => {
            name => 'assembly',
            message => 'The wise man fears the wrath of a gentle heart. Go to page 314159.',
            choices => {
                '(abyss): Run abyss to create a new assembly.' => \&Bio::Adventure::Assembly::Abyss,
                '(extract_trinotate): Extract the most likely hits from Trinotate.' => \&Bio::Adventure::Annotation::Extract_Trinotate,
                '(extend_kraken): Extend a kraken2 database with some new sequences.' => \&Bio::Adventure::Annotation::Extend_Kraken_DB,
                '(filterdepth): Filter contigs based on sequencing depth.' => \&Bio::Adventure::Assembly::Filter_Depth,
                '(kraken2): Taxonomically classify reads.' => \&Bio::Adventure::Annotation::Kraken,
                '(transdecoder):  Run transdecoder on a putative transcriptome.' => \&Bio::Adventure::Assembly::Transdecoder,
                '(trinotate): Perform de novo transcriptome annotation with trinotate.' => \&Bio::Adventure::Annotation::Trinotate,
                '(trinity): Perform de novo transcriptome assembly with trinity.' => \&Bio::Adventure::Assembly::Trinity,
                '(trinitypost): Perform post assembly analyses with trinity.' => \&Bio::Adventure::Assembly::Trinity_Post,
                '(velvet): Perform de novo genome assembly with velvet.' => \&Bio::Adventure::Assembly::Velvet,
                '(unicycler): Perform de novo assembly with unicycler.' => \&Bio::Adventure::Assembly::Unicycler,
                '(shovill): Perform the shovill pre/post processing with spades.' => \&Bio::Adventure::Assembly::Shovill,
                '(terminasereorder): Reorder an existing assembly to the most likely terminase.' => \&Bio::Adventure::Phage::Terminase_ORF_Reorder,
                '(prodigal): Look for ORFs in bacterial/viral sequence.' => \&Bio::Adventure::Annotation::Prodigal,
                '(glimmer): Look for ORFs in bacterial/viral sequence.' => \&Bio::Adventure::Annotation::Glimmer,
                '(watson_plus): Make sure the Watson strand has the most pluses.' => \&Bio::Adventure::Annotation::Watson_Plus,
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
                '(jellyfish): Perform a kmer count of some data.' => \&Bio::Adventure::Count::Jellyfish,
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
        Phage => {
            name => 'phage',
            message => '',
            choices => {
                '(filterkraken): Filter out host sequences using a kraken report.' => \&Bio::Adventure::Phage::Filter_Host_Kraken,
                '(classifyphage): Use ICTV data to classify a phage assembly.' => \&Bio::Adventure::Phage::Classify_Phage,
                '(phageterm): Invoke phageterm to search for terminal repeats.' => \&Bio::Adventure::Annotation::Phageterm,
                '(phastaf): Search for phage regions in arbitrary assemblies.' => \&Bio::Adventure::Phage::Phastaf,
                '(terminasereorder): Reorder a phage assembly to the putative terminase.' => \&Bio::Adventure::Phage::Terminase_ORF_Reorder,
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
                '(priboseq): Perform a preset pipeline of ribosome profiling tasks.' => \&Bio::Adventure::Pipeline::Riboseq,
                '(ptnseq): Perform a preset pipeline of TNSeq tasks.' => \&Bio::Adventure::Pipeline::TNSeq,
                '(pbt1): Perform a preset group of bowtie1 tasks.' => \&Bio::Adventure::Pipeline::Bowtie,
                '(pbt2): Use preset bowtie2 tasks.' => \&Bio::Adventure::Pipeline::Bowtie2,
                '(phisat): Use preset tophat tasks.' => \&Bio::Adventure::Pipeline::Hisat,
                '(pbwa): Try preset bwa tasks.' => \&Bio::Adventure::Pipeline::BWA,
                '(psalmon): Try preset salmon tasks.' => \&Bio::Adventure::Pipeline::Salmon,
                '(passemble): Preset assembly.' => \&Bio::Adventure::Pipeline::Assemble,
                '(phageassemble): Preset phage assembly.' => \&Bio::Adventure::Pipeline::Phage_Assemble,
                '(annotateassembly): Do some searches on an assembly.' => \&Bio::Adventure::Pipline::Annotation_Assembly,
            },
        },
        PopulationGenetics => {
            name => 'popgen',
            message => 'The front pattern DOES move - and no wonder! The woman behind shakes it!',
            choices => {
                '(angsd): Get the set of filtered variants given a list of bam files.' => \&Bio::Adventure::PopGen::Angsd_Filter,
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
        Visualize => {
            name => 'visualize',
            message => 'When Red wins, she stands alone.  Go to page: 96485.332',
            choices =>  {
                '(cgview): Invoke cgview to visualize a genome.' => \&Bio::Adventure::Visualization::CGView,
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

sub Get_TODOs {
    my %args = @_;
    my $todo_list = ();
    my $possible_todos = {
        "abricate+" => \$todo_list->{todo}{'Bio::Adventure::Resistance::Abricate'},
        "abyss+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Abyss'},
        "angsdfilter+" => \$todo_list->{todo}{'Bio::Adventure::PopGen::Angsd_Filter'},
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
        "cgview+" => \$todo_list->{todo}{'Bio::Adventure::Visualization::CGView'},
        "classifyphage+" => \$todo_list->{todo}{'Bio::Adventure::Phage::Classify_Phage'},
        "copyraw+" => \$todo_list->{todo}{'Bio::Adventure::Prepare::Copy_Raw'},
        "countstates+" => \$todo_list->{todo}{'Bio::Adventure::Riboseq::Count_States'},
        "concat+" => \$todo_list->{todo}{'Bio::Adventure::Align::Concatenate_Searches'},
        "cutadapt+" => \$todo_list->{todo}{'Bio::Adventure::Trim::Cutadapt'},
        "essentialitytas+" => \$todo_list->{todo}{'Bio::Adventure::TNSeq::Essentiality_TAs'},
        "extendkraken+" => \$todo_list->{todo}{'Bio::Adventure::Annotation::Extend_Kraken_DB'},
        "extracttrinotate+" => \$todo_list->{todo}{'Bio::Adventure::Annotation::Extract_Trinotate'},
        "splitalignfasta+" => \$todo_list->{todo}{'Bio::Adventure::Align_Fasta::Split_Align_Fasta'},
        "fastado+" => \$todo_list->{todo}{'Bio::Adventure::Align_Fasta::Do_Fasta'},
        "fastasplitalign+" => \$todo_list->{todo}{'Bio::Adventure::Align_Fasta::Split_Align_Fasta'},
        "fastamerge+" => \$todo_list->{todo}{'Bio::Adventure::Align_Fasta::Merge_Parse_Fasta'},
        "fastaparse+" => \$todo_list->{todo}{'Bio::Adventure::Align_Fasta::Parse_Fasta'},
        "fastqct+" => \$todo_list->{todo}{'Bio::Adventure::QA::Fastqc'},
        "fastqdump+" => \$todo_list->{todo}{'Bio::Adventure::Prepare::Fastq_Dump'},
        "filterdepth+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Filter_Depth'},
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
        "interproscan+" => \$todo_list->{todo}{'Bio::Adventure::Annotation::Interproscan'},
        "jellyfish+" => \$todo_list->{todo}{'Bio::Adventure::Count::Jellyfish'},
        "kallisto+" => \$todo_list->{todo}{'Bio::Adventure::Map::Kallisto'},
        "kraken+" => \$todo_list->{todo}{'Bio::Adventure::Annotation::Kraken'},
        "mergeannotations+" => \$todo_list->{todo}{'Bio::Adventure::Annotation::Merge_Annotations'},
        "mergeparse+" => \$todo_list->{todo}{'Bio::Adventure::Align_Blast::Merge_Parse_Blast'},
        "mergeprodigal+" => \$todo_list->{todo}{'Bio::Adventure::Annotation::Merge_Annot_Prodigal'},
        "mimap+" => \$todo_list->{todo}{'Bio::Adventure::MiRNA::Mi_Map'},
        "phageterm+" => \$todo_list->{todo}{'Bio::Adventure::Annotation::Phageterm'},
        "phastaf+" => \$todo_list->{todo}{'Bio::Adventure::Phage::Phastaf'},
        "prodigal+" => \$todo_list->{todo}{'Bio::Adventure::Annotation::Prodigal'},
        "parseblast+" => \$todo_list->{todo}{'Bio::Adventure::Align_Blast::Parse_Blast'},
        "posttrinity+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Trinity_Post'},
        "prokka+" => \$todo_list->{todo}{'Bio::Adventure::Annotation::Prokka'},
        "racer+" => \$todo_list->{todo}{'Bio::Adventure::Trim::Racer'},
        "readsample+" => \$todo_list->{todo}{'Bio::Adventure::Prepare::Read_Samples'},
        "resfinder+" => \$todo_list->{todo}{'Bio::Adventure::Annotation::Resfinder'},
        "rgi+" => \$todo_list->{todo}{'Bio::Adventure::Annotation::Rgi'},
        "rsem+" => \$todo_list->{todo}{'Bio::Adventure::Map::RSEM'},
        "runessentiality+" => \$todo_list->{todo}{'Bio::Adventure::TNSeq::Run_Essentiality'},
        "sam2bam+" => \$todo_list->{todo}{'Bio::Adventure::Convert::Sam2Bam'},
        "salmon+" => \$todo_list->{todo}{'Bio::Adventure::Map::Salmon'},
        "terminasereorder+" => \$todo_list->{todo}{'Bio::Adventure::Phage::Terminase_ORF_Reorder'},
        "shovill+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Shovill'},
        "snippy+" => \$todo_list->{todo}{'Bio::Adventure::SNP::Snippy'},
        "alignsnpsearch+" => \$todo_list->{todo}{'Bio::Adventure::SNP::Align_SNP_Search'},
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
        "trainprodigal+" => \$todo_list->{todo}{'Bio::Adventure::Annotation::Train_Prodigal'},
        "transdecoder+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Transdecoder'},
        "trimomatic+" => \$todo_list->{todo}{'Bio::Adventure::Trim::Trimomatic'},
        "trinity+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Trinity'},
        "trinitypost+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Trinity_Post'},
        "trinotate+" => \$todo_list->{todo}{'Bio::Adventure::Annotation::Trinotate'},
        "trnascan+" => \$todo_list->{todo}{'Bio::Adventure::Annotation::tRNAScan'},
        "unicycler+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Unicycler'},
        "variantgenome+" => \$todo_list->{todo}{'Bio::Adventure::SNP::Make_Genome'},
        "velvet+" => \$todo_list->{todo}{'Bio::Adventure::Assembly::Velvet'},
        "watsonplus+" => \$todo_list->{todo}{'Bio::Adventure::Annotation::Watson_Plus'},
        ## A set of pipelines combining the above.
        "pannotate+" => \$todo_list->{todo}{'Bio::Adventure::Pipeline::Annotate_Assembly'},
        "passemble+" => \$todo_list->{todo}{'Bio::Adventure::Pipeline::Assemble'},
        "pbt1+" => \$todo_list->{todo}{'Bio::Adventure::RNAseq_Pipeline_Bowtie'},
        "pbt2+" => \$todo_list->{todo}{'Bio::Adventure::RNAseq_Pipeline_Bowtie2'},
        "pbwa+" => \$todo_list->{todo}{'Bio::Adventure::Pipeline::BWA'},
        "phageassemble+" => \$todo_list->{todo}{'Bio::Adventure::Pipeline::Phage_Assemble'},
        "phisat+" => \$todo_list->{todo}{'Bio::Adventure::Pipeline::Hisat'},
        "pkallisto+" => \$todo_list->{todo}{'Bio::Adventure::RNAseq_Pipeline_Kallisto'},
        "priboseq+" => \$todo_list->{todo}{'Bio::Adventure::Riboseq_Pipeline'},
        "psalmon+" => \$todo_list->{todo}{'Bio::Adventure::Pipeline::Salmon'},
        "ptnseq+" => \$todo_list->{todo}{'Bio::Adventure::TNseq_Pipeline'},
        "ptophat+" => \$todo_list->{todo}{'Bio::Adventure::RNAseq_Pipeline_Tophat'},
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

    }
    $input =~ tr/[A-Z]/[a-z]/;
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
    for my $varname (keys %getvars_default_vars) {
        ## required and args are provided in the function call, so skip them.
        next if ($varname eq 'required' || $varname eq 'args');
        $returned_vars{$varname} = $getvars_default_vars{$varname};
    }
    ## Pick up options passed in the parent function call, e.g.
    ## my $bob = $class->Function(bob => 'jane');
    for my $varname (keys %function_override_vars) {
        ## required and args are provided in the function call, so skip them.
        next if ($varname eq 'required' || $varname eq 'args');
        $returned_vars{$varname} = $function_override_vars{$varname};
        ## Try to ensure that the shell reverts to the default 'bash'
        if (!defined($function_override_vars{language})) {
            $returned_vars{language} = 'bash';
            $returned_vars{shell} = '/usr/bin/env bash';
        }
    }
    ## Final loop to pick up options from the commandline or a TERM prompt.
    ## These supercede everything else.
    for my $varname (keys %getopt_override_vars) {
        next if ($varname eq 'required' || $varname eq 'args');
        next unless (defined($getopt_override_vars{$varname}));
        $returned_vars{$varname} = $getopt_override_vars{$varname};
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


=item C<Set_Vars>

  Handle the peculiar mix of instance options held in $class->{options},
  the set of arguments passed to %args, a list of required and potentially
  missing options, and whatever default values one wishes to set.

  $options = $class->Set_Vars(exclude => 'bob');

=cut
sub Set_Vars {
    my ($class, %args) = @_;
    my $options = $class->{options};
    my $ref = ref($options);
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


sub Module_Loader {
    my ($class, %args) = @_;
    my $action = 'load';
    if ($args{action}) {
        $action = $args{action};
    }
    my @mod_lst;
    my $mod_class = ref($args{modules});
    if ($mod_class eq 'SCALAR') {
        push(@mod_lst, $args{modules});
    } elsif ($mod_class eq 'ARRAY') {
        @mod_lst = @{$args{modules}};
    } elsif (!$mod_class) {
        push(@mod_lst, $args{modules});
    } else {
        print "I do not know this class: $mod_class, $args{modules}\n";
    }

    my $count = 0;
    my $mod;
    if ($action eq 'unload') {
        @mod_lst = reverse(@mod_lst);
    }
    foreach $mod (@mod_lst) {
        if ($args{verbose}) {
            print "(Un)loading $mod\n";
        }
        my ($stdout, $stderr, @result) = capture {
            try {
                my $test;
                if ($action eq 'load') {
                    $test = Env::Modulecmd::load($mod);
                } else {
                    $test = Env::Modulecmd::unload($mod);
                }
            } catch ($e) {
                print "There was an error loading ${mod}, ${e}\n";
            };
        };
        $count++;
    }
    return($count);
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
    my $genome_file = qq"$options->{basedir}/${fasta_name}.pdata";
    if (-r $genome_file) {
        $chromosomes = lock_retrieve($genome_file);
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
        } ## End reading each chromosome
        my $stored = lock_store($chromosomes, $genome_file);
    } ## End checking for a .pdata file
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
    my $genome_file = qq"$options->{basedir}/${gff_name}.pdata";
    print "Reading $options->{gff}, seeking features tagged (last column ID) $id_tag,
 type (3rd column): $feature_type.\n";
    if (-r $genome_file && !$options->{debug}) {
        $gff_out = lock_retrieve($genome_file);
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
        if (-f $genome_file) {
            unlink($genome_file);
        }
        my $stored = lock_store($gff_out, $genome_file);
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
    my $options = $class->Get_Vars(
        args => \%args);

    ## If we are invoking an indirect job, we need a way to serialize the options
    ## in order to get them passed to the eventual interpreter
    my $option_file = "";
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
        $option_file = File::Temp->new(
            TEMPLATE => 'optionsXXXX',
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
            use Data::Dumper;
            print Dumper \%saved_options;
        }
    }
    my $runner;
    if ($class->{sbatch_path}) {
        $runner = Bio::Adventure::Slurm->new();
    } elsif ($class->{qsub_path}) {
        $runner = Bio::Adventure::Torque->new();
    } elsif ($class->{bash_path}) {
        ## I should probably have something to handle gracefully bash jobs.
        $runner = Bio::Adventure::Local->new();
    } else {
        die("Could not find sbatch, qsub, nor bash.");
    }

    ## Add the current options to the runner:
    for my $k (keys %{$options}) {
        $runner->{$k} = $options->{$k};
    }
    my $result = $runner->Submit($class, %args);
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
