package HPGL::Job;
our %jobs;

=head 1 NAME
    HPGL::Job
    A super-simple state-keeper for torque jobs.
=cut
sub new {
    my ($class, %args) = @_;
    my $me = bless {}, $class;
    foreach my $key (keys %args) {
        $me->{$key} = $args{$key} if ($args{$key});
    }
    my $id = $me->{id};
    $jobs{$id} = $me;
    return($me);
}

package HPGL;
use common::sense;
use autodie qw":all";
use Getopt::Long;
use Cwd qw"cwd";
use File::Basename;
use File::Find;
use Term::ReadLine;
use Archive::Extract;
use PerlIO;
use Digest::MD5 qw"md5 md5_hex md5_base64";
use Net::Amazon::S3;
use Digest::MD5 qw(md5 md5_hex md5_base64);
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

    ## Default attributes for the HPGL class go here
    ## Note that they may be overwritten if specified below with Getoptions()
    my $bash = "/usr/bin/bash";
    my %options = (
        align_trimmed => 1,
        aligner => "bowtie",
        aws_access_key_id => 'PQOD56HDPQOXJPYQYHH9',
        aws_hostname => 'gembox.cbcb.umd.edu',
        aws_secret_access_key => 'kNl27xjeVSnChzEww9ziq1VkgUMNmNonNYjxWkGw',
        bam_small => "21-26",
        bam_mid => "27-32",
        bam_big => "33-40",
        basedir => cwd(),
        bt_args => {},
        btmulti => 0,
        debug => 1,
        depends => {},
        help => undef,
        hpgl => undef,
        hpglid => undef,
        htseq_identifier => "ID",
        htseq_stranded => "no",
        input => undef,
        jobs => [],
        jobids => {},
        len_min => 21,
        len_max => 40,
        libdir => "$ENV{HOME}/libraries",
        libtype => "genome",
        paired => 0,
        pbs => 1,
        species => undef,
        trimmer => "trimmomatic",
        type => "rnaseq",
    );
    $options{qsub_args} = "-j oe -V -m n";
    $options{qsub_queue} = "throughput";
    $options{qsub_queues} = ["throughput","workstation","long","large"];
    $options{qsub_shell} = "$bash";
    $options{qsub_mem} = 6;
    $options{qsub_wall} = "10:00:00";
    $options{qsub_cpus} = "4";
    $options{qsub_depends} = "depend=afterok:";
    $options{qsub_loghost} = "ibissub00.umiacs.umd.edu";
    $options{qsub_logdir} = qq"$options{qsub_loghost}:$options{basedir}/outputs";

    my $opt = GetOptions(
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
    foreach my $k (keys %args) {
        $options{$k} = $args{$k};
    }
    foreach my $k (keys %options) {
        $me->{$k} = $options{$k};
    }
    $me->Check_Options(["input",]);
    my $basename = $me->{input};
    my @suffixes = (".fastq",".gz",".xz", ".fasta", ".sam", ".bam", ".count");
    $basename = basename($basename, @suffixes);
    $basename = basename($basename, @suffixes);    
    $me->{basename} = $basename;
    $me->{bt_args} = {
            v0M1 => "--best -v 0 -M 1",
##            def => " ",
##            v1M1l20 => "--best -v 1 -M 1 -y -l 10",
            v1M1 => "--best -v 1 -M 1",
            v2M1 => "--best -v 2 -M 1",
#            v1 => "--best -v 1",
    };
    if (!defined($me->{hpglid})) {
        my $tmp = $me->{input};
        $tmp =~ s/^(hpgl\d+).*/$1/g;
        $me->{hpglid} = $tmp;
    }
    return($me);
}

=head2
    Biopieces_Graph()
    Reads in a fastq file and uses biopieces to make some graphs
    describing the sequences therein.
=cut
sub Biopieces_Graph {
    my $me = shift;
    my %args = @_;
    my $input = $me->{input};
    my $bp_depends_on;
    my $basename = $me->{basename};    
    $bp_depends_on = $args{depends} if ($args{depends});
    my @inputs = split(/\,/, $input);
    my $comment = qq!## This script uses biopieces to draw some simple graphs of the sequence.!;
    my $job_string = qq!
mkdir -p biopieces
xzcat -f ${input} | read_fastq -i - -e base_33 |\\
 plot_scores -T \'Quality Scores\' -t svg -o biopieces/${basename}_quality_scores.svg |\\
 plot_nucleotide_distribution -T 'NT. Distribution' -t svg -o biopieces/${basename}_ntdist.svg |\\
 plot_lendist -T 'Length Distribution' -k SEQ_LEN -t svg -o biopieces/${basename}_lendist.svg |\\
!;
    my $bp = $me->Qsub(job_name => "biop",
                       depends => $bp_depends_on,
                       job_string => $job_string,
                       comment => $comment,
                       input => $input,
		       prescript => $args{prescript},
		       postscript => $args{postscript},
	);
    return($bp);
}

=head2
    Bowtie()
=cut
sub Bowtie {
    my $me = shift;
    my %args = @_;
    my %bt_jobs = ();
    my $bt_input = $me->{input};
    my $bt_depends_on;
    $bt_depends_on = $args{depends} if ($args{depends});
    $me->Check_Options(["species"]);
    my $basename = $me->{basename};

    my $bt_type = "v0M1";
    $bt_type = $args{bt_type} if ($args{bt_type});
    my $bt_args = $me->{bt_args}->{$bt_type};
    my $jobname = qq"bt${bt_type}";
    $jobname = $args{jobname} if ($args{jobname});
    my $libtype = $me->{libtype};
    $libtype = $args{libtype} if ($args{libtype});
    my $count = 1;
    $count = $args{count} if (defined($args{count}));

    ## Check that the indexes exist
    my $bt_reflib = "$me->{libdir}/${libtype}/$me->{species}";
    my $bt_reftest = qq"${bt_reflib}.1.ebwt";
    if (!-r $bt_reftest) {
        my $index_job = $me->BT1_Index(depends => $bt_depends_on, libtype => $libtype);
        $bt_jobs{index} = $index_job;
        $bt_depends_on = $index_job->{pbs_id};
    }
    my $bowtie_input_flag = "-q";  ## fastq by default
    $bowtie_input_flag = "-f" if ($me->{input} =~ /\.fasta$/);
    my $basename = $me->{basename};
    my $species = $me->{species};
    my $error_file = qq"bowtie_out/${basename}-${bt_type}.err";
    my $trim_output_file = qq"outputs/${basename}-trimomatic.out";
    my $comment = qq!## This is a bowtie1 alignment of $bt_input against
## $bt_reflib using arguments: $bt_args.
## This jobs depended on: $bt_depends_on
!;
    my $aligned_filename = qq"bowtie_out/${basename}-${bt_type}_aligned_${species}.fastq";
    my $unaligned_filename = qq"bowtie_out/${basename}-${bt_type}_unaligned_${species}.fastq";
    my $sam_filename = qq"bowtie_out/${basename}-${bt_type}.sam";
    my $job_string = qq!mkdir -p bowtie_out && sleep 10 && bowtie $bt_reflib $bt_args \\
  -p 4 \\
  $bowtie_input_flag $bt_input \\
  --un ${unaligned_filename} \\
  --al ${aligned_filename} \\
  -S ${sam_filename} \\
  2>${error_file} \\
  1>bowtie_out/${basename}-${bt_type}.out
!;
    my $bt_job = $me->Qsub(job_name => qq"bt1_${bt_type}",
                           depends => $bt_depends_on,
                           job_string => $job_string,
                           input => $bt_input,
                           comment => $comment,
                           output => $sam_filename,
                           unaligned => $unaligned_filename,
                           aligned => $aligned_filename,
			   prescript => $args{prescript},
			   postscript => $args{postscript},
        );
    $bt_jobs{bowtie} = $bt_job;

    my $un_comp = $me->Recompress(depends => $bt_job->{pbs_id},
                                  job_name => "xzun",
                                  comment => qq"## Compressing the sequences which failed to align against $bt_reflib using options $bt_args\n",
                                  input => "bowtie_out/${basename}-${bt_type}_unaligned_${species}.fasta");
    $bt_jobs{unaligned_compression} = $un_comp;

    my $al_comp = $me->Recompress(input => "bowtie_out/${basename}-${bt_type}_aligned_${species}.fasta",
                                  comment => qq"## Compressing the sequences which successfully aligned against $bt_reflib using options $bt_args",
                                  job_name => "xzal",
                                  depends => $bt_job->{pbs_id},);
    $bt_jobs{aligned_compression} = $al_comp;

    my $stats = $me->Get_Stats(depends => $bt_job->{pbs_id},
                               job_name => "bt1stats",
                               bt_type => $bt_type,
                               count_table => qq"${basename}-${bt_type}.count.xz",
                               trim_input => ${trim_output_file},
                               bt_input => $error_file);
    
    my $sam_job = $me->Samtools(depends => $bt_job->{pbs_id},
                                job_name => "s2b",
                                input => $bt_job->{output},);
    
    $bt_jobs{samtools} = $sam_job;
    $me->{output} = $sam_job->{output};
    my $htmulti;
    if ($count) {
        if ($libtype eq 'rRNA') {
            $htmulti = $me->HTSeq(depends => $sam_job->{pbs_id},
                                  suffix => $bt_type,
                                  libtype => $libtype,
                                  input => $sam_job->{output},);
        } else {
            $htmulti = $me->HT_Multi(depends => $sam_job->{pbs_id},
                                     suffix => $bt_type,
                                     libtype => $libtype,
                                     input => $sam_job->{output},);
            $bt_jobs{htseq} = $htmulti;
        }
    }
    return(\%bt_jobs);
}

=head2
    BT_Multi()
=cut
sub BT_Multi {
    my $me = shift;
    my %args = @_;
    $me->Check_Options(["species"]);
    my $basename = $me->{basename};
    my $species = $me->{species};
    my $depends_on = $args{depends};
    my %bt_types = %{$me->{bt_args}};
    foreach my $type (keys %bt_types) {
        print "Starting $type\n";
        my $job = $me->Bowtie(bt_type => $type,
                              depends => $depends_on,
                              jobname => qq"bt${type}",
			      prescript => $args{prescript},
			      postscript => $args{postscript},
            );
    }
}

=head2
    Bowtie_RRNA()
    Perform an alignment against a home-curated set of ribosomal
    RNA/tRNA sequences.  The alignment requires a fastq input library
    and fasta library found in 'libraries/rRNA/$me->{species}.fasta'

  Example:
    my $rrna = $hpgl->Bowtie_RRNA();
    ## If you want to exclude the rRNA sequences from future alignments:
    my $rrna = $hpgl->Bowtie_RRNA(exclude => 1);
    
=cut
sub Bowtie_RRNA {
    my $me = shift;
    my %args = @_;
    $me->Check_Options(["species"]);
    my $basename = $me->{basename};
    $me->{basename} = qq"rRNA_${basename}";
    my $exclude = 0;
    $exclude = $args{exclude} if ($args{exclude});
    my $species = $me->{species};
    my $depends_on = $args{depends};
    if ($exclude) {
        my $in = $me->{input};
        $args{postscript} = qq"mv $in includerrna_${in} && mv bowtie_out/${basename}-rRNA_unaligned_${species}.fastq $in";
    }
    my $job = $me->Bowtie(depends => $depends_on,
                          job_name => qq"btrrna",
                          libtype => 'rRNA',
                          prescript => $args{prescript},
                          postscript => $args{postscript},
        );
    ## Return the basename back to normal so that future tasks don't
    ## get confuseled.
    $me->{basename} = $basename;
    return($job);
}

=head2
    BT1_Index()
=cut
sub BT1_Index {
    my $me = shift;
    my %args = @_;
    my $basename = $me->{basename};
    my $dep = "";
    $dep = $args{depends};
    my $libtype = $me->{libtype};
    $libtype = $args{libtype} if ($args{libtype});
    my $job_string = qq!bowtie-build $me->{libdir}/${libtype}/$me->{species}.fasta $me->{libdir}/genome/$me->{species}
!;
    my $comment = qq!## Generating bowtie1 indexes for species: $me->{species} in $me->{libdir}/${libtype}!;
    my $bt1_index = $me->Qsub(job_name => "bt1idx",
                              depends => $dep,
                              job_string => $job_string,
                              comment => $comment,
			      prescript => $args{prescript},
			      postscript => $args{postscript},
        );
    return($bt1_index);
}

=head2
    BT2_Index()
=cut
sub BT2_Index {
    my $me = shift;
    my %args = @_;
    my $basename = $me->{basename};
    my $dep = "";
    $dep = $args{depends};
    my $job_string = qq!
if [ \! -r "$me->{libdir}/genome/$me->{species}.fa" ]; then
  ln -s $me->{libdir}/genome/$me->{species}.fasta $me->{libdir}/genome/$me->{species}.fa
fi
bowtie2-build $me->{libdir}/genome/$me->{species}.fasta $me->{libdir}/genome/$me->{species}
!;
    my $comment = qq!## Generating bowtie2 indexes for species: $me->{species} in $me->{libdir}/genome!;
    my $jobid = $me->Qsub(job_name => "bt2idx",
                          depends => $dep,
                          job_string => $job_string,
                          comment => $comment,
			  prescript => $args{prescript},
			  postscript => $args{postscript},
        );
    return($jobid);
}

=head2
    BWA()
=cut
sub BWA {
    my $me = shift;
    my %args = @_;
    my %bwa_jobs = ();
    my $bwa_input = $me->{input};
    my $bwa_depends_on;
    $me->Check_Options(["species"]);
    my $basename = $me->{basename};

    ## Check that the indexes exist
    my $bwa_reflib = "$me->{libdir}/genome/$me->{species}";
    my $bwa_reftest = qq"${bwa_reflib}.bwa";
    if (!-r $bwa_reftest) {
        my $index_job = $me->BWA_Index(depends => $bwa_depends_on);
        $bwa_jobs{index} = $index_job;
        $bwa_depends_on = $index_job->{pbs_id};
    }
    my $basename = $me->{basename};
    my $species = $me->{species};
    my $comment = qq!## This is a BWA alignment of $bwa_input against
## $bwa_reflib.!;
    my $job_string = qq!mkdir -p bwa_out
cd bwa_out && bwa aln $me->{libdir}/genome/$me->{species}.fasta $bwa_input 2>bwa.err 1>${basename}.sai
!;
    my $bwa_job = $me->Qsub(job_name => "bwa",
			    depends => $bwa_depends_on,
			    job_string => $job_string,
			    input => $bwa_input,
			    comment => $comment,
			    output => qq"bwa/${basename}.sai",
			    prescript => $args{prescript},
			    postscript => $args{postscript},
        );
    $bwa_jobs{bowtie} = $bwa_job;

    my $sam_job = $me->Samtools(depends => $bwa_job->{pbs_id},
                                job_name => "s2b",
                                input => $bwa_job->{output},);
    $bwa_jobs{samtools} = $sam_job;

    my $htmulti = $me->HT_Multi(depends => $sam_job->{pbs_id},
                                input => $sam_job->{output},);
    $bwa_jobs{htseq} = $htmulti;
    
    my $c = 0;
    foreach my $ht_out (@{$htmulti}) {
        $c++;
        if ($ht_out) {
            my $count_compression = $me->Recompress(depends => $ht_out->{pbs_id},
                                                    job_name => "xz_hts${c}",
                                                    input => $ht_out->{output});
        }
    }
    return(\%bwa_jobs);
}

=head2
    BWA_Index()
=cut
sub BWA_Index {
    my $me = shift;
    my %args = @_;
    my $basename = $me->{basename};
    my $dep = "";
    $dep = $args{depends};
    my $job_string = qq!cd $me->{libdir}/genome/ && bwa index $me->{libdir}/genome/$me->{species}.fasta
!;
    my $comment = qq!## Generating bwa indexes for species: $me->{species} in $me->{libdir}/genome!;
    my $bwa_index = $me->Qsub(job_name => "bwaidx",
                              depends => $dep,
                              job_string => $job_string,
                              comment => $comment,
			      prescript => $args{prescript},
			      postscript => $args{postscript},
        );
    return($bwa_index);
}

=head2
    Cleanup_Seq()
=cut
sub Cleanup_Seq {
    my $me = shift;
    my %args = @_;
    my $depends = "";
    if ($args{depends}) {
        $depends = $args{depends};
    }
    my $job_string = qq!rm -rf status outputs scripts sequences && find . -name '*.fastq*' -exec rm {} ';'
!;
    my $clean = $me->Qsub(job_name => "th_clean",
                          qsub_queue => 'throughput',
                          depends => $depends,
                          job_string => $job_string,
			  prescript => $args{prescript},
			  postscript => $args{postscript},
        );
    return($clean);
}

=head2
    Cleanup_Tophat()
=cut
sub Cleanup_Tophat {
    my $me = shift;
    my %args = @_;
    my $depends = "";
    if ($args{depends}) {
        $depends = $args{depends};
    }
    my $job_string = qq!rm -rf tophat_out
!;
    my $clean = $me->Qsub(job_name => "th_clean",
                          qsub_queue => 'throughput',
                          depends => $depends,
                          job_string => $job_string,
			  prescript => $args{prescript},
			  postscript => $args{postscript},
        );
    return($clean);
}

=head2
    Cutadapt()
=cut
sub Cutadapt {
    my $me = shift;
    my %args = @_;
    my $input = $me->{input};
    my $basename = $me->{basename};
    my @inputs = split(/\,/, $input);
    my $output = qq"${basename}-trimmed.fastq";
    my $cutadapt_flags = " -a ACAGGTTGGATGATAAGTCCCCGGTCTGACACATCTCCCTAT -a AGATCGGAAGAGCACACGTCTGAAC -b AGATCGGAAGAGCACACGTCTGAAC ";
    my $minlen = 12;
    my $maxlen = 40;

    my $comment = qq!## This script makes use of biopieces and cutadapt to trim away adapters
## and separate the sequence file into a few pieces depending on size
## and adapter status.  It also performs some simple graphs of the data.!;
    my $job_string = qq!
mkdir -p cutadapt 
xzcat -f ${input} | cutadapt - ${cutadapt_flags} -e 0.1 -n 3 -m ${minlen} -M ${maxlen} \\
    --too-short-output=cutadapt/${basename}_tooshort.fastq \\
    --too-long-output=cutadapt/${basename}_toolong.fastq \\
    --untrimmed-output=cutadapt/${basename}_untrimmed.fastq \\
    2>cutadapt/cutadapt.err 1> ${output}!;
    my $cutadapt = $me->Qsub(job_name => "cutadapt",
                             qsub_wall => "8:00:00",
                             job_string => $job_string,
                             input => $input,
                             output => $output,
                             comment => $comment,
			     prescript => $args{prescript},
			     postscript => $args{postscript},
	);
    my $comp_short = $me->Recompress(job_name => "xzcutshort",
                                     depends => $cutadapt->{pbs_id},
                                     input => qq"cutadapt/${basename}_tooshort.fastq",
                                     comment => qq"## Compressing the tooshort sequences.",
	);
    my $comp_long = $me->Recompress(job_name => "xzcutlong",
                                    depends => $cutadapt->{pbs_id},
                                    input => qq"cutadapt/${basename}_toolong.fastq",
                                    comment => qq"## Compressing the toolong sequences.",
	);
    my $comp_un = $me->Recompress(job_name => "xzuncut",
                                  depends => $cutadapt->{pbs_id},
                                  input => qq"cutadapt/${basename}_untrimmed.fastq",
                                  comment => qq"## Compressing the toolong sequences.",
	);
    my $comp_original = $me->Recompress(job_name => "xzorig",
                                        depends => $cutadapt->{pbs_id},
                                        input => qq"$input",
                                        output => qq"sequences/${input}.xz",
                                        comment => qq"## Compressing the original sequence.",
	);
}

=head2
    Fastqc_Pairwise()
=cut
sub Fastqc_Pairwise {
    my $me = shift;
    my %args => @_;
    my $type = $args{type};
    $type = "unfiltered" unless ($type);
    my $input = $me->{input};
    my @input_list = split(/\:|\,/, $input);
    if (scalar(@input_list) <= 1) {
        my $ret = $me->Fastqc_Single(@_);
        return($ret);
    }
    my $r1 = $input_list[0];
    my $r2 = $input_list[1];
    my $basename = basename($r1, (".fastq"));
    $basename =~ s/\_R1$//g;
    $me->{basename} = $basename;
    my $job_string = qq!mkdir -p ${basename}-${type}_fastqc &&\\
  fastqc -o ${basename}-${type}_fastqc ${r1} ${r2} \\
  2>outputs/${basename}-${type}_fastqc.out 1>&2
!;
    my $comment = qq!## This FastQC run is against ${type} data and is used for
## an initial estimation of the overall sequencing quality.!;
    my $fqc_jobid = qq"${basename}_fqc";
    my $jobid = $me->Qsub(job_name => "fqc",
                          qsub_cpus => 8,
                          qsub_queue => "workstation",
                          job_string => $job_string,
                          comment => $comment,
			  prescript => $args{prescript},
			  postscript => $args{postscript},
        );
}

=head2
    Fastqc_Single()
=cut
sub Fastqc_Single {
    my $me = shift;
    my %args => @_;
    my $type = $args{type};
    $type = "unfiltered" unless ($type);
    my $input = $me->{input};
    my $basename = $me->{basename};
    my $job_string = qq!mkdir -p ${basename}-${type}_fastqc &&\\
  fastqc -o ${basename}-${type}_fastqc ${input} \\
  2>outputs/${basename}-${type}_fastqc.out 1>&2
!;
    my $comment = qq!## This FastQC run is against ${type} data and is used for
## an initial estimation of the overall sequencing quality.!;
    my $fqc_jobid = qq"${basename}_fqc";
    my $jobid = $me->Qsub(job_name => "fqc",
                          qsub_cpus => 8,
                          qsub_queue => "workstation",
                          job_string => $job_string,
                          comment => $comment,
			  prescript => $args{prescript},
			  postscript => $args{postscript},
        );
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

=head2
    Gff2Gtf()
=cut
sub Gff2Gtf {
    my $in = shift;
    use File::Basename;
    use Bio::FeatureIO;
 
    my ($name, $path, $suffix) = fileparse($in, qr/\.gff/);
    my $outFile = $path . $name . ".gtf";
 
    my $inGFF = new Bio::FeatureIO('-file' => "$in",
                                   '-format' => 'GFF',
                                   '-version' => 3);
    my $outGTF = new Bio::FeatureIO('-file' => ">$outFile",
                                    '-format' => 'GFF',
                                    '-version' => 2.5);
    
    while (my $feature = $inGFF->next_feature() ) {
        $outGTF->write_feature($feature);
    }
}

=head2
    Help()
=cut
sub Help {
    print "There are many help functions, but this is my own.\n";
}

=head2
    HT_Multi()
=cut
sub HT_Multi {
    my $me = shift;
    my %args = @_;
    my @jobs = ();
    $me->Check_Options(["species"]);
    my $basename = $me->{basename};
    my $species = $me->{species};
    my %gff_types = (
        misc => "transcript",
        rintron => "exon",
        nmd => "exon",
        linc => "gene",
        pseudo => "exon",
        antisense => "exon",
        exon => "exon",
        fiveputr => "five_prime_UTR",
        threeputr => "three_prime_UTR",
        sn => "exon",
        mi => "exon",
        sno => "exon",
        rrna => "gene",
        interCDS => "CDS",
        operons => "gene",
    );
    my $htseq_runs = 0;
    my $input = $args{input};
##    my $aligntype = $args{aligntype};
    my $suffix = $args{suffix};
    my $out_prefix = $input;
    $out_prefix =~ s/\.bam$//g;
    foreach my $type (keys %gff_types) {
        my $gff = qq"$me->{libdir}/genome/${species}_${type}.gff";
        my $output = qq"${out_prefix}_${type}.count";
        if (-r "$gff") {
            print "Found $gff, performing htseq with it.\n";
            my $ht = $me->HTSeq(htseq_type => $gff_types{$type},
                                htseq_gff => $gff,
                                jobname => "ht${type}${suffix}",
                                depends => $args{depends},
                                input => $input,
                                output => $output,
				prescript => $args{prescript},
				postscript => $args{postscript},
                                suffix => $args{suffix},
                );
            push(@jobs, $ht);
            $htseq_runs++;
        } else {
            print "Did not find $gff, skipping.\n";
        }
    } ## End foreach type

    ## Also perform a whole genome count
    my $gff = qq"$me->{libdir}/genome/${species}.gff";
    my $output = $input;
    $output =~ s/\.bam$/\.count/g;
    if (-r "$gff") {
        my $ht = $me->HTSeq(htseq_type => "all",
                            htseq_gff => $gff,
                            jobname => "htall${suffix}",
                            depends => $args{depends},
                            input => $input,
                            output => $output,
                            prescript => $args{prescript},
                            postscript => $args{postscript},
            );
        push(@jobs, $ht);
    }
    return(\@jobs);
}

=head2
    HTSeq()
=cut
sub HTSeq {
    my $me = shift;
    my %args = @_;
    $me->Check_Options(["species", "htseq_stranded", "htseq_identifier"]);
    my $basename = $me->{basename};
    my $species = $me->{species};
    my $stranded = $me->{htseq_stranded};
    my $identifier = $me->{htseq_identifier};
    my $htseq_jobname = qq"htseq";
    $htseq_jobname = $args{jobname} if ($args{jobname});
    my $depends = "";
    $depends = $args{depends} if ($args{depends});
    my $libtype = "genome";
    $libtype = $args{libtype} if ($args{libtype});
    my $type = "exon";
    $type = $args{htseq_type} if ($args{htseq_type});
    my $gff = qq"$me->{libdir}/${libtype}/${species}.gff";
    $gff = $args{htseq_gff} if ($args{htseq_gff});
    my $input = "${basename}.bam";
    $input = $args{input} if ($args{input});
    my $output = $input;
    if ($args{suffix}) {
        my $suffix = $args{suffix};
        $output =~ s/\.bam/_${suffix}\.count/g;
    } else {
        $output =~ s/\.bam/\.count/g;
    }
    $output = $args{output} if ($args{output});
    my $error = $input;
    $error =~ s/\.bam/\.error/g;
    if (!-r "$gff") {
        print "Unable to read $gff, please fix this and try again.\n";
    }
    $type = 'gene' if ($type eq 'all');
    my $job_string = qq!htseq-count -q -f bam -s $stranded -t $type -i $identifier $input $gff 1>$output 2>$error && pxz $output
!;
    my $comment = qq!## Counting the number of hits in $input for each feature found in $gff
## Is this stranded? $stranded.  The attribute type (3rd column of the gff file) is $type
## and the comment field (10th column) used to name the feature is $identifier!;
    my $htseq = $me->Qsub(job_name => $htseq_jobname,
                          qsub_mem => 6,
                          depends => $depends,
                          job_string => $job_string,
                          comment => $comment,
                          input => $input,
                          output => $output,
			  prescript => $args{prescript},
			  postscript => $args{postscript},
        );
    return($htseq);

}

=head2
    Recompress()
=cut
sub Recompress {
    my $me = shift;
    my %args = @_;
    my $input;
    if ($args{input}) {
        $input = $args{input};
    } else {
        $input = $me->{input};
    }
    my ($input1, $input2);
    if ($input =~ /\:|\,/) {
	($input1, $input2) = split(/\:|\,/, $input);
    } else {
	$input1 = $input;
    }
    my $job_name = "xz";
    $job_name = $args{job_name} if ($args{job_name});
    my $depends = "";
    $depends = $args{depends} if ($args{depends});
    my $comment = "";
    $comment = $args{comment} if ($args{comment});
    my $basename = $me->{basename};
    my $job_string = qq!pxz -f -9e $input1 $input2
!;
    if ($args{output}) {
        $job_string .= qq!mv ${input1}.xz $args{output}
!;
	if ($input2) {
	    $job_string .= qq!mv ${input2}.xz $args{output2}
!;
	}
    }

    my $trim_jobid = qq"${basename}_xz";
    my $compression = $me->Qsub(job_name => $job_name,
                                depends => $depends,
                                qsub_queue => "throughput",
                                qsub_cpus => 1,
                                qsub_mem => 4,
                                qsub_wall => "18:00:00",
                                job_string => $job_string,
                                input => $input,
                                comment => $comment,
				prescript => $args{prescript},
				postscript => $args{postscript},
        );
    return($compression);
}

=head2
    Uncompress()
=cut
sub Uncompress {
    my $me = shift;
    my %args = @_;
    my $input;
    if ($args{input}) {
        $input = $args{input};
    } else {
        $input = $me->{input};
    }
    my ($input1, $input2);
    if ($input =~ /\:|\,/) {
	($input1, $input2) = split(/\:|\,/, $input);
    } else {
	$input1 = $input;
    }
    my $job_name = "unxz";
    $job_name = $args{job_name} if ($args{job_name});
    my $depends = "";
    $depends = $args{depends} if ($args{depends});
    my $comment = "";
    $comment = $args{comment} if ($args{comment});
    my $basename = $me->{basename};
    my $job_string = qq!pxz -f -d $input1 $input2
!;

    my $trim_jobid = qq"${basename}_unxz";
    my $compression = $me->Qsub(job_name => $job_name,
                                depends => $depends,
                                qsub_queue => "throughput",
                                qsub_cpus => 1,
                                qsub_mem => 4,
                                qsub_wall => "18:00:00",
                                job_string => $job_string,
                                input => $input,
                                comment => $comment,
				prescript => $args{prescript},
				postscript => $args{postscript},
        );
    return($compression);
}

=head2
    Get_Stats()
=cut
sub Get_Stats {
    my $me = shift;
    my %args = @_;
    my $bt_input = $args{bt_input};
    my $trim_input = $args{trim_input};
    my $basename = $me->{basename};
    my $bt_type = "";
    $bt_type = $args{bt_type} if ($args{bt_type});
    my $depends = "";    
    $depends = $args{depends} if ($args{depends});
    my $job_name = "stats";
    $job_name = $args{job_name} if ($args{job_name});    
    my $jobid = qq"${basename}_stats";
    my $count_table = "";
    $count_table = $args{count_table} if ($args{count_table});
    my $comment = qq!
## This is a stupidly simple job to collect alignment statistics.
!;
    my $job_string = qq!
original_reads_tmp=\$(grep "^Input Reads" $trim_input | awk '{print \$3}' | sed 's/ //g')
original_reads=\${original_reads_tmp:-0}
reads_tmp=\$(grep "^# reads processed" $bt_input | awk -F: '{print \$2}' | sed 's/ //g')
reads=\${reads_tmp:-0}
one_align_tmp=\$(grep "^# reads with at least one reported" $bt_input | awk -F": " '{print \$2}' | sed 's/ .*//g')
one_align=\${one_align_tmp:-0}
failed_tmp=\$(grep "^# reads that failed to align" $bt_input | awk -F": " '{print \$2}' | sed 's/ .*//g')
failed=\${failed_tmp:-0}
sampled_tmp=\$(grep "^# reads with alignments sampled" $bt_input | awk -F": " '{print \$2}' | sed 's/ .*//g')
sampled=\${sampled_tmp:-0}
rpm_tmp=\$(perl -e "printf(1000000 / \${one_align})" 2>/dev/null)
rpm=\${rpm_tmp:-0}
stat_string=\$(printf "${basename},${bt_type},%s,%s,%s,%s,%s,%s,${count_table}" "\${original_reads}" "\${reads}" "\${one_align}" "\${failed}" "\${sampled}" "\$rpm")
echo "\$stat_string" >> ../bowtie_stats.csv
!;
    my $stats = $me->Qsub(job_name => $job_name,
                          depends => $depends,
                          qsub_queue => "throughput",
                          qsub_cpus => 1,
                          qsub_mem => 1,
                          qsub_wall => "00:10:00",
                          job_string => $job_string,
                          input => $bt_input,
                          comment => $comment,
        );
    return($stats);
    
}

=head2
    Sam2Bam()
=cut
sub Sam2Bam {
    my $me = shift;
    $me->Check_Options(["species"]);
    my %args = @_;
    ## A little logic to see if the invocation of this is for an extent .sam file
    ## Or calling it on an existing .fastq(.gz), in which case one must assume
    ## It is being called on one or more files in the bowtie_out/ directory
    ## And which start with $basename, include something like -trimmed-v0M1.sam
    my $basename = $me->{basename};
    my @input_list = ();
    if ($args{input}) {
        push(@input_list, $args{input});
    } elsif (-r $me->{input} and $me->{input} =~ /\.sam$/) {
        push(@input_list, $me->{input});
        my $sam = $me->Samtools(sam => \@input_list);
    } elsif (-r $me->{input} and $me->{input} =~ /\.fastq$/) {
        if (-d "bowtie_out") {
            find({ wanted => sub { push(@input_list, "bowtie_out/$_") if ($_ =~ /\.sam$/); }, follow => 1 }, 'bowtie_out/');
            my $sam = $me->Samtools(sam => \@input_list);
        } else {
            foreach my $k (%{$me->{bt_args}}) {
                my $output_string = "bowtie_out/${basename}-${k}.sam";
                push(@input_list, $output_string);
            }
            my $bt = $me->Bowtie();
            my $sam = $me->Samtools(depends => $bt->{pbs_id}, sam => \@input_list);
        }
    } else {
        die("I don't know what to do without a .fastq file nor a .sam file.\n");
    }
}

=head2
    Samtools()
=cut
sub Samtools {
    my $me = shift;
    my %args = @_;
    my $basename = $me->{basename};
    my $depends = "";
    if ($args{depends}) {
        $depends = $args{depends};
    }
    my $input = $args{input};
    my $output = $input;
    $output =~ s/\.sam$/\.bam/g;
    my $sorted = $input;
    $sorted =~ s/\.sam$//g;
    $sorted = qq"${sorted}-sorted";
    print "Converting to a compressed/sorted bam file.\n";
    my $job_string = qq!samtools view -u -t $me->{libdir}/genome/$me->{species}.fasta -S $input 1>$output
samtools sort -l 9 $output $sorted
rm $output && rm $input && mv ${sorted}.bam $output && samtools index $output
bamtools stats -in $output 2>${output}.stats 1>&2
!;
    my $comment = qq!## Converting the text sam to a compressed, sorted, indexed bamfile.
## Also printing alignment statistics to ${output}.stats!;
    my $samtools = $me->Qsub(job_name => "sam",
                             depends => $depends,
                             job_string => $job_string,
                             input => $input,
                             output => $output,
                             comment => $comment,
			     prescript => $args{prescript},
			     postscript => $args{postscript},
        );
    return($samtools);
}

=head2
    TopHat()
=cut
sub TopHat {
    my $me = shift;
    my %args = @_;
    $me->Check_Options(["species"]);
    my $depends = "";
    if ($args{depends}) {
        $depends = $args{depends};
    }
    my $tophat_args = ' -g 1 ';
    if ($args{tophat_args}) {
        $tophat_args = $args{tophat_args};
    }
    my $tophat_dir = 'tophat_out';
    if ($args{tophat_dir}) {
        $tophat_dir = $args{tophat_dir};
    }
    my $basename = $me->{basename};
    my $bt_reflib = "$me->{libdir}/genome/$me->{species}";
    my $bt_reftest = qq"${bt_reflib}.1.bt2";
    if (!-r $bt_reftest) {
        $me->BT2_Index();
    }
    if (!-r qq"$me->{libdir}/genome/$me->{species}.gtf") {
        print "Missing the gtf file for $me->{species}\n";
        Gff2Gtf("$me->{libdir}/genome/$me->{species}.gff");
    }
    my $inputs = $me->{input};
    my @in = split(/\:/, $inputs);
    $inputs =~ s/\:/ /g;

    my $job_string = qq!mkdir -p ${tophat_dir} && tophat ${tophat_args} --b2-very-sensitive -p 1 -o ${tophat_dir} \\
-G $me->{libdir}/genome/$me->{species}.gtf --no-novel-juncs \\
  $me->{libdir}/genome/$me->{species} \\
  $inputs && samtools index ${tophat_dir}/accepted_hits.bam
!;
    my $comment = qq!## I still have no clue what I am doing when I use tophat...
## However, I know that -g 1 will allow only 1 hit in the case of multihits, but randomly place it
## From the manual:  "If there are more alignments with the same score than this
## number, TopHat will randomly report only this many alignments"
## -N 1 will discard anything with >1 mismatch (default is 2)
## -r adjusts the allowable mean distance between the paired reads
## --mate-std-dev sets the deviation of -r
## --microexon-search will tell it to search short introns for reads >=50
!;
    my $tophat = $me->Qsub(job_name => "th",
                           qsub_cpus => 4,
                           qsub_queue => "long",
                           qsub_wall => "144:00:00",
                           qsub_mem => 10,
                           depends => $depends,
                           job_string => $job_string,
                           comment => $comment,
			   prescript => $args{prescript},
			   postscript => $args{postscript},
        );

    my $accepted = "accepted_hits.bam";
    $accepted = $args{accepted_hits} if ($args{accepted_hits});
    my $count_table = "accepted_hits.count";
    $count_table = $args{count_table} if ($args{count_table});

    my $htmulti = $me->HT_Multi(depends => $tophat->{pbs_id},
                                input => $accepted,
                                job_name => qq"$in[0]_count",
                                ##output => $count_table,
        );

    return({tophat => $tophat, htseq => $htmulti,});
}

=head2
    Trimomatic()
=cut
sub Trimomatic {
    my $me = shift;
    my %args = @_;
    my $trim;
    if ($me->{input} =~ /\:|\,/) {
	$trim = $me->Trimomatic_Pairwise(@_);
    } else {
	$trim = $me->Trimomatic_Single(@_);
    }
    return($trim);
}

=head2
    Trimomatic_Pairwise()
=cut
sub Trimomatic_Pairwise {
    my $me = shift;
    my %args = @_;
    my $input = $me->{input};
    $input = $args{input} if ($args{input} and !$input);
    my @input_list = split(/\:|\,/, $input);
    if (scalar(@input_list) <= 1) {
        my $ret = $me->Trimomatic_Single(input => $input);
        return($ret);
    }
    my $r1 = $input_list[0];
    my $r2 = $input_list[1];
    my $basename = basename($r1, (".gz"));
    my $r1b = basename($r1, (".gz"));    
    my $r2b = basename($r2, (".gz"));
    $basename = basename($r1, (".fastq"));
    $r1b = basename($r1, (".fastq"));    
    $r1b = basename($r2, (".fastq"));

    my $r1o = qq"${r1b}-trimmed.fastq";
    my $r2o = qq"${r2b}-trimmed.fastq";
    my $r1op = qq"${r1b}-trimmed_paired.fastq";
    my $r1ou = qq"${r1b}-trimmed_unpaired.fastq";
    my $r2op = qq"${r2b}-trimmed_paired.fastq";
    my $r2ou = qq"${r2b}-trimmed_unpaired.fastq";

    $basename =~ s/\_R1$//g;
    $me->{basename} = $basename;
    my $output = qq"${r1o}:${r2o}";
    my $comment = qq!## This call to trimmomatic removes illumina and epicentre adapters from ${input}.
## It also performs a sliding window removal of anything with quality <25;
## cutadapt provides an alternative to this tool.
## The original sequence data is recompressed and saved in the sequences/ directory.!;
    my $job_string = qq!
mkdir -p sequences
## In case a trimming needs to be redone...
if [[ \! -r "${r1}" ]]; then
  if [[ -r "sequences/${r1}.xz" ]]; then
    mv sequences/${r1}.xz . && pxz -d ${r1}.xz && mv sequences/${r2}.xz . && pxz -d ${r2}.xz
  else
    echo "Missing files. Did not find ${r1} nor sequences/${r1}.xz"
    exit 1
  fi
fi
trimmomatic PE -threads 1 -phred33 ${r1} ${r2} ${r1op} ${r1ou} ${r2op} ${r2ou} ILLUMINACLIP:$me->{libdir}/adapters.fa:2:20:4 SLIDINGWINDOW:4:25 2>outputs/${basename}-trimomatic.err 1>outputs/${basename}-trimomatic.out
excepted=\$(grep "Exception" outputs/${basename}-trimomatic.err)
## The following is in case the illumina clipping fails, which it does if this has already been run I think.
if [[ "\${excepted}" \!= "" ]]; then
  trimmomatic PE -threads 1 -phred33 ${r1} ${r2} ${r1op} ${r1ou} ${r2op} ${r2ou} SLIDINGWINDOW:4:25 1>>outputs/${basename}-trimomatic.out 2>>output/${basename}-trimomatic.err
fi
sleep 10
mv ${r1op} ${r1o} && mv ${r2op} ${r2o}
!;
    ## Example output from trimomatic:
    ## Input Reads: 35327530 Surviving: 34441992 (97.49%) Dropped: 885538 (2.51%)
    ## Perhaps I can pass this along to Get_Stats()
    my $trim = $me->Qsub(job_name => "trim",
			 qsub_wall => "4:00:00",
			 job_string => $job_string,
			 input => $input,
			 output => $output,
			 comment => $comment,
			 prescript => $args{prescript},
			 postscript => $args{postscript},
        );
    my $comp = $me->Recompress(job_name => "",
			       depends => $trim->{pbs_id},
			       input => $input,
			       output => qq"sequences/${r1}.xz",
			       output2 => qq"sequences/${r2}.xz",
			       comment => qq"## The original sequence file is in sequences/",
        );
    ## Set the input for future tools to the output from trimming.
    $me->{input} = $output;
    return($trim);
}

=head2
    Trimomatic_Single()
=cut
sub Trimomatic_Single {
    my $me = shift;
    my %args = @_;
    my $input = $me->{input};
    my $basename = $me->{basename};
    $basename = basename($basename, (".gz"));
    $basename = basename($basename, (".fastq"));    
    my $output = qq"${basename}-trimmed.fastq";
    my $comment = qq!## This call to trimomatic removes illumina and epicentre adapters from ${input}.
## It also performs a sliding window removal of anything with quality <25;
## cutadapt provides an alternative to this tool.
## The original sequence data is recompressed and saved in the sequences/ directory.!;
    my $job_string = qq!
mkdir -p sequences
## In case a trimming needs to be redone...
if [[ \! -r "${input}" ]]; then
  if [[ -r "sequences/${input}.xz" ]]; then
    mv sequences/${input}.xz . && pxz -d ${input}.xz
  else
    echo "Missing files. Did not find ${input} nor sequences/${input}.xz"
    exit 1
  fi
fi
trimmomatic SE -phred33 ${input} ${output} ILLUMINACLIP:$me->{libdir}/adapters.fa:2:20:4 SLIDINGWINDOW:4:25 1>outputs/${basename}-trimomatic.out 2>output/${basename}-trimomatic.err
!;
    my $trim = $me->Qsub(job_name => "trim",
			 qsub_wall => "4:00:00",
			 job_string => $job_string,
			 input => $input,
			 output => $output,
			 comment => $comment,
			 prescript => $args{prescript},
			 postscript => $args{postscript},
	);
    my $comp = $me->Recompress(job_name => "",
			       depends => $trim->{pbs_id},
			       input => $input,
			       output => qq"sequences/${input}.xz",
			       comment => qq"## The original sequence file is in sequences/",
        );
    ## Set the input for future tools to the output from this trimming operation.
    $me->{input} = $output;
    return($trim);
}

=head2
    Finishedp_Rawseq()
=cut
sub Finishedp_Rawseq {
    my $me = shift;
    my %args = @_;
    if (-r qq"$me->{hpgl}_forward.fastq" && -r qq"$me->{hpgl}_reverse.fastq") {
	return(1);
    } else {
	return(0);
    }
}

=head2
    Finishedp_Fastqc()
=cut
sub Finishedp_Fastqc {
    my $me = shift;
    my %args = @_;
    if (-d qq"$me->{hpgl}-unfiltered_fastqc" or -d qq"$me->{hpgl}-filtered_fastqc") {
	return(1);
    } else {
	return(0);
    }
}

=head2
    Finishedp_Trimomatic()
=cut
sub Finishedp_Trimomatic {
    my $me = shift;
    my %args = @_;
    if (-r qq"$me->{basename}-trimmed.fastq" or (-r qq"$me->{basename}_forward-trimmed.fastq" and -r qq"$me->{basename}_reverse-trimmed.fastq")) {
	return(1);
    } else {
	return(0);
    }
}

=head2
    Finishedp_Graphs()
=cut
sub Finishedp_Graphs {
    my $me = shift;
    my %args = @_;
    if (-r qq"$me->{basename}_ntdist.svg") {
	return(1);
    } else {
	return(0);
    }
}

=head2
    Finishedp_Tophat()
=cut
sub Finishedp_Tophat {
    my $me = shift;
    my %args = @_;
    my $subdir = ".";
    $subdir = $args{species} if ($args{species});
    if (-r qq"$subdir/tophat_out/accepted_hits.bam") {
	return(1);
    } else {
	return(0);
    }
}

=head2
    Finishedp_HTseq()
=cut
sub Finishedp_Htseq {
    my $me = shift;
    my %args = @_;
    my $count_table = qq"$me->{hpgl}.count";
    if ($args{species}) {
	$count_table = qq"$me->{hpgl}_$args{species}.count";
    }
    if (-r $count_table) {
	return(1);
    } else {
	return(0);
    }
}

=head2
    Qsub()
=cut
sub Qsub {
    my $me = shift;
    my %args = @_;
    my $queue = $args{queue};
    if ($args{depends}) {
        print "This job depends on $args{depends}\n";
    }
    
    my $job_name;
    my $name_suffix = substr($me->{basename}, 0, 8);
    if ($args{job_name}) {
        my $name_prefix = $args{job_name};
        $job_name = qq"${name_prefix}-${name_suffix}";
    } else {
        $job_name = $name_suffix;
    }
    my $qsub_args = $me->{qsub_args};
    my $qsub_queue = $me->{qsub_queue};
    $qsub_queue = $args{qsub_queue} if ($args{qsub_queue});
    my $qsub_shell = $me->{qsub_shell};
    $qsub_shell = $args{qsub_shell} if ($args{qsub_shell});
    my $qsub_mem = $me->{qsub_mem};
    $qsub_mem = $args{qsub_mem} if ($args{qsub_mem});
    my $qsub_wall = $me->{qsub_wall};
    $qsub_wall = $args{qsub_wall} if ($args{qsub_wall});
    my $qsub_cpus = $me->{qsub_cpus};
    $qsub_cpus = $args{qsub_cpus} if ($args{qsub_cpus});
    my $qsub_logdir = $me->{qsub_logdir};
    $qsub_logdir = $args{qsub_logdir} if ($args{qsub_logdir});
    my $qsub_log = qq"${qsub_logdir}/${job_name}.qsubout";
    my $depends_string = $me->{qsub_depends};
    $depends_string .= $args{depends};
    my $job_output = "";
    $job_output = $args{output} if ($args{output});
    my $job_input = "";
    $job_input = $args{input} if ($args{input});
    my $script_file = qq"$me->{basedir}/scripts/${job_name}.sh";
    my $mycwd = cwd();
    my $script_start = qq!#PBS -V -S $qsub_shell -q $qsub_queue -d $me->{basedir}
#PBS -N $job_name -l mem=${qsub_mem}gb -l walltime=$qsub_wall -l ncpus=$qsub_cpus
#PBS -o $qsub_log $qsub_args
echo "####Started $script_file at \$(date)" >> outputs/log.txt
source /cbcb/lab/nelsayed/scripts/dotbashrc
cd $me->{basedir} || exit
!;
    my $script_end = qq!## The following lines give status codes and some logging
echo \$? > status/${job_name}.status
echo "####Finished $script_file at \$(date), it took \$(( \$SECONDS / 60 )) minutes." >> outputs/log.txt
echo "####This job consisted of the following:" >> outputs/log.txt
cat "\$0" >> outputs/log.txt
!;
    system("cd $mycwd && mkdir -p status outputs scripts sequences");
    print "The job is:
$args{job_string}" if ($me->{debug});

    my $total_script_string = "";
    $total_script_string .= "$script_start\n";
    $total_script_string .= "$args{comment}\n" if ($args{comment});
    $total_script_string .= "$args{prescript}\n" if ($args{prescript});
    $total_script_string .= "$args{job_string}\n" if ($args{job_string});
    $total_script_string .= "if [ \$? == \"0\" ]; then\n";
    if ($args{postscript}) {
        $total_script_string .= qq!if [ \$? == "0" ]; then
   $args{postscript}
fi
!;
    }
    $total_script_string .= "$script_end\n";

    open(SCRIPT, ">$script_file");
    print SCRIPT $total_script_string;
    close(SCRIPT);

    my $qsub = qq"bash $script_file";
    if ($me->{pbs}) {
        $qsub = qq"qsub -W $depends_string $script_file";
    }
    my $job_id;
    ##sleep(1);
    my $qsub_pid = open(my $fh, "-|", $qsub);
    while(my $line = <$fh>) {
      chomp($line);
      $job_id = $line;
    } 
    close($fh);
    my $jobid_name = qq"$me->{basename}-$args{job_name}";
    print "Starting a new job: $jobid_name\nEOJ\n\n";
    my $job = new HPGL::Job(id => $jobid_name,
                            submitter => $qsub,
                            mem => $qsub_mem,
                            walltime => $qsub_wall,
                            cpus => $qsub_cpus,
                            jobname => $job_name,
                            log => $qsub_log,
                            depends_string => $depends_string,
                            queue => $qsub_queue,
                            qsub_args => $qsub_args,
                            basedir => $me->{basedir},
                            pbs_id => $job_id,
                            script_file => $script_file,
                            script_start => $script_start,
                            script_body => $args{job_string},
                            output => $job_output,
                            input => $job_input,
        );
    sleep(1);
    return($job);
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
    Dump_Reads()
=cut
#sub Dump_Reads {
#my %options;
#my $casava_passed = 0;
#my $line_count = 0;
#my $connection;
#main();
#
#sub main {
#    %options = Parse_Options();
#    $connection = Connect();
#    Download_Files();
#    my $seq = $line_count / 4;
#    print "$seq sequences were downloaded, of which $casava_passed passed casava's filtering.\n";    
#}
#
#sub Connect {
#    my $connection = new Net::Amazon::S3({aws_access_key_id => $options{access_key},
#                                    aws_secret_access_key => $options{secret_key},
#                                    host => 'gembox.cbcb.umd.edu',
#                                    secure => 1,
#                                    retry => 0,});
#    die("The Amazon::S3 connection failed.") if (!defined($connection));
#    return($connection);
#}
#
#sub Download_Files {
#    my $fwd = "$options{bucket}_forward.fastq";
#    my $rev = "$options{bucket}_reverse.fastq";
#    if ($options{orientation} eq 'both') {
#        Perform_Download('forward', $fwd);
#        Perform_Download('reverse', $rev);
#    } elsif ($options{orientation} eq 'forward') {
#        Perform_Download('forward', $fwd);
#    } elsif ($options{orientation} eq 'reverse') {
#        Perform_Download('reverse', $rev);
#    } else {
#        die("I can only download forward reads, reverse reads, or both.");
#    }
#}
#
#sub Perform_Download {
#    my $orientation = shift;
#    my $dest = shift;
#    if (-r $dest) {
#        print "$dest already exists, deleting it.\n" if ($options{verbose});
#        unlink($dest);
#    }
#
#    my $response = $connection->buckets;
#    my @buckets = @{$response->{buckets}};
#    my $c = 0;
#    my $found_bucket = 0;    
#  BUCKETS: foreach my $bucket (@buckets) {
#      my $bucket_name = $bucket->bucket;
#      if ($bucket_name eq $options{bucket}) {
#          $found_bucket++;
#          Download_Bucket($bucket, $dest, $orientation);
#
#          if ($options{multi}) {
#            INNER: foreach my $buck (@buckets) {
#              LETTERS: foreach my $letter ("a".."z") {
#                  my $test = join("", $options{bucket}, $letter);
#                  my $new_name = $buck->bucket();
#                  print "Now testing for $new_name and $test\n";
#                  if ($new_name eq $test) {
#                      Download_Bucket($buck, $dest, $orientation);
#                      next LETTERS;
#                  } else {
#                      next INNER;
#                  }
#              } ## End rechecking the list of buckets for xxxa -> xxxz
#            } ## End making the names xxxa -> xxxz
#          } ## End if multi-buckets are on
#      } else { ## End if we find hpgl0xxx
#          next BUCKETS;
#      }
#  }
#    if ($found_bucket == 0) {
#        die("Did not find the bucket: $options{bucket}");
#    }
#}
#
#sub Download_Bucket {
#    my $bucket = shift;
#    my $dest = shift;
#    my $orientation = shift;
#    open(OUT, ">>$dest");
#    my $test = $bucket->list_all or die $connection->err . ": " . $connection->errstr;
#    my @keys = @{$bucket->list_all->{keys} || []};
#    ## The following for loop creates a list of files to download from ceph.
#    ## The logic therein is very specific to Najib's data and should be changed for other
#    ## applications.
#    my @files_to_download = ();
#    foreach my $key (@keys) {
#        my $keyname = $key->{key};
#        next if ($keyname =~ /\.md5$/);
#        if ($orientation eq 'forward') {
#            ## The difference between forward/reverse is R1 vs R2
#            if ($keyname =~ /.*_[A|T|G|C]+_L\d+_R1_\d+\.fastq/) {
#                push(@files_to_download, $keyname);
#            }
#        } else {
#            if ($keyname =~ /.*_[A|T|G|C]+_L\d+_R2_\d+\.fastq/) {                
#                push(@files_to_download, $keyname);
#            }
#        }
#    }
#    ## Now that we have a list of files to download, do it!
#    foreach my $file (@files_to_download) {
#        my $tmp_dest = basename($file);
#        if (-r $tmp_dest) {
#            print "$tmp_dest already exists, deleting it.\n" if ($options{verbose});
#            unlink($tmp_dest);
#        }
#        print "Downloading $file.\n";
#        $bucket->get_key_filename($file, undef, $tmp_dest);
#        my $md5_file = "$file.md5";
#        my $md5_dest = "$tmp_dest.md5";
#        if (-r $md5_dest) {
#            print "$md5_dest already exists, deleting it.\n" if ($options{verbose});
#            unlink($md5_dest);
#        }
#        $bucket->get_key_filename($md5_file, undef, $md5_dest);
#        print "Checking md5 Sum\n";                
#        my $m = Check_MD5($tmp_dest, $md5_dest);
#        my @suffixes = (".gz",".bz2",".xz");
#        my $inter_dest = basename($file, @suffixes);
#        if ($file =~ /\.fastq$/) {  ## Then it wasn't compressed
#            print "The file was not compressed.\n";
#        } else {
#            if (-r $inter_dest) {
#                print "$inter_dest already exists, deleting it.\n" if ($options{verbose});
#                unlink($inter_dest);
#            }
#            print "Extracting $tmp_dest.\n";
#            my $ae = new Archive::Extract(archive => $tmp_dest);
#            my $ok = $ae->extract(to => $inter_dest);
#        }
#        print "Concatenating $inter_dest onto $dest.\n" if ($options{verbose});
#        open(IN, "<:mmap", $inter_dest);
#        while(<IN>) {
#            $line_count++;
#            $casava_passed++ if ($_ =~ m/^@.* [^:]*:N:[^:]*:/);
#            print OUT;
#        }
#        close(IN);
#        print "Cleaning up $inter_dest and $tmp_dest.\n" if ($options{verbose});
#        unlink($inter_dest) if (-r $inter_dest);
#        unlink($tmp_dest) if (-r $tmp_dest);
#    }
#    close(OUT);    
#}
#
#sub List {
#    if (!$options{bucket}) {
#        List_Buckets();
#    } else {
#        my $response = $connection->buckets;
#        my @buckets = @{$response->{buckets}};        
#        my $found_bucket = 0;
#        foreach my $bucket (@buckets) {
#            my $name = $bucket->bucket();
#            if ($name eq $options{bucket}) {
#                $found_bucket++;
#            } else {
#                next;
#            }
#            my @keys = @{$bucket->list_all->{keys} || []};
#            my @files_to_download = ();
#            foreach my $key (@keys) {
#                print "$key->{key}\t$key->{size}\t$key->{last_modified}\n";
#            }
#        }
#        if ($found_bucket == 0) {
#            print "Did not find $options{bucket}\n";
#            List_Buckets();
#        }
#    }
#}
#
#sub List_Buckets {
#    print "Available buckets are:\n";
#    my $response = $connection->buckets;
#    my @buckets = @{$response->{buckets}};    
#    foreach my $bucket (@buckets) {
#        my $name = $bucket->bucket();
#        print "$name\n";
#    }
#}
#
#sub Check_MD5 {
#    my $downloaded_file = shift;
#    my $md5_file = shift;
#    open(MD5_IN, "<$md5_file");
#    my $original_md5 = "";
#    while(my $line = <MD5_IN>) {
#        chomp($line);
#        my ($md5, $filename) = split(/\s+/, $line);
#        $original_md5 = $md5;
#    }
#    close(MD5_IN);
#    unlink($md5_file);
#    my $checksum = new Digest::MD5;
#    open(MD5_CHECK, "<$downloaded_file");
#    my $chk = \*MD5_CHECK;
#    $chk = \*MD5_CHECK;  ## I don't understand how this works...
#    $checksum->addfile($chk);
#    my $new_md5_checksum = $checksum->hexdigest;
#    my ($new_md5, $new_filename) = split(/\s+/, $new_md5_checksum);
#    if ($new_md5 eq $original_md5) {
#        return(1);
#    } else {
#        print "The checksums no longer agree.  This file has changed since it was uploaded,
#or the download failed and should be reperformed.\n";
#        return(undef);
#    }
#}
#
#sub Parse_Options {
#    my %options = (
#        bucket => undef,
#        orientation => 'both',
#        debug => 0,
#        md5 => 0,
#        sample => 0,
#        help => 0,
#        list => 0,
#        file => undef,
#        verbose => 1,
#        multi => 0,
#        access_key => 'PQOD56HDPQOXJPYQYHH9',
#        secret_key => 'kNl27xjeVSnChzEww9ziq1VkgUMNmNonNYjxWkGw',
#        );
#    my $opt = GetOptions(
#        "bucket|b:s" => \$options{bucket},
#        "orientation|o:s" => \$options{orientation},
#        "debug|d" => \$options{debug},
#        "md5|m" => \$options{md5},
#        "sample|s" => \$options{sample},
#        "help|h" => \$options{help},
#        "list|l" => \$options{list},
#        "file|f:s" => \$options{file},
#        "verbose|v:s" => \$options{verbose},
#        "multi|m:s" => \$options{multi},
#        );
#    
#    if ($options{help}) {
#        print "This script is intended to make downloading and viewing reads from ceph easier.
#An example invocation might be:
#> dump_reads.pl -b hpgl0223 -o r
#That will dump the reverse (-o r) reads from the 'bucket' hpgl0223 (forward is default)
#> dump_reads.pl -l
#This will list the available buckets.
#";
#        exit(0);
#    }
#    if ($options{orientation} =~ /^r|^R/) {
#        $options{orientation} = 'reverse';
#    } elsif ($options{orientation} =~ /^f|^F/) {
#        $options{orientation} = 'forward';
#    } else {
#        $options{orientation} = 'both';
#    }
#    if ($options{list}) {
#        List();
#        exit(0);
#    }
#
#    ## First ensure that we have all the information required to download sequence files.
#    my $term = new Term::ReadLine('>');
#    my $attribs = $term->Attribs;
#    $attribs->{completion_suppress_append} = 1;
#    my $OUT = $term->OUT || \*STDOUT;
#    if (!defined($options{bucket})) {
#        $options{bucket} = $term->readline("Please provide a hpgl to download (lowercase):");
#        $options{bucket} =~ s/\s+$//g;
#    }
#    return(%options);
#}
#
#}

1;
