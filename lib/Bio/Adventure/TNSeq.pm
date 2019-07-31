package Bio::Adventure::TNSeq;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Cwd qw(abs_path getcwd);
use File::Basename;
use FileHandle;
use File::Path qw"make_path";
use PerlIO;
use PerlIO::gzip;
use String::Approx qw"amatch";

=head1 Name

Bio::Adventure::TNSeq - some functions for working with TNSeq data.

=head1 SYNOPSIS

The functions here were copy/pasted from some other code pieces I
wrote while working on mgas 5448 tnseq data.  They have not been
incorporated into the general Bio::Adventure.pm scheme yet.

=head1 METHODS

=head2 C<Sort_Indexes>

$hpgl->Do_Sort_Indexes(inputdir => '.');
will assume there is a $(pwd)/indexes.txt with a list of
sequencing index to sample names, it will then read all the fastq
files in the inputdir and sort them into piles by index.  Most
tasks/tools of course do not/should not need this, but the tnseq
libraries used by the McIver lab are not quite standard.

=cut
sub Do_Sort_Indexes {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['index_file', 'input'],
    );
    my $input = $options->{input};
    $input = '.' unless (defined($input));
    my $outdir = $options->{output_dir};
    $outdir = 'outputs' unless (defined($outdir));
    $options = $class->Set_Vars(outdir => $outdir);

    if ($options->{help}) {
        print qq"This script is intended to split TNSeq reads as indexed by Yoann Le Breton.
It therefore makes some assumptions about the beginnings of reads and requires two arguments:
 --inputdir or -d <dir>:  A directory containing the raw sequences (compressed or not)
 --indexes or -i <file>:  A file containing the sequence barcodes and sample names, one per line.
";
        exit(0);
    }

    my $index_hash = Bio::Adventure::TNSeq::Read_Indexes(
        $class,
        index_file => $options->{index_file},
    );

    my $reads;
    if (-f $options->{input}) {
        $reads = Bio::Adventure::TNSeq::Sort_TNSeq_File_Approx(
            $class,
            file => $options->{input},
            index_hash => $index_hash,
        );
    } elsif (-d $options->{input}) {
        $reads = Bio::Adventure::TNSeq::Sort_TNSeq_Dir(
            $class,
            dir => $options->{input},
            index_hash => $options->{index_hash},
        );
    } else {
        die("I need a directory containing some sequence.") unless($options->{input});
    }

    my $log = new FileHandle(">$options->{outdir}/tnseq_sorting_stats.txt");
    my $message = qq"
The total number of reads read in $options->{input}: $reads->{total}->{read}.
The total number of reads sorted in $options->{input}: $reads->{total}->{written}.
";
    print $log $message;
    foreach my $k (keys %{$reads}) {
        next if ($k eq 'total');
        $message = qq"Index $k had $reads->{$k}->{written} reads written from $options->{input} to $reads->{$k}->{filename}.\n";
        print $log $message;
        $reads->{$k}->{handle}->close();
    }
    $log->close();
    return($reads->{total});
}

=head2 C<TA_Check>

Check an alignment of tnseq data for terminal TAs to ensure that the putative
insertion positions are actually mariner specific, e.g. they are 'TA'.

=cut
sub TA_Check {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
    );
    my $job_basename = $class->Get_Job_Name();
    my $jstring = qq"
use Bio::Adventure;
use Bio::Adventure::TNSeq;
my \$ret = Bio::Adventure::TNSeq::Do_TA_Check(
  \$h,
  input => '$options->{input}',
";
    foreach my $option (sort keys %args) {
        if ($args{$option}) {
            $jstring .= qq"  ${option} => '$args{$option}',\n";
        }
    }
    $jstring .= qq");";
    my $sort_job = $class->Submit(
        comment => "# Check for tailing TAs!",
        cpus => 1,
        depends => $options->{depends},
        jname => qq"tacheck_${job_basename}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        language => 'perl',
        mem => 8,
        queue => "throughput",
        walltime => "10:00:00",
    );
    return($sort_job);
}

=head2 C<Do_TA_Check>

Actually invoke the check for TAs.

=cut
sub Do_TA_Check {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $input = $options->{input};
    my $input_base = basename($input, ('.gz', '.xz', '.bz2'));
    $input_base = basename($input_base, ('.fastq', '.fasta'));
    my $with_ta = FileHandle->new("| gzip -f -9 > ${input_base}_ta.fastq.gz");
    my $without_ta = FileHandle->new("| gzip -f -9 > ${input_base}_nota.fastq.gz");
    my $inputted = FileHandle->new("less ${input} |"); ## PerlIO::gzip is failing on some files.
    my $in = new Bio::SeqIO(-fh => $inputted, -format => 'Fastq');
    my $count = 0;
    my $with_ta_count = 0;
    my $without_ta_count = 0;
  READS: while (my $in_seq = $in->next_dataset()) {
        $count++;
        my $id = $in_seq->{'-descriptor'};
        my $sequence = $in_seq->{'-seq'};
        my $qual = $in_seq->{'-raw_quality'};
        if ($sequence =~ /TA$/) {
            $with_ta_count++;
            my $fastq_string = qq"\@${id}
${sequence}
+
${qual}
";
            print $with_ta $fastq_string;
        } else {
            $without_ta_count++;
            my $fastq_string = qq"\@${id}
${sequence}
+
${qual}
";
            print $without_ta $fastq_string;
        }
    }
    my $counters = {
        total => $count,
        with => $with_ta_count,
        without => $without_ta_count,
    };
    $input->close();
    $inputted->close();
    return($counters);
}

=head2 C<Sort_Indexes>

Sort a pile of tnseq data for the set of known indexes.

=cut
sub Sort_Indexes {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    ## If options are required, feed them back into %args here?
    my $jstring = qq"
use Bio::Adventure;
my \$ret = Bio::Adventure::TNSeq::Do_Sort_Indexes(
  \$h,\n";
    foreach my $option (keys %args) {
        next unless defined($args{$options});
        $jstring .= qq"  ${option} => '$args{$option}',\n";
    }
    $jstring .= qq");";
    my $sort_job = $class->Submit(
        comment => "# Sort those indexes!",
        cpus => 1,
        language => 'perl',
        jname => "sort_indexes",
        jstring => $jstring,
        mem => 8,
        queue => "workstation",
        walltime => "60:00:00",
    );
    return($sort_job);
}

=head2 C<Sort_TNSeq_File_Approx>

hpgl->Sort_TNSeq_File(file => 'filename.fastq.gz', indexes =>*Index Hash*)
is called by Sort_Indexes() when a single file is to be demultiplexed.
It assumes the input file is gzipped and will use the PerlIO::gzip
filter to uncompress the data.  It will then extract the sequences
per the rules set up in the index hash.

=cut
sub Sort_TNSeq_File_Approx {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        index_hash => {},
    );
    my $data = $options->{index_hash};

    return(undef) unless ($options->{input} =~ /\.fastq/);
    my $out = FileHandle->new(">$options->{outdir}/tnseq_sorting_out.txt");
    ## Does FileHandle work here?
    my $inputted = FileHandle->new("less $options->{input} 2>/dev/null |");
    my $in = Bio::SeqIO->new(-fh => $inputted, -format => 'Fastq');
    my $count = 0;
  READS: while (my $in_seq = $in->next_dataset()) {
        $count++;
        $data->{total}->{read}++;
        print $out "Read $data->{total}->{read} entries.\n" if ($options->{debug});
        my $id = $in_seq->{'-descriptor'};
        my $sequence = $in_seq->{'-seq'};
        my $qual = $in_seq->{'-raw_quality'};
        my $found = 0;
        my @index_list = ();
        foreach my $index (keys %{$data}) {
            push(@index_list, $index) unless ($index eq 'total' or $index eq 'unknown' or $index eq 'ambiguous');
        }
        my $found_id = '';
        my $match_substring = substr($sequence, 0, 9);
        foreach my $index (@index_list) {
            my @matches = amatch($index, [ 'S1' ], ($match_substring));
            if (scalar(@matches) == 1) {
                $found++;
                $found_id = $index;
            }
        }                       ## End each index
        if ($found == 1) {
            my $new_sequence = substr($sequence, 9,);
            my $new_quality = substr($qual, 9,);
            $data->{total}->{written}++;
            if (($data->{total}->{written} % 10000) == 0) {
                print $out "Wrote $data->{total}->{written} reads from $options->{input}\n";
            }
            my $fastq_string = qq"\@${id}
${new_sequence}
+
${new_quality}
";
            my $handle = $data->{$found_id}->{handle};
            print $handle $fastq_string;
            $data->{$found_id}->{written}++;
        } elsif ($found > 1) {
            my $handle = $data->{ambiguous}->{handle};
            my $fastq_string = qq"\@${id}
${sequence}
+
${qual}
";
            print $handle $fastq_string;
            $data->{ambiguous}->{written}++;
        } elsif ($found < 1) {
            my $handle = $data->{unknown}->{handle};
            my $fastq_string = qq"\@${id}
${sequence}
+
${qual}
";
            print $handle $fastq_string;
            $data->{unknown}->{written}++;
        }
    } ## End reading the fastq file
    $out->close();
    return($data);
}

=head2 C<Sort_TNSeq_Dir>

Called by Sort_Indexes() when a directory of fastq files are to
be demultiplexed. It uses File::Find to find all the .fastq.gz files
in the target directory (recursively) which do not have the name of
the output directory in the filename.
It then passes each file to Sort_TNSeq_File()

=cut
sub Sort_TNSeq_Dir {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        index_hash => {},
    );
    my $cwd_dir = getcwd();
    my $searchdir = qq"$args{dir}";
    my $files = 0;
    unless ($searchdir =~ /^\//) {
        $searchdir = qq"${cwd_dir}/${searchdir}";
    }
    ## print "Searching: $searchdir  for files to read.\n";
    my @directory = ($searchdir);
    my @file_list = ();
    find(sub { push(@file_list, $File::Find::name) if ($File::Find::name =~ /\.fastq\.gz/ and
                                                           $file::Find::name !~ /$options->{outdir}/); }, @directory);
    my @approxes = ();
    foreach my $file (@file_list) {
        $files = $files++;
        next if ($file =~ /$options->{outdir}/);
        ## This might be incorrect, CHECKME!
        $file =~ s/\/\.\//\//;
        my $approx = Bio::Adventure::TNSeq::Sort_TNSeq_File_Approx(
            $class,
            file => $file,
            index_hash => $options->{index_hash},
        );
        push(@approxes, $approx);
    }
    return(@approxes);
}

=head2 C<Read_Indexes>

Reads a file containing the indexes and sets up the data structure
used to keep count of the reads demultiplexed along with the set
of file handles to which to write the data.

=cut
sub Read_Indexes {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    make_path($class->{outdir}) if (!-d $options->{outdir});
    unlink("$options->{outdir}/unknown.fastq.gz") if (-r "$options->{outdir}/unknown.fastq.gz");
    unlink("$options->{outdir}/ambiguous.fastq.gz") if (-r "$options->{outdir}/ambiguous.fastq.gz");
    my $unknown = FileHandle->new("| gzip -9 >> $options->{outdir}/unknown.fastq.gz");
    my $ambiguous = FileHandle->new("| gzip -9 >> $options->{outdir}/ambiguous.fastq.gz");
    my $indexes = {
        total => {read => 0, written => 0,
              },
        unknown => {name => 'unknown',
                    filename => 'uknown.fastq',
                    handle => $unknown,
                    written => 0,
                },
        ambiguous => {name => 'ambiguous',
                      filename => 'ambiguous.fastq',
                      handle => $ambiguous,
                      written => 0,
                  },
    };

    my $index_file = FileHandle->new("<$args{index_file}");
    while (my $line = <$index_file>) {
        chomp $line;
        next if ($line =~ /^#/);
        next unless ($line =~ /^A|T|G|C|a|t|g|c/);
        my ($ind, $sample) = split(/\s+|\,|;/, $line);
        ## Add a quick piece of logic in case I flip the indexes.txt file
        if ($ind =~ /^h|H/ && $sample =~ /^A|T|G|C/) {
            ($sample, $ind) = split(/\s+|\,|;/, $line);
        }
        my $output_filename = qq"$options->{outdir}/${sample}.fastq.gz";
        unlink($output_filename) if (-r $output_filename);
        my $handle = FileHandle->new("| gzip -9 >> $output_filename");
        $indexes->{$ind} = {
            name => $sample,
            filename => $output_filename,
            handle => $handle,
            written => 0,
        };
    }
    $index_file->close();
    return($indexes);
}

=head2 C<Essentiality_TAs>

This will read the $args{bam} bam file and look for reads which
start/end with TA in the orientation expected for a successful
insertion by the mariner retrotransposable element.

essentiality_tas
essentiality_tas -f genome.fasta -o output.txt -n genome_name -g genome.gff -b bamfile.bam

        Other options:  The following are possible
        --name  |-n : The name to print in the text file (required)
        --bam   |-b : The bam alignment (required)
        --fasta |-f : The fasta genome  (required)
        --output|-o : The output file   (required)
        --debug |-d : Print debugging information?
        --gff   |-g : The gff annotation (required)
        --kept_tas  : Keep the intermediate tas file (useful for coverage data)
        --type  |-t : Type of entry in the gff file to use (gene by default)
        --help      : Print this help information

=cut
sub Essentiality_TAs {
    my ($class, %args) = @_;
    my $output_directory = qq"outputs/essentiality";
    make_path($output_directory);
    ## A smarter way of handling counting inter-cds regions would be to use the Read_GFF() functions to pull the intercds regions rather than the hackish @inter_starts @inter_ends.
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input'],
        feature_type => 'exon',
    );
    my $input = $options->{input};
    my $checked_input = $class->Check_Input(files => $input);
    die("Essentiality_TAs requires a bam input, you gave: ${input}") unless ($input =~ /\.bam$/);
    my $input_name = basename($input, ('.bam'));
    ## This is a bit of a mess, I need to clean up how these arguments are passed.
    my $output = $options->{output};
    $output = qq"${input}_out" unless (defined($options->{output}));
    my $out = FileHandle->new(">>${output}.log");

    my $genome_fasta = qq"$options->{libdir}/$options->{libtype}/$options->{species}.fasta";
    my $genome_gff = $genome_fasta;
    $genome_gff =~ s/\.fasta/\.gff/;
    my $genome = Bio::Adventure::TNSeq::Allocate_Genome(
        $class,
        fasta => $genome_fasta,
    );
    my @data = @{$genome->{data}};
    my $length = scalar(@data);
    ## These 5 arrays will hold the lists of start/ends of each region of interest
    ## This is a somewhat dirty method for this, but I just wanted something easy
    my $cds_gff = $class->Read_Genome_GFF(
        feature_type => $options->{feature_type},
        gff => $genome_gff,
        gff_tag => $options->{gff_tag},
    );
    ## use Data::Dumper;
    ## print Dumper $cds_gff;
    my $chr_name = $cds_gff->{stats}->{chromosomes}->[0]; ## Grab the name of the first chromosome
    my @names = @{$cds_gff->{stats}->{feature_names}};
    my @cds_starts = @{$cds_gff->{stats}->{cds_starts}};
    my @cds_ends = @{$cds_gff->{stats}->{cds_ends}};
    my @inter_starts = @{$cds_gff->{stats}->{inter_starts}};
    my @inter_ends = @{$cds_gff->{stats}->{inter_ends}};
    my $num_starts = scalar(@cds_starts);
    print "About to start counting against chromosome: ${chr_name} with ${num_starts} features.\n";

    print $out "Counting TAs in ${input}, this takes a while.\n";
    my $datum = Bio::Adventure::TNSeq::Count_TAs(
        $class,
        data => \@data,
    );
    @data = @{$datum};

    my $tas_file = qq"${output_directory}/${input_name}_tas.txt";
    print $out "Printing list of #TAs observed by position to ${tas_file}.\n";
    my $ta_count = FileHandle->new(">${tas_file}");
    print $ta_count qq"Position\tn_fwd\tn_rev\tn_both\n";
    foreach my $c (0 .. $#data) {
        my $fwd = 0;
        my $rev = 0;
        my $mis = 0;
        my $num_found = 0;
        if ($data[$c]->{fwd}) {
            $fwd = $data[$c]->{fwd};
            $num_found++;
        }
        if ($data[$c]->{rev}) {
            $rev = $data[$c]->{rev};
            $num_found++;
        }
        if ($data[$c]->{mismatch}) {
            $mis = $data[$c]->{mismatch};
            $num_found++;
        }
        if ($num_found > 0) {
            print $ta_count "${c}\t${fwd}\t${rev}\t${mis}\n";
        }
    }
    $ta_count->close();

    my $inter_tas = FileHandle->new(">${output_directory}/${input_name}_interCDS_tas.txt");
    my $cds_tas = FileHandle->new(">${output_directory}/${input_name}_gene_tas.txt");
    my $wig = FileHandle->new(">${output_directory}/${input_name}.wig");

    my $file_start = qq"#type=COPY_NUMBER
#Chromosome\tStart\tEnd\tFeature reads\t TAs
";
    print $inter_tas $file_start;
    print $cds_tas $file_start;
    print $wig "#Start\tReads
variableStep chrom=${chr_name}\n";

    ##print $out "Printing #TAs observed by (inter)CDS region to ${input_name}_(gene|interCDS)_tas.txt\n";
    my $subtotal = 0;
    foreach my $c (0 .. $#data) {
        if ($data[$c]->{ta}) {
            my $fwd = $data[$c]->{fwd};
            my $rev = $data[$c]->{rev};
            my $mis = $data[$c]->{mismatch};
            my $total = $fwd + $rev;
            print $wig "${c}\t${total}\n";
            my $d = $c + 2;
            $subtotal = $subtotal + $total;
            my $feature = "unfound";
            ## $c and $d outline the position of the TA
            ## So I want to figure out if the TA is inside an ORF or intercds, and its name...
          NAMES: for my $i (0 .. $#names) {
                if ($c >= $inter_starts[$i] and $c <= $inter_ends[$i]) {
                    $feature = qq"inter_$names[$i]";
                    last NAMES;
                }
                if ($c >= $cds_starts[$i] and $c <= $cds_ends[$i]) {
                    $feature = qq"cds_$names[$i]";
                    last NAMES;
                }
            }                   ## End searching for feature names.
            ## Now look for inter vs cds features
            if ($feature =~ /^inter/) {
                print $inter_tas "${chr_name}\t${c}\t${d}\t${feature}\t${total}\t1\n";
            } else {
                print $cds_tas "${chr_name}\t${c}\t${d}\t${feature}\t${total}\t1\n";
            }
        }
    }
    $wig->close();
    $inter_tas->close();
    $cds_tas->close();
    $out->close();
    ## This should return something interesting!
    return($subtotal);
}

=head2 C<Run_Essentiality>

Invoke the DeJesus tn_hmm tool to search for essential genes.

=cut
sub Run_Essentiality {
    my ($class, %args) = @_;
    my @param_list = ('1','2','4','8','16','32');
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        runs => 1000,
    );
    my $input = $options->{input};
    print "Remember, tn_hmm requires a .wig input file while mh_ess wants a something_gene_tas.txt.\n";
    my $output = basename($input, ('.txt'));
    ## Set up the tn_hmm job -- these inputs are likely wrong
    my $input_wig = basename($input, ('.txt')) . '.wig';
    $input_wig =~ s/_gene//g;
    $input_wig =~ s/_tas//g;
    $input_wig = qq"outputs/essentiality/${input_wig}";
    my $gff = qq"$options->{libdir}/genome/$options->{species}.gff";
    my $output_file = qq"tn_hmm-${output}.csv";
    my $error_file =  qq"tn_hmm-${output}.err";
    my $comment = qq"## Performing tn_hmm on ${input_wig} using ${gff}\n";
    my $jstring = qq!## This only works with python 2.7.9
eval \$(modulecmd bash purge)
eval \$(modulecmd bash add Python2/2.7.9)
eval \$(modulecmd bash add essentiality)
if [ "Python 2.7.9" \!= "\$(python --version 2>&1)" ]; then
  echo "tn-hmm only works with python 2.7.9, do some module shenanigans."
  exit 1
fi

tn-hmm.py -f ${input_wig} -gff ${gff} \\
  1>outputs/essentiality/${output_file} \\
  2>outputs/essentiality/${error_file}

process_genes.py -f outputs/essentiality/${output_file} \\
  1>outputs/essentiality/genes_${output_file} \\
  2>outputs/essentiality/genes_${error_file}

process_segments.py -f outputs/essentiality/${output_file} \\
  1>outputs/essentiality/segments_${output_file} \\
  2>outputs/essentiality/segments_${error_file}
!;
    ## tn-hmm requires a wig file and gff
    my $tn_hmm = $class->Submit(jname => "tn_hmm",
                                jstring => $jstring,
                                jprefix => "15",
                                comment => $comment,
                                output => $output_file,
                            );

    foreach my $param (@param_list) {
        ## The inputs for this are:
        ##   1. the tas file for either cds or intercds provided by count_tas
        my $output_file = qq"mh_ess-${output}_m${param}.csv";
        my $error_file = qq"mh_ess-${output}_m${param}.err";
        $jstring = qq!## This only works with python 2.7.9
eval \$(modulecmd bash purge)
eval \$(modulecmd bash add Python2/2.7.9)
eval \$(modulecmd bash add essentiality)
## The input for mh-ess is the gene_tas.txt file.
if [ "Python 2.7.9" \!= "\$(python --version 2>&1)" ]; then
  echo "mh-ess only works with python 2.7.9, do some module shenanigans."
  exit 1
fi

gumbelMH.py -f ${input} -m ${param} -s $options->{runs} \\
  1>outputs/essentiality/${output_file} \\
  2>outputs/essentiality/${error_file}
!;
        my $comment = qq!# Run mh-ess using the min hit: $param!;
        my $mh_ess = $class->Submit(
            comment => $comment,
            depends => $options->{depends},
            jname => "mh_ess-${param}",
            jprefix => '16',
            jstring => $jstring,
            output => $output_file,
            prescript => $options->{prescript},
            postscript => $options->{postscript},
        );
    }
    ## This should return something interesting, but it doesn't
    return($tn_hmm);
}

=head2 C<Count_TAs>

Count the number of TAs in a genome.

=cut
sub Count_TAs {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species'],
    );
    my $length = $options->{length};
    my $input = $options->{input};
    my $data = $options->{data};
    eval {use Bio::DB::Sam; 1};
    my $fasta = qq"$options->{libdir}/$options->{libtype}/$options->{species}.fasta";
    my $sam = Bio::DB::Sam->new(-bam => ${input},
                                -fasta => ${fasta},);
    my @targets = $sam->seq_ids;
    my $num = scalar(@targets);
    my $bam = Bio::DB::Bam->open(${input});
    my $header = $bam->header;
    my $target_count = $header->n_targets;
    my $target_names = $header->target_name;
    my $align_count = 0;
    my $million_aligns = 0;
    my $alignstats = qx"samtools idxstats $input";
    my @alignfun = split(/\n/, $alignstats);
    my @aligns = split(/\t/, $alignfun[0]);
    my @unaligns = split(/\t/, $alignfun[1]);
    my $number_reads = $aligns[2] + $unaligns[3];
    my $output_name = qq"${input}.out";
    my $out = FileHandle->new(">$output_name");
    print $out "There are $number_reads alignments in ${input} made of $aligns[2] aligned reads and $unaligns[3] unaligned reads.\n";
  BAMLOOP: while (my $align = $bam->read1) {
        $align_count++;
        ##if ($class->{debug}) {  ## Stop after a relatively small number of reads when debugging.
        ##    last BAMLOOP if ($align_count > 200);
        ##}
        if (($align_count % 1000000) == 0) {
            $million_aligns++;
            print $out "Finished $million_aligns million alignments out of ${number_reads}.\n";
        }
        my $seqid = $target_names->[$align->tid];
        ## my $start = $align->pos + 1;
        ## my $end = $align->calend;
        my $start = $align->pos;
        my $end = $align->calend - 1;
        my $cigar = $align->cigar_str;
        my $strand = $align->strand;
        my $seq = $align->query->dna;
        my $qual= $align->qual;

        my @seq_array = split(//, $seq);
        ## At this point there are a few possibilities:
        ## 1.  We have only a forward hit and it is correct
        ## 2.  We have only a reverse hit and it is correct
        ## 3.  We have both and they are correct
        ###  The rest of the possibilities come up because we allow 1 mismatch in the alignment
        ## 4.  It is a forward hit but off by one base due to mismatch
        ## 5.  It is a reverse hit but off by one base due to mismatch
        ## 6.  It looks like both, but only one is correct due to mismatch
        ## 7.  No hit
        my $forward_hit = 0;
        my $reverse_hit = 0;
        if ($seq_array[0] eq 'T' and $seq_array[1] eq 'A') {
            $forward_hit = 1;
        }
        if ($seq_array[-1] eq 'A' and $seq_array[-2] eq 'T') {
            $reverse_hit = 1;
        }
        my $new_end = $end - 1;
        if ($cigar eq '') {
            ## print "This did not align.\n";
            next BAMLOOP;
        }
        if ($new_end < 0) {
            ##print "The alignment appears to end before the start of the genome.\n"; ## if ($class->{debug});
            ##print "start:$start end:$end cigar:$cigar seq:$seq $strand $qual\n"; ## if ($class->{debug});
            next BAMLOOP;
        }
        ## Now do the same check on the reference sequence (genome)
        my $ref_forward_hit = 0;
        my $ref_reverse_hit = 0;
        ## $data was passed to this function and
        if (!defined($data->[$start]->{nt})) {
            print $out "Undefined start:${start}, setting it to 'N'.\n";
            $data->[$start]->{nt} = 'N';
        }
        if (!defined($data->[$start + 1]->{nt})) {
            print $out "Undefined start:${start}, setting it to 'N'.\n";
            $data->[$start + 1]->{nt} = 'N';
        }
        if (!defined($data->[$new_end]->{nt})) {
            print $out "Undefined new end:${new_end}, setting it to 'N'.\n";
            $data->[$new_end]->{nt} = 'N';
        }
        if ($data->[$start]->{nt} eq 'T' and $data->[$start + 1]->{nt} eq 'A') {
            $ref_forward_hit = 1;
        }
        if ($data->[$new_end]->{nt} eq 'T' and $data->[$end]->{nt} eq 'A') {
            $ref_reverse_hit = 1;
        }

        if ($forward_hit == 0 and $reverse_hit == 0) {
            if ($ref_forward_hit == 1) {
                print $out "$start Mismatch reference forward hit $qual $cigar.\n" if ($class->{debug});
                $data->[$start]->{fwd}++;
                $data->[$start]->{mismatch}++;
            } elsif ($ref_reverse_hit == 1) {
                print $out "$new_end Mismatch reference reverse hit $qual $cigar.\n" if ($class->{debug});
                $data->[$new_end]->{rev}++;
                $data->[$new_end]->{mismatch}++;
            } else {
                print $out "$start $end No easily observable hit $qual $cigar.\n" if ($class->{debug});
                $data->[$start]->{mismatch}++;
                $data->[$new_end]->{mismatch}++;
            }
        } elsif ($forward_hit == 1 and $reverse_hit == 0) {
            if ($ref_forward_hit == 1) {
                print $out "$start fwd_only Query: ${seq} Ref:$data->[$start]->{nt}$data->[$start + 1]->{nt} $qual $cigar\n" if ($class->{debug});
                $data->[$start]->{fwd}++;
            } else {
                print $out "$start Mismatch forward hit Query: ${seq} Ref:$data->[$start]->{nt}$data->[$start + 1]->{nt} $qual $cigar\n" if ($class->{debug});
                $data->[$start]->{mismatch}++;
            }
        } elsif ($forward_hit == 0 and $reverse_hit == 1) {
            if ($ref_reverse_hit == 1) {
                print $out "$new_end rev_only Query: ${seq} Ref:$data->[$new_end]->{nt}$data->[$end]->{nt} $qual $cigar\n" if ($class->{debug});
                $data->[$new_end]->{rev}++;
            } else {
                print $out "$new_end Mismatch reverse hit Query: ${seq} Ref:$data->[$start]->{nt}$data->[$start + 1]->{nt} $qual $cigar\n" if ($class->{debug});
                $data->[$new_end]->{mismatch}++;
            }
        } elsif ($forward_hit == 1 and $reverse_hit == 1) {
            if ($ref_forward_hit == 1 and $ref_reverse_hit == 1) {
                print $out "$start $new_end fwd and reverse hits: Query: ${seq} Ref:$data->[$start]->{nt}$data->[$start + 1]->{nt} $data->[$new_end]->{nt}$data->[$end]->{nt} $qual $cigar\n" if ($class->{debug});
                $data->[$start]->{fwd}++;
                ## $data->[$new_end]->{rev}++  ## Removing this because double count
            } elsif ($ref_forward_hit == 1) {
                print $out "$start fwd hit confirmed: Query: ${seq} Ref:$data->[$start]->{nt}$data->[$start + 1]->{nt} $data->[$new_end]->{nt}$data->[$end]->{nt} $qual $cigar\n" if ($class->{debug});
                $data->[$start]->{fwd}++;
            } elsif ($ref_reverse_hit == 1) {
                print $out "$new_end rev hit confirmed: Query: ${seq} Ref:$data->[$start]->{nt}$data->[$start + 1]->{nt} $data->[$new_end]->{nt}$data->[$end]->{nt} $qual $cigar\n" if ($class->{debug});
                $data->[$new_end]->{rev}++;
            } else {
                print $out "$start $new_end Neither: Query: ${seq} Ref:$data->[$start]->{nt}$data->[$start + 1]->{nt} $data->[$new_end]->{nt}$data->[$end]->{nt} $qual $cigar\n" if ($class->{debug});
                $data->[$start]->{mismatch}++;
                ## $data->[$new_end]->{mismatch}++;  ## No double count here either
            }
        } else {
            die("This else statement should never happen.\n");
        }
    }                           ## End reading each bam entry
    $out->close();
    return($data);
}                               ## End Count_Tas

=head2 C<Allocate_Genome>

When counting a genome, first make a data structure to hold the counts.

=cut
sub Allocate_Genome {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $genome = {length => 0,
                  data => [],
                  sequences => 0,
              };
    my $data = [];
    my $infile = Bio::SeqIO->new(-file => $args{fasta}, -format => 'Fasta');
    while (my $in_seq = $infile->next_seq()) {
        $genome->{sequences}++;
        $genome->{length} = $genome->{length} + $in_seq->length;
        my $seq = $in_seq->seq;
        my @seq_array = split(//, $seq);
        my $count = 0;
        ## print "Allocating array of length $genome->{length}\n";
      LOOP: foreach my $base (@seq_array) {
            $data->[$count]->{nt} = $base;
            my $next;

            if (defined($seq_array[$count + 1])) {
                $next = $seq_array[$count + 1];
            } else {     ## This should only happen at the end of the genome.
                ## print "Setting $count + 1 to 'N', end of genome.\n";
                $seq_array[$count + 1] = 'N';
                $next = 'N';
                last LOOP;
            }
            if ($base eq 'T' and $next eq 'A') {
                $data->[$count]->{ta} = 1;
            } else {
                $data->[$count]->{ta} = 0;
            }
            $data->[$count]->{fwd} = 0;
            $data->[$count]->{rev} = 0;
            $data->[$count]->{mismatch} = 0;
            $count++;
        }
        ## print "Finished allocating array.\n";
    }
    $genome->{data} = $data;
    return($genome);
}

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

L<essentiality>

=cut

1;

## EOF;
