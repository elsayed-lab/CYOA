# NAME

    Bio::Adventure - A perl library to make preprocessing high-throughput data easier.

# SYNOPSIS

      Bio::Adventure.pm Bio::Adventure::Job.pm Some perl libraries to submit jobs to the
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

# DESCRIPTION

    This library should write out PBS compatible job files for torque and
    submit them to appropriate queues on the cluster.  It should also
    collect the outputs and clean up the mess.

## Methods

- `Help`

        Help() always gives 0.
        Before it returns, it will hopefully print some useful information
        regarding ways to invoke Bio::Adventure.pm.

- `Get_Input`

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

- `Get_Vars`

        Handle the peculiar mix of instance options held in $class->{options},
        the set of arguments passed to %args, a list of required and potentially
        missing options, and whatever default values one wishes to set.

        my $options = $class->Get_Vars(args => \%args, required => ['input', 'genome'], bob => 'jane');

- `Set_Vars`

        Handle the peculiar mix of instance options held in $class->{options},
        the set of arguments passed to %args, a list of required and potentially
        missing options, and whatever default values one wishes to set.

        $options = $class->Set_Vars(exclude => 'bob');

- `Last_Stat`

        Last_Stat() reads the final line of an input file and returns it
        as a string.

        This is useful because many of these tools append to .csv files
        with summary information (alignment statistics, sequence sizes,
        etc) and the tests keep a set of expected outputs for the most
        recent run.

- `Read_Genome_Fasta`

        Read a fasta file and return the chromosomes.

- `Read_Genome_GFF`

        Read a GFF file and extract the annotation information from it.

# AUTHOR - atb

Email abelew@gmail.com

# SEE ALSO

    L<Bio::Seq> L<autodie> L<local::lib> L<Bio::Adventure::RNASeq_Aligners>
    L<Bio::Adventure::RNASeq_Assembly> L<Bio::Adventure::RNASeq_Count> L<Bio::Adventure::RNASeq_Aligners> L<Bio::Adventure::RNASeq_QA>
    L<Bio::Adventure::RNASeq_Trim> L<Bio::Adventure::Align_Blast> L<Bio::Adventure::Align_Fasta> L<Bio::Adventure::Align>
    L<Bio::Adventure::Cleanup> L<Bio::Adventure::Compress> L<Bio::Adventure::Convert> L<Bio::Adventure::PBS>
    L<Bio::Adventure::Prepare> L<Bio::Adventure::Riboseq> L<Bio::Adventure::SeqMisc> L<Bio::Adventure::SNP> L<Bio::Adventure::Status>
    L<Bio::Adventure::TNSeq>
