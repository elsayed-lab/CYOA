package Bio::Adventure::Count;
## LICENSE: gplv2
## ABSTRACT:  Kitty!
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use Bio::DB::Sam;
use Bio::Tools::GFF;
use Cwd;
use File::Basename;
use File::Path qw"make_path";
use File::Which qw"which";
use String::Approx qw"amatch";

=head1 NAME

 Bio::Adventure::Count - Perform Sequence alignments counting with HTSeq

=head1 SYNOPSIS

 These functions handle the counting of reads, primarily via htseq.

=head1 METHODS

=cut
sub Guess_Strand {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species', 'gff_type', 'gff_tag',],
        genes => 'all',
        jmem => 12,
        jprefix => '41',
        output => 'strand_counts.txt',);
    my $job_name = $class->Get_Job_Name();
    my $output = $options->{output};
    my $comment = '## This submits a job to figure out the strandedness of a library.';
    my $jstring = qq!
use Bio::Adventure;
use Bio::Adventure::Count;
my \$result = Bio::Adventure::Count::Guess_Strand_Worker(\$h,
  input => '$options->{input}',
  genes => '$options->{genes}',
  gff_tag => '$options->{gff_tag}',
  gff_type => '$options->{gff_type}',
  output => '${output}',
  species => '$options->{species}',);
!;
    my $guess = $class->Submit(
        comment => $comment,
        genes => $options->{genes},
        gff_tag => $options->{gff_tag},
        gff_type => $options->{gff_type},
        input => $options->{input},
        jcpus => 1,
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
        jname => qq"guess_strand_${job_name}",
        jprefix => $options->{jprefix},
        jqueue => 'workstation',
        jstring => $jstring,
        language => 'perl',
        output => $output,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        species => $options->{species},);
    return($guess);
}

sub Guess_Strand_Worker {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species', 'gff_type', 'gff_tag',],
        coverage => 0.2,
        jcpus => 1,
        output_log => '',
        output => '',);
    my $log = FileHandle->new(">$options->{output}.log");
    my $gff = qq"$options->{libpath}/$options->{libtype}/$options->{species}.gff";
    my $fasta = qq"$options->{libpath}/$options->{libtype}/$options->{species}.fasta";
    my $features = $class->Read_Genome_GFF(
        gff => $gff,
        gff_type => $options->{gff_type},
        gff_tag => $options->{gff_tag},);
    my $sam = Bio::DB::Sam->new(-bam => $options->{input},
                                -fasta => ${fasta},);
    my @targets = $sam->seq_ids;
    my $num = scalar(@targets);
    my $bam = Bio::DB::Bam->open($options->{input});
    my $header = $bam->header;
    my $target_count = $header->n_targets;
    my $target_names = $header->target_name;
    my $align_count = 0;
    my $million_aligns = 0;
    my $alignstats = qx"samtools idxstats $options->{input}";
    my @alignfun = split(/\n/, $alignstats);
    my @aligns = split(/\t/, $alignfun[0]);
    my @unaligns = split(/\t/, $alignfun[1]);
    my $number_reads = $aligns[2] + $unaligns[3];
    my $output_name = qq"$options->{input}.out";
    print $log "There are $number_reads alignments in $options->{input} made of $aligns[2] aligned reads and $unaligns[3] unaligned reads.\n";
    my $forward_hit = 0;
    my $reverse_hit = 0;
    my %result = ();
  BAMLOOP: while (my $align = $bam->read1) {
      my $cigar = $align->cigar_str;
      if ($cigar eq '') {
          ## print "This did not align.\n";
          next BAMLOOP;
      }
      $align_count++;
      ##if ($class->{debug}) {  ## Stop after a relatively small number of reads when debugging.
      ##    last BAMLOOP if ($align_count > 200);
      ##}
      if (($align_count % 1000000) == 0) {
          $million_aligns++;
          print $log "Finished $million_aligns million alignments out of ${number_reads}.\n";
      }
      my $seqid = $target_names->[$align->tid];
      ## my $start = $align->pos + 1;
      ## my $end = $align->calend;
      my $start = $align->pos;
      my $end = $align->calend - 1;
      my $strand = $align->strand;
      my $seq = $align->query->dna;
      my $qual = $align->qual;
      my $pairedp = $align->paired;
      my $unmappedp = $align->unmapped;
      my $mate_unmappedp = $align->munmapped;
      my $reversedp = $align->reversed;
      my $mate_reversedp = $align->mreversed;
      my $mate_id = $align->mtid;
      my $mate_start = $align->mpos;
      my $properp = $align->proper_pair;
      if ($unmappedp && $mate_unmappedp) {
          next BAMLOOP;
      }

      my @seq_array = split(//, $seq);
      my @tags = $align->get_all_tags;
      ## A reminder when playing with tag values and get_all_tags()
      ## These are coming from the 2nd sam column and are the result of a decimal->binary conversion:
      ## E.g. a paired, proper, first of pair, reverse in binary is: 1100101 -> 1+2+16+64 -> 83
      ## 0x1(1) or first binary character: paired or unpaired
      ## 0x2(2) or second binary character: proper or not proper pair
      ## 0x4(4): this read is unmapped
      ## 0x8(8): the mate read is unmapped
      ## 0x10(16): this read is reverse
      ## 0x20(32): the mate is reverse
      ## 0x40(64): this is the first of a pair
      ## 0x80(128): this is the second of a pair
      ## 0x100(256): this is not the primary read
      ## 0x200(512): this failed the sequencer checks
      ## 0x400(1024): this is an optical duplicate
      ## 0x800(2048): this is a supplemental alignment.
      my $first = $align->get_tag_values('FIRST_MATE');

      ## Look for a feature at this position...
      my %chr_features = %{$features->{$seqid}};
      for my $feat_id (keys %chr_features) {
          my %feat = %{$chr_features{$feat_id}};

          ##  ------------------->>>>>-------<<<<<--------------------------
          ##                  >>>>>
          if ($start <= $feat{start} && $end >= $feat{start} && $end < $feat{end}) {
              use Data::Dumper;
              print Dumper %feat;
          }
      }
  } ## End reading each bam entry
    $log->close();
}



=head2 C<HT_Multi>

 Invoke htseq multiple times with options for counting different transcript types.

=item C<Arguments>

 species(required): Defines the gff/fasta files to query in {libpath}.
 input(required): The sorted/indexed .bam alignment file to count.
 htseq_stranded(required): Boolean of the stranded state of the library.
 gff_type(''): Type of feature to count.  If left alone, this will count up
   each category of feature and choose the most represented.
   (deprecated in favor of gff_type)
 gff_type('gene): Ibid.  I just have not converted everything to use this.
 gff_tag('ID'): GFF tag used to identify the genes/transcripts/etc.  In
   metazoans this is usually 'ID', in bacteria it is often
   'locus_tag', in parasites 'gene_id'.  There are lots of choices.
 libtype('genome'): Used to differentiate between genomic counts vs. rRNA
   vs. contaminants.
 modules('htseq'): List of environment modules.
 paired('1'): Is this a paired library?

=cut
sub HT_Multi {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species', 'input',],
        modules => ['htseq'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('htseq-count');
    die('Could not find htseq in your PATH.') unless($check);

    my %ro_opts = %{$options};
    my $species = $options->{species};
    my $htseq_input = $options->{input};
    my $stranded = $options->{stranded};
    if ($stranded eq '1') {
        $stranded = 'yes';
    } elsif ($stranded eq '0') {
        $stranded = 'no';
    }

    my $gff_type = $options->{gff_type};
    my $gff_tag = $options->{gff_tag};
    my @jobs = ();
    my $script_suffix = qq"";
    if ($args{suffix}) {
        $script_suffix = $args{suffix};
    }
    my $jprefix = "";
    $jprefix = $args{jprefix} if ($args{jprefix});
    my @gff_types = ('antisense', 'exon', 'fiveputr', 'interCDS',
                     'linc', 'mi', 'misc', 'nmd', 'operons', 'pseudo',
                     'rintron', 'rrna', 'sn', 'sno', 'threeputr');
    my $htseq_runs = 0;
    ## Top level directory containing the input files.
    my $top_dir = basename(getcwd());
    ## The sam/bam input basename
    my $output_name = basename($htseq_input, ('.bam', '.sam'));

    foreach my $gff_type (@gff_types) {
        my $gff = qq"$options->{libpath}/genome/${species}_${gff_type}.gff";
        my $gtf = $gff;
        $gtf =~ s/\.gff/\.gtf/g;
        my $htseq_jobname = qq"hts_${gff_type}_${output_name}_$options->{species}_s${stranded}_${gff_type}_${gff_tag}";
        if (-r "$gff") {
            print "Found $gff, performing htseq with it.\n";
            my $ht = $class->Bio::Adventure::Count::HTSeq(
                htseq_gff => $gff,
                input => $htseq_input,
                gff_type => $gff_type,
                gff_tag => $options->{gff_tag},
                jdepends => $options->{jdepends},
                jname => $htseq_jobname,
                jprefix => $jprefix,
                jqueue => 'throughput',
                postscript => $options->{postscript},
                prescript => $options->{prescript},
                suffix => $options->{suffix},);
            push(@jobs, $ht);
            $htseq_runs++;
        } elsif (-r "$gtf") {
            print "Found ${gtf}, performing htseq with it.\n";
            my $ht = $class->Bio::Adventure::Count::HTSeq(
                gff_type => $gff_type,
                htseq_gff => $gtf,
                gff_tag => $options->{gff_tag},
                input => $htseq_input,
                jdepends => $options->{jdepends},
                jname => $htseq_jobname,
                jprefix => $options->{jprefix},
                jqueue => 'throughput',
                postscript => $options->{postscript},
                prescript => $options->{prescript},
                suffix => $options->{suffix},);
            push(@jobs, $ht);
            $htseq_runs++;
        }
    } ## End foreach type
    ## Also perform a whole genome count
    my $gff = qq"$options->{libpath}/genome/${species}.gff";
    my $gtf = $gff;
    $gtf =~ s/\.gff/\.gtf/g;
    my $htall_jobname = qq"htall_${output_name}_$options->{species}_s${stranded}_$ro_opts{gff_type}_$ro_opts{gff_tag}";
    if (-r "$gff") {
        print "Found ${gff}, performing htseq_all with it.\n";
        my $ht = $class->Bio::Adventure::Count::HTSeq(
            htseq_gff => $gff,
            gff_tag => $ro_opts{gff_tag},
            gff_type => $ro_opts{gff_type},
            input => $htseq_input,
            jdepends => $options->{jdepends},
            jname => $htall_jobname,
            jprefix => $jprefix,
            jqueue => 'throughput',
            postscript => $options->{postscript},
            prescript => $options->{prescript},
            suffix => $options->{suffix},);
        push(@jobs, $ht);
    } elsif (-r "${gtf}") {
        print "Found ${gtf}, performing htseq_all with it.\n";
        my $ht = $class->Bio::Adventure::Count::HTSeq(
            htseq_gff => $gff,
            gff_type => 'none',
            input => $htseq_input,
            jdepends => $options->{jdepends},
            jname => $htall_jobname,
            jprefix => $jprefix,
            jqueue => 'throughput',
            postscript => $args{postscript},
            prescript => $args{prescript},
            suffix => $args{suffix},);
        push(@jobs, $ht);
    } else {
        print "Did not find ${gff} nor ${gtf}, not running htseq_all.\n";
    }
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return(\@jobs);
}

=head2 C<HT_Types>

 Guess about most appropriate flags for htseq-count.

 This function reads the first 100k lines of a gff file and use that
 to guess at the most likely type of feature when invoking htseq.

=item C<Arguments>

 gff_type('gene'): When set, this will just count that type.
 gff_tag('ID'): Ditto, but the GFF tag for IDs.

=cut
sub HT_Types {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,);
    my $my_type = $options->{gff_type};
    my $my_id = $options->{gff_tag};
    print "Calling htseq with options for type: ${my_type} and tag: ${my_id}.\n";
    $my_type = '' unless($my_type);
    $my_type = '' unless($my_id);
    my $gff_out = {};
    my $max_lines = 100000;
    my $count = 0;
    my %found_types = ();
    my %found_canonical = ();
    my $reader = new FileHandle("<$options->{annotation}");
  LOOP: while(my $line = <$reader>) {
        chomp $line;
        next LOOP if ($line =~ /^\#/);
        $count++;
        ## chr1    unprocessed_pseudogene  transcript      3054233 3054733 .       +       .       gene_id "ENSMUSG00000090025"; transcript_id "ENSMUST00000160944"; gene_name "Gm16088"; gene_source "havana"; gene_biotype "pseudogene";transcript_name "Gm16088-001"; transcript_source "havana"; tag "cds_end_NF"; tag "cds_start_NF"; tag "mRNA_end_NF"; tag "mRNA_start_NF";

        my ($chromosome, $source, $type, $start, $end, $score, $strand, $frame, $attributes) = split(/\t/, $line);
        my @attribs = split(/;\s*/, $attributes);
        my (@pairs, @names, @values) = ();
        my $splitter = '\s+';
        if ($attribs[0] =~ /=/) {
            $splitter = '=';
        }
        for my $attrib (@attribs) {
            my ($name, $value) = split(/$splitter/, $attrib);
            push(@names, $name);
            push(@values, $value);
        }
        my $canonical = $names[0];
        if ($found_types{$type}) {
            $found_types{$type}++;
        } else {
            $found_types{$type} = 1;
        }
        if ($found_canonical{$canonical}) {
            $found_canonical{$canonical}++;
        } else {
            $found_canonical{$canonical} = 1;
        }
        if ($count > $max_lines) {
            last LOOP;
        }
    } ## End loop
    $reader->close();
    my $found_my_type = 0;
    my $max_type = '';
    my $max = 0;
    foreach my $type (keys %found_types) {
        if ($found_types{$type} > $max) {
            $max_type = $type;
            $max = $found_types{$type};
        }
        if ($found_types{$my_type}) {
            print "The specified type: ${my_type} is in the gff file, comprising $found_types{$my_type} of the first 40,000.\n";
            $found_my_type = 1;
        }
    } ## End the loop

    my $max_can = 0;
    my $max_canonical = 0;
    foreach my $can (keys %found_canonical) {
        if (!defined($found_canonical{$can})) {
            $found_canonical{$can} = 0;
        }
        if ($found_canonical{$can} > $max_canonical) {
            $max_can = $can;
            $max_canonical = $found_canonical{$can};
        }
    } ## End the loop
    my $returned_canonical = $max_can;

    my $returned_type = "";
    if ($found_my_type == 1) {
        $returned_type = $my_type;
    } else {
        print "Did not find your specified type.  Changing it to: ${max_type} which had ${max} entries.\n";
        $returned_type = $max_type;
    }
    my $ret = [$returned_type, $returned_canonical];
    return($ret);
}

=head2 C<HTSeq>

 Run htseq-count on a sorted, indexed bam file.
 10.1093/bioinformatics/btu638

 Invoke htseq-count.  This should be able to automagically pick up
 other types of countable features and send the htseq-count results to
 a separate count file.

=item C<Arguments>

 input: Sorted/indexed bam file to count.
 species: Defines the gff/fasta files to read.
 stranded: Is this library stranded?
 htseq_args(--order=name --idattr=gene_id --minaqual=10 --type exon
  --stranded=yes --mode=union): Define arbitrary htseq arguments.
 gff_type(''): Redundant with gff_type, used to choose a specific
  type to count; when left blank, HT_Types() will make a guess.
 gff_type('gene'): Deprecated, Ibid.
 gff_tag('ID'): GFF tag used to identify the counted features.
 jname(''): Job name base for the cluster.
 jprefix(''): Prefix for the job name and output directory.
 libtype('genome'): Choose the library to count,
  genomic/rRNA/contaminants/etc.
 mapper('hisat2'): Specify the mapper used so that the output file
  will make it easier to tell where it came from.
 modules('htseq'): List of environment modules to load.
 paired('1'): Is this library paired?

=cut
sub HTSeq {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species', 'htseq_args',],
        jname => '',
        jprefix => '',
        modules => ['htseq'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $stranded = $options->{stranded};
    if ($stranded eq '1') {
        $stranded = 'yes';
    } elsif ($stranded eq '0') {
        $stranded = 'no';
    }
    print "Note, this is running with stranded: $options->{stranded}\n";
    my $gff_tag = $options->{gff_tag};
    my $htseq_input = $options->{input};
    my $gff_type = 'all';
    if (defined($options->{gff_type}) && $options->{gff_type} ne '') {
        $gff_type = $options->{gff_type};
    }
    ## Top level directory containing the input files.
    my $top_dir = basename(getcwd());
    ## The sam/bam input basename
    ## Keep in mind this input bam file likely has the genome as part of its name
    my $output_name = basename($htseq_input, ('.bam', '.sam'));

    ## And directory containing it.
    my $output_dir = dirname($htseq_input);
    my $output = qq"${output_dir}/${output_name}";
    my $gff = qq"$options->{libpath}/$options->{libtype}/$options->{species}.gff";
    $gff = $args{htseq_gff} if ($args{htseq_gff});
    my $gtf = $gff;
    $gtf =~ s/\.gff/\.gtf/g;
    my $htseq_args = "";

    ## Set the '-t FEATURETYPE --type' argument used by htseq-count
    ## This may be provided by a series of defaults in %HT_Multi::gff_types, overridden by an argument
    ## Or finally auto-detected by HT_Types().
    ## This is imperfect to say the nicest thing possible, I need to consider more appropriate ways of handling this.
    my $gff_type_arg = '';
    my $gff_tag_arg = '';
    my $annotation = $gtf;
    if (!-r "${gtf}") {
        $annotation = $gff;
    }
    if (!defined($gff_type) or $gff_type eq '' or $gff_type eq 'auto' or
            !defined($gff_tag) or $gff_tag eq '' or $gff_tag eq 'auto') {
        my $gff_type_pair = $class->Bio::Adventure::Count::HT_Types(
            annotation => $annotation,
            type => $gff_type,);
        $gff_type = $gff_type_pair->[0];
        $gff_type_arg = qq" --type ${gff_type}";
        $gff_tag_arg = qq" --idattr ${gff_tag}"
    } elsif (ref($gff_type) eq 'ARRAY') {
        $gff_type_arg = qq" --type $gff_type->[0]";
        $gff_tag_arg = qq" --idattr ${gff_tag}";
    } elsif ($gff_type eq 'none') {
        $gff_type_arg = qq'';
        $gff_tag_arg = qq'';
    } else {
        $gff_type_arg = qq" --type ${gff_type}";
        $gff_tag_arg = qq" --idattr ${gff_tag}";
    }

    $output .= qq"_${gff_type}_s${stranded}_${gff_type}_${gff_tag}.count";
    if (!-r "${gff}" and !-r "${gtf}") {
        die("Unable to read ${gff} nor ${gtf}, please fix this and try again.\n");
    }
    my $error = basename($output, ('.count'));
    $error = qq"${output_dir}/${error}.stderr";

    my $htseq_jobname = qq"hts_${top_dir}_${gff_type}_$options->{mapper}_$options->{species}_s${stranded}_${gff_type}_${gff_tag}";
    my $htseq_invocation = qq!htseq-count \\
  -q -f bam -s ${stranded} ${gff_type_arg} ${gff_tag_arg} \\!;
    my $jstring = qq!
${htseq_invocation}
  ${htseq_input} \\
  ${annotation} \\
  2>${error} \\
  1>${output}
xz -f -9e ${output} \\
  2>${error}.xz \\
  1>${output}.xz
!;
    $output = qq"${output}.xz";
    my $comment = qq!## Counting the number of hits in ${htseq_input} for each feature found in ${annotation}
## Is this stranded? ${stranded}.  The defaults of htseq are:
## $options->{htseq_args}
!;
    my $htseq = $class->Submit(
        comment => $comment,
        gff_type => $options->{gff_type},
        gff_tag => $options->{gff_tag},
        input => $htseq_input,
        jdepends => $options->{jdepends},
        jmem => 6,
        jname => $htseq_jobname,
        jprefix => $options->{jprefix},
        jqueue => 'throughput',
        jstring => $jstring,
        modules => $options->{modules},
        output => $output,
        postscript => $args{postscript},
        prescript => $args{prescript},);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($htseq);
}

=head2 C<Jellyfish>

 Run jellyfish with multiple values of K on a fast(a|q) file(s).
 10.1093/bioinformatics/btr011

 Use jellyfish to count up all of the kmers in a set of sequence(s).
 This should also send the wacky fasta format fo counted sequences to a
 tsv of kmers and numbers, which should be a rather more tractable
 format with which to play.

=item C<Arguments>

 input(required): Fast(a|q) file(s) to index.
 length('9,11,13,15'): Text list of k values to give to jellyfish.
 jprefix(18): Prefix for the job name.
 modules('jellyfish'): Module list.

=cut
sub Jellyfish {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        length => '9,11,13,15',
        jprefix => 18,
        modules => ['jellyfish'],);
    my @kmer_array = split(/\,|:/, $options->{length});
    my $count = 0;
    my $ret;
    if (scalar(@kmer_array) > 1) {
      KMERARR: for my $k (@kmer_array) {
          my $job = $class->Bio::Adventure::Count::Jellyfish(
              %{$options},
              length => $k);
          if ($count == 0) {
              $ret = $job;
          } else {
              $ret->{$count} = $job;
          }
          $count++;
      }
        return($ret);
    }

    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('jellyfish');
    die('Could not find jellyfish in your PATH.') unless($check);
    my $job_name = $class->Get_Job_Name();
    my $inputs = $class->Get_Paths($options->{input});
    my $cwd_name = basename(cwd());
    my $output_dir = qq"outputs/$options->{jprefix}jellyfish_${cwd_name}";

    my $input_string = qq"<(less $options->{input})";
    if ($options->{input} =~ /\:|\;|\,|\s+/) {
        ## Then multiple files were provided.
        my @file_list = split(/\:|\;|\,|\s+/, $options->{input});
        $input_string = '';
        for my $f (@file_list) {
            $input_string .= qq"<(less $f) ";
        }
    }

    my $jelly_base = qq"${output_dir}/${cwd_name}_$options->{length}";
    my $count_file = qq"${jelly_base}.count";
    my $info_file = qq"${jelly_base}.info";
    my $histogram_file = qq"${jelly_base}.hist";
    my $count_fasta = qq"${jelly_base}_by_count.fasta";
    my $matrix_file = qq"${jelly_base}_matrix.csv";
    my $comment = '## Invoke jellyfish on some sequence!';
    my $jstring = qq!mkdir -p ${output_dir}
jellyfish count -m $options->{length} \\
  -o ${count_file} \\
  -s 50000 -t 4 \\
  ${input_string} \\
  2>${jelly_base}.stderr 1>${jelly_base}.stdout
jellyfish info ${count_file} > ${info_file} \\
  2>>${jelly_base}.stderr
jellyfish histo ${count_file} > ${histogram_file} \\
  2>>${jelly_base}.stderr
jellyfish dump ${count_file} > ${count_fasta} \\
  2>>${jelly_base}.stderr
!;

    my $jelly = $class->Submit(
        cpus => $options->{cpus},
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => qq"jelly_${job_name}_$options->{length}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => 12,
        modules => $options->{modules},
        count_file => $count_file,
        info_file => $info_file,
        histogram_file => $histogram_file,
        count_fasta => $count_fasta,
        output => $matrix_file,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');

    $comment = qq"## This should create a matrix with rows as kmers and elements
## comprised of the number of occurrences.
";
    my $new_prefix = $options->{jprefix} + 1;
    $jstring = qq!
use Bio::Adventure;
use Bio::Adventure::Phage;
my \$result = \$h->Bio::Adventure::Count::Jellyfish_Matrix(
  input => '${count_fasta}',
  jdepends => '$jelly->{job_id}',
  jname => 'jelly_matrix',
  jprefix => '${new_prefix}',
  output => '${matrix_file}',);
!;
    my $matrix_job = $class->Submit(
        comment => $comment,
        input => $count_fasta,
        jdepends => $jelly->{job_id},
        jname => qq"jelly_matrix_$options->{length}",
        jprefix => $new_prefix,
        jstring => $jstring,
        output => $matrix_file,
        language => 'perl',);
    $jelly->{matrix_job} = $matrix_job;

    my $compress_files = qq"${count_file}:${info_file}:${histogram_file}:${count_fasta}:${matrix_file}";
    my $comp = $class->Bio::Adventure::Compress::Recompress(
        comment => '## Compress the jellyfish output files.',
        jdepends => $matrix_job->{job_id},
        jname => qq"xzjelly_$options->{length}",
        jprefix => $options->{jprefix} + 1,
        input => $compress_files,);
    $jelly->{compression} = $comp;
    $jelly->{output} = qq"${matrix_file}.xz";
    $jelly->{count_file} = qq"${count_file}.xz";
    $jelly->{info_file} = qq"${info_file}.xz";
    $jelly->{histogram_file} = qq"${histogram_file}.xz";
    $jelly->{count_fasta} = qq"${count_fasta}.xz";
    return($jelly);
}

=head2 C<Jellyfish_Matrix>

 Convert the jellyfish fasta format to a matrix-compatible tsv.

 This function is responsible for actually converting the fasta output
 from jellyfish into tsv.  It is pretty quick and dirty.

=item C<Arguments>

 input(required): The peculiar fasta file produced by jellyfish.
 output('fasta_matrix.csv'): Filename to which to convert the
  jellyfish output.
 jprefix(19): Prefix for the jobname/output directory.

=cut
sub Jellyfish_Matrix {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        output => 'fasta_matrix.csv',
        jprefix => 19,);
    my $job_name = $class->Get_Job_Name();
    my $inputs = $class->Get_Paths($options->{input});
    my $in = FileHandle->new("less $options->{input} |");
    my $counter = 1;
    my $counts = {};
    my $nmer_count = 0;
    my $nmer_identity = '';
    while (my $line = <$in>) {
        chomp($line);
        if ($counter == 1) {  ## Then it is the number line, also why in the flying hell did they do that
            $nmer_count = $line;
            $nmer_count =~ s/^>//g;
            $counter++;
        } elsif ($counter == 2) {
            $counter--;
            $nmer_identity = $line;
            $counts->{$nmer_identity} = $nmer_count;
        } else {
            die('Should not get here.');
        }
    }
    $in->close();

    my $out = FileHandle->new(">$options->{output}");
    foreach my $k (sort keys %{$counts}) {
        print $out "${k}\t$counts->{$k}\n";
    }
    $out->close();
    return($nmer_count);
}

=head2 C<Mash>

 Calculate distances using mash

=cut
sub Mash {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jprefix => 21,
        length => 9,
        sketch => 9,
        modules => ['mash'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('mash');
    die('Could not find mash in your PATH.') unless($check);
    my $job_name = $class->Get_Job_Name();
    my $inputs = $class->Get_Paths($options->{input});
    ## Unlike most of my jobs, the input argument here is a directory, so just grab it
    my $cwd_name = basename($options->{input});
    my $output_dir = qq"outputs/$options->{jprefix}mash_${cwd_name}/dist";
    my $sketch_dir = qq"outputs/$options->{jprefix}mash_${cwd_name}/sketch";
    my $comment = qq"## Playing with mash";
    my $jstring = qq!mkdir -p ${output_dir}
mkdir -p ${sketch_dir}
for fasta in $options->{input}/*; do
  name="\${fasta%.*}"
  mash sketch \\
    -s $options->{sketch} -k $options->{length} \\
    $options->{input}/\${fasta} \\
    -o ${sketch_dir}/\${name}
done

for outer in ${sketch_dir}/*; do
  name="\${outer%.*}"
  for inner in ${sketch_dir}/*; do
    mash dist \\
      ${sketch_dir}/\${outer} \\
      ${sketch_dir}/\${inner} >> \\
      ${output_dir}/pairwise_distances_s$options->{sketch}_k$options->{length}.txt
!;

    my $mash = $class->Submit(
        comment => $comment,
        jcpus => $options->{jcpus},
        jdepends => $options->{jdepends},
        jmem => 12,
        jname => qq"mash_${job_name}_$options->{length}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        modules => $options->{modules},);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($mash);
}

=head2 C<Mi_Map>

 Map reads to mature/immature miRNA.

 This function has not been used in a very long time and likely will
 require some work if I wish to use it again.

=item C<Arguments>

 mirbase_data(required): File containing miRNAs from the mirbase.
 mature_fasta(required): Fasta file of the mature miRNA species.
 mi_genome(required): The set of all miRNAs expected.
 bamfile(required): Input bam file to search.

=cut
sub Mi_Map {
    my ($class, %args) = @_;
    eval 'use Bio::DB::Sam; 1;';
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['mirbase_data', 'mature_fasta', 'mi_genome', 'bamfile',]);
    my $bam_base = basename($options->{bamfile}, ('.bam', '.sam'));
    my $pre_map = qq"${bam_base}_mirnadb.txt";
    my $final_map = qq"${bam_base}_mature.count";

    ## Step 1:  Read the mature miRNAs to get the global IDs and sequences of the longer transcripts.
    ## This provides the IDs as 'mm_miR_xxx_3p' and short sequences ~21 nt.
    print "Starting to read miRNA sequences.\n";
    my $sequence_db = $class->Bio::Adventure::Count::Read_Mi(seqfile => $options->{mature_fasta});

    ## Step 2:  Collate those IDs against the set of miRNA_transcripts->miRNA_mature IDs using
    ## the tab delimited file downloaded from mirbase.org
    print "Starting to read miRNA mappings.\n";
    my $sequence_mappings = $class->Bio::Adventure::Count::Read_Mappings_Mi(
        mappings => $options->{mirbase_data},
        output => $pre_map,
        seqdb => $sequence_db,);
    ## At this point, we have brought together mature sequence/mature IDs to parent transcript IDs
    ## When we read the provided bam alignment, we will simultaneously read the immature miRNA database
    ## and use them to make a bridge from the mature sequence/IDs to immature miRNA IDs.
    my $final_hits = $class->Bio::Adventure::Count::Read_Bam_Mi(
        mappings => $sequence_mappings,
        mi_genome => $options->{mi_genome},
        bamfile => $options->{bamfile},);
    my $printed = $class->Bio::Adventure::Count::Final_Print_Mi(
        data => $final_hits,
        output => $final_map,);

    my $job = $printed;
    $job->{final_hits} = $final_hits;
    return($job);
} ## End of Mi_Map

=head2 C<Read_Mi>

 Read an miRNA database.

 This takes the fasta file from the mirbase and extracts the IDs and
 sequences.

=item C<Arguments>

 seqfile: The fasta file in question.

=cut
sub Read_Mi {
    my ($class, %args) = @_;
    my $fasta = Bio::SeqIO->new(-file => $args{seqfile}, -format => 'Fasta');
    my %sequences = ();
    while (my $mi_seq = $fasta->next_seq()) {
        next unless(defined($mi_seq->id));
        my $id = $mi_seq->id;
        my $length = $mi_seq->length;
        my $seq = $mi_seq->seq;
        $sequences{$id}->{sequence} = $seq;
    }
    return(\%sequences);
}

=head2 C<Read_Mappings_Mi>

 Read an miRNA database and get the connections between the various IDs, mature
 sequences, and immature sequences.

=item C<Arguments>

 output: Output file to write.
 seqdb: Hash of IDs to sequences from Read_Mi().
 mappings: Hash of IDs to precursors/etc.

=cut
sub Read_Mappings_Mi {
    my ($class, %args) = @_;
    my $output = $args{output};
    my $seqdb = $args{seqdb};
    my $inmap = new FileHandle($args{mappings});
    my $mimap = {};
    while (my $line = <$inmap>) {
        chomp $line;
        $line =~ s/"//g;
        my ($hit_id, $ensembl, $mirbase_id, $fivep_mature, $fivep_id, $threep_mature, $threep_id) = split(/\s+/, $line);
        if (defined($fivep_id) and $fivep_id ne "") {
            $mimap->{$fivep_id}->{mirbase_id} = $mirbase_id;
            $mimap->{$fivep_id}->{hit_id} = $hit_id;
            $mimap->{$fivep_id}->{ensembl} = $ensembl;
            $mimap->{$fivep_id}->{mimat} = $fivep_id;
        }
        if (defined($threep_id) and $threep_id ne "") {
            $mimap->{$threep_id}->{mirbase_id} = $mirbase_id;
            $mimap->{$threep_id}->{hit_id} = $hit_id;
            $mimap->{$threep_id}->{ensembl} = $ensembl;
            $mimap->{$threep_id}->{mimat} = $threep_id;
        }
    }
    $inmap->close();
    my $out = FileHandle->new(">${output}");
    ## Now re-key the database of miRNAs so that there are (potentially) multiple elements to each mirbase ID
    ## This is a bit redundant, but hopefully clearer.  The idea is to first gather all miR data, then make a list of specific
    ## ID/sequences for each mature miRNA child of the mirbase_id parent entry.  Eg:
    ## mirbase_id parent ->  [ mature_id1, mature_id2, ... ]
    ## where each mature_id is a { sequence, count, id, ensembl, ... }
    my $newdb = {};
  LOOP: foreach my $id (keys %{$seqdb}) {
        next LOOP unless ($id =~ /^mmu/);
        my @hit_list = ();
        if (defined($mimap->{$id})) {
            my $new_id = $id;
            my $sequence = $seqdb->{$id}->{sequence};
            if ($mimap->{$id}->{hit_id}) {
                $new_id = $mimap->{$id}->{hit_id};
            }
            my $element = {
                count => 0,
                ensembl => $mimap->{$id}->{ensembl},
                hit_id => $mimap->{$id}->{hit_id},
                id => $id,
                mimat => $mimap->{$id}->{mimat},
                mirbase_id => $mimap->{$id}->{mirbase_id},
                sequence => $sequence,
            };
            my $ensembl_id = $element->{ensembl};
            if (defined($newdb->{$ensembl_id})) {
                ## Then there should be an existing element, append to it.
                @hit_list = @{$newdb->{$ensembl_id}};
                push(@hit_list, $element);
                $newdb->{$ensembl_id} = \@hit_list;
                print $out "${id}; $element->{mirbase_id}; $element->{ensembl}; $element->{hit_id}; $element->{mimat}; $element->{sequence}\n";
            } else {
                @hit_list = ();
                push(@hit_list, $element);
                $newdb->{$ensembl_id} = \@hit_list;
            }
        } else {
            ## print "Did not find $id\n";
        }
    }
    $out->close();
    return($newdb);
}

=head2 C<Read_Bam_Mi>

 Read a bam file and cross reference it against an miRNA database.

=item C<Arguments>

 mappings: Set of database entries to query against.
 bamfile: Input alignments to search.
 mi_genome: Fasta file containing the genomic context of the miRNAs.

=cut
sub Read_Bam_Mi {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
    );
    my $mappings = $args{mappings};
    my $sam = Bio::DB::Sam->new(-bam => $args{bamfile}, -fasta=> $args{mi_genome},);
    my $bam = $sam->bam;
    my $header = $bam->header;
    my $target_count = $header->n_targets;
    my $target_names = $header->target_name;
    my $align_count = 0;
    print "Started reading bamfile.\n";
    ## I probably don't need to acquire most of this information, as I am only really taking
    ## the read_seq and read_seqid
  BAM: while (my $align = $bam->read1) {
        if ($options->{debug}) {
            last BAM if ($align_count > 4000);
        }
        my $read_seqid = $target_names->[$align->tid];
        ## my $read_start = $align->pos + 1;
        ## my $read_end = $align->calend;
        ## my $read_strand = $align->strand;
        ## my $read_cigar = $align->cigar_str;
        ## my @read_scores = $align->qscore;        # per-base quality scores
        ## my $read_match_qual= $align->qual;       # quality of the match
        ## my $read_length = $align->query->end;
        my $read_seq = $align->query->seq->seq; ## maybe?
        $align_count++;

        ## Check each element in the mappings to see if they are contained in the read's sequence and vice versa.
        ## If so, increment the count for that element and move on.
        $read_seqid =~ s/^chr[\d+]_//g;
        if ($mappings->{$read_seqid}) {
            my @element_list = @{$mappings->{$read_seqid}};
            my $found_element = 0;
            my $element_length = scalar(@element_list);
            ## print "Found map, checking against ${element_length} mature RNAs.\n";
            foreach my $c (0 .. $#element_list) {
                my $element_datum = $element_list[$c];
                my $element_seq = $element_list[$c]->{sequence};
                ## print "Comparing: $read_seq vs. $element_seq\n";
                my @read_vs_element = amatch($read_seq, [ 'S1' ], ($element_seq));
                my @element_vs_read = amatch($element_seq, [ 'S1' ], ($read_seq));
                if (scalar(@read_vs_element) > 0 or scalar(@element_vs_read) > 0) {
                    $mappings->{$read_seqid}->[$c]->{count}++;
                    $found_element = 1;
                }
                ## if ($read_seq =~ /$element_seq/ or $element_seq =~ /$read_seq/) {
                ##     $mappings->{$read_seqid}->[$c]->{count}++;
                ##     print "Using an exact match: ";
                ##     print "Incremented $mappings->{$read_seqid}->[$c]->{mirbase_id} to $mappings->{$read_seqid}->[$c]->{count}\n";
                ##     $found_element = 1;
                ## }
            }
        }
    } ## Finish iterating through every read
    return($mappings);
}

=head2 C<Final_Print_Mi>

 Print out the final counts of miRNA mappings.

=item C<Arguments>

 data: Result from cross referencing the bam/genome/miRNAs.
 output: Output filename for writing a tsv of the counts.

=cut
sub Final_Print_Mi {
    my ($class, %args) = @_;
    my $final = $args{data};
    my $output = FileHandle->new(">$args{output}");
    my $hits = 0;
    foreach my $immature (keys %{$final}) {
        foreach my $mature (@{$final->{$immature}}) {
            $hits = $hits++;
            print $output "$mature->{mimat}\t$mature->{count}\n";
        }
    }
    $output->close();
    return($hits);
} ## End of Final_Print

=head2 C<Count_Alignments>

 Compare alignments across species.

 Count alignments across a host and parasite.  This function should
 give a sense of how many reads were single-mapped to the host and
 parasite along with multi-hits across both.  This was first used to
 distinguish between T. cruzi and human hits, ergo the 'Tc' as the
 default parasite pattern.

=item C<Arguments>

 input(required): Input bam file to count against.
 species(required): Set the host or parasite species to count against.
 para_patterm('^Tc'): Pattern used to differentiate parasite-mapped
  reads.
 host_pattern(''): Pattern used to differentiate host-mapped reads.
 libpath: Change the dirname of the fasta libraries.
 libtype('genome'): Change the subdirectory of the fasta libraries.

=cut
sub Count_Alignments {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        para_pattern => '^Tc',
        host_pattern => '',);
    my $result = {
        mapped => 0,
        unmapped => 0,
        multi_count => 0,
        single_count => 0,
        unmapped_count => 0,
        single_para => 0, ## Single-hit parasite
        single_host => 0, ## Single-hit host
        multi_host => 0, ## Multi-hit host -- keep in mind htseq will not count these.
        multi_para => 0, ## Ditto parasite
        single_both => 0, ## These are the dangerzone reads, 1 hit on both. -- false positives
        single_para_multi_host => 0,
        single_host_multi_para => 0,
        multi_both => 0, ## These have multi on both and are not a problem.
        zero_both => 0, ## This should stay zero.
        wtf => 0,
    };

    my %group = ( para => 0, host => 0,);
    my %null_group = ( para => 0, host => 0,);
    my $fasta = qq"$options->{libpath}/$options->{libtype}/$options->{species}.fasta";
    my $sam = Bio::DB::Sam->new(-bam => $options->{input},
                                -fasta => $fasta,);
    my @targets = $sam->seq_ids;
    my $num = scalar(@targets);
    my $bam = Bio::DB::Bam->open($options->{input});
    my $header = $bam->header;
    my $target_count = $header->n_targets;
    my $target_names = $header->target_name;
    my $align_count = 0;
    my $million_aligns = 0;
    my $alignstats = qx"samtools idxstats $options->{input}";
    my @alignfun = split(/\n/, $alignstats);
    my @aligns = split(/\t/, $alignfun[0]);
    my @unaligns = split(/\t/, $alignfun[1]);
    my $number_reads = $aligns[2] + $unaligns[3];
    my $output_name = qq"$options->{input}.out";
    my $out = FileHandle->new(">${output_name}");
    print $out "There are ${number_reads} alignments in $options->{input} made of $aligns[2] aligned reads and $unaligns[3] unaligned reads.\n";
    my $last_readid = "";
  BAMLOOP: while (my $align = $bam->read1) {
        $align_count++;
        if (($align_count % 1000000) == 0) {
            $million_aligns++;
            print $out "Finished $million_aligns million alignments out of ${number_reads}.\n";
        }
        my $seqid = $target_names->[$align->tid];
        my $readid = $align->qname;
        my $start = $align->pos;
        my $end = $align->calend - 1;
        my $cigar = $align->cigar_str;
        my $strand = $align->strand;
        my $seq = $align->query->dna;
        my $qual= $align->qual;
        if ($cigar eq '') {
            $result->{unmapped}++;
            next BAMLOOP;
        } else {
            ## Everything which follows is for a mapped read.
            $result->{mapped}++;
            my $type;
            if ($options->{para_pattern}) {
                if ($seqid =~ m/$options->{para_pattern}/) {
                    $type = 'para';
                } else {
                    $type = 'host';
                }
            } elsif ($options->{host_pattern}) {
                if ($seqid =~ m/$options->{host_pattern}/) {
                    $type = 'host';
                } else {
                    $type = 'para';
                }
            }
            if ($readid ne $last_readid) {
                ## Count up what is currently in the group, then reset it.
                my $reads_in_group = $group{host} + $group{para};
                if ($reads_in_group > 1) {
                    $result->{multi_count}++;
                } elsif ($reads_in_group == 1) {
                    $result->{single_count}++;
                } else {
                    $result->{unmapped_count}++;
                }
                if ($group{host} == 0) {
                    if ($group{para} == 0) {
                        $result->{zero_both}++;
                    } elsif ($group{para} == 1) {
                        $result->{single_para} = $result->{single_para} + $reads_in_group;
                    } elsif ($group{para} > 1) {
                        $result->{multi_para} = $result->{multi_para} + $reads_in_group;
                    } else {
                        $result->{wtf}++;
                    }
                } elsif ($group{host} == 1) {
                    if ($group{para} == 0) {
                        $result->{single_host} = $result->{single_host} + $reads_in_group;
                    } elsif ($group{para} == 1) {
                        $result->{single_both} = $result->{single_both} + $reads_in_group;
                    } elsif ($group{para} > 1) {
                        $result->{single_host_multi_para} = $result->{single_host_multi_para} + $reads_in_group;
                    } else {
                        $result->{wtf}++;
                    }
                } elsif ($group{host} > 1) {
                    if ($group{para} == 0) {
                        $result->{multi_host} = $result->{multi_host} + $reads_in_group;
                    } elsif ($group{para} == 1) {
                        $result->{single_para_multi_host} = $result->{single_para_multi_host} + $reads_in_group;
                    } elsif ($group{host} > 1) {
                        $result->{multi_both} = $result->{multi_both} + $reads_in_group;
                    } else {
                        $result->{wtf}++;
                    }
                } else {
                    $result->{wtf}++;
                }
                print "Map:$result->{mapped} Un:$result->{unmapped} SingleC:$result->{single_count} MultC:$result->{multi_count} sp:$result->{single_para} sh:$result->{single_host} mh:$result->{multi_host} mp:$result->{multi_para} SB:$result->{single_both} spmh:$result->{single_para_multi_host} shmp:$result->{single_host_multi_para} bm:$result->{multi_both} bz:$result->{zero_both} wtf:$result->{wtf}\n" if ($options->{debug});
                %group = %null_group;
                ## Now set the first read for the new group to the type of this read.
                $group{$type}++;
            } else {
                $group{$type}++;
            }
        }
        $last_readid = $readid;
    } ## End reading each bam entry
    print $out "Mapped: $result->{mapped}
Unmapped: $result->{unmapped}
Multi-mapped: $result->{multi}
Single-parasite: $result->{single_para}
Single-host: $result->{single_host}
Multi-parasite, no host: $result->{multi_para}
Multi-host, no parasite: $result->{multi_host}
DANGER Single-both: $result->{single_both}
DANGER Single-parasite, multi-host: $result->{single_para_multi_host}
DANGER Single-host, multi-parasite: $result->{single_host_multi_para}
Multi-both: $result->{both_multi}\n";
    $out->close();
    return($result);
}

=head2 C<SLSearch>

 Search a pile of reads for the trypanosome spliced leader.

 Use some simple pattern matching on a pile of reads to look for
 sequences of interest.  By default this looks for a portion of the
 spliced leader sequence from Leishmania major.  Having written and
 used this, I realized that it is ... dumb.  I should have used
 jellyfish and simply counted up the hits across a range of the SL.
 Doing this with jellyfish would also let me count up each window of
 the SL and plot a histogram of the occurrences, thus showing the
 optimal starting point to discriminate the SL as opposed to my ad hoc
 choice of position generated via repeated invocations of 'grep | wc'.

 With that in mind, this counts up the number of times a string of
 interest appears in the input files in the forward and RC directions,
 and prints that in an easy-to-read format.

=item C<Arguments>

 input(required): Set of input files to read.
 search('AGTTTCTGTACTTTATTGG'): SL substring to search.
 jmem(24): Memory to allocate for this task.
 jprefix(50): Default jobname/output prefix.

=cut
sub SLSearch {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jmem => 24,
        jprefix => '50',
        search => 'AGTTTCTGTACTTTATTGG',);
    my $output_dir =  qq"outputs/$options->{jprefix}SL_search";
    my $output_made = make_path($output_dir);
    my $comment = '## Search for SL sub-sequences.';
    my $jstring = qq?
use Bio::Adventure::Count;
my \$result = \$h->Bio::Adventure::Count::SLSearch_Worker(
  input => '$options->{input}',
  output => '${output_dir}',
  jprefix => '$options->{jprefix}',
  jname => 'slsearch',);
?;
    my $slsearch = $class->Submit(
        comment => $comment,
        input => $options->{input},
        jcpus => 1,
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
        jname => 'slsearch',
        jprefix => $options->{jprefix},
        jstring => $jstring,
        language => 'perl',
        output => $output_dir,);
    $class->{language} = 'bash';
    $class->{shell} = '/usr/bin/env bash';
    return($slsearch);
}

=head2 C<SLSearch_Worker>

 This function does the actual work for SLSearch().

=cut
sub SLSearch_Worker {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'output'],
        jmem => 24,
        jprefix => '50',
        search => 'AGTTTCTGTACTTTAT',);
    ## Ideally, I would like to run this function on the mapping
    ## directories and search the aligned/unaligned sequences
    ## separately.
    ## Unfortunately, in our cleaning, we deleted those files.
    ## Therefore it (for the moment) will instead just read over the
    ## trimmed files.
    ## I think I will write the following: the read IDs for every
    ## read which matches the forward or RC versions of the SL and
    ## a summary txt file of the result.
    my @input_lst = ();
    if ($options->{input} =~ /\:|\;|\,|\s+/) {
        my @tmp_lst = split(/\:|\;|\,|\s+/, $options->{input});
        for my $i (@tmp_lst) {
            push (@input_lst, $i) if (-r $i);
        }
    } else {
        if (-r $options->{input}) {
            push(@input_lst, $options->{input});
        }
    }
    if (scalar(@input_lst) == 0) {
        return(undef);
    }

    my $fwd_search = $options->{search};
    my $rc_search = reverse($fwd_search);
    $rc_search =~ tr/ATGCUatgcu/TACGAtacga/;
    my $search_length = length($fwd_search);

    my $log_fh = FileHandle->new(">$options->{output}/slsearch_log.txt");
    print $log_fh qq"Starting search for a portion of the SL sequence: $options->{search}
in the file(s): $options->{input}.\n";
    ## Now start the main loop, open file handles for a main log
    ## and for the per-input outputs.  Create a couple of global counters.

    my %global_search_result = (
        found => 0,  ## The total number of observed SL.
        fwd_found => 0,  ## The number found in the forward orientation.
        rc_found => 0,  ## The number of revcomp found.
        searched => 0); ## The total number of sequences searched.
    ## One may reasonably ask why all_found is not just the sum of the next two.
    ## I am thinking to catch the pathological case where a single input sequence
    ## has both the forward and RC sequences.

    for my $i (@input_lst) {
        my $ind_name = basename($i, ('.gz', '.bz2', '.xz'));
        my $format = 'Fastq';
        $format = 'Fasta' if ($i =~ /\.fasta/);

        $ind_name = basename($ind_name, ('.fastq', '.fasta'));
        my $ind_result = FileHandle->new(">$options->{output}/${ind_name}.tsv");
        ## I want to save the position information too so that I can see if I my search string
        ## is too close to the beginning of the read.
        print $ind_result "ReadID\tposition\torientation\n";
        my %ind_search_result = (
            found => 0,
            fwd_found => 0,
            rc_found => 0,
            searched => 0);
        my $reader = FileHandle->new("less ${i} |");
        my $seqio = Bio::SeqIO->new(-format => $format, -fh => $reader);
      FSA: while (my $seq = $seqio->next_seq) {
          $ind_search_result{searched}++;
          $global_search_result{searched}++;
          my $seq_id = $seq->id;
          my $read_seq = $seq->seq;
          my $fwd_end = undef;
          if ($read_seq =~ m/$fwd_search/g) {
              $fwd_end = pos($read_seq);
          }
          my $rc_end = undef;
          if ($read_seq =~ m/$rc_search/g) {
              $rc_end = pos($read_seq);
          }
          ## Get out if we do not find the SL portion.
          if (!defined($fwd_end) && !defined($rc_end)) {
              next FSA;
          }

          if ($fwd_end) {
              my $fwd_start = $fwd_end - ($search_length - 1);
              print $ind_result "${seq_id}\t${fwd_start}\tFWD\n";
              $ind_search_result{found}++;
              $ind_search_result{fwd_found}++;
              $global_search_result{found}++;
              $global_search_result{fwd_found}++;
          }

          if ($rc_end) {
              my $rc_start = $rc_end - ($search_length - 1);
              print $ind_result "${seq_id}\t${rc_start}\tRC\n";
              $ind_search_result{found}++;
              $ind_search_result{rc_found}++;
              $global_search_result{found}++;
              $global_search_result{rc_found}++;
          }
      } ## End reading the input fastq/fasta file.

        print $log_fh qq"
${ind_name} results:
  sequences searched: $ind_search_result{searched}
  subsequences observed: $ind_search_result{found}
  forward observed: $ind_search_result{fwd_found}
  reverse-complement observed: $ind_search_result{rc_found}\n";
        $ind_result->close();
        $reader->close();
    } ## End of the input loop
        print $log_fh qq"
Total results:
  sequences searched: $global_search_result{searched}
  subsequences observed: $global_search_result{found}
  forward observed: $global_search_result{fwd_found}
  reverse-complement observed: $global_search_result{rc_found}\n";
    $log_fh->close();
}

=head1 AUTHOR - atb

Email <abelew@gmail.com>

=head1 SEE ALSO

L<htseq-count> L<Bio::DB::Sam> L<Bio::SeqIO>

=cut

1;
