package Bio::Adventure::RNASeq_Count;
## LICENSE: gplv2
## ABSTRACT:  Kitty!
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use Bio::Tools::GFF;
use File::Basename;
use File::Path qw"make_path";
use File::Which qw"which";
use String::Approx qw"amatch";

=head1 NAME

    Bio::Adventure::RNASeq_Count - Perform Sequence alignments counting with HTSeq

=head1 SYNOPSIS

    use Bio::Adventure;
    my $hpgl = new Bio::Adventure;
    $hpgl->HT_Multi();

=head2 Methods

=over 2

=item C<HT_Multi>

    some stuff here

=cut
sub HT_Multi {
    my ($class, %args) = @_;
    my $check = which('htseq-count');
    die("Could not find htseq in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        jname => '',
        required => ["species", "htseq_input"],
    );
    my $species = $options->{species};
    my $htseq_input = $options->{htseq_input};
    my $job_basename = $options->{jname};
    my @jobs = ();

    my $script_suffix = qq"";
    if ($args{suffix}) {
        $script_suffix = $args{suffix};
    }
    my $jprefix = "";
    $jprefix = $args{jprefix} if ($args{jprefix});
    my %gff_types = (
        antisense => ["exon", "ID"],
        exon => ["exon", "ID"],
        fiveputr => ["five_prime_UTR", "ID"],
        interCDS => ["CDS", "ID"],
        linc => ["gene", "ID"],
        mi => ["exon", "ID"],
        misc => ["transcript", "ID"],
        nmd => ["exon", "ID"],
        operons => ["gene", "ID"],
        pseudo => ["exon", "ID"],
        rintron => ["exon", "ID"],
        rrna => ["gene", "ID"],
        sn => ["exon", "ID"],
        sno => ["exon", "ID"],
        threeputr => ["three_prime_UTR", "ID"],
    );
    my $htseq_runs = 0;
    ## my $aligntype = $args{aligntype};
    ## my $out_prefix = $input;
    my $out_basename = basename($htseq_input, (".bam", ".sam"));
    my $out_suffixdir = dirname($htseq_input); ## The directory name we will change from bowtie/tophat/whatever to htseq
    foreach my $type (keys %gff_types) {
        my $gff = qq"$options->{libdir}/genome/${species}_${type}.gff";
        my $gtf = $gff;
        $gtf =~ s/\.gff/\.gtf/g;
        my $output = qq"${out_suffixdir}/${out_basename}_${type}.count";
        my $htseq_jobname = qq"ht${type}_$options->{species}";
        if ($options->{jname}) {
            $htseq_jobname = qq"ht${type}_$options->{jname}";
        }
        if (-r "$gff") {
            print "Found $gff, performing htseq with it.\n";
            my $ht = Bio::Adventure::RNASeq_Count::HTSeq(
                $class,
                htseq_gff => $gff,
                htseq_input => $htseq_input,
                htseq_type => $gff_types{$type},
                depends => $options->{depends},
                jname => $htseq_jobname,
                job_output => $output,
                jprefix => $jprefix,
                postscript => $options->{postscript},
                prescript => $options->{prescript},
                queue => 'throughput',
                suffix => $options->{suffix},
            );
            push(@jobs, $ht);
            $htseq_runs++;
        } elsif (-r "$gtf") {
            print "Found $gtf, performing htseq with it.\n";
            my $ht = Bio::Adventure::RNASeq_Count::HTSeq(
                $class,
                htseq_gff => $gtf,
                htseq_input => $htseq_input,
                htseq_type => $gff_types{$type},
                depends => $options->{depends},
                jname => $htseq_jobname,
                job_output => $output,
                jprefix => $options->{jprefix},
                postscript => $options->{postscript},
                prescript => $options->{prescript},
                queue => 'throughput',
                suffix => $options->{suffix},
            );
            push(@jobs, $ht);
            $htseq_runs++;
        }
        ##else {
        ##print "Did not find ${gff} nor ${gtf}, skipping htseq_${type}.\n";
        ##}
    }                           ## End foreach type
    ## Also perform a whole genome count
    my $gff = qq"$options->{libdir}/genome/${species}.gff";
    my $gtf = $gff;
    $gtf =~ s/\.gff/\.gtf/g;
    my $output = $htseq_input;
    $output =~ s/\.bam$/\.count/g;
    my $htall_jobname = qq"htall_$options->{species}";
    if ($options->{jname}) {
        $htall_jobname = qq"htall_$options->{jname}";
    }
    if (-r "$gff") {
        print "Found $gff, performing htseq_all with it.\n";
        my $ht = Bio::Adventure::RNASeq_Count::HTSeq(
            $class,
            htseq_gff => $gff,
            htseq_id => $options->{htseq_id},
            htseq_input => $htseq_input,
            htseq_type => $options->{htseq_type},
            depends => $options->{depends},
            jname => $htall_jobname,
            job_output => $output,
            jprefix => $jprefix,
            postscript => $options->{postscript},
            prescript => $options->{prescript},
            queue => 'throughput',
            suffix => $options->{suffix},
        );
        push(@jobs, $ht);
    } elsif (-r "$gtf") {
        print "Found $gtf, performing htseq_all with it.\n";
        my $ht = Bio::Adventure::RNASeq_Count::HTSeq(
            $class,
            htseq_gff => $gff,
            htseq_input => $htseq_input,
            htseq_type => "none",
            depends => $options->{depends},
            jname => $htall_jobname,
            job_output => $output,
            jprefix => $jprefix,
            postscript => $args{postscript},
            prescript => $args{prescript},
            queue => 'throughput',
            suffix => $args{suffix},
        );
        push(@jobs, $ht);
    } else {
        print "Did not find ${gff} nor ${gtf}, not running htseq_all.\n";
    }
    return(\@jobs);
}

sub HT_Types {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
    );
    my $my_type = $options->{type};
    $my_type = "" unless($my_type);
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
    my $max_type = "";
    my $max = 0;
    my $max_canonical = 0;
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

    my $max_can = "";
    foreach my $can (keys %found_canonical) {
        if (!defined($found_canonical{$can})) {
            $found_canonical{$can} = 0;
        }
        if ($found_canonical{$can} > $max_canonical) {
            $max_can = $can;
            $max_canonical = $found_types{$can};
        }
    }                           ## End the loop
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

=head2

    HTSeq()

=cut
sub HTSeq {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ["species", "htseq_stranded", "htseq_args",],
        htseq_input => $class->{options}->{input},
        depends => '',
        jname => '',
        jprefix => '',
        libtype => 'genome',
    );
    my $basename = $options->{basename};
    my $stranded = $options->{htseq_stranded};
    my $htseq_type = $options->{htseq_type};
    my $htseq_id = $options->{htseq_id};
    my $output = $options->{job_output};
    my $htseq_jobname = qq"hts_$options->{species}";
    if ($options->{jname}) {
        $htseq_jobname = $options->{jname};
    }
    my $htseq_input = $options->{htseq_input};
    my $error = $htseq_input;
    my $gff = qq"$options->{libdir}/$options->{libtype}/$options->{species}.gff";
    $gff = $args{htseq_gff} if ($args{htseq_gff});
    my $gtf = $gff;
    $gtf =~ s/\.gff/\.gtf/g;
    my $htseq_args = "";
    if ($options->{suffix}) {
        my $suffix = $options->{suffix};
        $output =~ s/\.bam/_${suffix}\.count/g;
    } else {
        $output =~ s/\.bam/\.count/g;
    }
    $error =~ s/\.bam/_htseq\.err/g;
    if (!-r "${gff}" and !-r "${gtf}") {
        die("Unable to read ${gff} nor ${gtf}, please fix this and try again.\n");
    }
    my $annotation = $gtf;
    if (!-r "${gtf}") {
        $annotation = $gff;
    }
    ## Set the '-t FEATURETYPE --type' argument used by htseq-count
    ## This may be provided by a series of defaults in %HT_Multi::gff_types, overridden by an argument
    ## Or finally auto-detected by HT_Types().
    ## This is imperfect to say the nicest thing possible, I need to consider more appropriate ways of handling this.
    my $htseq_type_arg = "";
    my $htseq_id_arg = "";
    if (!defined($htseq_type) or $htseq_type eq '' or $htseq_type eq 'auto' or
            !defined($htseq_id) or $htseq_id eq '' or $htseq_id eq 'auto') {
        $htseq_type = Bio::Adventure::RNASeq_Count::HT_Types(
            $class,
            annotation => $annotation,
            type => $htseq_type,
        );
        $htseq_type_arg = qq" --type $htseq_type->[0]";
        $htseq_id_arg = qq" --idattr $htseq_type->[1]"
    } elsif (ref($htseq_type) eq "ARRAY") {
        $htseq_type_arg = qq" --type $htseq_type->[0]";
        $htseq_id_arg = qq" --idattr $htseq_type->[1]"
    } elsif ($htseq_type eq 'none') {
        $htseq_type_arg = qq"";
        $htseq_id_arg = qq"";
    } else {
        $htseq_type_arg = qq" --type ${htseq_type}";
        $htseq_id_arg = qq" --idattr ${htseq_id}";
    }

    ## Much like samtools, htseq versions on travis are old.
    ## Start with the default, non-stupid version.
    my $htseq_version = qx"htseq-count -h | grep version";
    my $htseq_invocation = qq!htseq-count -q -f bam -s ${stranded} ${htseq_id_arg} ${htseq_type_arg} \\!;
    if ($htseq_version =~ /0\.5/) {
        ## Versions older than 0.6 are stupid.
        $htseq_invocation = qq!samtools view ${htseq_input} | htseq-count -q -s ${stranded} ${htseq_id_arg} ${htseq_type_arg} \\!;
        $htseq_input = '-';
    }
    my $jstring = qq!
${htseq_invocation}
  ${htseq_input} \\
  ${annotation} \\
  1>${output} 2>${error} && \\
    xz -f -9e ${output}
!;
    my $comment = qq!## Counting the number of hits in ${htseq_input} for each feature found in ${annotation}
## Is this stranded? ${stranded}.  The defaults of htseq are:
## $options->{htseq_args}->{default}
!;
    my $htseq = $class->Submit(
        comment => $comment,
        input => $htseq_input,
        depends => $options->{depends},
        jname => $htseq_jobname,
        job_output => $output,
        jprefix => $options->{jprefix},
        jstring => $jstring,
        mem => 6,
        postscript => $args{postscript},
        prescript => $args{prescript},
        queue => "throughput",
    );
    return($htseq);
}

sub Mi_Map {
    my ($class, %args) = @_;

    eval "use Bio::DB::Sam; 1;";
    my $options = $class->Get_Vars(
        args => \%args,
        required => ["mirbase_data", "mature_fasta",
                     "mi_genome", "bamfile",]
    );
    my $bam_base = basename($options->{bamfile}, (".bam", ".sam"));
    my $pre_map = qq"${bam_base}_mirnadb.txt";
    my $final_map = qq"${bam_base}_mature.count";

    ## Step 1:  Read the mature miRNAs to get the global IDs and sequences of the longer transcripts.
    ## This provides the IDs as 'mm_miR_xxx_3p' and short sequences ~21 nt.
    print "Starting to read miRNA sequences.\n";
    my $sequence_db = Bio::Adventure::RNASeq_Count::Read_Mi(
        $class,
        seqfile => $options->{mature_fasta},
    );

    ## Step 2:  Collate those IDs against the set of miRNA_transcripts->miRNA_mature IDs using
    ## the tab delimited file downloaded from mirbase.org
    print "Starting to read miRNA mappings.\n";
    my $sequence_mappings = Bio::Adventure::RNASeq_Count::Read_Mappings_Mi(
        $class,
        mappings => $options->{mirbase_data},
        output => $pre_map,
        seqdb => $sequence_db,
    );
    ## At this point, we have brought together mature sequence/mature IDs to parent transcript IDs
    ## When we read the provided bam alignment, we will simultaneously read the immature miRNA database
    ## and use them to make a bridge from the mature sequence/IDs to immature miRNA IDs.
    my $final_hits = Bio::Adventure::RNASeq_Count::Read_Bam_Mi(
        $class,
        mappings => $sequence_mappings,
        mi_genome => $options->{mi_genome},
        bamfile => $options->{bamfile},
    );
    my $printed = Bio::Adventure::RNASeq_Count::Final_Print_Mi(
        $class,
        data => $final_hits,
        output => $final_map,
    );

} ## End of Mi_Map

## Quick and dirty sequence reader
sub Read_Mi {
    my ($class, %args) = @_;
    my $fasta = Bio::SeqIO->new(-file => $args{seqfile}, -format => "Fasta");
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
            ##print "Found $id across the mappings.\n";
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
                ## $newdb->{$new_id} = \@hit_list;
                $newdb->{$ensembl_id} = \@hit_list;
                print $out "$id; $element->{mirbase_id}; $element->{ensembl}; $element->{hit_id}; $element->{mimat}; $element->{sequence}\n";
            } else {
                @hit_list = ();
                push(@hit_list, $element);
                ## $newdb->{$new_id} = \@hit_list;
                $newdb->{$ensembl_id} = \@hit_list;
            }
        } else {
            ## print "Did not find $id\n";
        }
    }
    $out->close();
    return($newdb);
}

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
                    ##print "Found hit with amatch against $element_seq: ";
                    ##print "We need to get the slot in position: {ensembl_id} -> [c] -> {id}\n";
                    ##print "is read_seqid ensembl? $read_seqid  \n";
                    ##print "Incremented $mappings->{$read_seqid}->[$c]->{id} to $mappings->{$read_seqid}->[$c]->{count}\n";
                    $found_element = 1;
                }
                ## if ($read_seq =~ /$element_seq/ or $element_seq =~ /$read_seq/) {
                ##     $mappings->{$read_seqid}->[$c]->{count}++;
                ##     print "Using an exact match: ";
                ##     print "Incremented $mappings->{$read_seqid}->[$c]->{mirbase_id} to $mappings->{$read_seqid}->[$c]->{count}\n";
                ##     $found_element = 1;
                ## }
            }
            ##if ($found_element == 0) {
            ##    print "No map for $read_seq\n";
            ##}
        }
    }                           ## Finish iterating through every read
    return($mappings);
}

sub Final_Print_Mi {
    my ($class, %args) = @_;
    my $final = $args{data};
    my $output = FileHandle->new(">$args{output}");
    my $hits = 0;
    foreach my $immature (keys %{$final}) {
        foreach my $mature (@{$final->{$immature}}) {
            $hits = $hits++;
            print $output "$mature->{mimat} $mature->{count}\n";
        }
    }
    $output->close();
    return($hits);
}                               ## End of Final_Print

=back

=head1 AUTHOR - atb

Email <abelew@gmail.com>

=head1 SEE ALSO

    L<htseq-count>

=cut

1;
