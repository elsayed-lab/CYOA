package Bio::Adventure::Convert;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use Bio::FeatureIO;
use Bio::Tools::GFF;
use File::Basename;
use File::Which qw"which";
use List::MoreUtils qw"uniq";
use TryCatch;

=head1 NAME

Bio::Adventure::Convert - Perform conversions between various formats.


=head1 SYNOPSIS

The functions here handle the various likely conversions one may perform.
sam to bam, gff to fasta, genbank to fasta, etc.

=head1 METHODS

=head2 C<Gff2Fasta>

$hpgl->Gff2Fasta(genome => 'something.fasta', gff => 'something.gff')
will read the genome's fasta and annotation files and return a set
of fasta files which include the coding sequence regions which are
members of args{feature_type} in the gff file.

It writes one file of amino acid sequence and one of nucleotides.
Upon completion, it returns the number of entries written.

=over

=item I<tag> - Which gff types to add to the fasta file?

=back

=cut
sub Gff2Fasta {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args,
                                   required => ['gff', 'genome'],
                                   tag => 'gene_id',
                               );
    my $genome = $options->{genome};
    my $gff = $options->{gff};
    my $tag = $options->{tag};
    my $genome_basename = basename($genome, ('.fasta'));
    my $chromosomes = $class->Read_Genome_Fasta(genome => $genome);
    my $gff_handle = FileHandle->new("less ${gff} |");
    my $out_fasta_amino = FileHandle->new(qq">${genome_basename}_cds_aa.fasta");
    my $out_fasta_nt = FileHandle->new(qq">${genome_basename}_cds_nt.fasta");
    my @tag_list = ('ID','gene_id','locus_tag','transcript_id');
    my $feature_type = Bio::Adventure::RNASeq_Count::HT_Types($class, annotation => $gff);
    $feature_type = $feature_type->[0];
    my $annotation_in = Bio::Tools::GFF->new(-fh => $gff_handle, -gff_version => 3);
    my $features_written = 0;
    my $features_read = 0;

  LOOP: while(my $feature = $annotation_in->next_feature()) {
        $features_read++;
        next LOOP unless ($feature->{_primary_tag} eq $feature_type);
        my $location = $feature->{_location};
        my $start = $location->start();
        my $end = $location->end();
        my $strand = $location->strand();
        my @something = $feature->each_tag_value();
        my $chosen_tag = $options->{tag};
        my @ids;
        try {
            @ids = $feature->each_tag_value($chosen_tag);
        } catch {
            next LOOP;
        }
        my $gff_chr = $feature->{_gsf_seq_id};
        my $gff_string = $annotation_in->{$gff_chr};
        if (!defined($chromosomes->{$gff_chr})) {
            print STDERR "Something is wrong with $gff_chr\n";
            next LOOP;
        }
        my $id = "";
        foreach my $i (@ids) {
            $id .= "$i ";
        }
        $id =~ s/(\s+$)//g;
        my $genome_obj = $chromosomes->{$gff_chr}->{obj};
        my $cds = $genome_obj->subseq($start, $end);
        if ($strand == -1) {
            $cds = reverse($cds);
            $cds =~ tr/ATGCatgc/TACGtacg/;
        }
        my $seq_obj = Bio::Seq->new();
        $seq_obj->seq($cds);
        my $aa_cds = $seq_obj->translate->seq();
        $aa_cds = join("\n", ($aa_cds =~ m/.{1,80}/g));
        $cds = join("\n", ($cds =~ m/.{1,80}/g));
        print $out_fasta_amino ">${gff_chr}_${id}
${aa_cds}
";
        print $out_fasta_nt ">${gff_chr}_${id}
${cds}
";
        $features_written++;
    }                           ## End LOOP
    $gff_handle->close();
    $out_fasta_amino->close();
    $out_fasta_nt->close();
    return($features_written);
}

=head2 C<Read_GFF>

use Bio::Tools::GFF to gather information from a gff file.

=cut
sub Read_GFF {
    my ($class, %args) = @_;
    my $chromosomes = $args{chromosomes};
    my $gff = FileHandle->new("<$args{gff}");

    my $annotation_in = Bio::Tools::GFF->new(-fh => \$gff, -gff_version => 3);
    my $gff_out = {};
    print "Starting to read gff: $args{gff}\n";
  LOOP: while(my $feature = $annotation_in->next_feature()) {
        next LOOP unless ($feature->{_primary_tag} eq $args{gff_type});
        my $location = $feature->{_location};
        my $start = $location->start();
        my $end = $location->end();
        my $strand = $location->strand();
        my @ids = $feature->each_tag_value("ID");
        my $id = "";
        my $gff_chr = $feature->{_gsf_seq_id};
        my $gff_string = $annotation_in->gff_string($feature);
        if (!defined($chromosomes->{$gff_chr})) {
            print STDERR "Something is wrong with $gff_chr\n";
            next LOOP;
        }
        foreach my $i (@ids) {
            $i =~ s/^cds_//g;
            $i =~ s/\-\d+$//g;
            $id .= "$i ";
        }
        $id =~ s/\s+$//g;
        my @gff_information = split(/\t+/, $gff_string);
        my $description_string = $gff_information[8];
        my $orf_chromosome = $gff_chr;
        my $annot = {id => $id,
                     start => $start, ## Genomic coordinate of the start codon
                     end => $end,     ## And stop codon
                     strand => $strand,
                     description_string => $description_string,
                     chromosome => $gff_chr,
                 };
        $gff_out->{$gff_chr}->{$id} = $annot;
    }                           ## End looking at every gene in the gff file
    $gff->close();
    return($gff_out);
}

=head2 C<Gff2Gtf>

$hpgl->Gff2Gtf(gff => 'mmusculus.gff')
reads a given gff file and writes a gtf file from the features
found therein.  It returns the number of features written.

note that I only use gtf files for tophat, thus they must have a tag 'transcript_id'!
This is woefully untrue for the tritrypdb gff files.  Thus I need to have a regex in this
to make sure that web_id or somesuch is changed to transcript_id.

=cut
sub Gff2Gtf {
    my ($class, %args) = @_;
    my $input = $args{gff};

    my ($name, $path, $suffix) = fileparse($input, qr/\.gff/);
    my $out_file = $path . $name . ".gtf";

    my $in_gff = Bio::FeatureIO->new('-file' => "$input",
                                     '-format' => 'GFF',
                                     '-version' => 3);
    my $out_gtf = FileHandle->new(">${out_file}");
    local $| = 1;
    ##my $out_gtf = new FileHandle();
    ##$out_gtf->open(">$out_file");

    my $features_written = 0;
    my $in_features = 0;
  FEATURES: while (my $feature = $in_gff->next_feature()) {
        $in_features++;
        ## the output handle is reset for every file
        ## According to UCSC, GTF file contains 9 column:
        ## <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes]
        ## An example working gtf line: (except the double spaces are tabs...)
        ## TcChr20-S  TriTrypDB  CDS  14765  15403  .  -  0    gene_id "cds_TcCLB.397937.10-1"; transcript_id "cds_TcCLB.397937.10-1";

        my @tags = $feature->get_all_tags();
        my $seqid = $feature->seq_id();
        my $location = $feature->location();
        my $start = $feature->start();
        my $end = $feature->end();
        my $strand = $feature->strand();
        if ($strand == 1) {
            $strand = '+';
        } else {
            $strand = '-';
        }
        my $source = $feature->source_tag();
        my $primary_id = $feature->primary_tag();
        my $phase = '.';
        if ($primary_id eq 'CDS') {
            $phase = $feature->phase();
        }
        next FEATURES if ($primary_id eq 'supercontig' or $primary_id eq 'region' or $primary_id eq 'chromosome');
        ## my $string = qq"${seqid}\t${source}\tCDS\t$start\t$end\t.\t$strand\t0\t";
        my $string = qq"${seqid}\t${source}\t${primary_id}\t${start}\t${end}\t.\t${strand}\t${phase}\t";

        my $last_column = "";
        my $annot = $feature->{annotation};
        my $stringified;
        my $transcript_string = "";
        my $geneid_string = "";
        foreach my $key ($annot->get_all_annotation_keys()) {
            my @values = $annot->get_Annotations($key);
            foreach my $value (@values) {
                $stringified = $value->{value};
                $stringified = $value->{term}->{name} unless($stringified);
                $last_column .= qq!${key} "${stringified}"; !;
            }
            if ($key eq 'ID') {
                $transcript_string = qq!transcript_id "${stringified}"; !;
            }
            if ($key eq 'seq_id') {
                $geneid_string = qq!gene_id "${stringified}"; !;
            }
        }
        my $seq_string = $string . ${transcript_string} . ${geneid_string} . $last_column;
        $seq_string =~ s/\s+$//g;
        print $out_gtf $seq_string;
        $features_written++;
    }                           ## End iterating over every FEATURES
    $out_gtf->close();
    ##close(OUT_GTF);
    return($features_written);
}

=head2 C<Sam2Bam>

$hpgl->Sam2Bam();
Used to invoke samtools to take the sam output from bowtie/bwa and
convert it to an compressed-sorted-indexed bam alignment.

This function just calls $class->Samtools(), but has a little logic
to see if the invocation of this is for an extent .sam file or
calling it on an existing .fastq(.gz), in which case one must
assume it is being called on one or more files in the bowtie_out/
directory and which start with $basename, include something like
-trimmed-v0M1.sam.

=cut
sub Sam2Bam {
    my ($class, %args) = @_;
    my $check = which('samtools');
    die("Could not find samtools in your PATH.") unless($check);
    my $options = $class->Get_Vars(args => \%args,
                                   required => ["species", "input"]);
    my $basename = $options->{basename};
    my @input_list = ();
    my $depends = "";
    if ($options->{input}) {
        push(@input_list, $options->{input});
    } elsif (-r $options->{input} and $options->{input} =~ /\.sam$/) {
        push(@input_list, $options->{input});
    } elsif (-r $options->{input} and $options->{input} =~ /\.fastq$/) {
        if (-d "bowtie_out") {
            find({ wanted => sub { push(@input_list, "bowtie_out/$_") if ($_ =~ /\.sam$/); }, follow => 1 }, 'bowtie_out/');

        } else {
            foreach my $k (%{$options->{bt_args}}) {
                my $output_string = "bowtie_out/${basename}-${k}.sam";
                push(@input_list, $output_string);
            }
            my $bt = Bio::Adventure::RNASeq_Map::Bowtie($class, %args);
            $depends = $bt->{pbs_id};
        }
    } else {
        die("I don't know what to do without a .fastq file nor a .sam file.\n");
    }
    my $sam = Bio::Adventure::Convert::Samtools($class, %args,
                                                depends => $depends,
                                                sam => \@input_list);
    return($sam);
}

=head2 C<Samtools>

$hpgl->Samtools() calls (in order): samtools view, samtools sort,
and samtools index.  Upon completion, it invokes bamtools stats to
see what the alignments looked like.

It explicitly does not pipe one samtools invocation into the next,
not for any real reason but because when I first wrote it, it
seemed like the sorting was taking too long if I did not already
have the alignments in a bam file.

=cut
sub Samtools {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args,
                                   required => ['input', 'species'],
                                   jname => 'sam',
                                   jprefix => '',
                               );

    my $job_basename = $options->{job_basename};
    my $input = $options->{input};

    my $output = $input;
    $output =~ s/\.sam$/\.bam/g;
    my $sorted = $input;
    $sorted =~ s/\.sam$//g;
    $sorted = qq"${sorted}-sorted";
    my $paired = $sorted;
    $paired =~ s/\-sorted/\-paired/g;
    ## Add a samtools version check because *sigh*
    my $samtools_version = qx"samtools 2>&1 | grep Version";
    ## Start out assuming we will use the new samtools syntax.
    my $samtools_first = qq"samtools view -u -t $options->{libdir}/genome/$options->{species}.fasta \\
  -S ${input} -o ${output} \\
  2>${output}.err 1>${output}.out && \\";
    my $samtools_second = qq"  samtools sort -l 9 ${output} -o ${sorted}.bam \\
  2>${sorted}.err 1>${sorted}.out && \\";
    ## If there is a 0.1 in the version string, then use the old syntax.
    if ($samtools_version =~ /0\.1/) {
        $samtools_first = qq"samtools view -u -t $options->{libdir}/genome/$options->{species}.fasta \\
  -S ${input} 1>${output} && \\";
        $samtools_second = qq"  samtools sort -l 9 ${output} ${sorted} \\
  2>${sorted}.err 1>${sorted}.out && \\";
    }
    my $jstring = qq!
if \$(test \! -r ${input}); then
    echo "Could not find the samtools input file."
    exit 1
fi
${samtools_first}
${samtools_second}
  rm ${output} && \\
  rm ${input} && \\
  mv ${sorted}.bam ${output} && \\
  samtools index ${output}
## The following will fail if this is single-ended.
samtools view -b -f 2 -o ${paired}.bam ${output} && \\
  samtools index ${paired}.bam
bamtools stats -in ${output} 2>${output}.stats 1>&2 && \\
  bamtools stats -in ${paired}.bam 2>${paired}.stats 1>&2
!;
    my $comment = qq!## Converting the text sam to a compressed, sorted, indexed bamfile.
## Also printing alignment statistics to ${output}.stats
## This job depended on: $options->{depends}!;
    my $samtools = $class->Submit(comment => $comment,
                                  depends => $options->{depends},
                                  input => $input,
                                  jname => $options->{jname},
                                  job_output => qq"${output}",
                                  job_paired => qq"${paired}.bam",
                                  jprefix => $options->{jprefix},
                                  jstring => $jstring,
                                  postscript => $options->{postscript},
                                  prescript => $options->{prescript},
                              );
    return($samtools);
}

=head2 C<Gb2Gff>

$hpgl->Gb2Gff() takes a genbank genome file and splits it into:
A genomic fasta file, CDS fasta, peptide fasta, gff file of all
entries, CDS, and interCDS regions.

=cut
sub Gb2Gff {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args, required => ['input']);
    my $input = $options->{input};

    my @suffix = (".gb", ".genbank");
    my $base = basename($input, @suffix);

    my $in = FileHandle->new("less ${input} |");
    my $seqio = Bio::SeqIO->new(-format => 'genbank', -fh => $in);
    my $seq_count = 0;
    my $total_nt = 0;
    my $feature_count = 0;
    my $fasta = Bio::SeqIO->new(-file => qq">${base}.fasta", -format => 'fasta', -flush => 0);
    my $gffout = Bio::Tools::GFF->new(-file => ">${base}_all.gff", -gff_version => 3);
    my $gene_gff = Bio::Tools::GFF->new(-file => ">${base}.gff", -gff_version => 3);
    my $cds_gff = Bio::Tools::GFF->new(-file => ">${base}_cds.gff", -gff_version => 3);
    my $inter_gffout = Bio::Tools::GFF->new(-file => ">${base}_interCDS.gff", -gff_version => 3);
    my $cds_fasta = Bio::SeqIO->new(-file => qq">${base}_cds.fasta", -format => 'Fasta');
    my $pep_fasta = Bio::SeqIO->new(-file => qq">${base}_pep.fasta", -format => 'Fasta');
    while (my $seq = $seqio->next_seq) {
        $seq_count++;
        $total_nt = $total_nt + $seq->length();
        $fasta->write_seq($seq);
        print "Wrote $seq_count features.\n";
        my @feature_list = ();
        ## defined a default name

        ## $feature is an object of type: Bio::SeqFeatureI
        ## Some of the things you can call on it are:
        ## display_name(), primary_tag(), has_tag(something), get_tag_values(), get_all_tags(), attach_seq(Bio::Seq),
        ## seq(), entire_seq(), seq_id(), gff_string(), spliced_seq(), primary_id(), phase(),
        foreach my $feature ($seq->top_SeqFeatures()) {
            $gffout->write_feature($feature);
            my $feat_count = 0;
        }
        for my $feat_object ($seq->get_SeqFeatures) {
            $feature_count++;
            my $feat_count = 0;
            if ($feat_object->primary_tag eq "CDS") {
                $cds_gff->write_feature($feat_object);
                $feat_count++;
                my $id_string = "";
                my $id_hash = {};
                foreach my $thing (keys %{$feat_object->{_gsf_tag_hash}}) {
                    my @arr = @{$feat_object->{_gsf_tag_hash}->{$thing}};
                    $id_hash->{$thing} = $arr[0];
                }
                my $id = qq"$id_hash->{protein_id}";
                my $desc = "";
                if (defined($id_hash->{gene})) {
                    $desc .= "$id_hash->{gene} ; ";
                }
                if (defined($id_hash->{db_xref})) {
                    $desc .= "$id_hash->{db_xref} ; ";
                }
                if (defined($id_hash->{product})) {
                    $desc .= "$id_hash->{product}";
                }
                $desc .= "\n";
                my $start = $feat_object->start;
                my $end = $feat_object->end;
                my $len = $feat_object->length;
                my $seq = $feat_object->spliced_seq->seq;
                my $pep = $feat_object->spliced_seq->translate->seq;
                my $size = {
                    id => $id,
                    start => $start,
                    end => $end,
                    len => $len,
                    feat => $feat_object,
                };
                push(@feature_list, $size);
                my $seq_object = Bio::PrimarySeq->new(-id => $id, -seq => $seq, description => $desc);
                my $pep_object = Bio::PrimarySeq->new(-id => $id, -seq => $pep, description => $desc);
                $cds_fasta->write_seq($seq_object);
                my $ttseq = $seq_object->seq;
                $pep_fasta->write_seq($pep_object);
            } elsif ($feat_object->primary_tag eq "gene") {
                $gene_gff->write_feature($feat_object);
            }
        } ## End looking at every feature and putting them into the @feature_list
        ## Now make a hash from it and fill in the inter-cds data

        my %inter_features = ();
        my ($p_st, $p_en, $c_st, $c_en, $inter_start, $inter_end) = 0;
      INTER: for my $c (0 .. $#feature_list) {
            my $p_size = {};
            if ($c == 0) {
                $p_size = $feature_list[$#feature_list];
            } else {
                $p_size = $feature_list[$c - 1];
            }

            my $c_size = $feature_list[$c];
            next INTER if (!defined($p_size->{start}));
            next INTER if (!defined($p_size->{end}));
            next INTER if (!defined($c_size->{start}));
            next INTER if (!defined($c_size->{end}));

            my $p_start = $p_size->{start};
            my $p_end = $p_size->{end};
            my $c_start = $c_size->{start};
            my $c_end = $c_size->{end};
            my $interstart = 0;
            my $interend = 0;

            if ($p_start > $p_end) {
                $interstart = $p_start;
            } else {
                $interstart = $p_end;
            }
            if ($c_start < $c_end) {
                $interend = $c_start;
            } else {
                $interend = $c_end;
            }
            my $tmp_start = $interstart;
            my $tmp_end = $interend;
            if ($tmp_start > $tmp_end) {
                $interstart = $tmp_end;
                $interend = $tmp_start;
            }
            ## I need a little logic to catch the first feature
            ## As it is currently written, it grabs the position of the last feature
            ## and uses that number rather than the beginning of the chromosome.
            ## Therefore we end up with a first feature which contains the entire chromosome.
            if ($c <= 2) {   ## Arbitrarily set it to check the first 2 features
                if ($interend >= 100000) { ## If the first feature take >= 1/2 the chromosome
                    $interend = $interstart;
                    $interstart = 1;
                }
            }

            my $new_feature = $c_size->{feat};
            my $inter_location = Bio::Location::Atomic->new(-start => $interstart, -end => $interend, -strand => 1);
            $new_feature->{_location} = $inter_location;
            $inter_gffout->write_feature($new_feature);
        }
    }                           ## End while every sequence
    my $ret_stats = {num_sequences => $seq_count,
                     total_nt => $total_nt,
                     num_features => $feature_count,
                 };
    close($in);
    return($ret_stats);
}

=back

=head1 AUTHOR - atb

Email abelew@gmail.com

=head1 SEE ALSO

    L<samtools> L<Bio::FeatureIO>

=cut

1;
