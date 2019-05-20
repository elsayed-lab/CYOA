package Bio::Adventure::SNP;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use File::Basename;
use Math::Round qw":all";

sub Align_SNP_Search {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        vcf_cutoff => 4,
        vcf_minpct => 0.8,
    );

    my $genome = qq"$options->{libdir}/$options->{libtype}/$options->{species}.fasta";
    my $query = $options->{input};
    my $query_home = dirname(${query});
    my $query_base = basename(${query}, (".bam"));
    $query = qq"${query_home}/${query_base}";
    print "About to start Bowtie2 search against of $query against $options->{species}.\n";
    my $bt2_job = Bio::Adventure::RNASeq_Map::Bowtie2(
        $class,
        htseq_type => "exon",
        input => $query,
        species => $options->{species},
    );
    my $bamfile = $bt2_job->{samtools}->{job_output};
    print "About to start SNP search of $bamfile against $options->{species}\n";
    my $search = Bio::Adventure::SNP::SNP_Search(
        $class,
        input => $bamfile,
        depends => $bt2_job->{samtools}->{job_id},
        species => $options->{species},
    );
    return($search);
}

sub SNP_Search {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        varfilter => 0,
        vcf_cutoff => 10,
        vcf_minpct => 0.8,
    );
    my $genome = qq"$options->{libdir}/$options->{libtype}/$options->{species}.fasta";
    my $query = $options->{input};
    my $query_home = dirname(${query});
    my $query_base = basename(${query}, (".bam"));
    $query = qq"${query_home}/${query_base}";

    my $varfilter = $options->{varfilter};
    my $vcf_cutoff = $options->{vcf_cutoff};
    my $vcf_minpct = $options->{vcf_minpct};

    my $vcfutils_dir = qq"outputs/vcfutils_$options->{species}";
    my $pileup_input = qq"${vcfutils_dir}/${query_base}.bam";
    my $pileup_error = qq"${vcfutils_dir}/${query_base}_pileup.err";
    my $pileup_output = qq"${vcfutils_dir}/${query_base}_pileup.vcf";
    my $filter_error = qq"${vcfutils_dir}/${query_base}_varfilter.err";
    my $filter_output = qq"${vcfutils_dir}/${query_base}_filtered.vcf";
    my $call_output = qq"${vcfutils_dir}/${query_base}_summary.vcf";
    my $call_error = qq"${vcfutils_dir}/${query_base}_summary.err";
    my $final_output = qq"${vcfutils_dir}/${query_base}.bcf";
    my $final_error = qq"${vcfutils_dir}/${query_base}_bcf.err";
    my $jstring = qq!mkdir -p ${vcfutils_dir}
echo "Started samtools sort at \$(date)" >> outputs/vcfutils_$options->{species}.out
!;
    unless (-r "${pileup_input}") {
        if ($query =~ m/sorted/) {
            $jstring .= qq!if test \! -e \$(pwd)/${pileup_input}; then
  ln -sf \$(pwd)/${query}.bam \$(pwd)/${pileup_input}
fi
!;
        } else {
            $jstring .= qq!samtools sort -l 9 -@ 4 ${query}.bam -o ${pileup_input} 2>outputs/samtools_sort.out 1>&2
if [ "\$?" -ne "0" ]; then
    echo "samtools sort failed."
exit 1
fi
!;
        }  ## End checking if the pileup input does not have sorted.
    } ## Found the input for samtools mpileup
    $jstring .= qq!
if [ \! -r "${genome}.fai" ]; then
    samtools faidx ${genome}
fi
samtools mpileup -uvf ${genome} 2>samtools_mpileup.err \\
    ${pileup_input} |\\
  bcftools call -c - 2>bcftools_call.err |\\
  bcftools view -l 9 -o ${final_output} -O b - \\
    2>${call_error}
if [ "\$?" -ne "0" ]; then
    echo "mpileup/bcftools failed."
    exit 1
fi
bcftools index ${final_output} 2>bcftools_index.err
echo "Successfully finished." >> outputs/vcfutils_$options->{species}.out
!;

    my $comment_string = qq!## Use samtools, bcftools, and vcfutils to get some idea about how many variant positions are in the data.!;
    my $pileup = $class->Submit(
        comment => $comment_string,
        depends => $options->{depends},
        jname => "bcf_${query_base}_$options->{species}",
        jprefix => "80",
        job_type => "snpsearch",
        jstring => $jstring,
    );
    $comment_string = qq!## This little job should make unique IDs for every detected
## SNP and a ratio of snp/total for all snp positions with > 20 reads.
## Further customization may follow.
!;
    print "TESTME: $pileup->{job_id}\n";
    my $parse = Bio::Adventure::SNP::SNP_Ratio(
        $class,
        input => ${final_output},
        vcf_cutoff => $vcf_cutoff,
        vcf_minpct => $vcf_minpct,
        species => $options->{species},
        depends => $pileup->{job_id},
    );
    return([$pileup, $parse]);
}

sub SNP_Ratio {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        vcf_cutoff => 4,
        vcf_minpct => 0.8,
    );
    my $print_input = $options->{input};
    my $print_output = $print_input;
    $print_output =~ s/\.bcf//g;
    $print_output = qq"${print_output}_$options->{species}";
    my $genome = qq"$options->{libdir}/$options->{libtype}/$options->{species}.fasta";

    my $jname = qq"parsesnp";
    my $comment_string = qq!
## Parse the SNP data and generate a modified $options->{species} genome.
##  This should read the file:
## ${print_input}
##  and provide 4 new files:
## ${print_output}_count.txt
## ${print_output}_all.txt
## ${print_output}_pct.txt
## and a modified genome: ${print_output}_modified.fasta
!;

    my $jstring = qq"
use Bio::Adventure::SNP;
Bio::Adventure::SNP::Make_SNP_Ratio(
  \$h,
  input => '$print_input',
  output => '$print_output',
  species => '$options->{species}',
  vcf_cutoff => '$options->{vcf_cutoff}',
  vcf_minpct => '$options->{vcf_minpct}',
);
";

    my $parse_job = $class->Submit(
        comment => $comment_string,
        depends => $options->{depends},
        jname => "${jname}_$options->{species}",
        jprefix => "81",
        jstring => $jstring,
        language => 'perl',
        queue => 'throughput',
    );
    return($parse_job);
}

sub Make_SNP_Ratio {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        vcf_cutoff => 4,
        vcf_minpct => 0.8,
    );
    my $input = $options->{input};
    my $output_base = $options->{output};
    my $species = $options->{species};
    my $genome = "$options->{libdir}/$options->{libtype}/${species}.fasta";
    my $gff = "$options->{libdir}/$options->{libtype}/${species}.gff";
    my $in = FileHandle->new("bcftools view ${input} |");
    my $pct_out = FileHandle->new(">${output_base}_pct.txt");
    print $pct_out "# chr_pos_ref_alt\tpct\n";
    my $count_out = FileHandle->new(">${output_base}_count.txt");
    print $count_out "# chr_pos_ref_alt\tdiff_count\n";
    my $all_out = FileHandle->new(">${output_base}_all.txt");
    print $all_out "# chr_pos_ref_alt\tdiff_count\tall_count\tpct\n";

    my $count = 0;     ## Use this to count the positions changed in the genome.

    ## Ok, so I want to simplify the vcf output so that I can create a pseudo
    ## count table of every potential SNP position I want to create a 2 column
    ## file where the first column is a unique ID containing 'chr_pos_ref_new'
    ## and the second column is some sort of score which I will be able to use
    ## for a PCA when I merge a bunch of these together.  Candidates include:
    ## quality score, DP*AF1 (depth * max likelihood estimate) Maybe also consider
    ## DP4 column which is comma separated: #forward-ref, #reverse-ref, #forward-alt, #reverse-alt
    ## If that, then sum all 4 and only take those greater than x (20?), then take
    ## forward-alt+reverse-alt/(sum of all) to get simple SNP ratio

    my $num_variants = 0; ## Use this to count the variant positions of reasonably high confidence.
  READER: while (my $line = <$in>) {
        next READER if ($line =~ /^#/);
        ## When using samtools mpileup with the -t LIST, we get some extra
        ## fields at the end which are called tags and tag_info here.
        my ($chr, $pos, $id, $ref, $alt, $qual, $filt, $info, $tags, $tag_info) = split(/\s+/, $line);
        ## Only accept the simplest non-indel mutations.
        next READER unless ($alt eq 'A' or $alt eq 'a' or $alt eq 'T' or $alt eq 't' or
                                $alt eq 'G' or $alt eq 'g' or $alt eq 'C' or $alt eq 'c');

        ## I do not know if indels are here anymore, we will see.
        my @info_list = split(/;/, $info);
        my $snp_id = qq"chr_${chr}_pos_${pos}_ref_${ref}_alt_${alt}";
        my $snp_pct = 0.0;
        my $all_sum = 0;
        my $diff_sum = 0;
        foreach my $element (@info_list) {
            my ($element_type, $element_value) = split(/=/, $element);
            if ($element_type eq 'DP4') {
                my ($same_forward, $same_reverse, $alt_forward, $alt_reverse) = split(/,/, $element_value);
                $diff_sum = $alt_forward + $alt_reverse;
                $all_sum = $same_forward + $same_reverse + $diff_sum;
                $snp_pct = ($alt_forward + $alt_reverse) / $all_sum;
                $snp_pct = nearest(0.01, $snp_pct);
            }
        }
        if ($snp_pct >= $options->{vcf_minpct} && $all_sum >= $options->{vcf_cutoff}) {
            $num_variants = $num_variants + 1;
            print $all_out qq"${snp_id}\t${diff_sum}\t${all_sum}\t${snp_pct}\n";
            print $pct_out qq"${snp_id}\t${snp_pct}\n";
            print $count_out qq"${snp_id}\t${diff_sum}\n";
        }
    }
    $in->close();
    $count_out->close();
    $all_out->close();
    $pct_out->close();

    if ($num_variants > 0) {
        my $input_genome = Bio::Adventure::Read_Genome_Fasta($class, %args, fasta => $genome,);
        my $input_data = FileHandle->new("<${output_base}_pct.txt");
        my $annotations = Bio::Adventure::Read_Genome_GFF($class, gff => $gff, feature_type => 'gene', %args);
        my $output_by_gene = FileHandle->new(">${output_base}_hits_by_gene.txt");
        my $vars_by_gene = {};
      READER: while (my $line = <$input_data>) {
            chomp $line;
            next READER if ($line =~ /^\#/);
            $count = $count++;
            next READER if ($line =~ /^#/);
            my ($position_data, $count) = split(/\t/, $line);
            ## LmxM.28_1042130_G_A
            my ($chr, $pos, $orig, $alt, $rest);
            ($chr, $rest) = split(/_pos_/, $position_data);
            $chr =~ s/^chr_//g;
            ($pos, $rest) = split(/_ref_/, $rest);
            ($orig, $alt) = split(/_alt_/, $rest);
            my $tmp_seq = $input_genome->{$chr}->{sequence};
            my $initial_search_seq = $tmp_seq;
            my $found_nt = substr($initial_search_seq, $pos - 2, 5);
            my $replace_seq = $tmp_seq;
            my $replaced = substr($replace_seq, $pos - 1, 1, $alt);
            ##my $final_search_seq = $tmp_seq;
            ##my $new_nt =  substr($final_search_seq, $pos - 2, 5);
            ##my $length_test = length($replace_seq);
            ##print "Original region at pos: $pos are: $found_nt, changing $orig to $alt.  Now they are $new_nt length: $length_test\n";
            $input_genome->{$chr}->{sequence} = $replace_seq;

            if (!defined($annotations->{$chr})) {
                print "$chr was not defined in the annotations database.\n";
            } else {
                my %tmp = %{$annotations->{$chr}};
                foreach my $gene_id (keys %tmp) {
                    my $gene = $tmp{$gene_id};
                    if (ref($gene) eq 'HASH') {
                        $vars_by_gene->{$gene->{id}}->{length} = $gene->{end} - $gene->{start};
                        if ($gene->{start} <= $pos &&
                                $gene->{end} >= $pos) {
                            if (!defined($vars_by_gene->{$gene->{id}}->{count})) {
                                $vars_by_gene->{$gene->{id}}->{count} = 1;
                            } else {
                                $vars_by_gene->{$gene->{id}}->{count} = $vars_by_gene->{$gene->{id}}->{count} + 1;
                            }
                            print $output_by_gene "$gene->{id}\t$chr\t$pos\t${orig}_${alt}\n";
                        }
                    }
                }
            }
        }
        $input_data->close();
        $output_by_gene->close();
        my $var_by_gene = FileHandle->new(">${output_base}_var_by_genelen.txt");
        foreach my $geneid (keys %{$vars_by_gene}) {
            my $gene_ratio = 0;
            my $good = 1;
            if (! $vars_by_gene->{$geneid}->{count}) {
                $good = 0;
            }
            if (! $vars_by_gene->{$geneid}->{length}) {
                $good = 0;
            }
            if (! $vars_by_gene->{$geneid}->{length}) {
                $good = 0;
            }
            if ($good) {
                $gene_ratio = ($vars_by_gene->{$geneid}->{count} / $vars_by_gene->{$geneid}->{length}) * 1000.0;
            }
            print $var_by_gene "$geneid\t$gene_ratio\n";
        }
        $var_by_gene->close();

        my $output_genome = FileHandle->new(">${output_base}_modified.fasta");
        foreach my $ch (sort keys %{$input_genome}) {
            ## my $formatted = $text->format($input_genome->{$ch}->{sequence});
            print $output_genome ">${ch}\n";
            ## Take from: https://www.biostars.org/p/70448/
            foreach my $seq_line (unpack('(a[80])*', $input_genome->{$ch}->{sequence})) {
                print $output_genome "$seq_line\n";
            }
            ##${formatted}
            ##";
        }
        $output_genome->close();
    } else {
        print "No differences were detected by bcftools call.\n";
    }
    return($count);
}

1;
