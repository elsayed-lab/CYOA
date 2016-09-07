package CYOA;
use Modern::Perl;
use autodie;

sub SNP_Search {
    my $me = shift;
    my %args = @_;
    ## The following comes from my make_snp_genomes.sh
    ##export START=$(pwd)
    ##export REF="/cbcb/nelsayed-scratch/atb/libraries/genome/lpanamensis.fasta"
    ##for hpgl in $(find . -name 'hpgl*'); do
    ##export HPGL=${hpgl}
    ##export CMD="samtools mpileup -uf ${REF} ${START}/${HPGL}/outputs/tophat_lpanamensis/accepted_paired_sorted.bam | bcftools view - cg - | vcfutils.pl vcf2fq | read_fastq -i - | write_fasta -x -o ${HPGL}_snp_genome.fasta"
    ##echo "${CMD}"
    ##cat <<"EOF" | qsub -V -d $(pwd) -q long -
    ##samtools mpileup -uf ${REF} ${START}/${HPGL}/outputs/tophat_lpanamensis/accepted_paired_sorted.bam | bcftools view -cg - | vcfutils_hacked.pl vcf2fq > ${HPGL}_snp_genome.fasta
    ##EOF
    ##done

    ## While the following comes from:
    ## http://samtools.sourceforge.net/mpileup.shtml
    ##samtools mpileup -uf ref.fa aln1.bam aln2.bam | bcftools view -bvcg - > var.raw.bcf
    ##bcftools view var.raw.bcf | vcfutils.pl varFilter -D100 > var.flt.vcf

    $me->Check_Options(args => \%args, needed => ["genome","query","name"]);
    my $genome = $me->{genome};
    my $query = $me->{query};
    my $name = $me->{name};
    my $new = qq"${name}_";
    $new .= basename($query, (".bam"));

    my $job_string = qq!
samtools sort -l 9 -@ 4 ${query} outputs/${new} 2>outputs/samtools_sort.out 1>&2
if [ \!-r "${genome}.fai" ]; then
    samtools faidx ${genome}
fi
samtools mpileup -u -f ${genome} outputs/${new}.bam 2>outputs/${new}_pileup.err |\
  bcftools view -bvcg - 1>outputs/${name}_raw.bcf 2>outputs/${new}_bcf.err &&\
bcftools view outputs/${name}_raw.bcf | vcfutils_hacked.pl varFilter -D100 1>outputs/${name}_filtered.vcf 2>outputs/${new}_filter.err &&\
rm outputs/${new}
!;

    my $pileup_job = $me->Qsub(job_name => "${name}_vcf",
                               job_string => $job_string,
                               qsub_mem => 10,
                               qsub_queue => "workstation",
        );

    my $print_input = qq"outputs/${name}_raw.bcf";
    my $print_output = qq"outputs/${name}_snp_ratio.txt";
    my $comment_string = qq!## This little job should make unique IDs for every detected
## SNP and a ratio of snp/total for all snp positions with > 20 reads.
## Further customization may follow.
!;
    my $parse_job = $me->SNP_Ratio(input => $print_input, output => $print_output, depends => $pileup_job->{pbs_id});
    return($pileup_job, $parse_job);
}

sub SNP_Ratio {
    my $me = shift;
    my %args = @_;
    my $depends = $args{depends} if ($args{depends});
    $me->Check_Options(args => \%args, needed => ["input","output"]);
    my $print_input = $args{input} unless ($me->{input});
    my $print_output = $args{output} unless ($me->{output});
    my $job_name = qq"parsesnp";
    my $comment_string = "# Parse the SNP data.
";
    my $job_string = qq?
use CYOA;
my \$h = new CYOA;
\$h->SNP_Ratio(input => '$print_input', output => '$print_output');
?;
    my $parse_job = $me->Qsub_Perl(job_name => "${job_name}_ratio",
                              depends => $depends,
                              job_string => $job_string,
                              comment => $comment_string,
        );
    my $un_comp = $me->Recompress(depends => $parse_job->{pbs_id},
                                  job_name => "xz_snpparse",
                                  comment => qq"## Compressing the snp scorer output",
                                  input => "${print_output}");
}

sub Make_SNP_Ratio {
    my $me = shift;
    my %args = @_;
    my $minimum = 20;
    my $ratio_cutoff = 0.3;
    my $in = new FileHandle;
    $in->open("bcftools view $args{input} |");
    my $out = new FileHandle;
    $out->open(">$args{output}");

## Ok, so I want to simplify the vcf output so that I can create a pseudo count table of every potential SNP position
## I want to create a 2 column file where the first column is a unique ID containing 'chr_pos_ref_new' and the second
## column is some sort of score which I will be able to use for a PCA when I merge a bunch of these together.
## Candidates include: quality score, DP*AF1 (depth * max likelihood estimate)
## Maybe also consider DP4 column which is comma separated: #forward-ref, #reverse-ref, #forward-alt, #reverse-alt
##   If that, then sum all 4 and only take those greater than 20, then take forward-alt+reverse-alt/(sum of all) to get simple ratio SNP

  OUT: while (my $line = <$in>) {
      next OUT if ($line =~ /^#/);
      #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  hpgl0242_accepted_paired.bam
      my ($chr, $pos, $id, $ref, $alt, $qual, $filt, $info) = split(/\s+/, $line);
      next OUT if (length($alt) > 1);  ## I don't feel like dealing with indels
     my @info_list = split(/;/, $info);
      my $snp_id = qq"${chr}_${pos}_${ref}_${alt}";
      my $snp_ratio = 0.0;
      my $snp_sum = 0;
    INFO: foreach my $element (@info_list) {
        my ($element_type, $element_value) = split(/=/, $element);
        if ($element_type eq 'DP4') {
            my ($same_forward, $same_reverse, $alt_forward, $alt_reverse) = split(/,/, $element_value);
            $snp_sum = $same_forward + $same_reverse + $alt_forward + $alt_reverse;
            $snp_ratio = ($alt_forward + $alt_reverse) / $snp_sum;
            last INFO;
        } else {
            next INFO;
        }
    } ## End info list loop
      if ($snp_sum >= $minimum && $snp_ratio >= $ratio_cutoff) {
          print $out qq"${snp_id}\t${snp_ratio}\n";
      }
  }
    $in->close();
    $out->close();
}

1;
