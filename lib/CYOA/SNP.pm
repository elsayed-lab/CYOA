package CYOA;
use Modern::Perl;
use autodie;

sub SNP_Search {
    my $me = shift;
    my %args = @_;
    $me->Check_Options(args => \%args, needed => ["genome","query"]);
    my $genome = $me->{genome};
    my $query = $me->{query};
    $query = basename(${query}, (".bam"));
    my $new = qq"${query}_sorted";

    my $min_hits = 10;
    if ($args{min_hits}) {
        $min_hits = $args{min_hits};
    }

    my $job_string = qq!
samtools sort -l 9 -@ 4 ${query}.bam outputs/${new} 2>outputs/samtools_sort.out 1>&2
if [ \!-r "${genome}.fai" ]; then
    samtools faidx ${genome}
fi
mkdir -p outputs/vcfutils
samtools mpileup -u -f ${genome} outputs/${new}.bam 2>outputs/vcfutils/${new}_pileup.err |\\
  bcftools view -bvcg - 1>outputs/vcfutils/${new}_raw.bcf 2>outputs/vcfutils/${new}_bcf.err &&\\
  bcftools view outputs/vcfutils/${new}_raw.bcf |\\
  vcfutils_hacked.pl varFilter -d ${min_hits} 1>outputs/vcfutils/${new}_vcf_summary.txt 2>outputs/vcfutils/${new}_filter.err &&\\
rm outputs/${new}
!;
    my $comment_string = qq!## Use samtools, bcftools, and vcfutils to get some idea about how many variant positions are in the data.!;
    my $pileup_job = $me->Qsub(
        comment => $comment_string,
        job_name => "${query}_vcf",
        job_prefix => "80",
        job_string => $job_string,
        qsub_mem => 10,
        qsub_queue => "workstation",
        );

    my $print_input = qq"outputs/vcfutils/${new}_raw.bcf";
    my $print_output = qq"outputs/vcfutils/${new}_snp";
    $comment_string = qq!## This little job should make unique IDs for every detected
## SNP and a ratio of snp/total for all snp positions with > 20 reads.
## Further customization may follow.
!;
    my $parse_job = $me->SNP_Ratio(input => $print_input,
                                   output => $print_output,
                                   depends => $pileup_job->{pbs_id},);
    return($pileup_job, $parse_job);
}

sub SNP_Ratio {
    my $me = shift;
    my %args = @_;
    my $depends = $args{depends} if ($args{depends});
    $me->Check_Options(args => \%args, needed => ["input","output"]);
    my $print_input = $me->{input};
    $print_input = $args{input} if ($args{input});
    my $print_output = $me->{output};
    $print_output = $args{output} if ($args{output});
    my $job_name = qq"parsesnp";
    my $comment_string = "# Parse the SNP data.
";
    my $job_string = qq?
use CYOA;
my \$h = new CYOA;
\$h->Make_SNP_Ratio(input => '$print_input', output => '$print_output');
?;
    my $parse_job = $me->Qsub_Perl(
        comment => $comment_string,
        depends => $depends,
        job_name => "${job_name}_ratio",
        job_prefix => "81",
        job_string => $job_string,
        qsub_queue => "throughput",
        );
    my $un_comp = $me->Recompress(
        comment => qq"## Compressing the snp scorer output",
        depends => $parse_job->{pbs_id},
        input => "${print_output}",
        job_name => "xz_snpparse",
        job_prefix => "82",
        );
}

sub Make_SNP_Ratio {
    my $me = shift;
    my %args = @_;
    my $ratio_cutoff = 0.75;
    my $in = new FileHandle;
    $in->open("bcftools view $args{input} |");
    my $ratio_out = new FileHandle;
    $ratio_out->open(">$args{output}_ratio.txt");
    my $count_out = new FileHandle;
    $count_out->open(">$args{output}_count.txt");

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
      my $diff_sum = 0;
    INFO: foreach my $element (@info_list) {
        my ($element_type, $element_value) = split(/=/, $element);
        if ($element_type eq 'DP4') {
            my ($same_forward, $same_reverse, $alt_forward, $alt_reverse) = split(/,/, $element_value);
            $diff_sum = $alt_forward + $alt_reverse;
            $snp_sum = $same_forward + $same_reverse + $diff_sum;
            $snp_ratio = ($alt_forward + $alt_reverse) / $snp_sum;
            last INFO;
        } else {
            next INFO;
        }
    } ## End info list loop
      if ($snp_ratio >= $ratio_cutoff) {
          print $count_out qq"${snp_id}\t${snp_sum}\n";
          print $ratio_out qq"${snp_id}\t${snp_ratio}\n";
      }
  }
    $in->close();
    $count_out->close();
    $ratio_out->close();
}

1;
