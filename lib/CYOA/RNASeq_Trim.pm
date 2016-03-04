package CYOA;

=head1 NAME
    CYOA::RNASeq_Trim - Use trimomatic/cutadapt/etc to trim libraries

=head1 SYNOPSIS

    use CYOA;
    my $hpgl = new CYOA;
    $hpgl->Cutadapt();

=head2 Methods

=item C<Cutadapt>

    $hpgl->Cutadapt(); will use biopieces/cutadapt to attempt to
    remove sequence adapters from a library.  This is most common used
    by me for ribosome profiling libraries.

=cut
sub Cutadapt {
    my $me = shift;
    my %args = @_;
    $me->Check_Options(['input']);

    my $type = $me->{type};
    $type = $args{type} if (defined($args{type}));
    $type = 'tnseq' unless(defined($type));
    my $input = $me->{input};
    my $basename = basename($input, @{$me->{suffixes}});
    $basename = basename($basename, @{$me->{suffixes}});
    my $adapter_flags = "";
    my $minlength = 7;
    my $maxlength = 42;
    if ($type eq 'tnseq') {
        $adapter_flags = qq" -a ACAGGTTGGATGATAAGTCCCCGGTCTGACACATC -a ACAGTCCCCGGTCTGACACATCTCCCTAT -a ACAGTCCNCGGTCTGACACATCTCCCTAT ";
        $maxlength = 20;
    } elsif ($type eq 'riboseq') {
        $adapter_flags = qq" -a ACAGGTTGGATGATAAGTCCCCGGTCTGACACATCTCCCTAT -a AGATCGGAAGAGCACACGTCTGAAC -b AGATCGGAAGAGCACACGTCTGAAC ";
        $minlength = 16;
        if ($input =~ /\.csfasta/) {
            $adapter_flags = qq" -a CGCCTTGGCCGTACAGCAGCATATTGGATAAGAGAATGAGGAACCCGGGGCAG -a GCGGAACCGGCATGTCGTCGGGCATAACCCTCTCTTACTCCTTGGGCCCCGTC ";
        }
    } else {
        $adapter_flags = qq" -a ACAGGTTGGATGATAAGTCCCCGGTCTGACACATCTCCCTAT ";
    }

    my $comment = qq!## This script makes use of biopieces and cutadapt to trim away adapters
## and separate the sequence file into a few pieces depending on size
## and adapter status.  It also performs some simple graphs of the data.!;
    my $out_dir = qq"$me->{basedir}/outputs/cutadapt";
    my $type_flag = '';
    my $out_suffix = 'fastq';
    my $input_flags = qq"less ${input} | cutadapt - ";
    if ($input =~ /\.csfasta/) {
        $type_flag = '-c -t --strip-f3';
        $me->Check_Options(["qual"]);
        ## $input_flags = qq"less ${input} | sed 's/^T//g' | cutadapt - $me->{qual} ";
        $input_flags = qq"less ${input} | cutadapt - $me->{qual} "; ## If we are keeping quality files
        ## $input_flats = qq"less ${input} | cutadapt - ";
        $out_suffix = 'fastq';
    }
    my $output = qq"${basename}-trimmed.${out_suffix}";
    my $job_string = qq!
mkdir -p ${out_dir} && \\
 ${input_flags} ${type_flag} ${adapter_flags} -e 0.1 -n 3 -m ${minlength} -M ${maxlength} \\
  --too-short-output=${out_dir}/${basename}_tooshort.${out_suffix} \\
  --too-long-output=${out_dir}/${basename}_toolong.${out_suffix} \\
  --untrimmed-output=${out_dir}/${basename}_untrimmed.${out_suffix} \\
  2>outputs/cutadapt.err 1>${output}
!;
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
                                     input => qq"${out_dir}/${basename}_tooshort.fastq",
                                     comment => qq"## Compressing the tooshort sequences.",
	);
    my $comp_long = $me->Recompress(job_name => "xzcutlong",
                                    depends => $cutadapt->{pbs_id},
                                    input => qq"${out_dir}/${basename}_toolong.fastq",
                                    comment => qq"## Compressing the toolong sequences.",
	);
    my $comp_un = $me->Recompress(job_name => "xzuncut",
                                  depends => $cutadapt->{pbs_id},
                                  input => qq"${out_dir}/${basename}_untrimmed.fastq",
                                  comment => qq"## Compressing the toolong sequences.",
	);
    my $comp_original = $me->Recompress(job_name => "xzorig",
                                        depends => $cutadapt->{pbs_id},
                                        input => qq"$input",
                                        output => qq"sequences/${input}.xz",
                                        comment => qq"## Compressing the original sequence.",
	);
    return($cutadapt);
}

=item C<Trimomatic>

    $hpgl->Trimomatic(); calls the java tool trimomatic to remove
    adapters/low quality sequences.  If $args{input} has a ':' or ','
    then this will assume the input is comprised of two pairwise files
    and will call 'Trimomatic_Pairwise()', otherwise
    'Trimomatic_Single()'.

=cut
sub Trimomatic {
    my $me = shift;
    my %args = @_;
    my $trim;
    if ($me->{input} =~ /:|\,/) {
	$trim = $me->Trimomatic_Pairwise(%args);
    } else {
	$trim = $me->Trimomatic_Single(%args);
    }
    return($trim);
}

=item C<Trimomatic_Pairwise>

    $hpgl->Trimomatic_Pairwise(); invokes trimomatic with parameters
    suitable for pairwise sequence libraries.

=cut
sub Trimomatic_Pairwise {
    my $me = shift;
    my %args = @_;
    my $input = $me->{input};
    $input = $args{input} if ($args{input} and !$input);
    my @input_list = split(/:|\,/, $input);
    if (scalar(@input_list) <= 1) {
        my $ret = $me->Trimomatic_Single(input => $input);
        return($ret);
    }
    my $r1 = $input_list[0];
    my $r2 = $input_list[1];
    my @suff = (".fastq",".gz",".xz");

    my $basename = basename($r1, @suff);
    $basename = basename($basename, @suff);
    $basename =~ s/_R1$//g;
    $basename =~ s/_forward$//g;

    my $r1b = basename($r1, @suff);
    $r1b = basename($r1b, @suff);
    my $r2b = basename($r2, @suff);
    $r2b = basename($r2b, @suff);

    my $r1o = qq"${r1b}-trimmed.fastq.gz";
    my $r1op = qq"${r1b}-trimmed_paired.fastq.gz";
    my $r1ou = qq"${r1b}-trimmed_unpaired.fastq.gz";

    my $r2o = qq"${r2b}-trimmed.fastq.gz";
    my $r2op = qq"${r2b}-trimmed_paired.fastq.gz";
    my $r2ou = qq"${r2b}-trimmed_unpaired.fastq.gz";

    $me->{basename} = $basename;
    my $output = qq"${r1o}:${r2o}";
    my $comment = qq!## This call to trimomatic removes illumina and epicentre adapters from ${input}.
## It also performs a sliding window removal of anything with quality <25;
## cutadapt provides an alternative to this tool.
## The original sequence data is recompressed and saved in the sequences/ directory.!;
    my $job_string = qq!
## Trimomatic_Pairwise: In case a trimming needs to be redone...
if [[ \! -r "${r1}" ]]; then
  if [[ -r "sequences/${r1b}.fastq.xz" ]]; then
    mv sequences/${r1b}.fastq.xz . && pxz -d ${r1b}.fastq.xz && pigz ${r1b}.fastq && mv sequences/${r2b}.fastq.xz . && pxz -d ${r2b}.fastq.xz && pigz ${r2b}.fastq
  else
    echo "Missing files. Did not find ${r1} nor sequences/${r1b}.fastq.xz"
    exit 1
  fi
fi
trimomatic PE -threads 1 -phred33 ${r1} ${r2} ${r1op} ${r1ou} ${r2op} ${r2ou} ILLUMINACLIP:$me->{libdir}/adapters.fa:2:20:4 SLIDINGWINDOW:4:25 1>outputs/${basename}-trimomatic.out 2>&1
excepted=\$(grep "Exception" outputs/${basename}-trimomatic.out)
## The following is in case the illumina clipping fails, which it does if this has already been run I think.
if [[ "\${excepted}" \!= "" ]]; then
  trimomatic PE -threads 1 -phred33 ${r1} ${r2} ${r1op} ${r1ou} ${r2op} ${r2ou} SLIDINGWINDOW:4:25 1>>outputs/${basename}-trimomatic.out 2>&1
fi
sleep 10
mv ${r1op} ${r1o} && mv ${r2op} ${r2o}
!;
    ## Example output from trimomatic:
    ## Input Read Pairs: 10000 Both Surviving: 9061 (90.61%) Forward Only Surviving: 457 (4.57%) Reverse Only Surviving: 194 (1.94%) Dropped: 288 (2.88%)
    ## Perhaps I can pass this along to Get_Stats()
    my $trim = $me->Qsub(job_name => "trim",
			 qsub_wall => "12:00:00",
			 job_string => $job_string,
			 input => $input,
			 output => $output,
			 comment => $comment,
			 prescript => $args{prescript},
			 postscript => $args{postscript},
        );
    ## Set the input for future tools to the output from trimming.
    $me->{input} = $output;
    $trim_stats = $me->Trimomatic_Stats(basename => $basename, depends => $trim->{pbs_id}, pairwise => 1);
    $trim->{stats} = $trim_stats;
    return($trim);
}

=item C<Trimomatic_Single>

    $hpgl->Trimomatic_Single(); invokes trimomatic with parameters
    suitable for single-read sequence libraries.

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
## Trimomatic_Single: In case a trimming needs to be redone...
if [[ \! -r "${input}" ]]; then
  if [[ -r "sequences/${input}.xz" ]]; then
    mv sequences/${input}.xz . && pxz -d ${input}.xz
  else
    echo "Missing files. Did not find ${input} nor sequences/${input}.xz"
    exit 1
  fi
fi
trimomatic SE -phred33 ${input} ${output} ILLUMINACLIP:$me->{libdir}/adapters.fa:2:20:4 SLIDINGWINDOW:4:25 1>outputs/${basename}-trimomatic.out 2>&1
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
    $trim_stats = $me->Trimomatic_Stats(basename => $basename, depends => $trim->{pbs_id});
    ## Set the input for future tools to the output from this trimming operation.
    $me->{input} = $output;
    return($trim);
}

=item C<Trimomatic_Stats>

    Collect the trimming statistics from the output file
    'trimomatic.out' and report them in a .csv file by library.

=cut
sub Trimomatic_Stats {
    my $me = shift;
    my %args = @_;
    my $basename = $args{basename};
    my $input_file = "outputs/${basename}-trimomatic.out";
    my $depends = "";
    $depends = $args{depends} if ($args{depends});
    my $job_name = 'trimst';
    $job_name = $args{job_name} if ($args{job_name});
    my $comment = qq!## This is a stupidly simple job to collect trimomatic statistics!;
    my $stat_output = qq"outputs/trimomatic_stats.csv";
    my $job_string = qq!
if [ \! -r ${stat_output} ]; then
  echo "total_reads,surviving_reads,dropped_reads" > ${stat_output}
fi
total_reads_tmp=\$(grep "^Input Reads" ${input_file} | awk '{print \$3}')
total_reads=\${total_reads_tmp:-0}
surviving_reads_tmp=\$(grep "^Input Reads" ${input_file} | awk '{print \$5}')
surviving_reads=\${surviving_reads_tmp:-0}
dropped_reads_tmp=\$(grep "^Input Reads" ${input_file} | awk '{print \$8}')
dropped_reads=\${dropped_reads_tmp:-0}

stat_string=\$(printf "${basename},%s,%s,%s" "\${total_reads}" "\${surviving_reads}" "\${dropped_reads}")
echo "\$stat_string" >> ${stat_output}
!;
    if ($args{pairwise}) {
        ## The output looks a bit different for pairwise input:
        ## Input Read Pairs: 10000 Both Surviving: 9061 (90.61%) Forward Only Surviving: 457 (4.57%) Reverse Only Surviving: 194 (1.94%) Dropped: 288 (2.88%)
        $job_string = qq!
if [ \! -r ${stat_output} ]; then
  echo "total_reads,surviving_both,surviving_forward,surviving_reverse,dropped_reads" > ${stat_output}
fi
total_reads_tmp=\$(grep "^Input Read Pairs" ${input_file} | awk '{print \$4}')
total_reads=\${total_reads_tmp:-0}
surviving_both_tmp=\$(grep "^Input Read Pairs" ${input_file} | awk '{print \$7}')
surviving_both=\${surviving_both_tmp:-0}
surviving_forward_tmp=\$(grep "^Input Read Pairs" ${input_file} | awk '{print \$12}')
surviving_forward=\${surviving_forward_tmp:-0}
surviving_reverse_tmp=\$(grep "^Input Read Pairs" ${input_file} | awk '{print \$17}')
surviving_reverse=\${surviving_reverse_tmp:-0}
dropped_reads_tmp=\$(grep "^Input Read Pairs" ${input_file} | awk '{print \$20}')
dropped_reads=\${dropped_reads_tmp:-0}

stat_string=\$(printf "${basename},%s,%s,%s,%s,%s" "\${total_reads}" "\${surviving_both}" "\${surviving_forward}" "\${surviving_reverse}" "\${dropped_reads}")
echo "\$stat_string" >> ${stat_output}
!;
    }
    my $stats = $me->Qsub(job_name => $job_name,
                          depends => $depends,
                          qsub_queue => "throughput",
                          qsub_cpus => 1,
                          qsub_mem => 1,
                          qsub_wall => "00:10:00",
                          job_string => $job_string,
                          input => $input_file,
                          comment => $comment,
        );
    return($stats);
}

=back

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

    L<Bio::Tools::Run::StandAloneBlast>

=cut

1;

