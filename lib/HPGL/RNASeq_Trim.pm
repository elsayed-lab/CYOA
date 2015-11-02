package HPGL;

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
xzcat -f ${input} | cutadapt - ${cutadapt_flags} -e 0.1 -n 3 -m ${minlen} -M ${maxlen} \\
    --too-short-output=cutadapt/${basename}_tooshort.fastq \\
    --too-long-output=cutadapt/${basename}_toolong.fastq \\
    --untrimmed-output=cutadapt/${basename}_untrimmed.fastq \\
    2>outputs/cutadapt.err 1> ${output}!;
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
    my $comment = qq!## This call to trimomatic removes illumina and epicentre adapters from ${input}.
## It also performs a sliding window removal of anything with quality <25;
## cutadapt provides an alternative to this tool.
## The original sequence data is recompressed and saved in the sequences/ directory.!;
    my $job_string = qq!
## In case a trimming needs to be redone...
if [[ \! -r "${r1}" ]]; then
  if [[ -r "sequences/${r1}.xz" ]]; then
    mv sequences/${r1}.xz . && pxz -d ${r1}.xz && mv sequences/${r2}.xz . && pxz -d ${r2}.xz
  else
    echo "Missing files. Did not find ${r1} nor sequences/${r1}.xz"
    exit 1
  fi
fi
trimomatic PE -threads 1 -phred33 ${r1} ${r2} ${r1op} ${r1ou} ${r2op} ${r2ou} ILLUMINACLIP:$me->{libdir}/adapters.fa:2:20:4 SLIDINGWINDOW:4:25 2>outputs/${basename}-trimomatic.err 1>outputs/${basename}-trimomatic.out
excepted=\$(grep "Exception" outputs/${basename}-trimomatic.err)
## The following is in case the illumina clipping fails, which it does if this has already been run I think.
if [[ "\${excepted}" \!= "" ]]; then
  trimomatic PE -threads 1 -phred33 ${r1} ${r2} ${r1op} ${r1ou} ${r2op} ${r2ou} SLIDINGWINDOW:4:25 1>>outputs/${basename}-trimomatic.out 2>>output/${basename}-trimomatic.err
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
## In case a trimming needs to be redone...
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
    my $comp = $me->Recompress(job_name => "",
			       depends => $trim->{pbs_id},
			       input => $input,
			       output => qq"sequences/",
			       comment => qq"## The original sequence file is in sequences/",
        );
    $trim_stats = $me->Trimomatic_Stats(basename => $basename, depends => $trim->{pbs_id});
    ## Set the input for future tools to the output from this trimming operation.
    $me->{input} = $output;
    return($trim);
}

sub Trimomatic_Stats {
    my $me = shift;
    my %args = @_;
    my $basename = $args{basename};
    my $input_file = "outputs/${basename}-trimomatic.out";
    my $depends = "";
    $depends = $args{depends} if ($args{depends});
    my $job_name = 'trimst';
    $job_name = $args{job_name} if ($args{job_name});
    my $comment = qq!
## This is a stupidly simple job to collect alignment statistics.
!;
    my $job_string = qq!
total_reads_tmp=\$(grep "^Input Reads" $input_file | awk '{print \$3}')
total_reads=\${total_reads_tmp:-0}
surviving_reads_tmp=\$(grep "^Input Reads" $input_file | awk '{print \$5}')
surviving_reads=\${surviving_reads_tmp:-0}
dropped_reads_tmp=\$(grep "^Input Reads" $input_file | awk '{print \$8}')
dropped_reads=\${dropped_reads_tmp:-0}

stat_string=\$(printf "${basename},%s,%s,%s" "\${total_reads}" "\${surviving_reads}" "\${dropped_reads}")
echo "\$stat_string" >> outputs/trimomatic_stats.csv
!;
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

1;
