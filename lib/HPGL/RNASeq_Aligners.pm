package HPGL;

=head2
    Bowtie()
=cut
sub Bowtie {
    my $me = shift;
    my %args = @_;
    my %bt_jobs = ();
    my $bt_input = $me->{input};
    my $bt_depends_on;
    $bt_depends_on = $args{depends} if ($args{depends});
    $me->Check_Options(["species"]);
    my $basename = $me->{basename};

    my $bt_type = "v0M1";
    $bt_type = $args{bt_type} if ($args{bt_type});
    my $bt_args = $me->{bt_args}->{$bt_type};
    my $jobname = qq"bt${bt_type}";
    $jobname = $args{jobname} if ($args{jobname});
    my $libtype = $me->{libtype};
    $libtype = $args{libtype} if ($args{libtype});
    my $count = 1;
    $count = $args{count} if (defined($args{count}));

    ## Check that the indexes exist
    my $bt_reflib = "$me->{libdir}/${libtype}/$me->{species}";
    my $bt_reftest = qq"${bt_reflib}.1.ebwt";
    if (!-r $bt_reftest) {
        my $index_job = $me->BT1_Index(depends => $bt_depends_on, libtype => $libtype);
        $bt_jobs{index} = $index_job;
        $bt_depends_on = $index_job->{pbs_id};
    }
    my $bowtie_input_flag = "-q";  ## fastq by default
    $bowtie_input_flag = "-f" if ($me->{input} =~ /\.fasta$/);
    my $basename = $me->{basename};
    my $species = $me->{species};
    my $error_file = qq"bowtie_out/${basename}-${bt_type}.err";
    my $trim_output_file = qq"outputs/${basename}-trimomatic.out";
    my $comment = qq!## This is a bowtie1 alignment of $bt_input against
## $bt_reflib using arguments: $bt_args.
## This jobs depended on: $bt_depends_on
!;
    my $aligned_filename = qq"bowtie_out/${basename}-${bt_type}_aligned_${species}.fastq";
    my $unaligned_filename = qq"bowtie_out/${basename}-${bt_type}_unaligned_${species}.fastq";
    my $sam_filename = qq"bowtie_out/${basename}-${bt_type}.sam";
    my $job_string = qq!mkdir -p bowtie_out && sleep 10 && bowtie $bt_reflib $bt_args \\
  -p 4 \\
  $bowtie_input_flag $bt_input \\
  --un ${unaligned_filename} \\
  --al ${aligned_filename} \\
  -S ${sam_filename} \\
  2>${error_file} \\
  1>bowtie_out/${basename}-${bt_type}.out
!;
    my $bt_job = $me->Qsub(job_name => qq"bt1_${bt_type}",
                           depends => $bt_depends_on,
                           job_string => $job_string,
                           input => $bt_input,
                           comment => $comment,
                           output => $sam_filename,
                           unaligned => $unaligned_filename,
                           aligned => $aligned_filename,
			   prescript => $args{prescript},
			   postscript => $args{postscript},
        );
    $bt_jobs{bowtie} = $bt_job;

    my $un_comp = $me->Recompress(depends => $bt_job->{pbs_id},
                                  job_name => "xzun",
                                  comment => qq"## Compressing the sequences which failed to align against $bt_reflib using options $bt_args\n",
                                  input => "bowtie_out/${basename}-${bt_type}_unaligned_${species}.fasta");
    $bt_jobs{unaligned_compression} = $un_comp;

    my $al_comp = $me->Recompress(input => "bowtie_out/${basename}-${bt_type}_aligned_${species}.fasta",
                                  comment => qq"## Compressing the sequences which successfully aligned against $bt_reflib using options $bt_args",
                                  job_name => "xzal",
                                  depends => $bt_job->{pbs_id},);
    $bt_jobs{aligned_compression} = $al_comp;

    my $stats = $me->BT1_Stats(depends => $bt_job->{pbs_id},
                               job_name => "bt1stats",
                               bt_type => $bt_type,
                               count_table => qq"${basename}-${bt_type}.count.xz",
                               trim_input => ${trim_output_file},
                               bt_input => $error_file);

    my $sam_job = $me->Samtools(depends => $bt_job->{pbs_id},
                                job_name => "s2b",
                                input => $bt_job->{output},);

    $bt_jobs{samtools} = $sam_job;
    $me->{output} = $sam_job->{output};
    my $htmulti;
    if ($count) {
        if ($libtype eq 'rRNA') {
            $htmulti = $me->HTSeq(depends => $sam_job->{pbs_id},
                                  suffix => $bt_type,
                                  libtype => $libtype,
                                  input => $sam_job->{output},);
        } else {
            $htmulti = $me->HT_Multi(depends => $sam_job->{pbs_id},
                                     suffix => $bt_type,
                                     libtype => $libtype,
                                     input => $sam_job->{output},);
            $bt_jobs{htseq} = $htmulti;
        }
    }
    return(\%bt_jobs);
}

=head2
    BT_Multi()
=cut
sub BT_Multi {
    my $me = shift;
    my %args = @_;
    $me->Check_Options(["species"]);
    my $basename = $me->{basename};
    my $species = $me->{species};
    my $depends_on = $args{depends};
    my %bt_types = %{$me->{bt_args}};
    foreach my $type (keys %bt_types) {
        print "Starting $type\n";
        my $job = $me->Bowtie(bt_type => $type,
                              depends => $depends_on,
                              jobname => qq"bt${type}",
			      prescript => $args{prescript},
			      postscript => $args{postscript},
            );
    }
}

=head2
    Bowtie_RRNA()
    Perform an alignment against a home-curated set of ribosomal
    RNA/tRNA sequences.  The alignment requires a fastq input library
    and fasta library found in 'libraries/rRNA/$me->{species}.fasta'

  Example:
    my $rrna = $hpgl->Bowtie_RRNA();
    ## If you want to exclude the rRNA sequences from future alignments:
    my $rrna = $hpgl->Bowtie_RRNA(exclude => 1);

=cut
sub Bowtie_RRNA {
    my $me = shift;
    my %args = @_;
    $me->Check_Options(["species"]);
    my $basename = $me->{basename};
    $me->{basename} = qq"rRNA_${basename}";
    my $exclude = 0;
    $exclude = $args{exclude} if ($args{exclude});
    my $species = $me->{species};
    my $depends_on = $args{depends};
    if ($exclude) {
        my $in = $me->{input};
        $args{postscript} = qq"mv $in includerrna_${in} && mv bowtie_out/${basename}-rRNA_unaligned_${species}.fastq $in";
    }
    my $job = $me->Bowtie(depends => $depends_on,
                          job_name => qq"btrrna",
                          libtype => 'rRNA',
                          prescript => $args{prescript},
                          postscript => $args{postscript},
        );
    ## Return the basename back to normal so that future tasks don't
    ## get confuseled.
    $me->{basename} = $basename;
    return($job);
}

=head2
    BT1_Index()
=cut
sub BT1_Index {
    my $me = shift;
    my %args = @_;
    my $basename = $me->{basename};
    my $dep = "";
    $dep = $args{depends};
    my $libtype = $me->{libtype};
    $libtype = $args{libtype} if ($args{libtype});
    my $job_string = qq!bowtie-build $me->{libdir}/${libtype}/$me->{species}.fasta $me->{libdir}/genome/$me->{species}
!;
    my $comment = qq!## Generating bowtie1 indexes for species: $me->{species} in $me->{libdir}/${libtype}!;
    my $bt1_index = $me->Qsub(job_name => "bt1idx",
                              depends => $dep,
                              job_string => $job_string,
                              comment => $comment,
			      prescript => $args{prescript},
			      postscript => $args{postscript},
        );
    return($bt1_index);
}

=head2
    BT2_Index()
=cut
sub BT2_Index {
    my $me = shift;
    my %args = @_;
    my $basename = $me->{basename};
    my $dep = "";
    $dep = $args{depends};
    my $job_string = qq!
if [ \! -r "$me->{libdir}/genome/$me->{species}.fa" ]; then
  ln -s $me->{libdir}/genome/$me->{species}.fasta $me->{libdir}/genome/$me->{species}.fa
fi
bowtie2-build $me->{libdir}/genome/$me->{species}.fasta $me->{libdir}/genome/$me->{species}
!;
    my $comment = qq!## Generating bowtie2 indexes for species: $me->{species} in $me->{libdir}/genome!;
    my $jobid = $me->Qsub(job_name => "bt2idx",
                          depends => $dep,
                          job_string => $job_string,
                          comment => $comment,
			  prescript => $args{prescript},
			  postscript => $args{postscript},
        );
    return($jobid);
}

=head2
    BWA()
=cut
sub BWA {
    my $me = shift;
    my %args = @_;
    my %bwa_jobs = ();
    my $bwa_input = $me->{input};
    my $bwa_depends_on;
    $me->Check_Options(["species"]);
    my $basename = $me->{basename};

    ## Check that the indexes exist
    my $bwa_reflib = "$me->{libdir}/genome/$me->{species}";
    my $bwa_reftest = qq"${bwa_reflib}.bwa";
    if (!-r $bwa_reftest) {
        my $index_job = $me->BWA_Index(depends => $bwa_depends_on);
        $bwa_jobs{index} = $index_job;
        $bwa_depends_on = $index_job->{pbs_id};
    }
    my $basename = $me->{basename};
    my $species = $me->{species};
    my $comment = qq!## This is a BWA alignment of $bwa_input against
## $bwa_reflib.!;
    my $job_string = qq!mkdir -p bwa_out
cd bwa_out && bwa aln $me->{libdir}/genome/$me->{species}.fasta $bwa_input 2>bwa.err 1>${basename}.sai
!;
    my $bwa_job = $me->Qsub(job_name => "bwa",
			    depends => $bwa_depends_on,
			    job_string => $job_string,
			    input => $bwa_input,
			    comment => $comment,
			    output => qq"bwa/${basename}.sai",
			    prescript => $args{prescript},
			    postscript => $args{postscript},
        );
    $bwa_jobs{bowtie} = $bwa_job;

    my $sam_job = $me->Samtools(depends => $bwa_job->{pbs_id},
                                job_name => "s2b",
                                input => $bwa_job->{output},);
    $bwa_jobs{samtools} = $sam_job;

    my $htmulti = $me->HT_Multi(depends => $sam_job->{pbs_id},
                                input => $sam_job->{output},);
    $bwa_jobs{htseq} = $htmulti;

    my $c = 0;
    foreach my $ht_out (@{$htmulti}) {
        $c++;
        if ($ht_out) {
            my $count_compression = $me->Recompress(depends => $ht_out->{pbs_id},
                                                    job_name => "xz_hts${c}",
                                                    input => $ht_out->{output});
        }
    }
    return(\%bwa_jobs);
}

=head2
    BWA_Index()
=cut
sub BWA_Index {
    my $me = shift;
    my %args = @_;
    my $basename = $me->{basename};
    my $dep = "";
    $dep = $args{depends};
    my $job_string = qq!cd $me->{libdir}/genome/ && bwa index $me->{libdir}/genome/$me->{species}.fasta
!;
    my $comment = qq!## Generating bwa indexes for species: $me->{species} in $me->{libdir}/genome!;
    my $bwa_index = $me->Qsub(job_name => "bwaidx",
                              depends => $dep,
                              job_string => $job_string,
                              comment => $comment,
			      prescript => $args{prescript},
			      postscript => $args{postscript},
        );
    return($bwa_index);
}

=head2
    Bowtie_Stats()
=cut
sub Bowtie_Stats {
    my $me = shift;
    my %args = @_;
    my $bt_input = $args{bt_input};
    my $trim_input = $args{trim_input};
    my $basename = $me->{basename};
    my $bt_type = "";
    $bt_type = $args{bt_type} if ($args{bt_type});
    my $depends = "";
    $depends = $args{depends} if ($args{depends});
    my $job_name = "stats";
    $job_name = $args{job_name} if ($args{job_name});
    my $jobid = qq"${basename}_stats";
    my $count_table = "";
    $count_table = $args{count_table} if ($args{count_table});
    my $comment = qq!
## This is a stupidly simple job to collect alignment statistics.
!;
    my $job_string = qq!
original_reads_tmp=\$(grep "^Input Reads" $trim_input | awk '{print \$3}' | sed 's/ //g')
original_reads=\${original_reads_tmp:-0}
reads_tmp=\$(grep "^# reads processed" $bt_input | awk -F: '{print \$2}' | sed 's/ //g')
reads=\${reads_tmp:-0}
one_align_tmp=\$(grep "^# reads with at least one reported" $bt_input | awk -F": " '{print \$2}' | sed 's/ .*//g')
one_align=\${one_align_tmp:-0}
failed_tmp=\$(grep "^# reads that failed to align" $bt_input | awk -F": " '{print \$2}' | sed 's/ .*//g')
failed=\${failed_tmp:-0}
sampled_tmp=\$(grep "^# reads with alignments sampled" $bt_input | awk -F": " '{print \$2}' | sed 's/ .*//g')
sampled=\${sampled_tmp:-0}
rpm_tmp=\$(perl -e "printf(1000000 / \${one_align})" 2>/dev/null)
rpm=\${rpm_tmp:-0}
stat_string=\$(printf "${basename},${bt_type},%s,%s,%s,%s,%s,%s,${count_table}" "\${original_reads}" "\${reads}" "\${one_align}" "\${failed}" "\${sampled}" "\$rpm")
echo "\$stat_string" >> ../bowtie_stats.csv
!;
    my $stats = $me->Qsub(job_name => $job_name,
                          depends => $depends,
                          qsub_queue => "throughput",
                          qsub_cpus => 1,
                          qsub_mem => 1,
                          qsub_wall => "00:10:00",
                          job_string => $job_string,
                          input => $bt_input,
                          comment => $comment,
        );
    return($stats);
}


=head2
    TopHat()
=cut
sub TopHat {
    my $me = shift;
    my %args = @_;
    $me->Check_Options(["species"]);
    my $depends = "";
    if ($args{depends}) {
        $depends = $args{depends};
    }
    my $tophat_args = ' -g 1 ';
    if ($args{tophat_args}) {
        $tophat_args = $args{tophat_args};
    }
    my $tophat_dir = 'tophat_out';
    if ($args{tophat_dir}) {
        $tophat_dir = $args{tophat_dir};
    }
    my $basename = $me->{basename};
    my $bt_reflib = "$me->{libdir}/genome/$me->{species}";
    my $bt_reftest = qq"${bt_reflib}.1.bt2";
    if (!-r $bt_reftest) {
        $me->BT2_Index();
    }
    if (!-r qq"$me->{libdir}/genome/$me->{species}.gtf") {
        print "Missing the gtf file for $me->{species}\n";
        Gff2Gtf("$me->{libdir}/genome/$me->{species}.gff");
    }
    my $inputs = $me->{input};
    my @in = split(/\:/, $inputs);
    $inputs =~ s/\:/ /g;

    my $job_string = qq!mkdir -p ${tophat_dir} && tophat ${tophat_args} --b2-very-sensitive -p 1 -o ${tophat_dir} \\
-G $me->{libdir}/genome/$me->{species}.gtf --no-novel-juncs \\
  $me->{libdir}/genome/$me->{species} \\
  $inputs && samtools index ${tophat_dir}/accepted_hits.bam
!;
    my $comment = qq!## I still have no clue what I am doing when I use tophat...
## However, I know that -g 1 will allow only 1 hit in the case of multihits, but randomly place it
## From the manual:  "If there are more alignments with the same score than this
## number, TopHat will randomly report only this many alignments"
## -N 1 will discard anything with >1 mismatch (default is 2)
## -r adjusts the allowable mean distance between the paired reads
## --mate-std-dev sets the deviation of -r
## --microexon-search will tell it to search short introns for reads >=50
!;
    my $tophat = $me->Qsub(job_name => "th",
                           qsub_cpus => 4,
                           qsub_queue => "long",
                           qsub_wall => "144:00:00",
                           qsub_mem => 10,
                           depends => $depends,
                           job_string => $job_string,
                           comment => $comment,
			   prescript => $args{prescript},
			   postscript => $args{postscript},
        );

    my $accepted = "accepted_hits.bam";
    $accepted = $args{accepted_hits} if ($args{accepted_hits});
    my $count_table = "accepted_hits.count";
    $count_table = $args{count_table} if ($args{count_table});

    my $htmulti = $me->HT_Multi(depends => $tophat->{pbs_id},
                                input => $accepted,
                                job_name => qq"$in[0]_count",
                                ##output => $count_table,
        );

    return({tophat => $tophat, htseq => $htmulti,});
}

=head2
    BT1_Stats()
=cut
sub BT1_Stats {
    my $me = shift;
    my %args = @_;
    my $bt_input = $args{bt_input};
    my $trim_input = $args{trim_input};
    my $basename = $me->{basename};
    my $bt_type = "";
    $bt_type = $args{bt_type} if ($args{bt_type});
    my $depends = "";
    $depends = $args{depends} if ($args{depends});
    my $job_name = "stats";
    $job_name = $args{job_name} if ($args{job_name});
    my $jobid = qq"${basename}_stats";
    my $count_table = "";
    $count_table = $args{count_table} if ($args{count_table});
    my $comment = qq!
## This is a stupidly simple job to collect alignment statistics.
!;
    my $job_string = qq!
original_reads_tmp=\$(grep "^Input Reads" $trim_input | awk '{print \$3}' | sed 's/ //g')
original_reads=\${original_reads_tmp:-0}
reads_tmp=\$(grep "^# reads processed" $bt_input | awk -F: '{print \$2}' | sed 's/ //g')
reads=\${reads_tmp:-0}
one_align_tmp=\$(grep "^# reads with at least one reported" $bt_input | awk -F": " '{print \$2}' | sed 's/ .*//g')
one_align=\${one_align_tmp:-0}
failed_tmp=\$(grep "^# reads that failed to align" $bt_input | awk -F": " '{print \$2}' | sed 's/ .*//g')
failed=\${failed_tmp:-0}
sampled_tmp=\$(grep "^# reads with alignments sampled" $bt_input | awk -F": " '{print \$2}' | sed 's/ .*//g')
sampled=\${sampled_tmp:-0}
rpm_tmp=\$(perl -e "printf(1000000 / \${one_align})" 2>/dev/null)
rpm=\${rpm_tmp:-0}
stat_string=\$(printf "${basename},${bt_type},%s,%s,%s,%s,%s,%s,${count_table}" "\${original_reads}" "\${reads}" "\${one_align}" "\${failed}" "\${sampled}" "\$rpm")
echo "\$stat_string" >> ../bowtie_stats.csv
!;
    my $stats = $me->Qsub(job_name => $job_name,
                          depends => $depends,
                          qsub_queue => "throughput",
                          qsub_cpus => 1,
                          qsub_mem => 1,
                          qsub_wall => "00:10:00",
                          job_string => $job_string,
                          input => $bt_input,
                          comment => $comment,
        );
    return($stats);
}

1;
