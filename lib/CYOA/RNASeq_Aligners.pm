package CYOA;
use common::sense;
use autodie;

=head1 NAME

    CYOA::RNASeq_Aligners - Perform highthroughput sequence alignments with tools like bowtie/tophat/etc

=head1 SYNOPSIS

    use CYOA;
    my $hpgl = new CYOA;
    $hpgl->Bowtie();

=head2 Methods

=over 4

=item C<Bowtie>

    $hpgl->Bowtie() performs a bowtie alignment.  Unless instructed
    otherwise, it will do so with 0 mismatches and with multi-matches
    randomly placed 1 time among the possibilities. (options -v 0 -M 1)


    It checks to see if a bowtie1 compatible index is in
    $libdir/$libtype/indexes/$species, if not it attempts to create
    them.

    It will continue on to convert the bowtie sam output to a
    compressed, sorted, indexed bam file, and pass that to htseq-count
    using a gff file of the same species.

=cut
sub Bowtie {
    my $me = shift;
    my %args = @_;
    $me->Check_Options(args => \%args, needed => ["species","input"]);
    my %bt_jobs = ();
    my $bt_input = $me->{input};
    my $bt_depends_on;
    $bt_depends_on = $args{depends} if ($args{depends});
    my $basename = $me->{basename};
    my $species = $me->{species};
    my $bt_type = "v0M1";
    $bt_type = $args{bt_type} if ($args{bt_type});
    my $bt_args = $me->{bt_args}->{$bt_type};
    my $jobname = qq"bt${bt_type}_${species}";
    $jobname = $args{jobname} if ($args{jobname});
    my $libtype = $me->{libtype};
    $libtype = $args{libtype} if ($args{libtype});
    my $count = 1;
    $count = $args{count} if (defined($args{count}));
    my $bt_dir = qq"outputs/bowtie_$me->{species}";
    if ($args{bt_dir}) {
        $bt_dir = $args{bt_dir};
    }

    if ($bt_input =~ /\.gz$|\.bz2$|\.xz$/ ) {
        print "The input needs to be uncompressed, doing that now.\n" if ($me->{debug});
        ## Previously this step was done on the cluster, but I submitted too many and made the file servers mad
        ## In addition, these jobs are often done as part of a group and I don't want a race condition
        ##my $uncomp = $me->Uncompress(input => $bt_input, depends => $bt_depends_on);
        print "Uncompressing input files before submitting to the cluster.\n";
        if (${bt_input} =~ /\.gz/) {
            qx"nice -n 10 gunzip -f ${bt_input}";
        } elsif (${bt_input} =~ /\.bz2/) {
            qx"nice -n 10 bunzip2 -f ${bt_input}";
        } elsif (${bt_input} =~ /\.xz/) {
            qx"nice -n 10 xz -d -f ${bt_input}";
        }
        $bt_input =  basename($bt_input, ('.gz','.bz2','.xz'));
        $me->{input} = $bt_input;
        ##bt_depends_on = $uncomp->{pbs_id};
    }

    ## Check that the indexes exist
    my $bt_reflib = "$me->{libdir}/${libtype}/indexes/$me->{species}";
    my $bt_reftest = qq"${bt_reflib}.1.ebwt";
    if (!-r $bt_reftest && !$me->{bt1_indexjobs}) {
        $me->{bt1_indexjobs} = 1;
        my $index_job = $me->BT1_Index(depends => $bt_depends_on, libtype => $libtype);
        $bt_jobs{index} = $index_job;
        $bt_depends_on = $index_job->{pbs_id};
    }
    my $bowtie_input_flag = "-q";  ## fastq by default
    $bowtie_input_flag = "-f" if ($me->{input} =~ /\.fasta$/);

    my $cpus = $me->{cpus};
    my $error_file = qq"${bt_dir}/${basename}-${bt_type}.err";
    my $comment = qq!## This is a bowtie1 alignment of $bt_input against
## $bt_reflib using arguments: $bt_args.
## This jobs depended on: $bt_depends_on
!;
    my $aligned_filename = qq"${bt_dir}/${basename}-${bt_type}_aligned_${species}.fastq";
    my $unaligned_filename = qq"${bt_dir}/${basename}-${bt_type}_unaligned_${species}.fastq";
    my $sam_filename = qq"${bt_dir}/${basename}-${bt_type}.sam";
    my $job_string = qq!mkdir -p ${bt_dir} && sleep 10 && bowtie $bt_reflib $bt_args \\
  -p ${cpus} $bowtie_input_flag $bt_input \\
  --un ${unaligned_filename} \\
  --al ${aligned_filename} \\
  -S ${sam_filename} \\
  2>${error_file} \\
  1>${bt_dir}/${basename}-${bt_type}.out
!;
##    my $bt_job = $me->Qsub(job_name => qq"bt1_${bt_type}",
    my $bt_job = $me->Qsub(job_name => $jobname,
                           depends => $bt_depends_on,
                           job_string => $job_string,
                           input => $bt_input,
                           comment => $comment,
                           output => $sam_filename,
                           unaligned => $unaligned_filename,
                           aligned => $aligned_filename,
			   prescript => $args{prescript},
			   postscript => $args{postscript},
                           qsub_queue => "workstation",
        );
    $bt_jobs{bowtie} = $bt_job;

    my $un_comp = $me->Recompress(depends => $bt_job->{pbs_id},
                                  job_name => "xzun",
                                  comment => qq"## Compressing the sequences which failed to align against $bt_reflib using options $bt_args\n",
                                  input => "${bt_dir}/${basename}-${bt_type}_unaligned_${species}.fastq");
    $bt_jobs{unaligned_compression} = $un_comp;

    my $al_comp = $me->Recompress(input => "${bt_dir}/${basename}-${bt_type}_aligned_${species}.fastq",
                                  comment => qq"## Compressing the sequences which successfully aligned against $bt_reflib using options $bt_args",
                                  job_name => "xzal",
                                  depends => $bt_job->{pbs_id},);
    $bt_jobs{aligned_compression} = $al_comp;

    ## BT1_Stats also reads the trimomatic output, which perhaps it should not.
    my $trim_output_file = qq"outputs/${basename}-trimomatic.out";
    my $stats = $me->BT1_Stats(depends => $bt_job->{pbs_id},
                               job_name => "${jobname}_stats",
                               bt_type => $bt_type,
                               count_table => qq"${basename}-${bt_type}.count.xz",
                               trim_input => ${trim_output_file},
                               bt_input => $error_file);

    my $sam_job = $me->Samtools(depends => $bt_job->{pbs_id},
                                job_name => "${jobname}_s2b",
                                input => $bt_job->{output},);

    $bt_jobs{samtools} = $sam_job;
    $me->{output} = $sam_job->{output};
    my $htmulti;
    if ($count) {
        if ($libtype eq 'rRNA') {
            $htmulti = $me->HTSeq(depends => $sam_job->{pbs_id},
                                  suffix => $bt_type,
                                  htseq_type => $me->{htseq_type},
                                  libtype => $libtype,
                                  qsub_queue => "workstation",
                                  input => $sam_job->{output},);
        } else {
            $htmulti = $me->HT_Multi(depends => $sam_job->{pbs_id},
                                     suffix => $bt_type,
                                     htseq_type => $me->{htseq_type},
                                     libtype => $libtype,
                                     qsub_queue => "workstation",
                                     input => $sam_job->{output},);
            $bt_jobs{htseq} = $htmulti;
        }
    }
    return(\%bt_jobs);
}

=item C<Bowtie2>

    $hpgl->Bowtie2() performs a bowtie2 alignment.  Unless instructed
    otherwise, it will do so with 0 mismatches and with multi-matches
    randomly placed 1 time among the possibilities. (options -v 0 -M 1)


    It checks to see if a bowtie1 compatible index is in
    $libdir/$libtype/indexes/$species, if not it attempts to create
    them.

    It will continue on to convert the bowtie sam output to a
    compressed, sorted, indexed bam file, and pass that to htseq-count
    using a gff file of the same species.

=cut
sub Bowtie2 {
    my $me = shift;
    my %args = @_;
    $me->Check_Options(args => \%args, needed => ["species","input", "htseq_type"]);
    my %bt_jobs = ();
    my $libtype = 'genome';
    my $bt_input = $me->{input};
    my $bt_depends_on;
    $bt_depends_on = $args{depends} if ($args{depends});
    my $basename = $me->{basename};
    my $species = $me->{species};
    my $bt2_args = $me->{bt2_args};
    my $jobname = qq"bt2_${species}";
    $jobname = $args{jobname} if ($args{jobname});
    my $count = 1;
    $count = $args{count} if (defined($args{count}));

    my $bt_dir = qq"outputs/bowtie2_$me->{species}";
    if ($args{bt_dir}) {
        $bt_dir = $args{bt_dir};
    }

    if ($bt_input =~ /\.bz2$|\.xz$/ ) {
        print "The input needs to be uncompressed, doing that now.\n" if ($me->{debug});
        my $uncomp = $me->Uncompress(input => $bt_input, depends => $bt_depends_on);
        $bt_input =  basename($bt_input, ('.bz2','.xz'));
        $me->{input} = $bt_input;
        $bt_depends_on = $uncomp->{pbs_id};
    }

    if ($bt_input =~ /\s+/) {
        my @pair_listing = split(/\s+/, $bt_input);
        $bt_input = qq" -1 $pair_listing[0] -2 $pair_listing[1] ";
    }

    ## Check that the indexes exist
    my $bt_reflib = "$me->{libdir}/${libtype}/indexes/$me->{species}";
    my $bt_reftest = qq"${bt_reflib}.1.bt2";
    if (!-r $bt_reftest) {
        my $index_job = $me->BT2_Index(depends => $bt_depends_on, libtype => $libtype);
        $bt_jobs{index} = $index_job;
        $bt_depends_on = $index_job->{pbs_id};
    }
    my $bowtie_input_flag = "-q ";  ## fastq by default
    $bowtie_input_flag = "-f " if ($me->{input} =~ /\.fasta$/);

    my $cpus = $me->{cpus};
    my $error_file = qq"${bt_dir}/${basename}.err";
    my $comment = qq!## This is a bowtie2 alignment of $bt_input against
## $bt_reflib using arguments: $bt2_args.
## This jobs depended on: $bt_depends_on
!;
    my $aligned_filename = qq"${bt_dir}/${basename}_aligned_${species}.fastq";
    my $unaligned_filename = qq"${bt_dir}/${basename}_unaligned_${species}.fastq";
    my $sam_filename = qq"${bt_dir}/${basename}.sam";
    my $job_string = qq!mkdir -p ${bt_dir} && sleep 10 && bowtie2 -x $bt_reflib $bt2_args \\
  -p ${cpus} \\
  $bowtie_input_flag $bt_input \\
  --un ${unaligned_filename} \\
  --al ${aligned_filename} \\
  -S ${sam_filename} \\
  2>${error_file} \\
  1>${bt_dir}/${basename}.out
!;
##    my $bt_job = $me->Qsub(job_name => qq"bt1_${bt_type}",
    my $bt2_job = $me->Qsub(job_name => $jobname,
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
    $bt_jobs{bowtie} = $bt2_job;

    my $un_comp = $me->Recompress(depends => $bt2_job->{pbs_id},
                                  job_name => "xzun",
                                  comment => qq"## Compressing the sequences which failed to align against $bt_reflib using options $bt2_args\n",
                                  input => "${bt_dir}/${basename}_unaligned_${species}.fastq");
    $bt_jobs{unaligned_compression} = $un_comp;

    my $al_comp = $me->Recompress(input => "${bt_dir}/${basename}_aligned_${species}.fastq",
                                  comment => qq"## Compressing the sequences which successfully aligned against $bt_reflib using options $bt2_args",
                                  job_name => "xzal",
                                  depends => $bt2_job->{pbs_id},);
    $bt_jobs{aligned_compression} = $al_comp;

    ## BT1_Stats also reads the trimomatic output, which perhaps it should not.
    my $trim_output_file = qq"outputs/${basename}-trimomatic.out";
    my $stats = $me->BT2_Stats(depends => $bt2_job->{pbs_id},
                               job_name => "${jobname}_stats",
                               count_table => qq"${basename}.count.xz",
                               trim_input => ${trim_output_file},
                               bt_input => $error_file);

    my $sam_job = $me->Samtools(depends => $bt2_job->{pbs_id},
                                job_name => "${jobname}_s2b",
                                input => $bt2_job->{output},);
    $bt_jobs{samtools} = $sam_job;
    $me->{output} = $sam_job->{output};
    my $htmulti;
    if ($count) {
        if ($libtype eq 'rRNA') {
            $htmulti = $me->HTSeq(depends => $sam_job->{pbs_id},
                                  libtype => $libtype,
                                  input => $sam_job->{output},);
        } else {
            $htmulti = $me->HT_Multi(depends => $sam_job->{pbs_id},
                                     libtype => $libtype,
                                     input => $sam_job->{output},);
            $bt_jobs{htseq} = $htmulti;
        }
    }
    return(\%bt_jobs);
}

=item C<BT_Multi>

    $hpgl->BT_Multi() attempts to run multiple bowtie1 runs for a
    given species.  One run is performed for each of a few parameters
    which are kept in the variable '$hpgl->{bt_types}' and generally
    include: 0 mismatch, 1 mismatch, 2 mismatches, 1 randomly placed
    hit, 0 randomly placed hits, or default options.

=cut
sub BT_Multi {
    my $me = shift;
    my %args = @_;
    if ($args{interactive}) {
        print "Run with: cyoa --task rnaseq --method btmulti --input file1:file2 --species critter\n";
    }
    $me->Check_Options(args => \%args, needed => ["species","input", "htseq_type"]);
    my $basename = $me->{basename};
    my $species = $me->{species};
    my $depends_on = $args{depends};
    my %bt_types = %{$me->{bt_args}};
    foreach my $type (keys %bt_types) {
        print "Starting $type\n";
        my $job = $me->Bowtie(bt_type => $type,
                              depends => $depends_on,
                              jobname => qq"bt${type}_${species}",
			      prescript => $args{prescript},
			      postscript => $args{postscript},
            );
    }
}

=item C<Bowtie_RRNA>

    $hpgl->Bowtie_RRNA() performs an alignment against a home-curated
    set of ribosomal RNA/tRNA sequences.  The alignment requires a
    fastq input library and fasta library found in
    'libraries/rRNA/$me->{species}.fasta'

  Example:
    my $rrna = $hpgl->Bowtie_RRNA();
    ## If you want to exclude the rRNA sequences from future alignments:
    my $rrna = $hpgl->Bowtie_RRNA(exclude => 1);

=cut
sub Bowtie_RRNA {
    my $me = shift;
    my %args = @_;
    $me->Check_Options(args => \%args, needed => ["species"]);
    my $basename = $me->{basename};
    $me->{basename} = qq"rRNA_${basename}";
    my $exclude = 0;
    $exclude = $args{exclude} if ($args{exclude});
    my $species = $me->{species};
    my $depends_on = $args{depends};
    my ${bt_dir} = qq"outputs/bowtie_$me->{species}";
    if ($exclude) {
        my $in = $me->{input};
        $args{postscript} = qq"mv $in includerrna_${in} && mv ${bt_dir}/${basename}-rRNA_unaligned_${species}.fastq $in";
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

=item C<BT1_Index>

    $hpgl->BT1_Index() creates a bowtie1 index using
    $hpgl->{species}.fasta and leaves it in the indexes/ directory.

=cut
sub BT1_Index {
    my $me = shift;
    my %args = @_;
    $me->Check_Options(args => \%args, needed => ["species"]);
    my $basename = $me->{basename};
    my $dep = "";
    $dep = $args{depends} if (defined($args{depends}));
    my $libtype = $me->{libtype};
    $libtype = $args{libtype} if ($args{libtype});
    my $job_string = qq!bowtie-build $me->{libdir}/${libtype}/$me->{species}.fasta $me->{libdir}/${libtype}/indexes/$me->{species}
!;
    my $comment = qq!## Generating bowtie1 indexes for species: $me->{species} in $me->{libdir}/${libtype}/indexes!;
    my $bt1_index = $me->Qsub(job_name => "bt1idx",
                              depends => $dep,
                              job_string => $job_string,
                              comment => $comment,
			      prescript => $args{prescript},
			      postscript => $args{postscript},
        );
    return($bt1_index);
}

=item C<BT2_Index>

    $hpgl->BT2_Index() creates a bowtie2 index using
    $hpgl->{species}.fasta and leaves it in the indexes/ directory.

=cut
sub BT2_Index {
    my $me = shift;
    my %args = @_;
    my $basename = $me->{basename};
    my $dep = "";
    $dep = $args{depends};
    my $libtype = $me->{libtype};
    my $libdir = File::Spec->rel2abs($me->{libdir});
    my $job_string = qq!
if [ \! -r "${libdir}/genome/$me->{species}.fa" ]; then
  ln -s ${libdir}/genome/$me->{species}.fasta ${libdir}/genome/$me->{species}.fa
fi
if [ \! -r "${libdir}/genome/indexes/$me->{species}.fa" ]; then
  ln -s ${libdir}/genome/$me->{species}.fasta ${libdir}/genome/indexes/$me->{species}.fa
fi

bowtie2-build $me->{libdir}/genome/$me->{species}.fasta $me->{libdir}/${libtype}/indexes/$me->{species}
!;
    my $comment = qq!## Generating bowtie2 indexes for species: $me->{species} in $me->{libdir}/${libtype}/indexes!;
    my $indexer = $me->Qsub(job_name => "bt2idx",
                            depends => $dep,
                            job_string => $job_string,
                            comment => $comment,
                            prescript => $args{prescript},
                            postscript => $args{postscript},
        );
    return($indexer);
}

=item C<BWA>

    $hpgl->BWA() performs a bwa alignment using both the sam(s|p)e and
    aln algorithms.  It then converts the output (when appropriate) to
    sorted/indexed bam and passes them to htseq.

=cut
sub BWA {
    my $me = shift;
    my %args = @_;
    my %bwa_jobs = ();
    my $bwa_depends_on;
    $me->Check_Options(args => \%args, needed => ["species"]);
    my $basename = $me->{basename};
    my $inputs = $me->{input};
    my @in = split(/\:/, $inputs);
    $inputs =~ s/\:/ /g;
    ## Check that the indexes exist
    my $libtype = $me->{libtype};
    my $bwa_reflib = "$me->{libdir}/${libtype}/indexes/$me->{species}.fasta";
    my $bwa_reftest = qq"${bwa_reflib}.bwa";
    if (!-r $bwa_reftest) {
        my $index_job = $me->BWA_Index(depends => $bwa_depends_on);
        $bwa_jobs{index} = $index_job;
        $bwa_depends_on = $index_job->{pbs_id};
    }
    my $reporter_string = qq"bwa samse ${bwa_reflib} outputs/bwa/${basename}_aln-forward.sai $in[0] \\
1>outputs/bwa/${basename}_aln.sam 2>outputs/bwa/${basename}.samerr";
    my $aln_string = qq"bwa aln ${bwa_reflib} $in[0] 1>outputs/bwa/${basename}_aln-forward.sai 2>outputs/bwa/${basename}_aln-forward.err";
    if (scalar(@in) == 2) {
        $aln_string = qq"${aln_string}
bwa aln ${bwa_reflib} $in[1] 1>outputs/bwa/${basename}_aln-reverse.sai 2>outputs/bwa/${basename}_aln-reverse.err";
        $reporter_string = qq"bwa sampe ${bwa_reflib} \\
outputs/bwa/${basename}_aln-forward.sai outputs/bwa/${basename}_aln-reverse.sai \\
$in[0] $in[1] \\
1>outputs/bwa/${basename}_aln.sam 2>outputs/bwa/${basename}.samerr";
    }

    my $species = $me->{species};
    my $comment = qq!## This is a BWA alignment of $inputs against
## $bwa_reflib.
## It will perform a separate bwa mem run and bwa aln run.!;
    my $job_string = qq!mkdir -p outputs/bwa
bwa mem -a ${bwa_reflib} ${inputs} 2>outputs/bwa/bwa.err 1>outputs/bwa/${basename}_mem.sam
${aln_string}
${reporter_string}
!;
    my $bwa_job = $me->Qsub(job_name => "bwa_${species}",
			    depends => $bwa_depends_on,
			    job_string => $job_string,
			    input => $inputs,
			    comment => $comment,
			    output => qq"outputs/bwa/${basename}_mem.sam",
			    prescript => $args{prescript},
			    postscript => $args{postscript},
        );
    $bwa_jobs{bwa} = $bwa_job;

    my $mem_sam = $me->Samtools(depends => $bwa_job->{pbs_id},
                                job_name => "s2b_mem",
                                input => $bwa_job->{output},);
    $bwa_jobs{samtools_mem} = $mem_sam;
    my $mem_htmulti = $me->HT_Multi(job_name => "ht_mem", depends => $mem_sam->{pbs_id}, input => $mem_sam->{output},);
    $bwa_jobs{htseq_mem} = $mem_htmulti;

    my $aln_sam = $me->Samtools(depends => $bwa_job->{pbs_id},
                                job_name => 's2b_aln',
                                input => qq"outputs/bwa/${basename}_aln.sam",);
    $bwa_jobs{samtools_aln} = $aln_sam;
    my $aln_htmulti = $me->HT_Multi(jobname => "ht_aln", depends => $aln_sam->{pbs_id}, input => $aln_sam->{output},);
    $bwa_jobs{htseq_aln} = $aln_htmulti;

    my $trim_output_file = qq"outputs/${basename}-trimomatic.out";
    my $bwa_stats = $me->BWA_Stats(depends => $aln_sam->{pbs_id},
                                   job_name => 'bwastats',
                                   aln_output => $aln_sam->{output},
                                   mem_output => $mem_sam->{output},
                                   trim_input => $trim_output_file,);
##    my $c = 0;
##    foreach my $ht_out (@{$aln_htmulti}) {
##        $c++;
##        if ($ht_out) {
##            my $count_compression = $me->Recompress(depends => $ht_out->{pbs_id},
##                                                    job_name => "xz_hts${c}",
##                                                    input => $ht_out->{output});
##        }
##    }
    return(\%bwa_jobs);
}

=item C<BWA_Stats>

    $hpgl->BWA_Stats() collects some alignment statistics from bwa.

=cut
sub BWA_Stats {
    my $me = shift;
    my %args = @_;
    my $basename = $me->{basename};
    my $aln_input = $args{aln_output};
    $aln_input = qq"${aln_input}.stats";
    my $mem_input = $args{mem_output};
    $mem_input = qq"${mem_input}.stats";
    my $output = qq"${basename}.stats";
    my $trim_input = $args{trim_input};

    my $depends = "";
    $depends = $args{depends} if ($args{depends});
    my $job_name = "bwa_stats";
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
reads_tmp=\$(grep "^Total reads: " $aln_input | awk '{print \$3}' | sed 's/ //g')
reads=\${reads_tmp:-0}
aln_aligned_tmp=\$(grep "^Mapped reads" $aln_input | awk '{print \$3}' | sed 's/ .*//g')
aln_aligned=\${aln_aligned_tmp:-0}
mem_aligned_tmp=\$(grep "^Mapped reads" $mem_input | awk '{print \$3}' | sed 's/ .*//g')
mem_aligned=\${mem_aligned_tmp:-0}
rpm_tmp=\$(perl -e "printf(1000000 / \${aligned})" 2>/dev/null)
rpm=\${rpm_tmp:-0}
stat_string=\$(printf "${basename},%s,%s,%s,%s,%s,${count_table}" "\${original_reads}" "\${reads}" "\${aln_aligned}" "\${mem_aligned}" "\$rpm")
echo "\$stat_string" >> outputs/bwa_stats.csv
!;
    my $stats = $me->Qsub(job_name => $job_name,
                          depends => $depends,
                          qsub_queue => "throughput",
                          qsub_cpus => 1,
                          qsub_mem => 1,
                          qsub_wall => "00:10:00",
                          job_string => $job_string,
                          input => $aln_input,
                          comment => $comment,
        );
    return($stats);
}

=item C<BWA_Index>

    $hpgl->BWA_Index() creates bwa indexes.

=cut
sub BWA_Index {
    my $me = shift;
    my %args = @_;
    my $basename = $me->{basename};
    my $dep = "";
    $dep = $args{depends};
    my $libtype = $me->{libtype};
    my $job_string = qq!
start=\$(pwd)
cd \$start/$me->{libdir}/${libtype}/indexes && ln -s \$start/$me->{libdir}/genome/$me->{species}.fasta . && bwa index $me->{species}.fasta
cd \$start
!;
    my $comment = qq!## Generating bwa indexes for species: $me->{species} in $me->{libdir}/${libtype}/indexes!;
    my $bwa_index = $me->Qsub(job_name => "bwaidx",
                              depends => $dep,
                              job_string => $job_string,
                              comment => $comment,
			      prescript => $args{prescript},
			      postscript => $args{postscript},
        );
    return($bwa_index);
}

=item C<BT2_Stats>

    $hpgl->BT2_Stats() collects alignment statistics from bowtie 2.

=cut
sub BT2_Stats {
    my $me = shift;
    my %args = @_;
    my $bt_input = $args{bt_input};
    my $trim_input = $args{trim_input};
    my $basename = $me->{basename};
    my $bt_type = "";
    $bt_type = $args{bt_type} if ($args{bt_type});
    my $depends = "";
    $depends = $args{depends} if ($args{depends});
    my $job_name = "bt2_stats";
    $job_name = $args{job_name} if ($args{job_name});
    my $jobid = qq"${basename}_stats";
    my $count_table = "";
    $count_table = $args{count_table} if ($args{count_table});
    my $comment = qq!
## This is a stupidly simple job to collect alignment statistics.
!;
    my $job_string = qq!
original_reads_tmp=\$(grep "^Input Reads" "${trim_input}" 2>/dev/null | awk '{print \$3}' | sed 's/ //g')
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
echo "\$stat_string" >> outputs/bowtie_stats.csv
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

=item C<Kallisto_Index

    $hpgl->Kallisto_Index() uses kallisto and an annotated_CDS fasta sequence library to do the expected.

=cut
sub Kallisto_Index {
    my $me = shift;
    my %args = @_;
    $me->Check_Options(args => \%args, needed => ["species","input"]);
    my $basename = $me->{basename};
    my $dep = "";
    $dep = $args{depends};
    my $libtype = $me->{libtype};
    my $input = File::Spec->rel2abs($me->{input});
    unless (-r qq"${input}") {
        die("The indexing operation for kallisto will fail because ${input} does not exist.")
    }

    my $job_string = qq!
kallisto index -i $me->{libdir}/${libtype}/indexes/$me->{species}.idx ${input}
!;
    my $comment = qq!## Generating kallisto indexes for species: $me->{species} in $me->{libdir}/${libtype}/indexes!;
    my $jobid = $me->Qsub(job_name => "kalidx",
                          depends => $dep,
                          job_string => $job_string,
                          comment => $comment,
			  prescript => $args{prescript},
			  postscript => $args{postscript},
        );
    return($jobid);
}

=item C<Kallisto>

    $hpgl->Kallisto() should perform a kallisto alignment, I have not
    yet properly tested this  TODO!

=cut
sub Kallisto {
    my $me = shift;
    my %args = @_;
    $me->Check_Options(args => \%args, needed => ["species","input"]);
    my %ka_jobs = ();
    my $ka_depends_on;
    $ka_depends_on = $args{depends} if ($args{depends});
    my $libtype = 'genome';
    $libtype = $args{libtype} if ($args{libtype});
    my $basename = $me->{basename};
    my $species = $me->{species};
    my $ka_input = $me->{input};

    my $jobname = qq"kall_${species}";
    $jobname = $args{jobname} if ($args{jobname});
    my $ka_args = qq"";

    if ($ka_input =~ /\.bz2$|\.xz$/ ) {
        my $uncomp = $me->Uncompress(input => $ka_input, depends => $ka_depends_on);
        $ka_input =  basename($ka_input, ('.gz','.bz2','.xz'));
        $me->{input} = $ka_input;
        $ka_depends_on = $uncomp->{pbs_id};
    }
    my @input_test = split(/\:/, $ka_input);
    if (scalar(@input_test) == 1) {
        $ka_args .= " --single -l 40 -s 10 ";
    } elsif (scalar(@input_test) == 2) {
        ## Do nothing
        print "There are 2 inputs\n";
    } else {
        print "There are not 1 nor 2 inputs.\n";
    }
    $ka_input =~ s/\:/ /g;

    ## Check that the indexes exist
    my $ka_reflib = "$me->{libdir}/${libtype}/indexes/$me->{species}.idx";
    if (!-r $ka_reflib) {
        my $index_job = $me->Kallisto_Index(depends => $ka_depends_on, libtype => $libtype);
        $ka_jobs{index} = $index_job;
        $ka_depends_on = $index_job->{pbs_id};
    }

    my $outdir = qq"outputs/kallisto_${species}";
    my $error_file = qq"${outdir}/kallisto_${species}.stderr";
    my $output_sam = qq"${outdir}/kallisto_${species}.sam";
    my $output_bam = qq"${outdir}/kallisto_${species}.bam";
    my $output_stats = qq"${outdir}/kallisto_${species}.stats";
    my $sorted_bam = qq"${outdir}/kallisto_${species}-sorted";
    my $comment = qq!## This is a kallisto pseudoalignment of ${ka_input} against
## ${ka_reflib}.
## This jobs depended on: ${ka_depends_on}
## Other candidates for making a pretty count table include:
##  perl -F'\\t' -a -n -e 'print "\$F[0] \$F[3]\\n"' ${outdir}/abundance.tsv > ${outdir}/abundance.count
##   awk '{printf("%s %s\\n", \$1, \$4)}' ${outdir}/abundance.tsv > ${outdir}/abundance.count
## The sam->bam conversion is copy/pasted from Sam2Bam() but I figured why start another job
## because kallisto is so fast
!;
    my $job_string = qq!mkdir -p ${outdir} && sleep 10 && \\
kallisto quant ${ka_args} --plaintext --pseudobam -t 4 -b 100 -o ${outdir} -i ${ka_reflib} \\
  ${ka_input} 2>${error_file} 1>${output_sam} && \\
  cut -d "	" -f 1,4 ${outdir}/abundance.tsv > ${outdir}/abundance.count && \\
  samtools view -u -t ${ka_reflib} -S ${output_sam} 1>${output_bam} && \\
  samtools sort -l 9 ${output_bam} ${sorted_bam} && \\
  rm ${output_bam} && mv ${sorted_bam}.bam ${output_bam} && samtools index ${output_bam} && \\
  bamtools stats -in ${output_bam} 2>${output_stats} 1>&2
!;
    my $ka_job = $me->Qsub(job_name => qq"${jobname}",
                           depends => $ka_depends_on,
                           job_string => $job_string,
                           input => $ka_input,
                           comment => $comment,
                           qsub_queue => "workstation",
			   prescript => $args{prescript},
			   postscript => $args{postscript},
        );
    $ka_jobs{kallisto} = $ka_job;

    return(\%ka_jobs);
}

=item C<Tophat>

    $hpgl->Tophat() runs... guess... tophat!  It also sorts/indexes
    the accepted_hits file, collects some statistics, and passes the
    hits to htseq-count.

=cut
sub Tophat {
    my $me = shift;
    my %args = @_;
    $me->Check_Options(args => \%args, needed => ["species","input","htseq_type"]);
    my $depends = "";
    my $tophat_cpus = 4;
    if ($args{depends}) {
        $depends = $args{depends};
    }
    my $inputs = $me->{input};
    my @in = split(/:/, $inputs);
    $inputs =~ s/:/ /g;
    my $number_inputs = scalar(@in);
    my $tophat_args = ' -g 1 --microexon-search --b2-very-sensitive ';
    if ($args{tophat_args}) {
        $tophat_args = $args{tophat_args};
    }
    ## $tophat_args .= ' --no-mixed --no-discordant ' if (scalar(@in) > 1);
    ## $tophat_args .= ' ' if (scalar(@in) > 1);

    my $tophat_queue = $me->{qsub_queue};
    my $tophat_walltime = '18:00:00';
    my $tophat_mem = 8;
    if ($me->{species} eq 'hsapiens' or $me->{species} eq 'mmusculus') {
        $tophat_queue = 'workstation';
        $tophat_walltime =  "144:00:00";
        $tophat_mem = 20;
    }

    my $tophat_dir = qq"outputs/tophat_$me->{species}";
    if ($args{tophat_dir}) {
        $tophat_dir = $args{tophat_dir};
    }
    my $basename = $me->{basename};
    my $libtype = $me->{libtype};
    my $bt_reflib = "$me->{libdir}/${libtype}/indexes/$me->{species}";
    my $bt_reftest = qq"${bt_reflib}.1.bt2";
    my $index_job = undef;
    if (!-r $bt_reftest) {
        $index_job = $me->BT2_Index(depends => $depends);
        ## Use a colon separated append to make tophat depend on multiple jobs
        ## $depends .= qq":$index_job->{jobid}";
        ## Or just replace the dependency string with this job's and insert it into the stack
        $depends = $index_job->{pbs_id};
    }
    my $gtf_file = qq"$me->{libdir}/genome/$me->{species}.gtf";
    if (!-r $gtf_file) {
        print "Missing the gtf file for $me->{species}\n";
        print "Using the gff file.\n";
        $gtf_file =~ s/\.gtf/\.gff/;
        ##my $written = $me->Gff2Gtf(gff => "$me->{libdir}/genome/$me->{species}.gff");
        ##print STDERR "Gff2Gtf wrote $written features\n";
    }

    my $spliced = 0;
    if ($me->{spliced}) {
        $spliced = 1;
    }
    if ($args{spliced}) {
        $spliced = 1;
    }
    my $jobname = qq"th_$me->{species}";
    my $job_string = qq!mkdir -p ${tophat_dir} && tophat ${tophat_args} \\
  -G ${gtf_file} \\
  -p ${tophat_cpus} -o ${tophat_dir} \\
!;
    if ($spliced) {
        $job_string .= qq! --no-novel-juncs \\
!;
    }
    $job_string .= qq!  $me->{libdir}/genome/indexes/$me->{species} \\
  ${inputs} 2>outputs/tophat.out 1>&2 && \\
 samtools sort -l 9 -n ${tophat_dir}/accepted_hits.bam ${tophat_dir}/accepted_sorted && \\
 mv ${tophat_dir}/accepted_sorted.bam ${tophat_dir}/accepted_hits.bam && \\
 samtools index ${tophat_dir}/accepted_hits.bam && \\
 samtools sort -l 9 -n ${tophat_dir}/unmapped.bam ${tophat_dir}/unmapped_sorted && \\
 mv ${tophat_dir}/unmapped_sorted.bam ${tophat_dir}/unmapped.bam && \\
 samtools index ${tophat_dir}/unmapped.bam
!;
    if ($number_inputs > 1) {
        $job_string .= qq!
if [ -r "${tophat_dir}/accepted_hits.bam" ]; then
  samtools view -b -f 2 ${tophat_dir}/accepted_hits.bam > ${tophat_dir}/accepted_paired.bam && samtools index ${tophat_dir}/accepted_paired.bam
fi
!;
    }
    my $comment = qq!## I still have no clue what I am doing when I use tophat...
## However, I know that -g 1 will allow only 1 hit in the case of multihits, but randomly place it
## From the manual:  "If there are more alignments with the same score than this
## number, TopHat will randomly report only this many alignments"
## -N 1 will discard anything with >1 mismatch (default is 2)
## -r adjusts the allowable mean distance between the paired reads
## --mate-std-dev sets the deviation of -r
## --microexon-search will tell it to search short exons for reads >=50
!;
    my $tophat = $me->Qsub(job_name => ${jobname},
                           qsub_cpus => ${tophat_cpus},
                           qsub_queue => ${tophat_queue},
                           qsub_wall => ${tophat_walltime},
                           qsub_mem => ${tophat_mem},
                           depends => ${depends},
                           job_string => ${job_string},
                           comment => ${comment},
			   prescript => $args{prescript},
			   postscript => $args{postscript},
        );

    ## Set the input for htseq
    my $accepted = "${tophat_dir}/accepted_hits.bam";
    $accepted = $args{accepted_hits} if ($args{accepted_hits});
    my $count_table = "accepted_hits.count";
    $count_table = $args{count_table} if ($args{count_table});
    my $htmulti = $me->HT_Multi(depends => $tophat->{pbs_id},
                                input => $accepted,
                                htseq_type => $me->{htseq_type},
                                job_name => qq"hts_$me->{species}",);
    $tophat->{htseq} = $htmulti;
    ## Perform a separate htseq run using only the successfully paired hits
    if ($number_inputs > 1) {
        my $ht_paired = $me->HT_Multi(depends => $tophat->{pbs_id},
                                      htseq_type => $me->{htseq_type},
                                      input => qq"${tophat_dir}/accepted_paired.bam",
                                      job_name => qq"htsp_$me->{species}",);
    }
    ## Tophat_Stats also reads the trimomatic output, which perhaps it should not.
    my $trim_output_file = qq"outputs/${basename}-trimomatic.out";
    my $unaccepted = $accepted;
    $unaccepted =~ s/accepted_hits/unmapped/g;
    my $input_read_info = $accepted;
    $input_read_info =~ s/accepted_hits\.bam/prep_reads\.info/g;
    my $stats = $me->Tophat_Stats(depends => $tophat->{pbs_id},
                                  job_name => "tpstats_$me->{species}",
                                  count_table => qq"${count_table}.xz",
                                  trim_input => ${trim_output_file},
                                  accepted_input => $accepted,
                                  unaccepted_input => $unaccepted,
                                  prep_input => $input_read_info,);
    $tophat->{stats} = $stats;

    return($tophat);
}

=item C<BT1_Stats>

    $hpgl->BT1_Stats() collects some alignment statistics from
    bowtie1.

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
    my $stat_output = qq"outputs/bowtie_stats.csv";
    my $comment = qq!
## This is a stupidly simple job to collect alignment statistics.
!;
    my $job_string = qq!
if [ \! -r "${stat_output}" ]; then
  echo "name,type,original_reads,reads,one_hits,failed,samples,rpm,count_table" > $stat_output
fi
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
echo "\$stat_string" >> ${stat_output}
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

=item C<Tophat_Stats>

    $hpgl->Tophat_Stats() collects alignment statistics from the
    accepted_hits.bam/unaligned.bam files generated by a tophat run.

=cut
sub Tophat_Stats {
    my $me = shift;
    my %args = @_;
    my $accepted_input = $args{accepted_input};
    my $accepted_output = qq"${accepted_input}.stats";
    my $unaccepted_input = $args{unaccepted_input};
    my $unaccepted_output = qq"${unaccepted_input}.stats";
    my $read_info = $args{prep_input};
    my $trim_input = $args{trim_input};
    my $basename = $me->{basename};
    my $depends = "";
    $depends = $args{depends} if ($args{depends});
    my $job_name = "stats";
    $job_name = $args{job_name} if ($args{job_name});
    my $jobid = qq"${basename}_stats";
    my $count_table = "";
    $count_table = $args{count_table} if ($args{count_table});
    my $comment = qq!
## This is a stupidly simple job to collect tophat alignment statistics.
!;
    my $job_string = qq!
if [ \! -r outputs/tophat_stats.csv ]; then
  echo "basename,species,original_reads,aligned_reads,failed_reads,rpm,count_table" > outputs/tophat_stats.csv
fi
bamtools stats < ${accepted_input} 2>${accepted_output} 1>&2 && bamtools stats < ${unaccepted_input} 2>${unaccepted_output} 1>&2
original_reads_tmp=\$(grep "^Input Reads" $trim_input | awk '{print \$3}' | sed 's/ //g')
original_reads=\${original_reads_tmp:-0}
reads_tmp=\$(grep "^reads_in " $read_info | awk -F= '{print \$2}' | sed 's/ //g')
reads=\${reads_tmp:-0}
aligned_tmp=\$(grep "^Total reads" $accepted_output | awk '{print \$3}' | sed 's/ .*//g')
aligned=\${aligned_tmp:-0}
failed_tmp=\$(grep "^Total reads" $unaccepted_output | awk '{print \$3}' | sed 's/ .*//g')
failed=\${failed_tmp:-0}
rpm_tmp=\$(perl -e "printf(1000000 / \${aligned})" 2>/dev/null)
rpm=\${rpm_tmp:-0}
stat_string=\$(printf "${basename},$me->{species},%s,%s,%s,%s,%s,${count_table}" "\${original_reads}" "\${reads}" "\${aligned}" "\${failed}" "\$rpm")
echo "\$stat_string" >> outputs/tophat_stats.csv
!;
    my $stats = $me->Qsub(job_name => $job_name,
                          depends => $depends,
                          qsub_queue => "throughput",
                          qsub_cpus => 1,
                          qsub_mem => 1,
                          qsub_wall => "00:10:00",
                          job_string => $job_string,
                          input => $accepted_input,
                          comment => $comment,
        );
    return($stats);
}

=back

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

    L<bowtie> L<bowtie2> L<tophat> L<bwa> L<kallisto> L<samtools> L<htseq>

=cut

1;
