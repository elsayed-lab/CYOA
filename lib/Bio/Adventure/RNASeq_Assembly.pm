package Bio::Adventure::RNASeq_Assembly;
## LICENSE: gplv2
## ABSTRACT:  Kitty!
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use File::Basename;
use File::Spec;
use File::Which qw"which";
use Parse::CSV;

=head1 NAME

    Bio::Adventure::RNASeq_Assembly - Perform some denovo assembly tasks

=head1 SYNOPSIS

    use Bio::Adventure;
    my $hpgl = new Bio::Adventure;
    $hpgl->Trinity();

=head2 Methods

=over 4

=item C<Transdecoder>

    $hpgl->Transdecoder() submits a trinity denovo sequence assembly and
    runs its default post-processing tools.

=cut
sub Transdecoder {
    my ($class, %args) = @_;
    my $check = which('TransDecoder.LongOrfs');
    die("Could not find transdecoder in your PATH.") unless($check);
    my $transdecoder_exe_dir = dirname($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
    );
    my $transdecoder_input = File::Spec->rel2abs($options->{input});
    my $job_basename = $class->Get_Job_Name();
    my $output_dir = qq"outputs/trinity_${job_basename}";
    my $comment = qq!## This is a transdecoder submission script
!;
    my $jstring = qq!mkdir -p ${output_dir} && cd ${output_dir} && \\
  TransDecoder.LongOrfs -t ${transdecoder_input} \\
    2>${output_dir}/transdecoder_longorfs_${job_basename}.err \\
    1>${output_dir}/transdecoder_longorfs_${job_basename}.out
TransDecoder.Predict -t ${transdecoder_input} \\
  2>${output_dir}/transdecoder_predict_${job_basename}.err \\
  1>${output_dir}/transdecoder_predict_${job_basename}.out
${transdecoder_exe_dir}/util/cdna_alignment_orf_to_genome_orf.pl \\
  ${output_dir}/transcripts.fasta.transdecoder.gff3 \\
  ${output_dir}/transcripts.gff3 \\
  ${transdecoder_input} \\
  2>${output_dir}/cdna_alignment_${job_basename}.err \\
  1>${output_dir}/transcripts.fasta.transdecoder.genome.gff3
!;
    my $transdecoder_job = $class->Submit(
        comment => $comment,
        jname => "transdecoder_${job_basename}",
        jprefix => "47",
        jstring => $jstring,
        output => qq"outputs/transdecoder.sbatchout",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => "workstation",
    );
    return($transdecoder_job);
}

=item C<Trinotate>

    $hpgl->Trinotate() submits a trinity denovo sequence assembly and
    runs its default post-processing tools.

=cut
sub Trinotate {
    my ($class, %args) = @_;
    my $check = which('Trinotate');
    die("Could not find trinotate in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
    );
    my $job_basename = $class->Get_Job_Name();
    my $output_dir = qq"outputs/trinity_${job_basename}";
    my $trinotate_exe_dir = dirname($check);
    my $input_dir = dirname($options->{input});
    my $input_name = basename($options->{input});
    my $comment = qq!## This is a trinotate submission script
!;
    my $jstring = qq!cd ${input_dir}
${trinotate_exe_dir}/auto/autoTrinotate.pl \\
  --Trinotate_sqlite ${trinotate_exe_dir}/sample_data/Trinotate.boilerplate.sqlite \\
  --transcripts $input_name \\
  --gene_to_trans_map Trinity.fasta.gene_trans_map \\
  --conf ${trinotate_exe_dir}/auto/conf.txt \\
  --CPU 6 \\
  2>trinotate_${job_basename}.err \\
  1>trinotate_${job_basename}.out
!;
    my $trinotate_job = $class->Submit(
        cpus => 6,
        comment => $comment,
        depends => $options->{depends},
        jname => "trinotate_${job_basename}",
        jprefix => "48",
        jstring => $jstring,
        mem => 80,
        output => qq"trinotate.sbatchout",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => "large",
        walltime => "144:00:00",
    );
    return($trinotate_job);
}

=item C<Trinity>

    $hpgl->Trinity() submits a trinity denovo sequence assembly and
    runs its default post-processing tools.

=cut
sub Trinity {
    my ($class, %args) = @_;
    my $check = which('Trinity');
    die("Could not find trinity in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        contig_length => 600,
        required => ['input'],
    );
    my $job_basename = $class->Get_Job_Name();
    my %trin_jobs = ();
    my $trin_depends_on;
    my $output_dir = qq"outputs/trinity_${job_basename}";
    my $input_string = "";
    if ($options->{input} =~ /\:|\;|\,|\s+/) {
        my @in = split(/\:|\;|\,|\s+/, $options->{input});
        $input_string = qq"--left $in[0] --right $in[1] ";
    } else {
        $input_string = qq"--single $options->{input} ";
    }
    my $comment = qq!## This is a trinity submission script
!;
    my $jstring = qq!mkdir -p ${output_dir} && \\
  Trinity --seqType fq --min_contig_length $options->{contig_length} --normalize_reads \\
    --trimmomatic --max_memory 90G --CPU 6 \\
    --output ${output_dir} \\
    ${input_string} \\
    2>${output_dir}/trinity_${job_basename}.err \\
    1>${output_dir}/trinity_${job_basename}.out
!;
    my $trin_job = $class->Submit(
        cpus => 6,
        comment => $comment,
        depends => $trin_depends_on,
        jname => "trin_${job_basename}",
        jprefix => "45",
        jstring => $jstring,
        mem => 96,
        output => qq"outputs/trinity.sbatchout",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => "large",
        walltime => "144:00:00",
    );
    my $rsem_job = Bio::Adventure::RNASeq_Assembly::Trinity_Post(
        $class,
        %args,
        depends => $trin_job->{pbs_id},
        jname => "trin_rsem",
        input => $options->{input},
    );
    my $trinotate_job = Bio::Adventure::RNASeq_Assembly::Trinotate(
        $class,
        %args,
        depends => $trin_job->{pbs_id},
        jname => "trinotate",
        input => qq"${output_dir}/Trinity.fasta",
    );
    my $jobs = {
        trinity => $trin_job,
        trinity_post => $rsem_job,
        trinotate => $trinotate_job,
    };
    return($jobs);
}

sub Trinity_Post {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jname => "trin_rsem",
    );
    my $job_basename = $class->Get_Job_Name();
    my $trinity_out_dir = qq"outputs/trinity_${job_basename}";

    my $rsem_input = qq"${trinity_out_dir}/Trinity.fasta";
    my $trinity_path = which('Trinity');
    my $trinity_exe_dir = dirname($trinity_path);
    my $input_string = "";
    if ($options->{input} =~ /\:|\;|\,|\s+/) {
        my @in = split(/\:|\;|\,|\s+/, $options->{input});
        $input_string = qq"--left $in[0] --right $in[1] ";
    } else {
        $input_string = qq"--single $options->{input} ";
    }

    my $comment = qq!## This is a trinity post-processing submission script.
!;
    my $jstring = qq!
start=\$(pwd)
cd ${trinity_out_dir}
${trinity_exe_dir}/util/TrinityStats.pl Trinity.fasta \\
  2>${trinity_out_dir}/trinpost_stats.err \\
  1>${trinity_out_dir}/trinpost_stats.out &

${trinity_exe_dir}/util/align_and_estimate_abundance.pl \\
  --output_dir align_estimate.out \\
  --transcripts ${rsem_input} \\
  --seqType fq \\
  ${input_string} \\
  --est_method RSEM \\
  --aln_method bowtie \\
  --trinity_mode --prep_reference \\
  2>trinpost_align_estimate.err \\
  1>trinpost_align_estimate.out

${trinity_exe_dir}/util/abundance_estimates_to_matrix.pl \\
  --est_method RSEM \\
  --gene_trans_map Trinity.fasta.gene_trans_map \\
  align_estimate.out/RSEM.isoform.results \\
  2>trinpost_estimate_to_matrix.err \\
  1>trinpost_estimate_to_matrix.out

${trinity_exe_dir}/util/SAM_nameSorted_to_uniq_count_stats.pl \\
  bowtie_out/bowtie_out.nameSorted.bam \\
  2>trinpost_count_stats.err \\
  1>trinpost_count_stats.out

cd \${start}
!;
    my $trinpost_job = $class->Submit(
        comment => $comment,
        input => $options->{input},
        depends => $options->{depends},
        jname => "trinpost_$job_basename",
        jprefix => "46",
        jstring => $jstring,
        mem => 90,
        output => qq"outputs/trinitypost.sbatchout",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => "large",
        walltime => "144:00:00",
##        queue => "long",
##        walltime => "144:00:00",
    );
    return($trinpost_job);
}

sub Extract_Trinotate {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jname => "trin_rsem",
        output => 'interesting.fasta',
        evalue => 1e-10,
        identity => 70,
    );
    my $job_basename = $class->Get_Job_Name();
    my $trinity_out_dir = qq"outputs/trinity_${job_basename}";


    my $input = FileHandle->new("<$options->{input}");
    my $parser = Parse::CSV->new(
        handle => $input,
        sep_char => "\t",
        names => 1,
    );

    my $count = 0;
    my $all_annotations = [];
    my $out = FileHandle->new(">$options->{output}");
    while (my $object = $parser->fetch) {
        my $filled_in = Read_Write_Annotation($class, $object,
                                              db => $all_annotations,
                                              count => $count,
                                              evalue => $options->{evalue},
                                              identity => $options->{identity},
                                              fh => $out);
        $count = $count + 1;
    }
    $out->close();
    $input->close();
}


sub Read_Write_Annotation {
    my ($class, $object, %args) = @_;
    my $annotations = $args{db};
    my $element_number = $args{count};
    my $fh = $args{fh};
    ## $object is a hash reference containing the various csv fields.
    ## Extract them by name and dump them to the array of material
    my $element = {};
    my $filled_in = 0;
    my $ids = {};

    foreach my $key (keys %{$object}) {
        my @result_list = ();
        my @field_elements = ();
        my $field_text = $object->{$key};
        if ($field_text =~ /\`/) {
            @field_elements = split(/\`/, $field_text);
        } else {
            $field_elements[0] = $field_text;
        }

        for my $t (0..$#field_elements) {
            my $ret = {};
            if ($key eq 'sprot_Top_BLASTX_hit') {
                if ($field_text eq '.') {
                    ## Null case: if a . then just drop out.
                    $element->{swissprot} = undef;
                } else {
                    ## The swissprot elements look something like:
                    ## Name        ## ARIA_ARATH^
                    ## Name again? ## ARIA_ARATH^
                    ## Query range ## Q:7-573,H:513-701^
                    ## %identity   ## 75.66%ID^
                    ## E-value     ## E:3e-102^
                    ## Full name   ## RecName: Full=ARM REPEAT PROTEIN INTERACTING WITH ABF2;^
                    ## Ontology    ## Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliophyta; eudicotyledons; Gunneridae; Pentapetalae; rosids; malvids; Brassicales; Brassicaceae; Camelineae; Arabidopsis',
                    ## Start by pulling apart the data using the chosen (somewhat odd) separator '^'.
                    my ($name, $name_again, $query_coords, $identity, $e_value, $fullname, $ontology) = split(/\^/, $field_text);
                    ## Clean up the fields a little
                    $identity =~ s/\%ID//g;
                    $e_value =~ s/E://g;
                    $fullname =~ s/RecName: //g;
                    $fullname =~ s/Full=//g;
                    $fullname =~ s/\;//g;
                    my ($query, $h) = split(/\,/, $query_coords);
                    $query =~ s/Q\://g;
                    $h =~ s/H\://g;
                    $name =~ s/Full=//g;
                    my @ontology_list = split(/; /, $ontology);
                    my $species = qq"$ontology_list[$#ontology_list - 1] $ontology_list[$#ontology_list]";
                    ## Put them back into the $ret
                    $ret->{name} = $name;
                    $ret->{coords} = $query_coords;
                    $ret->{identity} = $identity;
                    $ret->{e_value} = $e_value;
                    $ret->{fullname} = $fullname;
                    $ret->{species} = $species;
                    $filled_in = $filled_in + 6;
                    ## And refill field_elements with this information.
                    $field_elements[$t] = $ret;
                    ## Finally, fill in the swissprot entry with this information.
                    $element->{swissprot} = \@field_elements;
                }
            } elsif ($key eq 'gene_ontology_pfam') {
                if ($field_text eq '.') {
                    $element->{pfam_go} = undef;
                } else {
                    my ($id, $ontology, $name) = split(/\^/, $field_text);
                    $ret->{id} = $id;
                    $ret->{ontology} = $ontology;
                    $ret->{name} = $name;
                    $filled_in = $filled_in + 3;
                    $field_elements[$t] = $ret;
                    $element->{pfam_go} = \@field_elements;
                }
            } elsif ($key eq 'gene_ontology_blast') {
                if ($field_text eq '.') {
                    $element->{blast_go} = undef;
                } else {
                    my ($id, $ontology, $name) = split(/\^/, $field_text);
                    $ret->{id} = $id;
                    $ret->{ontology} = $ontology;
                    $ret->{name} = $name;
                    $filled_in = $filled_in + 3;
                    $field_elements[$t] = $ret;
                    $element->{blast_go} = \@field_elements;
                }
            } elsif ($key eq 'RNAMMER') {
                if ($field_text eq '.') {
                    $element->{rnammer} = undef;
                } else {
                    my ($name, $region) = split(/\^/, $field_text);
                    $ret->{name} = $name;
                    $ret->{region} = $region;
                    $filled_in = $filled_in + 2;
                    $field_elements[$t] = $ret;
                    $element->{rnammer} = \@field_elements;
                }
            } elsif ($key eq 'eggnog') {
                if ($field_text eq '.') {
                    $element->{eggnog} = undef;
                } else {
                    my ($id, $name) = split(/\^/, $field_text);
                    $ret->{id} = $id;
                    $ret->{name} = $name;
                    $filled_in = $filled_in + 2;
                    $field_elements[$t] = $ret;
                    $element->{eggnog} = \@field_elements;
                }
            } elsif ($key eq 'transcript_id') {
                $filled_in = $filled_in + 1;
                $element->{transcript_id} = $field_text;
            } elsif ($key eq 'Kegg') {
                if ($field_text eq '.') {
                    $element->{kegg} = undef;
                } else {
                    my ($id, $name) = split(/\^/, $field_text);
                    $ret->{id} = $id;
                    $ret->{name} = $name;
                    $filled_in = $filled_in + 2;
                    $field_elements[$t] = $ret;
                    $element->{kegg} = \@field_elements;
                }
            } elsif ($key eq 'prot_coords') {
                if ($field_text eq '.') {
                    $element->{prot_coords} = undef;
                } else {
                    $filled_in = $filled_in + 1;
                    $element->{prot_coords} = $field_text;
                }
            } elsif ($key eq 'prot_id') {
                $element->{prot_id} = $field_text;
            } elsif ($key eq 'transcript') {
                $element->{sequence} = $field_text;
            } elsif ($key eq 'TmHMM') {
                ## An example field
                ## ExpAA=124.31^PredHel=6^Topology=i59-78o83-105i126-148o152-174i187-209o229-251i
                if ($field_text eq '.') {
                    $element->{tmhmm} = undef;
                } else {
                    my ($expaa, $pred_helixes, $topology) = split(/\^/, $field_text);
                    $expaa =~ s/ExpAA=//g;
                    $pred_helixes =~ s/PredHel=//g;
                    $topology =~ s/Topology=//g;
                    $ret->{expaa} = $expaa;
                    $ret->{pred_helixes} = $pred_helixes;
                    $ret->{topology} = $topology;
                    $filled_in = $filled_in + 3;
                    $field_elements[$t] = $ret;
                    $element->{tmhmm} = \@field_elements;
                }
            } elsif ($key eq 'SignalP') {
                ## sigP:1^18^0.613^YES
                if ($field_text eq '.') {
                    $element->{signalp} = undef;
                } else {
                    my ($first, $second, $third, $boolean) = split(/\^/, $field_text);
                    $first =~ s/sigP\://g;
                    $ret->{first} = $first;
                    $ret->{second} = $second;
                    $ret->{third} = $third;
                    $ret->{boolean} = $boolean;
                    $filled_in = $filled_in + 4;
                    $field_elements[$t] = $ret;
                    $element->{tmhmm} = \@field_elements;
                }
            }
            ## Set the nth annotation to its annotation.
            $annotations->[$element_number] = $element;
        } ## End iterating over every element in a multi-element field
    } ## End the foreach every key in the object
    $ids = Extract_Annotations($class, $element,
                               ids => $ids,
                               fh => $fh,
                               evalue => $args{evalue},
                               identity => $args{identity},
                           );
    return($filled_in);
}

sub Extract_Annotations {
    my ($class, $datum, %args) = @_;
    my $min_identity = $args{identity};
    my $eval_max = $args{evalue};
    my $ids = $args{ids};
    my $fh = $args{fh};
    my $options = $class->Get_Vars();

    my $id = $datum->{prot_id};
    $id = $datum->{transcript_id} if ($id eq '.');
    if (defined($options->{species})) {
        my $species = $options->{species};
        $species =~ s/\s+/_/g;
        $id =~ s/TRINITY/${species}/g;
    }
    $id =~ s/\:\:/; /g;

    if (defined($ids->{$id})) {
        ## Then this record has been processed, just increment it and move on.
        $ids->{$id}++;
    } elsif (!defined($datum->{swissprot} && !defined($datum->{rnammer}))) {
        $ids->{$id} = 1;
    } elsif (defined($datum->{rnammer})) {
        my $seq = $datum->{sequence};
        $seq = join("\n", ($seq =~ m/.{1,80}/g));
        $ids->{$id} = 1;
        my $header_string = qq"${id}; ${id}; undef; undef; $datum->{rnammer}->[0]->{name}, $datum->{rnammer}->[0]->{region}; undef; undef";
        print $fh qq">${header_string}
${seq}
";
    } elsif ($datum->{swissprot}->[0]->{identity} >= $min_identity &&
                 $datum->{swissprot}->[0]->{e_value} <= $eval_max) {
        $ids->{$id} = 1;
        my $seq = $datum->{sequence};
        $seq = join("\n", ($seq =~ m/.{1,80}/g));
        my $header_string = qq"${id}; $datum->{swissprot}->[0]->{fullname}; $datum->{swissprot}->[0]->{species}; $datum->{swissprot}->[0]->{e_value}; $datum->{swissprot}->[0]->{identity}";
        print $fh qq">${header_string}
${seq}
";
    }
    return($ids);
}

=back

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

    L<bowtie> L<bowtie2> L<tophat> L<bwa> L<kallisto> L<samtools> L<htseq>

=cut

1;
