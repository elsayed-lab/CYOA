package Bio::Adventure::Annotation;
## LICENSE: gplv2
## ABSTRACT:  Kitty!
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use File::Basename;
use Cwd qw(abs_path getcwd);
use File::Spec;
use File::Which qw"which";
use File::ShareDir qw":ALL";

=head1 NAME

Bio::Adventure::Annotation - Do some searches to help annotate genes.

=head1 SYNOPSIS

=cut

sub Aragorn {
    my ($class, %args) = @_;
    my $check = which('aragorn');
    die("Could not find aragorn in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        species => undef,
        arbitrary => ' -rp -fasta -w -m -t -mt ',
        );
    my $aragorn_args = $options->{arbitrary};
    my $job_basename = $class->Get_Job_Name();
    my %aragorn_jobs = ();
    my $aragorn_depends_on;
    my $output_dir = qq"outputs/aragorn";
    my $species_string = qq"";
    my $comment = qq!## This is a script to run aragorn.
!;
    my $jstring = qq!mkdir -p ${output_dir} && \\
  aragorn $options->{arbitrary} \\
    -o ${output_dir}/aragorn.txt \\
    $options->{input}
!;

    my $aragorn_job = $class->Submit(
        cpus => 6,
        comment => $comment,
        depends => $aragorn_depends_on,
        jname => "aragorn_${job_basename}",
        jprefix => "64",
        jstring => $jstring,
        mem => 24,
        output => qq"outputs/aragorn.sbatchout",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => "workstation",
    );
    my $jobs = {
        aragorn => $aragorn_job,
    };
    return($jobs);
}


sub Extend_Kraken_DB {
    my ($class, %args) = @_;
    my $check = which('kraken2');
    die("Could not find kraken2 in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        library => 'viral',
        );
    ## kraken2 --db ${DBNAME} --paired --classified-out cseqs#.fq seqs_1.fq seqs_2.fq
    my $job_basename = $class->Get_Job_Name();
    my %kraken_jobs = ();
    my $kraken_depends_on;
    my $output_dir = qq"outputs/extend_kraken";

    my $comment = qq!## This is a script to extend an existing kraken2 library with some new sequences.
!;
    my $jstring = qq!mkdir -p ${output_dir}
kraken2-build --download-taxonomy --db \${KRAKEN_DB_PATH}/$options->{library} \\
              2>${output_dir}/kraken2-build.out 1>&2
kraken2-build --add-to-library $options->{input} --db \${KRAKEN_DB_PATH}/$options->{library} \\
              2>>${output_dir}/kraken2-build.out 1>&2
kraken2-build --download-library $options->{library} \\
              --db \${KRAKEN_DB_PATH}/$options->{library} \\
              2>>${output_dir}/kraken2-build.out 1>&2
kraken2-build --build --db \${KRAKEN_DB_PATH}/$options->{library} \\
              2>>${output_dir}/kraken2-build.out 1>&2
!;
    my $kraken_job = $class->Submit(
        cpus => 6,
        comment => $comment,
        depends => $kraken_depends_on,
        jname => "kraken_${job_basename}",
        jprefix => "99",
        jstring => $jstring,
        mem => 96,
        output => qq"outputs/kraken_extend.sbatchout",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => "large",
        walltime => "144:00:00",
    );
    my $jobs = {
        kraken => $kraken_job,
    };
    return($jobs);
}


sub Glimmer {
    my ($class, %args) = @_;
    my $check = which('glimmer3');
    die("Could not find glimmer in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        );
    my $job_basename = $class->Get_Job_Name();
    my %glimmer_jobs = ();
    my $glimmer_depends_on;
    my $output_dir = qq"outputs/glimmer";

    my $comment = qq!## This is a script to run glimmer.
!;
##    my $jstring = qq!mkdir -p ${output_dir}
##long-orfs -n -t 1.15 $options->{input} ${output_dir}/first_run_longorfs.txt \\ 
##  2> ${output_dir}/first_run_longorfs.out 1>&2
##extract -t $options->{input} ${output_dir}/first_run_longorfs.txt > \\
##  ${output_dir}/first_run_training.txt
##build-icm -r ${output_dir}/first_run.icm < ${output_dir}/first_run_training.txt
##glimmer3 -o50 -g110 -t30 \\
##  $options->{input} \\
##  ${output_dir}/first_run.icm \\
##  ${output_dir}/first_run.out \\
##  2>${output_dir}first_run_glimmer3.out 1>&2
##
#### Use this first run output to go again
##
##extract -t $options->{input} train.coords > ${output_dir}/second_run_training.txt
##build-icm -r ${output_dir}/second_run.icm < ${output_dir}/second_run_training.txt
##upstream-coords.awk 25 0 train.coords | extract $options->{input} - > \\
##  ${output_dir}/second_run_upstream.txt
##elph ${output_dir}/second_run_upstream.txt LEN=6 | get-motif-counts.awk > \\
##  ${output_dir}/second_run_motif.txt
##startuse='start-codon-distrib -3 $options->{input} train.coords'
##glmmer3 -o50 -g110 -t30 -b ${output_dir}/second_run_motif.txt -P \${startuse} \\
##  $options->{input} ${output_dir}/second_run.icm \\
##  ${output_dir}/second_run.out
##!;
    my $jstring = qq!
cyoa_invoke_glimmer.pl --input $options->{input}
!;
    
    my $glimmer_job = $class->Submit(
        cpus => $options->{cpus},
        comment => $comment,
        depends => $glimmer_depends_on,
        jname => "glimmer_${job_basename}",
        jprefix => "65",
        jstring => $jstring,
        mem => 24,
        output => qq"outputs/glimmer.sbatchout",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => "workstation",
    );
    my $jobs = {
        glimmer => $glimmer_job,
    };
    return($jobs);
}


sub Interproscan {
    my ($class, %args) = @_;
    my $check = which('interproscan.sh');
    die("Could not find interproscan in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
    );
    my $job_basename = $class->Get_Job_Name();
    my $cwd_name = basename(getcwd());

    my $interproscan_exe_dir = dirname($check);
    my $input_path = abs_path($options->{input});
    my $input_dir = dirname($options->{input});
    my $input_name = basename(dirname($options->{input}));
    ## I am assuming I can get the assembly type from basename()
    my $dirname = basename($input_dir);
    my $output_dir = qq"outputs/interproscan_${input_name}";
    my $comment = qq!## This is a interproscan submission script
!;
    my $jstring = qq!mkdir -p ${output_dir}
start=\$(pwd)
cd ${output_dir}
interproscan.sh -i ${input_path}
cd \${start}
!;
    my $interproscan_job = $class->Submit(
        cpus => 6,
        comment => $comment,
        depends => $options->{depends},
        jname => "interproscan_${job_basename}",
        jprefix => "48",
        jstring => $jstring,
        mem => 80,
        output => qq"interproscan.sbatchout",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => "large",
        walltime => "144:00:00",
    );
    return($interproscan_job);
}

sub Kraken {
    my ($class, %args) = @_;
    my $check = which('kraken2');
    die("Could not find kraken2 in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        library => 'viral',
        );
    ## kraken2 --db ${DBNAME} --paired --classified-out cseqs#.fq seqs_1.fq seqs_2.fq
    my $job_basename = $class->Get_Job_Name();
    my %kraken_jobs = ();
    my $kraken_depends_on;
    my $output_dir = qq"outputs/kraken_${job_basename}";
    my $input_string = "";
    if ($options->{input} =~ /\:|\;|\,|\s+/) {
        my @in = split(/\:|\;|\,|\s+/, $options->{input});
        $input_string = qq" --paired $in[0]  $in[1] ";
    } else {
        $input_string = qq"$options->{input} ";
    }
    my $comment = qq!## This is a kraken2 submission script
!;
    my $jstring = qq!mkdir -p ${output_dir} && \\
  kraken2 --db $options->{library} \\
    --report ${output_dir}/kraken_report.txt --use-mpa-style \\
    --use-names ${input_string} \\
    --classified-out ${output_dir}/classified#.fastq.gz \\
    --unclassified-out ${output_dir}/unclassified#.fastq.gz \\
    2>${output_dir}/kraken.out 1>&2
!;
    my $kraken_job = $class->Submit(
        cpus => 6,
        comment => $comment,
        depends => $kraken_depends_on,
        jname => "kraken_${job_basename}",
        jprefix => "45",
        jstring => $jstring,
        mem => 96,
        output => qq"outputs/kraken.sbatchout",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => "large",
        walltime => "144:00:00",
    );
    my $jobs = {
        kraken => $kraken_job,
    };
    return($jobs);
}


sub Phageterm {
    my ($class, %args) = @_;
    my $check = which('PhageTerm.py');
    die("Could not find phageterm in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'library'],
        cpus => '8',
        );
    my $job_basename = $class->Get_Job_Name();
    my %phageterm_jobs = ();
    my $phageterm_depends_on;
    my $cwd_name = basename(getcwd());
    my $assembly_full = abs_path($options->{library});
    my $assembly_name = basename(dirname($options->{library}));
    my $output_dir = qq"outputs/phageterm_${assembly_name}";

    my $uncompress_string = qq"";
    my $input_string = qq"";
    my $delete_string = qq"";

    if ($options->{input} =~ /\:|\;|\,|\s+/) {
        my @in = split(/\:|\;|\,|\s+/, $options->{input});
        my $r1 = abs_path($in[0]);
        my $r2 = abs_path($in[1]);
        $uncompress_string = qq!
less ${r1} > r1.fastq && less ${r2} > r2.fastq
!;
        $input_string = qq! -f r1.fastq -p r2.fastq !;
        $delete_string = qq!rm r1.fastq && rm r2.fastq!;
    } else {
        my $r1 = abs_path($options->{input});
        $uncompress_string = qq!
less ${r1} > r1.fastq
!;
        $input_string = qq! -f r1.fastq !;
        $delete_string = qq!rm r1.fastq!;
    }
    
    my $comment = qq!## This is a script to run phageterm.
!;
    my $jstring = qq!start=\$(pwd)
mkdir -p ${output_dir}
cd ${output_dir}
${uncompress_string}
PhageTerm.py ${input_string} \\
  -r ${assembly_full} \\
  -c $options->{cpus} \\
  --report_title ${cwd_name}
${delete_string}
if [[ \! -f "${cwd_name}_direct-term-repeats.fasta" ]]; then 
  ## phageterm did not find a terminal repeat
  rm ${cwd_name}_sequence.fasta
  ln -s ${assembly_full} ${cwd_name}_sequence.fasta
fi
cd \${start}
!;   
    my $phageterm_job = $class->Submit(
        cpus => $options->{cpus},
        comment => $comment,
        depends => $phageterm_depends_on,
        jname => "phageterm_${job_basename}",
        jprefix => "64",
        jstring => $jstring,
        mem => 24,
        output => qq"outputs/phageterm.sbatchout",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => "workstation",
    );
    my $jobs = {
        phageterm => $phageterm_job,
    };
    return($jobs);
}


sub Prodigal {
    my ($class, %args) = @_;
    my $check = which('prodigal');
    die("Could not find prodigal in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        gcode => '11',
        );
    my $job_basename = $class->Get_Job_Name();
    my %prodigal_jobs = ();
    my $prodigal_depends_on;

    my $library_file = qq"$options->{libdir}/hmm/$options->{species}_gc$options->{gcode}.training";
    unless (-r $library_file) {
        die("Need a training hmm for this species and genetic code: ${library_file}.");
    }
    my $output_dir = qq"outputs/prodigal_${input_name}";
    my $comment = qq!## This is a script to run prodigal.
!;
    my $in_name = basename($options->{input}, ('.fasta'));
    my $jstring = qq!mkdir -p ${output_dir}
prodigal -i $options->{input} \\
  -a ${output_dir}/${in_name}_translated.fasta \\
  -d ${output_dir}/${in_name}_cds.fasta \\
  -s ${output_dir}/${in_name}_scores.txt \\
  -t ${library_file} \\
  2>${output_dir}/prodigal.err \\
  1>${output_dir}/prodigal.out
!;
    my $prodigal_job = $class->Submit(
        cpus => 1,
        comment => $comment,
        depends => $prodigal_depends_on,
        jname => "prodigal_${job_basename}",
        jprefix => "63",
        jstring => $jstring,
        mem => 24,
        output => qq"outputs/prodigal.sbatchout",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => "workstation",
    );
    my $jobs = {
        prodigal => $prodigal_job,
    };
    return($jobs);
}


sub Train_Prodigal {
    my ($class, %args) = @_;
    my $check = which('prodigal');
    die("Could not find prodigal in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        gcode => '11',
        );
    my $job_basename = $class->Get_Job_Name();
    my %prodigal_jobs = ();
    my $prodigal_depends_on;
    my $kingdom_string = '';
    my $output = qq"$options->{libdir}/hmm/$options->{species}_gc$options->{gcode}.training";
    my $comment = qq!## This is a script to train prodigal.
!;
    my $jstring = qq!mkdir -p ${output_dir}
prodigal -i $options->{input} \\
  -t ${output} \\
  2>${output_dir}/prodigal_training.stderr \\
  1>${output_dir}/prodigal_training.stdout
!;
    my $prodigal_job = $class->Submit(
        cpus => 1,
        comment => $comment,
        depends => $prodigal_depends_on,
        jname => "prodigal_training_${job_basename}",
        jprefix => "63",
        jstring => $jstring,
        mem => 24,
        output => qq"outputs/prodigal.sbatchout",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => "workstation",
    );
    my $jobs = {
        prodigal => $prodigal_job,
    };
    return($jobs);
}

sub Prokka {
    my ($class, %args) = @_;
    my $check = which('prokka');
    die("Could not find prokka in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input',],
        kingdom => '',
        gcode => '',
        arbitrary => ' --addgenes --force --compliant --genus unknown --species virus --strain phage ',
        );
    my $job_basename = $class->Get_Job_Name();
    my %prokka_jobs = ();
    my $prokka_depends_on;
    my $kingdom_string = '';
    if ($options->{kingdom} ne '') {
        $kingdom_string = qq" --kingdom $options->{kingdom} --gcode $options->{gcode} ";
    }
    
    my $cwd_name = basename(getcwd());
    my $input_name = basename(dirname($options->{input}));
    my $output_dir = qq"outputs/prokka_${input_name}";
    my $comment = qq!## This is a script to run prokka.
!;
    my $jstring = qq!mkdir -p ${output_dir}
prokka $options->{arbitrary} \\
  ${kingdom_string} \\
  --outdir ${output_dir} \\
  --prefix ${cwd_name} \\
  $options->{input} \\
  2>${output_dir}/prokka.stderr \\
  1>${output_dir}/prokka.stdout
!;
    my $prokka_job = $class->Submit(
        cpus => 6,
        comment => $comment,
        depends => $prokka_depends_on,
        jname => "prokka_${job_basename}",
        jprefix => "63",
        jstring => $jstring,
        mem => 24,
        output => qq"outputs/prokka.sbatchout",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => "workstation",
    );
    my $jobs = {
        prokka => $prokka_job,
    };
    return($jobs);
}

sub Resfinder {
    my ($class, %args) = @_;
    my $check = which('run_resfinder.py');
    die("Could not find resfinder in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        species => undef,
        arbitrary => ' -l 0.6 -t 0.8 --acquired ',
        );
    my $resfinder_args = $options->{arbitrary};
    my $job_basename = $class->Get_Job_Name();
    my %resfinder_jobs = ();
    my $resfinder_depends_on;
    my $assembly_name = basename(dirname($options->{input}));
    my $output_dir = qq"outputs/resfinder_${assembly_name}";
    my $species_string = qq"";
    if (defined($options->{species})) {
        $species_string = qq" --species $options->{species} ";
    }

    my $comment = qq!## This is a script to run resfinder.
!;
    my $jstring = qq!mkdir -p ${output_dir} && \\
  run_resfinder.py -ifa $options->{input} \\
    -o ${output_dir} \\
    ${resfinder_args} ${species_string}
!;
    my $resfinder_job = $class->Submit(
        cpus => 6,
        comment => $comment,
        depends => $resfinder_depends_on,
        jname => "resfinder_${job_basename}",
        jprefix => "63",
        jstring => $jstring,
        mem => 24,
        output => qq"outputs/resfinder.sbatchout",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => "workstation",
    );
    my $jobs = {
        resfinder => $resfinder_job,
    };
    return($jobs);
}


sub tRNAScan {
    my ($class, %args) = @_;
    my $check = which('trnascan');
    die("Could not find trnascan in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        species => undef,
        arbitrary => ' -G ',
        );
    my $trnascan_args = $options->{arbitrary};
    my $job_basename = $class->Get_Job_Name();
    my %trnascan_jobs = ();
    my $trnascan_depends_on;
    my $output_dir = qq"outputs/trnascan";
    my $species_string = qq"";
    my $comment = qq!## This is a script to run trnascan.
!;
    my $jstring = qq!mkdir -p ${output_dir} && \\
  tRNAscan-SE $options->{arbitrary} \\
    -o ${output_dir}/trnascan.txt \\
    $options->{input}
!;

    my $trnascan_job = $class->Submit(
        cpus => 6,
        comment => $comment,
        depends => $trnascan_depends_on,
        jname => "trnascan_${job_basename}",
        jprefix => "64",
        jstring => $jstring,
        mem => 24,
        output => qq"outputs/trnascan.sbatchout",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => "workstation",
    );
    my $jobs = {
        trnascan => $trnascan_job,
    };
    return($jobs);
}

=head2 C<Extract_Annotations>

Used by Extract_Trinotate to extract annotations from the strangely encoded trinotate outputs.

=over

=item I<evalue> Evalue cutoff for trinotate output.

=item I<identity> Minimum percent identity cutoff for trinotate output.

=item I<ids> IDs to extract.

=item I<fh> Filehandle to parse.

=back

=cut
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

=head2 C<Extract_Trinotate>

The trinotate output format is a bit... unwieldy.  This seeks to parse out the
useful information from it.

=over

=item I<input> * Input csv file from trinotate.

=item I<output> (interesting.fasta) Output for the parsed csv.

=item I<evalue> (1e-10) Evalue cutoff for trinotate output.

=item I<identity> (70) Minimum percent identity cutoff for trinotate output.

=back

=head3 C<Invocation>

> cyoa --task assembly --method extract --input trinotate_output.csv

=cut
sub Extract_Trinotate {
    my ($class, %args) = @_; my
        $options = $class->Get_Vars(args => \%args,
                                    required => ['input'],
                                    jname => "trin_rsem",
                                    output => 'interesting.fasta',
                                    evalue => 1e-10,
                                    identity => 70, );
    my $job_basename = $class->Get_Job_Name();
    my $trinity_out_dir =
        qq"outputs/trinity_${job_basename}";


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


=head2 C<Read_Write_Annotation>

Called by Extract_Trinotate() to help write out the trinotate csv information
into an easier-to-read format.

=cut
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

=head2 C<Transdecoder>

$hpgl->Transdecoder() submits a trinity denovo sequence assembly and runs its
default post-processing tools.

=over

=item I<input> * Output from trinity for post processing.

=back

=head3 C<Invocation>

> cyoa --task assembly --method transdecoder --input trinity.fasta

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

=head2 C<Trinotate>

Submit a trinity denovo sequence assembly to trinotate.

=over

=item I<input> * Input fasta from trinity.

=back

=head3 C<Invocation>

> cyoa --task assembly --method trinotate --input trinity.fasta

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
    my $cwd_name = basename(getcwd());

    my $trinotate_exe_dir = dirname($check);
    my $input_path = abs_path($options->{input});
    my $input_dir = dirname($options->{input});
    my $input_name = basename(dirname($options->{input}));
    my $input_basename = basename($options->{input});
    ## I am assuming I can get the assembly type from basename()
    my $dirname = basename($input_dir);
    my $output_dir = qq"outputs/trinotate_${input_name}";
    my $comment = qq!## This is a trinotate submission script
!;
    my $jstring = qq!mkdir -p ${output_dir}
start=\$(pwd)
cd ${output_dir}
ln -s ${input_path} .
if [ -f ${input_path}.gene_trans_map ]; then
  ln -s ${input_path}.gene_trans_map .
else
  ids=\$(grep "^>" ${input_path} | sed 's/>//g' | awk '{print \$1}')
  rm -f ${input_name}.gene_trans_map
  for i in \${ids}; do
    echo "\${i}	\${i}" >> ${input_name}.gene_trans_map
  done
fi
${trinotate_exe_dir}/auto/autoTrinotate.pl \\
  --Trinotate_sqlite ${trinotate_exe_dir}/sample_data/Trinotate.boilerplate.sqlite \\
  --transcripts ${input_basename} \\
  --gene_to_trans_map ${input_name}.gene_trans_map \\
  --conf ${trinotate_exe_dir}/auto/conf.txt \\
  --CPU 6 \\
  2>trinotate_${job_basename}.err \\
  1>trinotate_${job_basename}.out
mv Trinotate.xls Trinotate.tsv
cd \${start}
!;
    my $trinotate_job = $class->Submit(
        cpus => 6,
        comment => $comment,
        depends => $options->{depends},
        jname => "trinotate_${input_name}_${job_basename}",
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

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

=cut

1;
