package Bio::Adventure::Annotation;
## LICENSE: gplv2
## ABSTRACT:  Kitty!
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use Bio::SeqIO;
use Bio::Seq;
use Bio::SeqFeature::Generic;
use Bio::Tools::GuessSeqFormat;
use Cwd qw"abs_path getcwd cwd";
use Data::Table;
use Data::Table::Excel qw"tables2xlsx";
use File::Basename;
use File::Copy qw"cp";
use File::Spec;
use File::Path qw"make_path";
use File::Which qw"which";
use File::ShareDir qw":ALL";
use List::MoreUtils qw"any";
use Template;
use Text::CSV_XS::TSV;

=head2 C<Aragorn>

Use aragorn to search for tRNA genes in a sequence database.  By
default this will explictly search for tmRNA as well.

=cut
sub Aragorn {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        modules => ['aragorn'],
        jprefix => 21,
        species => undef,
        arbitrary => ' -rp -fasta -w -m -t -mt ',);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('aragorn');
    die("Could not find aragorn in your PATH.") unless($check);
    my $aragorn_args = $options->{arbitrary};
    my $job_name = $class->Get_Job_Name();
    my $output_dir = qq"outputs/$options->{jprefix}aragorn";
    my $species_string = qq"";
    my $comment = qq!## This is a script to run aragorn.
!;
    my $jstring = qq!mkdir -p ${output_dir} && \\
  aragorn $options->{arbitrary} \\
    -o ${output_dir}/aragorn.txt \\
    $options->{input}
!;

    my $aragorn = $class->Submit(
        cpus => 6,
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => "aragorn_${job_name}",
        jprefix => "64",
        jstring => $jstring,
        jmem => 24,
        modules => $options->{modules},
        output => qq"${output_dir}/aragorn.txt",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => 'workstation',);
    $loaded = $class->Module_Loader(modules => $options->{modules},
        action => 'unload',);
    return($aragorn);
}

=head2 C<Glimmer>

Use glimmer in two passes to search for ORFs in a sequence database.

=cut
sub Glimmer {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        modules => ['glimmer'],
        jprefix => '16',);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('glimmer3');
    die("Could not find glimmer in your PATH.") unless($check);
    my $job_name = $class->Get_Job_Name();
    my $output_dir = qq"outputs/$options->{jprefix}glimmer";

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
cyoa_invoke_glimmer.pl --input $options->{input} --jprefix $options->{jprefix}
!;

    ## FIXME: There are a bunch of potentially useful glimmer outputs which should be put here.

    my $glimmer = $class->Submit(
        cpus => $options->{cpus},
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => "glimmer_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => 24,
        modules => $options->{modules},
        output => qq"${output_dir}/${job_name}_glimmer.out",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => "workstation",);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($glimmer);
}

=head2 C<Glimmer_Single

Run glimmer in a single pass instead of the two passes as per the paper.

=cut
sub Glimmer_Single {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        modules => ['glimmer'],
        cutoff => 1.1,
        overlap => 20,
        minlength => 45,
        threshold => 30,
        jprefix => '16',);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('glimmer3');
    die("Could not find glimmer in your PATH.") unless($check);
    my $job_name = $class->Get_Job_Name();
    my $output_dir = qq"outputs/$options->{jprefix}glimmer";

    my $comment = qq!## This is a script to run glimmer.!;
    my $final_output = qq"${output_dir}/glimmer3.txt";
    my $final_error = qq"${output_dir}/glimmer3.err";
    my $jstring = qq!mkdir -p ${output_dir}
long-orfs -n -t $options->{cutoff} $options->{input} ${output_dir}/longorfs.txt \\
  2>${output_dir}/longorfs.err \\
  1>${output_dir}/longorfs.out
extract -t $options->{input} ${output_dir}/longorfs.txt \\
  1>${output_dir}/training.txt \\
  2>${output_dir}/training.err
build-icm -r ${output_dir}/single_run.icm < ${output_dir}/training.txt
glimmer3 -o$options->{overlap} -g$options->{minlength} -t$options->{threshold} \\
  $options->{input} \\
  ${output_dir}/single_run.icm \\
  ${output_dir}/glimmer3 \\
  2>${final_error} \\
  1>${final_output}
!;
    my $glimmer = $class->Submit(
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => "glimmer_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        modules => $options->{modules},
        output => qq"${output_dir}/glimmer3.predict",
        output_detail => qq"${output_dir}/glimmer3.detail",
        output_icm => qq"${output_dir}/single_run.icm",
        output_longorfs => qq"${output_dir}/longorfs.txt",
        output_training => qq"${output_dir}/training.txt",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => "workstation",);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($glimmer);
}

=head2 C<Interproscan>

Use interproscan to look for existing annotations which are similar to
a set of provided ORFs.

=cut
sub Interproscan {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jprefix => '21',
        modules => ['interproscan'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('interproscan.sh');
    die("Could not find interproscan in your PATH.") unless($check);

    my $job_name = $class->Get_Job_Name();
    my $inputs = $class->Get_Paths($options->{input});
    my $cwd_name = basename(cwd());

    my $interproscan_exe_dir = dirname($check);
    ## Hey, don't forget abs_path requires a file which already exists.
    my $input_filename = basename($options->{input});
    my $output_filename = qq"${input_filename}.tsv";
    my $input_dir = dirname($options->{input});
    my $input_dirname = basename($input_dir);
    my $input_path = abs_path($input_dir);
    $input_path = qq"${input_path}/${input_filename}";
    my $output_dir = qq"outputs/$options->{jprefix}interproscan_${input_dirname}";
    my $comment = qq!## This is a interproscan submission script
!;
    my $jstring = qq!mkdir -p ${output_dir}
start=\$(pwd)
cd ${output_dir}
interproscan.sh -i ${input_path} 2>interproscan.err \\
  1>interproscan.out
ln -sf ${output_filename} interproscan.tsv
cd \${start}
!;
    my $interproscan = $class->Submit(
        cpus => 8,
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => "interproscan_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => 80,
        modules => $options->{modules},
        output => qq"${output_dir}/interproscan.tsv",
        output_gff => qq"${output_dir}/${input_dirname}.faa.gff3",
        output_tsv => qq"${output_dir}/interproscan.tsv",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => "large",
        walltime => "144:00:00",);

    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($interproscan);
}

=head2 C<Kraken>

Use kraken2 to taxonomically classify reads.

=cut
sub Kraken {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        library => 'viral',
        jprefix => '11',
        modules => ['kraken'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('kraken2');
    die("Could not find kraken2 in your PATH.") unless($check);
    ## kraken2 --db ${DBNAME} --paired --classified-out cseqs#.fq seqs_1.fq seqs_2.fq
    my $job_name = $class->Get_Job_Name();
    my $input_directory = basename(cwd());
    my $output_dir = qq"outputs/$options->{jprefix}kraken_$options->{library}";
    make_path($output_dir);
    my $input_string = "";
    if ($options->{input} =~ /\:|\;|\,|\s+/) {
        my @in = split(/\:|\;|\,|\s+/, $options->{input});
        $input_string = qq" --paired <(less $in[0]) <(less $in[1]) ";
    } else {
        $input_string = qq"<(less $options->{input}) ";
    }
    my $comment = qq!## This is a kraken2 submission script
!;
    my $jstring = qq!kraken2 --db $ENV{KRAKEN2_DB_PATH}/$options->{library} \\
  --report ${output_dir}/kraken_report.txt --use-mpa-style \\
  --use-names ${input_string} \\
  --classified-out ${output_dir}/classified#.fastq.gz \\
  --unclassified-out ${output_dir}/unclassified#.fastq.gz \\
  2>${output_dir}/kraken.out 1>&2
!;
    my $kraken = $class->Submit(
        cpus => 6,
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => "kraken_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => 96,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => 'large',
        modules => $options->{modules},
        output => qq"${output_dir}/kraken_report.txt",);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($kraken);
}

=head2 C<Prodigal>

Invoke prodigal on an assembly to search for ORFs.  This will by
default look for an existing training file provided by
Train_Prodigal()'.

=cut
sub Prodigal {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        species => undef,
        gcode => '11',
        output_dir => undef,
        jprefix => '17',
        modules => ['prodigal'],
        prodigal_outname => undef,);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('prodigal');
    die("Could not find prodigal in your PATH.") unless($check);

    my $inputs = $class->Get_Paths($options->{input});
    my $job_name = $class->Get_Job_Name();
    $job_name = basename($job_name, ('.fsa'));
    my $train_string = '';
    my $library_file;
    if ($options->{species}) {
        $library_file = qq"$options->{libdir}/hmm/$options->{species}_gc$options->{gcode}.training";
        if (!-r $library_file) {
            print "Could not find the training file for this species.\n";
            print "Sleeping for a moment so that you can do something.\n";
            sleep(5);
            $train_string = '';
        }
        $train_string = qq" -t ${library_file} ";
    }

    my $in_name = basename($options->{input}, ('.fasta'));
    my $output_dir;
    if (defined($options->{output_dir})) {
        if ($options->{output_dir}) {
            $output_dir = $options->{output_dir};
        }
    } else {
        $output_dir = qq"outputs/$options->{jprefix}prodigal_${in_name}";
    }

    my ($cds_file, $translated_file, $scores_file, $gff_file, $gbk_file);
    if ($options->{prodigal_outname}) {
        $cds_file = qq"${output_dir}/$options->{prodigal_outname}_cds.fasta";
        $translated_file = qq"${output_dir}/$options->{prodigal_outname}_translated.fasta";
        $scores_file = qq"${output_dir}/$options->{prodigal_outname}_scores.txt";
        $gff_file = qq"${output_dir}/$options->{prodigal_outname}.gff";
        $gbk_file = qq"${output_dir}/$options->{prodigal_outname}.gb";
    } else {
        $cds_file = qq"${output_dir}/predicted_cds.fasta";
        $translated_file = qq"${output_dir}/predicted_translated.fasta";
        $scores_file = qq"${output_dir}/predicted_scores.txt";
        $gff_file = qq"${output_dir}/predicted_cds.gff";
        $gbk_file = qq"${output_dir}/predicted_cds.gb";
    }

    my $comment = qq!## This is a script to run prodigal.
!;
    my $jstring = qq!mkdir -p ${output_dir}
prodigal ${train_string} \\
  -i $options->{input} \\
  -a ${translated_file} \\
  -d ${cds_file} \\
  -s ${scores_file} \\
  -f gff -o ${gff_file} \\
  2>${output_dir}/prodigal_gff.err \\
  1>${output_dir}/prodigal_gff.out
prodigal ${train_string} \\
  -i $options->{input} \\
  -f gbk -o ${gbk_file} \\
  2>${output_dir}/prodigal_gbk.err \\
  1>${output_dir}/prodigal_gbk.out
!;
    my $prodigal = $class->Submit(
        cpus => 1,
        comment => $comment,
        jdepends => $options->{jdepends},
        jmem => 24,
        jname => "prodigal_${job_name}",
        jprefix => $options->{jprefix},
        jqueue => "workstation",
        jstring => $jstring,
        modules => $options->{modules},
        output => $gbk_file,
        output_cds => $cds_file,
        output_gff => $gff_file,
        output_scores => $scores_file,
        output_translated => $translated_file,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        training_input => $library_file,);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($prodigal);
}

=head2 C<Train_Prodigal>

Some assemblies I have been performing are on sets of sequence which
are too small for prodigal to train itself sufficiently; so this
function was written to provide an opportunity for one to collate a
larger sequence database for training.

=cut
sub Train_Prodigal {
    my ($class, %args) = @_;
    my $check = which('prodigal');
    die("Could not find prodigal in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        gcode => '11',
        modules => ['prodigal'],);
    my $job_name = $class->Get_Job_Name();
    my $kingdom_string = '';
    my $output_dir = qq"$options->{libdir}/hmm";
    my $output = qq"${output_dir}/$options->{species}_gc$options->{gcode}.training";
    my $comment = qq!## This is a script to train prodigal.
!;
    my $jstring = qq!mkdir -p ${output_dir}
prodigal -i $options->{input} \\
  -t ${output} \\
  2>${output_dir}/prodigal_training.stderr \\
  1>${output_dir}/prodigal_training.stdout
!;
    my $prodigal = $class->Submit(
        cpus => 1,
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => "prodigal_training_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => 24,
        modules => $options->{modules},
        output => $output,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => "workstation",);
    return($prodigal);
}

=head2 C<Prokka>

Prokka is an automagic annotation tool for bacterial assemblies.  It
seems useful for other relatively small genomes.

=cut
sub Prokka {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        arbitrary => '',
        coverage => 30,
        evalue => '1e-05',
        gcode => '11',
        genus => 'phage',
        jprefix => '19',
        kingdom => 'bacteria',
        locus_tag => 'unknownphage',
        modules => ['prokka'],
        species => 'virus',);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('prokka');
    die("Could not find prokka in your PATH.") unless($check);

    my $job_name = $class->Get_Job_Name();
    my $input_paths = $class->Get_Paths($options->{input});
    my $kingdom_string = '';
    if ($options->{kingdom} ne '') {
        $kingdom_string = qq" --kingdom $options->{kingdom} --gcode $options->{gcode} ";
    }

    my $training_file = qq"$options->{libdir}/hmm/$options->{species}_gc$options->{gcode}.training";
    my $cwd_name = basename(cwd());
    ## If this is being run a part of a pipeline, the basename(dirname(input)) will be something
    ## useful/interesting.  If not, it will not likely be of use but instead be '.'
    ## Let us check for that and set $input_name to something useful in that case.
    my $input_name;
    if (defined($input_paths->{dirname})) {
        $input_name = $input_paths->{dirname};
    } else {
        $input_name = $input_paths->{filebase_extension};
    }
    my $locus_tag;
    if ($options->{locus_tag}) {
        $locus_tag = $options->{locus_tag};
    } else {
        $locus_tag = basename($input_name, ('.fasta'));
    }
    my $output_dir = qq"outputs/$options->{jprefix}prokka_${input_name}";
    my $comment = qq!## This is a script to run prokka.
!;

##    prokka --addgenes -rfam --force --locustag EAb03 --genus Phage --species virus --cdsrnaolap --usegenus --kingdom Bacteria --gcode 11 --prodigaltf /home/trey/libraries/hmm/abaumannii_phage_samples_gc11.training --outdir outputs/19prokka_15termreorder_13unicycler --prefix EAb03 outputs/15termreorder_13unicycler/final_assembly_reordered.fasta --evalue '1e-05' --coverage 30 2>&1 | less

    my $jstring = qq!mkdir -p ${output_dir}
prokka --addgenes --rfam --force ${kingdom_string} \\
  --locustag ${locus_tag} --genus $options->{genus} \\
  --compliant --cdsrnaolap --usegenus \\
  --prodigaltf ${training_file} \\
  --outdir ${output_dir} \\
  --prefix ${cwd_name} \\
  $options->{input} \\
  2>${output_dir}/prokka.stderr \\
  1>${output_dir}/prokka.stdout
!;

    my $error_file =  qq"${output_dir}/${cwd_name}.err";
    my $peptide_file = qq"${output_dir}/${cwd_name}.faa";
    my $cds_file = qq"${output_dir}/${cwd_name}.ffn";
    my $assembly_copy = qq"${output_dir}/${cwd_name}.fna";
    my $assembly_renamed_contigs = qq"${output_dir}/${cwd_name}.fsa";
    my $genbank_file = qq"${output_dir}/${cwd_name}.gbk";
    my $gff_file = qq"${output_dir}/${cwd_name}.gff";
    my $sqn_file = qq"${output_dir}/${cwd_name}.sqn";
    my $tbl_file = qq"${output_dir}/${cwd_name}.tbl";
    my $tsv_file = qq"${output_dir}/${cwd_name}.tsv";

    my $prokka = $class->Submit(
        cpus => 6,
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => "prokka_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => 24,
        modules => $options->{modules},
        output => $cds_file,
        output_error => $error_file,
        output_peptide => $peptide_file,
        output_cds => $cds_file,
        output_assembly => $assembly_copy,
        output_assemblyrenamed => $assembly_renamed_contigs,
        output_fsa => $assembly_renamed_contigs,
        output_genbank => $genbank_file,
        output_gff => $gff_file,
        output_sqn => $sqn_file,
        output_tbl => $tbl_file,
        output_tsv => $tsv_file,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => "workstation",);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($prokka);
}

=head2 C<tRNAScan>

Alternative to aragorn.  Search for tRNAs!

=cut
sub tRNAScan {
    my ($class, %args) = @_;
    my $check = which('trnascan');
    die("Could not find trnascan in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        modules => ['trnascan'],
        species => undef,
        arbitrary => ' -G ',);
    my $trnascan_args = $options->{arbitrary};
    my $job_name = $class->Get_Job_Name();
    my $output_dir = qq"outputs/trnascan";
    my $species_string = qq"";
    my $comment = qq!## This is a script to run trnascan.
!;
    my $jstring = qq!mkdir -p ${output_dir} && \\
  tRNAscan-SE $options->{arbitrary} \\
    -o ${output_dir}/trnascan.txt \\
    $options->{input}
!;

    my $trnascan = $class->Submit(
        cpus => 6,
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => "trnascan_${job_name}",
        jprefix => "64",
        jstring => $jstring,
        jmem => 24,
        modules => $options->{modules},
        output => qq"${output_dir}/trnascan.txt",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => "workstation",);
    return($trnascan);
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
sub Trinotate_Extract_Annotations {
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
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jname => "trin_rsem",
        output => 'interesting.fasta',
        evalue => 1e-10,
        identity => 70,);
    my $job_name = $class->Get_Job_Name();
    my $trinity_out_dir = qq"outputs/trinity_${job_name}";
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
    return($count);
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
    $ids = Trinotate_Extract_Annotations($class, $element,
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
        modules => ['transdecoder'],);
    my $transdecoder_input = File::Spec->rel2abs($options->{input});
    my $job_name = $class->Get_Job_Name();
    my $output_dir = qq"outputs/trinity_${job_name}";
    my $comment = qq!## This is a transdecoder submission script
!;
    my $jstring = qq!mkdir -p ${output_dir} && cd ${output_dir} && \\
  TransDecoder.LongOrfs -t ${transdecoder_input} \\
    2>${output_dir}/transdecoder_longorfs_${job_name}.err \\
    1>${output_dir}/transdecoder_longorfs_${job_name}.out
TransDecoder.Predict -t ${transdecoder_input} \\
  2>${output_dir}/transdecoder_predict_${job_name}.err \\
  1>${output_dir}/transdecoder_predict_${job_name}.out
${transdecoder_exe_dir}/util/cdna_alignment_orf_to_genome_orf.pl \\
  ${output_dir}/transcripts.fasta.transdecoder.gff3 \\
  ${output_dir}/transcripts.gff3 \\
  ${transdecoder_input} \\
  2>${output_dir}/cdna_alignment_${job_name}.err \\
  1>${output_dir}/transcripts.fasta.transdecoder.genome.gff3
!;
    my $transdecoder = $class->Submit(
        comment => $comment,
        jname => "transdecoder_${job_name}",
        jprefix => "47",
        jstring => $jstring,
        modules => $options->{modules},
        output => qq"${output_dir}/transcripts.fasta.transdecoder.genome.gff3",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => "workstation",);
    return($transdecoder);
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
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jprefix => '20',
        modules => ['divsufsort', 'transdecoder', 'blast', 'blastdb', 'signalp', 'hmmer',
                    'tmhmm', 'rnammer', 'trinotate', ],
        trinotate => 'autoTrinotate.pl',
        config => 'conf.txt',
        );

    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('Trinotate');
    die("Could not find Trinotate in your PATH:
 $ENV{PATH}.") unless($check);

    my $job_name = $class->Get_Job_Name();
    my $cwd_name = basename(cwd());
    my $trinotate_exe_dir = dirname($check);
    ## Once again, abs_path only works on stuff which already exists.
    ## So create the output directory, and use that.
    my $input_paths = $class->Get_Paths($options->{input});
    my $input_full = $input_paths->{fullpath};
    my $output_name = basename($input_full, ('.fasta', '.fa', '.fna', '.fsa'));
    $output_name = qq"${output_name}.tsv";
    my $output_dir = qq"outputs/$options->{jprefix}trinotate_$input_paths->{dirname}";
    my $comment = qq!## This is a trinotate submission script
!;
    my $jstring = qq!mkdir -p ${output_dir}
start=\$(pwd)
cd ${output_dir}
ln -sf $input_paths->{fullpath} .
if [ -f $input_paths->{filename}.gene_trans_map ]; then
  ln -sf $input_paths->{filename}.gene_trans_map .
else
  ids=\$(grep "^>" $input_paths->{fullpath} | sed 's/>//g' | awk '{print \$1}')
  rm -f $input_paths->{filename}.gene_trans_map
  for i in \${ids}; do
    echo "\${i}	\${i}" >> $input_paths->{filename}.gene_trans_map
  done
fi

${trinotate_exe_dir}/auto/$options->{trinotate} \\
  --conf ${trinotate_exe_dir}/auto/$options->{config} \\
  --Trinotate_sqlite ${trinotate_exe_dir}/sample_data/Trinotate.boilerplate.sqlite \\
  --transcripts $input_paths->{filename} \\
  --gene_to_trans_map $input_paths->{filename}.gene_trans_map \\
  --CPU 6 \\
  2>trinotate_${job_name}.err \\
  1>trinotate_${job_name}.out
mv Trinotate.tsv ${output_name}
cd \${start}
!;
    my $trinotate = $class->Submit(
        cpus => 6,
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => qq"trinotate_$input_paths->{filename}_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => 12,
        modules => $options->{modules},
        output => qq"${output_dir}/${output_name}",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => 'large',
        jwalltime => '144:00:00',
        );
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($trinotate);
}

=head2 C<Watson_Plus>

Use prodigal to count up putative ORFs on the plus and minus strands.
If the number of minus strand ORFs is larger than plus, reverse
complement the sequence.

This function just puts the actual function which does the work onto
the dependency chain.

=cut
sub Watson_Plus {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jprefix => '44',);

    my $input_seq = $options->{input};
    my $job_name = 'watsonplus';
    my $output_dir = qq"outputs/$options->{jprefix}${job_name}";
    make_path($output_dir);
    my $watson_prodigal = $class->Bio::Adventure::Annotation::Prodigal(
        gcode => '11',
        input => $options->{input},
        jdepends => $options->{jdepends},
        jname => $job_name,
        jprefix => $options->{jprefix},
        modules => ['prodigal'],
        output_dir => $output_dir,
        species => 'phages',);
    $options->{jdepends} = $watson_prodigal->{job_id};
    my $input_gff = $watson_prodigal->{output_gff};
    my $output_file = basename($options->{input});
    $output_file = qq"${output_dir}/${output_file}";
    my $comment_string = qq"## This takes the prodigal output and checks to see which strand has
## more ORFs, if it is the crick strand, then the chromsome is reverse complemented.
";
    my $jstring = qq!
use Bio::Adventure::Annotation;
\$result = Bio::Adventure::Annotation::Watson_Rewrite(\$h,
  gff => '${input_gff}',
  input => '$options->{input}',
  jdepends => '$options->{jdepends}',
  jname => '$options->{jname}',
  jprefix => '$options->{jprefix}',
  output => '${output_file}',
  output_dir => '${output_dir}',);
!;
    my $rewrite = $class->Submit(
        comment => $comment_string,
        gff => $options->{gff},
        input => $options->{input},
        jdepends => $options->{jdepends},
        jname => $options->{jname},
        jprefix => $options->{jprefix},
        jstring => $jstring,
        language => 'perl',
        output => $output_file,
        output_dir => $output_dir,
        shell => '/usr/bin/env perl',);
    return($rewrite);
}

=head2 C<Watson_Rewrite>

Does the actual work of rewriting a genome to put the majority ORFs on
the plus strand.  Watson_Plus() calls this.

=cut
sub Watson_Rewrite {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],);

    my $input_gff = $options->{gff};
    my $input_seq = $options->{input};
    my $job_name = 'watsonplus';
    my $output_dir = qq"outputs/$options->{jprefix}${job_name}";
    make_path($output_dir);

    my $orf_count = {};
    my $annotation_in = FileHandle->new("<${input_gff}");
    my $write_file = basename($options->{input});
    $write_file = qq"${output_dir}/${write_file}";
    my $log_file = basename($write_file, ('.fasta'));
    $log_file = qq"${output_dir}/${log_file}.log";
    my $read_fasta = Bio::SeqIO->new(-file => qq"<$options->{input}", -format => 'Fasta');
    my $write_fasta = Bio::SeqIO->new(-file => qq">${write_file}", -format => 'Fasta');
    my $write_log = FileHandle->new(">${log_file}");
    print $write_log "Counting ORFs on the current watson and crick strand.
If the crick strand is larger, flipping them.
";
    ## A prodigal line looks like:
    ## EAb06   Prodigal_v2.6.3 CDS     3205    3414    8.7     -       0       ID=1_10;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.248;conf=88.17;score=8.74;cscore=-0.38;sscore=9.12;rscore=8.34;uscore=-0.93;tscore=2.86;
    while (my $line = <$annotation_in>) {
        next if ($line =~ /^#/);
        my ($name, $prod, $type, $start, $end, $score, $strand, $phase, $annot) = split(/\t/, $line);
        my @annot_list = split(/;/, $annot);
        my $id_annot = $annot_list[0];
        my ($contig, $orf);
        if ($id_annot =~ /ID=(.*)?_(\d+)$/) {
            $contig = $1;
            $orf = $2;
        } else {
            die("Could not get contig and orf.");
        }
        if (!defined($orf_count->{$contig})) {
            $orf_count->{$contig} = {
                plus => 0,
                minus => 0,
            };
        }
        if ($strand eq '+' || $strand eq '1') {
            $orf_count->{$contig}->{plus}++;
        } else {
            $orf_count->{$contig}->{minus}++;
        }
    } ## End iterating over every detected ORF
    $annotation_in->close();
    ## Now we should have some idea of which strand has the most ORFs for each contig.

  TESTLOOP: while (my $seq = $read_fasta->next_seq) {
      my $id = $seq->id;
      my $test = $orf_count->{$id};
      if ($test->{plus} >= $test->{minus}) {
          ## Then leave this contig alone.
          $write_fasta->write_seq($seq);
          print $write_log "$id was unchanged and has $test->{plus} plus and $test->{minus} minus ORFs.\n";
      } else {
          my $tmp = $seq->revcom();
          $write_fasta->write_seq($tmp);
          print $write_log "$id was reverse-complemented and now has $test->{minus} plus and $test->{plus} minus ORFs.\n";
      }
  } ## Done flipping sequence files.
    $write_log->close();
    return($orf_count);
}

1;

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

=cut
