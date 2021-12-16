package Bio::Adventure::Feature_Prediction;
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

One thing which fascinates me about tRNA genes: mitochondrial tRNA
genes have introns.  In addition, mtRNA genes are in the introns of
other genes; and in that context have regulatory effects.  So its like
an island on a lake with an island in a pond.

=cut
sub Aragorn {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        modules => ['aragorn'],
        jmem => 8,
        jprefix => 21,
        species => undef,
        arbitrary => ' -rp -fasta -w -m -t ',);
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
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => $options->{jmem},
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

The script template which follows was annoyingly brittle, so I instead
rewrote it as a short perl script which has some error handling and
captures each step in a separate function in the hopes that the
outputs will be clearer and easier to follow.

=cut
sub Glimmer {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        modules => ['glimmer'],
        jmem => 8,
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
        jmem => $options->{jmem},
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

Run glimmer in a single pass instead of the two passes as per the
paper.  I am somewhat against running glimmer like this, it does not
take much time to use glimmer's training functions; but some assembly
methods explicitly run it in a single-pass fashion.

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
        jmem => 8,
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
        jmem => $options->{jmem},
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

=head2 C<Phanotate>

Phanotate is a cousin to glimmer/prodigal, but with some assumptions
built in regarding the packing of viral/phage ORFs.  It therefore
results in a higher number of very closely packed ORFs, with some
increased risk of spurious calls.

This function invokes phanotate on a viral assembly to search for
ORFs.

=cut
sub Phanotate {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jmem => 6,
        jprefix => '17',
        modules => ['phanotate'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('phanotate.py');
    die("Could not find glimmer in your PATH.") unless($check);
    my $job_name = $class->Get_Job_Name();
    my $output_dir = qq"outputs/$options->{jprefix}phanotate";
    my $output_file = qq"${output_dir}/${job_name}_phanotate.tsv";
    my $input_paths = $class->Get_Paths($output_file);
    my $comment = qq"## This is a script to run phanotate.";
    my $jstring = qq!
phanotate.py \\
  --outfile ${output_file} \\
  $options->{input}
xz -9e -f ${output_file}
!;
    my $phanotate = $class->Submit(
        comment => $comment,
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
        jname => "phanotate_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        modules => $options->{modules},
        output => qq"${output_file}.xz",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => "workstation",);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($phanotate);
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
        jmem => 8,
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
        jmem => $options->{jmem},
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

=head2 C<Rho_Terminase_Predict>

Use the rho terminase prediction tool to hunt for rho.

=cut
sub RhoTermPredict {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input',],
        modules => ['rhotermpredict'],
        jmem => 12,
        jprefix => '51',);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $input_paths = $class->Get_Paths($options->{input});
    my $input_full = $input_paths->{fullpath};
    my $output_dir = qq"outputs/$options->{jprefix}rhotermpredict_$input_paths->{dirname}";
    my $output_file = qq"${output_dir}/predictions_coordinates_seqname.csv";
    my $info_file = qq"${output_dir}/info_about_predictions_seqname.csv";
    my $jstring = qq?mkdir -p ${output_dir}
start=\$(pwd)
cd ${output_dir}
cp $options->{input} .
echo $input_paths->{filename} | RhoTermPredict_algorithm.py

?;

    my $rhoterm = $class->Submit(
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
        jname => 'rhotermpredict',
        jprefix => $options->{jprefix},
        jstring => $jstring,
        modules => $options->{modules},
        output_file => $output_file,
        output_info => $info_file,);
    return($rhoterm);
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
        jmem => 8,
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
        jmem => $options->{jmem},
        modules => $options->{modules},
        output => $output,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => "workstation",);
    return($prodigal);
}

=head2 C<tRNAScan>

Alternative to aragorn.  Search for tRNAs!  Having spent a little time
poking at the two methods, it may prove worth while to employ them
both.

=cut
sub tRNAScan {
    my ($class, %args) = @_;
    my $check = which('trnascan');
    die("Could not find trnascan in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jmem => 6,
        species => undef,
        arbitrary => ' -G ',
        suffix => 'general',
        tool => 'trnascan',
        modules => ['infernal', 'trnascan'],);

    my $trnascan_args = $options->{arbitrary};
    my $job_name = $class->Get_Job_Name();
    my $output_dir = qq"outputs/trnascan";
    my $output_file = qq"${output_dir}/trnascan_$options->{suffix}.txt";
    my $species_string = qq"";
    my $comment = qq!## This is a script to run trnascan.
!;
    my $jstring = qq!mkdir -p ${output_dir}
$options->{tool} $options->{arbitrary} \\
  -o ${output_file} \\
  $options->{input} \\
  2>${output_dir}/trnascan.stderr 1>${output_dir}/trnascan.stdout
!;

    my $trnascan = $class->Submit(
        cpus => 6,
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => "trnascan_${job_name}",
        jprefix => "64",
        jstring => $jstring,
        jmem => $options->{jmem},
        modules => $options->{modules},
        output => $output_file,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => "workstation",);
    return($trnascan);
}

1;

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

=cut
