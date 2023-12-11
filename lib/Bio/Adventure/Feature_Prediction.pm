package Bio::Adventure::Feature_Prediction;
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

 Search for potential tRNA genes across the tree of life.
 aragorn: https://doi.org/10.1093/nar/gkh152

 Use aragorn to search for tRNA genes in a sequence database.  This
 invocation defaults to also searching for tmRNA, but not
 mitochondrial.

 One thing which fascinates me about tRNA genes: mitochondrial tRNA
 genes have introns.  In addition, mtRNA genes are in the introns of
 other genes; and in that context have regulatory effects.  So its like
 an island on a lake with an island in a pond.

=over

=item C<Arguments>

 input(required): Fasta of long sequences, presumably an assembly.
 arbitrary('-rp -fasta -w -m -t'): My favorite aragorn options.
 species(undef):  Currently unused, and TBH I do not remember what I
  intended for it.
 jmem(4): Expected memory usage.
 jprefix('21'): Prefix for the job name and output directory.
 modules('aragorn'): Environment module used.

=cut
sub Aragorn {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        arbitrary => ' -rp -fasta -w -m -t ',
        species => undef,
        jmem => 4,
        jprefix => 21,);
    my $aragorn_args = $options->{arbitrary};
    my $job_name = $class->Get_Job_Name();
    my $output_dir = qq"outputs/$options->{jprefix}aragorn";
    my $output_file = qq"${output_dir}/aragorn.txt";
    my $stdout = qq"${output_dir}/aragorn.stdout";
    my $stderr = qq"${output_dir}/aragorn.stderr";
    my $species_string = '';
    my $comment = '## This is a script to run aragorn.';
    my $jstring = qq!mkdir -p ${output_dir} && \\
  aragorn $options->{arbitrary} \\
    -o ${output_file} \\
    $options->{input} \\
    2>${stderr} \\
    1>${stdout}
!;
    my $aragorn = $class->Submit(
        comment => $comment,
        jcpu => 1,
        jdepends => $options->{jdepends},
        jname => qq"aragorn_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => $options->{jmem},
        output => $output_file,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        stderr => $stderr,
        stdout => $stdout);
    return($aragorn);
}

=back

=head2 C<Glimmer>

 Use glimmer in two passes to search for ORFs in a sequence database.
 Glimmer: https://doi.org/10.1093/nar/26.2.544

 The script template which follows was annoyingly brittle, so I instead
 rewrote it as a short perl script which has some error handling and
 captures each step in a separate function in the hopes that the
 outputs will be clearer and easier to follow.

=over

=item C<Arguments>

 input(required): Fasta of long sequence, likely an assembly.
 jmem(8): Expected memory usage.
 jprefix('16'): prefix for jobnames/output directories.
 modules('glimmer'): Environment module used.

=cut
sub Glimmer {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jmem => 8,
        jprefix => '16',);
    my $job_name = $class->Get_Job_Name();
    my $output_dir = qq"outputs/$options->{jprefix}glimmer";
    my $stdout = qq"${output_dir}/cyoa_glimmer.stdout";
    my $stderr = qq"${output_dir}/cyoa_glimmer.stderr";
    my $comment = '## This is a script to run glimmer.';
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
cyoa_invoke_glimmer.pl --input $options->{input} \\
  --jprefix $options->{jprefix} \\
  2>${stderr} 1>{$stdout}
!;

    ## FIXME: There are a bunch of potentially useful glimmer outputs which should be put here.
    my $glimmer = $class->Submit(
        comment => $comment,
        jcpu => 1,
        jdepends => $options->{jdepends},
        jname => qq"glimmer_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => $options->{jmem},
        output => qq"${output_dir}/${job_name}_glimmer.out",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        stderr => $stderr,
        stdout => $stdout,);
    return($glimmer);
}

=back

=head2 C<Glimmer_Single>

 Invoke glimmer in a single pass.

 This is an alternate invocation of glimmer which does not perform the
 training run. I am somewhat against running glimmer like this, it
 does not take much time to use glimmer's training functions; but some
 assembly methods explicitly run it in a single-pass fashion.

=over

=item C<Arguments>

 input(required): Fasta of long sequence, likely an assembly.
 cutoff(1.1): Used with the invocation of glimmer3.
 overlap(20): Allow overlaps of this length?
 minlength(45): Minimum number of nucleotides.
 threshold(30): the -t option to glimmer3.
 jmem(8): Expected memory usage.
 jprefix('16'): Output directory and job name prefix.
 modules('glimmer'): Environment module to use.

=cut
sub Glimmer_Single {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        cutoff => 1.1,
        overlap => 20,
        minlength => 45,
        threshold => 30,
        jmem => 8,
        jprefix => '16',);
    my $job_name = $class->Get_Job_Name();
    my $output_dir = qq"outputs/$options->{jprefix}glimmer";
    my $comment = '## This is a script to run glimmer.';
    my $final_output = qq"${output_dir}/glimmer3.txt";
    my $final_error = qq"${output_dir}/glimmer3.err";
    my $jstring = qq!mkdir -p ${output_dir}
long-orfs -n -t $options->{cutoff} $options->{input} ${output_dir}/longorfs.txt \\
  2>${output_dir}/longorfs.stderr \\
  1>${output_dir}/longorfs.stdout
extract -t $options->{input} ${output_dir}/longorfs.txt \\
  1>${output_dir}/training.txt \\
  2>${output_dir}/training.stderr
build-icm -r ${output_dir}/single_run.icm < ${output_dir}/training.txt \\
  2>${output_dir}/build-icm.stderr
glimmer3 -o$options->{overlap} -g$options->{minlength} -t$options->{threshold} \\
  $options->{input} \\
  ${output_dir}/single_run.icm \\
  ${output_dir}/glimmer3 \\
  2>${final_error} \\
  1>${final_output}
!;
    my $glimmer = $class->Submit(
        comment => $comment,
        jcpu => 1,
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
        jname => qq"glimmer_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        output => qq"${output_dir}/glimmer3.predict",
        output_detail => qq"${output_dir}/glimmer3.detail",
        output_icm => qq"${output_dir}/single_run.icm",
        output_longorfs => qq"${output_dir}/longorfs.txt",
        output_training => qq"${output_dir}/training.txt",
        stderr => $final_error,
        stdout => $final_output,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($glimmer);
}

=back

=head2 C<Phagepromoter>

 Use the phagepromoter tool to search for potential promoters.
 https://doi.org/10.1093/bioinformatics/btz580

 I pulled this tool off the galaxy toolshed and was disappointed in
 how pooly documented it is.  I have a pretty strong desire to rewrite
 it and make it able to train on other genera and/or make it generic;
 as well as document and clean it up.

=over

=item C<Arguments>

 input(required): Viral assembly fasta file, or genbank file.
 format(fasta): Format of the input(this is the first thing I would fix)
 both_strands(True): Analyze both strands?
 cutoff(0.5): Score cutoff, in my limited experience, 0.6 is a good score.
 family(Podoviridae): Virus family, I would like to make this smarter,
  at least have it read the viral taxonomy and use it.
 host(Pseudomonas): Host genus.  This could be smarter, too.
 phage_type(virulent): phagepromoter was trained on lytic and not strains.
 model(SVM2400): Choose a trained model to feed phagepromoter.
 jmem(6): Expected memory usage.
 jprefix(30): Prefix for the output directory and job name.
 modules(phagepromoter): Environment module to load.

=cut
sub Phagepromoter {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        format => 'fasta',
        both_strands => 'True',
        cutoff => 0.5,
        family => 'Podoviridae',
        host => 'Pseudomonas',
        phage_type => 'virulent',
        model => 'SVM2400',
        jmem => 6,
        jprefix => '30',);
    ## Here are the lines which define what is passed to it on ARGV : defaults
    ## gen_format = sys.argv[1] : genbank or fasta
    ## genome_file = sys.argv[2] : input file
    ## both = sys.argv[3] : True or False
    ## threshold = sys.argv[4] : 0.5
    ## family = sys.argv[5] : Podoviridae
    ## host = sys.argv[6] : Pseudomonas
    ## phage_type = sys.argv[7] : virulent
    ## model = sys.argv[8] : SVM2400
    my $job_name = $class->Get_Job_Name();
    my $output_dir = qq"outputs/$options->{jprefix}phagepromoter";
    my $stdout = qq"${output_dir}/phagepromoter.stdout";
    my $stderr = qq"${output_dir}/phagepromoter.stderr";
    my $output_file = qq"${output_dir}/${job_name}_phagepromoter.tsv";
    my $input_paths = $class->Get_Paths($options->{input});
    my $input_full = $input_paths->[0]->{fullpath};
    my $comment = '## This is a script to run phagepromoter.';
    my $output_fasta = qq"${output_dir}/output.fasta";
    my $jstring = qq!start=\$(pwd)
mkdir -p ${output_dir}
cd ${output_dir}
if phagepromoter.py $options->{format} \\
  ${input_full} \\
  $options->{both_strands} $options->{cutoff} \\
  $options->{family} $options->{host} \\
  $options->{phage_type} $options->{model} \\
  2>phagepromoter.stderr 1>phagepromoter.stdout; then

   echo "phagepromoter passed."
else

  echo "phagepromoter failed, probably because there is more than one contig."
fi
cd \${start}
!;
    my $phagepromoter = $class->Submit(
        comment => $comment,
        jcpu => 1,
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
        jname => qq"phagepromoter",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        output => qq"${output_file}.xz",
        output_fasta => $output_fasta,
        stderr => $stderr,
        stdout => $stdout,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($phagepromoter);
}

=back

=head2 C<Phanotate>

 Search for phage CDS in an assembly.
 phanotate: https://doi.org/10.1093/bioinformatics/btz265

 Phanotate is a cousin to glimmer/prodigal, but with some assumptions
 built in regarding the packing of viral/phage ORFs.  It therefore
 results in a higher number of very closely packed ORFs, with some
 increased risk of spurious calls.

 This function invokes phanotate on a viral assembly to search for
 ORFs.

=over

=item C<Arguments>

 input(required): Viral assembly fasta file.
 jmem(6): Expected memory usage.
 jprefix('17'): Prefix for the jobname and output directory.
 modules('trnascan', 'phanotate'): Environment module used.

=cut
sub Phanotate {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jmem => 6,
        jprefix => '17',);
    my $job_name = $class->Get_Job_Name();
    my $output_dir = qq"outputs/$options->{jprefix}phanotate";
    my $stderr = qq"${output_dir}/phanotate.stderr";
    my $stdout = qq"${output_dir}/phanotate.stdout";
    my $output_file = qq"${output_dir}/${job_name}_phanotate.tsv";
    my $comment = '## This is a script to run phanotate.';
    my $jstring = qq!
mkdir -p ${output_dir}
phanotate.py \\
  --outfile ${output_file} \\
  $options->{input} \\
  2>${stderr} \\
  1>${stdout}
xz -9e -f ${output_file}
!;
    my $phanotate = $class->Submit(
        comment => $comment,
        jcpu => 1,
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
        jname => qq"phanotate_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        output => qq"${output_file}.xz",
        stderr => $stderr,
        stdout => $stdout,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($phanotate);
}

=back

=head2 C<Prodigal>

 Search for CDS using prodigal.
 prodigal: 10.1186/1471-2105-11-119

 Invoke prodigal on an assembly to search for ORFs.  This will by
 default look for an existing training file provided by
 Train_Prodigal()'.

=over

=item C<Arguments>

 input(required): Fasta assembly of non-intron containing sequence.
 species(undef): Training species name, if left undefined prodigal
  will train on the assembly itself.
 gcode(11): Use this genomic code (bacterial).
 output_dir(undef): Put outputs here.
 prodigal_outname(undef): Prefix for the outputs.
 jmem(8): Expected memory usage.
 jprefix('17'): Prefix for the outputs and jobname.
 modules('prodigal'): Load this environment module.

=cut
sub Prodigal {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        species => undef,
        gcode => '11',
        output_dir => undef,
        prodigal_outname => undef,
        edge => 0,
        jmem => 8,
        jprefix => '17',);
    my $edge_string = ' -c ';
    $edge_string = '' if ($options->{edge});

    my $inputs = $class->Get_Paths($options->{input});
    my $job_name = $class->Get_Job_Name();
    $job_name = basename($job_name, ('.fsa'));
    my $train_string = '';
    my $library_file;
    if ($options->{species}) {
        $library_file = qq"$options->{libpath}/hmm/$options->{species}_gc$options->{gcode}.training";
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
    my $stdout = qq"${output_dir}/prodigal.stdout";
    my $stderr = qq"${output_dir}/prodigal.stderr";
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

    my $comment = '## This is a script to run prodigal.';
    my $jstring = qq!mkdir -p ${output_dir}
if prodigal ${train_string} ${edge_string} \\
  -i $options->{input} \\
  -a ${translated_file} \\
  -d ${cds_file} \\
  -s ${scores_file} \\
  -f gff -o ${gff_file} \\
  2>${stderr} \\
  1>${stdout} ; then
  prodigal ${edge_string} ${train_string} \\
    -i $options->{input} \\
    -f gbk -o ${gbk_file} \\
    2>${stderr} \\
    1>${stdout}
else
  echo "Prodigal failed, perhaps the contig was too short?"
fi
sleep 3
!;
    my $prodigal = $class->Submit(
        comment => $comment,
        jcpu => 1,
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
        jname => qq"prodigal_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        stderr => $stderr,
        stdout => $stdout,
        output => $gbk_file,
        output_cds => $cds_file,
        output_gff => $gff_file,
        output_scores => $scores_file,
        output_translated => $translated_file,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        training_input => $library_file,);
    return($prodigal);
}

=back

=head2 C<Rho_Terminase_Predict>

 Use the rho terminase prediction tool to hunt for rho.
 rhotermpredict: https://doi.org/10.1186/s12859-019-2704-x

 Polymerase termination is neat!  Look for rho termination factors
 with this tool.

=over

=item C<Arguments>

 input(required): Fasta of sequence to search.
 jmem(12): Expected memory usage.
 jprefix('51'): Prefix for the jobname/output directory.
 modules('rhotermpredict'): Use this environment module.

=cut
sub Rho_Predict {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        jmem => 12,
        jprefix => '51',
        required => ['input',],);
    my $input_paths = $class->Get_Paths($options->{input});
    my $input_full = $input_paths->[0]->{fullpath};
    my $input_file = $input_paths->[0]->{filename};
    my $cwd_name = basename(cwd());
    my $output_dir = qq"outputs/$options->{jprefix}rhotermpredict_${cwd_name}";
    my $stderr = qq"${output_dir}/rhotermpredict.stderr";
    my $stdout = qq"${output_dir}/rhotermpredict.stdout";
    my $output_file = qq"${output_dir}/predictions_coordinates_${cwd_name}_1.csv";
    my $info_file = qq"${output_dir}/info_about_predictions_${cwd_name}.csv";
    my $jstring = qq?mkdir -p ${output_dir}
start=\$(pwd)
cp $options->{input} ${output_dir}/
cd ${output_dir}
echo ${input_file} | RhoTermPredict_algorithm.py \\
  2>rhotermpredict.stderr \\
  1>rhotermpredict.stdout
cd \${start}
?;

    my $rhoterm = $class->Submit(
        jcpu => 1,
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
        jname => 'rhotermpredict',
        jprefix => $options->{jprefix},
        jstring => $jstring,
        stderr => $stderr,
        stdout => $stdout,
        output => $output_file,
        output_info => $info_file,);
    return($rhoterm);
}

=back

=head2 C<Train_Prodigal>

 Train prodigal to improve its predictions!

 Some assemblies I have been performing are on sets of sequence which
 are too small for prodigal to train itself sufficiently; so this
 function was written to provide an opportunity for one to collate a
 larger sequence database for training.

=over

=item C<Arguments>

 input(required): Fasta training input.
 species(required): Used to set the name of the output training file.
 gcode(11): Use this genomic code.
 jmem(8): Expected memory usage.
 modules('prodigal'): Use this environment module.

=cut
sub Train_Prodigal {
    my ($class, %args) = @_;
    my $check = which('prodigal');
    die('Could not find prodigal in your PATH.') unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        gcode => '11',
        jmem => 8,);
    my $job_name = $class->Get_Job_Name();
    my $kingdom_string = '';
    my $output_dir = qq"$options->{libpath}/hmm";
    my $made = make_path($output_dir);
    my $output = qq"${output_dir}/$options->{species}_gc$options->{gcode}.training";
    my $stdout = qq"${output_dir}/prodigal_training.stdout";
    my $stderr = qq"${output_dir}/prodigal_training.stderr";
    my $comment = '## This is a script to train prodigal.';
    my $jstring = qq!mkdir -p ${output_dir}
prodigal -i $options->{input} \\
  -t ${output} \\
  2>${stderr} \\
  1>${stdout}
!;
    my $prodigal = $class->Submit(
        jcpu => 1,
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => qq"prodigal_training_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => $options->{jmem},
        output => $output,
        stderr => $stderr,
        stdout => $stdout,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($prodigal);
}

=back

=head2 C<tRNAScan>

 Search for tRNA sequence-like elements with infernal.
 tRNAScan: 10.1093/nar/25.5.955

 Alternative to aragorn.  Search for tRNAs!  Having spent a little time
 poking at the two methods, it may prove worth while to employ them
 both.  By default, this invocation is less stringent than the
 previous aragorn invocation.

=over

=item C<Arguments>

 input(required): Fasta file in which to hunt.
 arbitrary('-G'): Any of the many arguments one may pass to trnascan.
 species(undef): Currently unused, I forgot why I wanted it.
 suffix('general'): jobname suffix -- I guess this should be jsuffix.
 tool('trnascan'): tRNAScan comes with a few executables, use this
  one.
 jmem(6): Expected memory usage.
 modules('infernal', 'trnascan'): Use these environment modules.

=cut
sub tRNAScan {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        arbitrary => ' -G ',
        species => undef,
        suffix => 'general',
        tool => 'trnascan',
        jmem => 6,
        jprefix => '28',);
    my $trnascan_args = $options->{arbitrary};
    my $job_name = $class->Get_Job_Name();
    my $output_dir = qq"outputs/$options->{jprefix}trnascan";
    my $output_file = qq"${output_dir}/trnascan_$options->{suffix}.txt";
    my $output_file_v2 = qq"${output_dir}/trnascan_$options->{suffix}_se.txt";
    my $species_string = '';
    my $comment = '## This is a script to run trnascan.';
    my $stdout = qq"${output_dir}/trnascan.stdout";
    my $stderr = qq"${output_dir}/trnascan.stderr";
    my $stdout_v2 = qq"${output_dir}/trnascan_se.stdout";
    my $stderr_v2 = qq"${output_dir}/trnascan_se.stderr";
    my $jstring = qq!mkdir -p ${output_dir}
## Note, trnascan often dies with a SEGFAULT
## I am wrapping it in a subshell, because this seems to happen
## at the end of the process.
set +o errexit
echo "Starting trnascan at $(date)" > ${output_dir}/trnascan_return.txt
first=\$($options->{tool} $options->{arbitrary} \\
  -o ${output_file} \\
  $options->{input} \\
  2>${stderr} \\
  1>${stdout})
echo "Finished first run with $? at $(date)" >> ${output_dir}/trnascan_return.txt

second=\$(tRNAscan-SE -Q -G $options->{input} \\
  -o ${output_file_v2} \\
  2>>${stderr_v2} \\
  1>>${stdout_v2})
echo "Finished second run with $? at $(date)" >> ${output_dir}/trnascan_return.txt
!;
    my $trnascan = $class->Submit(
        comment => $comment,
        jcpu => 1,
        jdepends => $options->{jdepends},
        jname => qq"trnascan_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => $options->{jmem},
        output => $output_file,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        stdout => $stdout,
        stderr => $stderr);
    return($trnascan);
}

1;

=back

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

=cut
