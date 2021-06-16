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
use File::Spec;
use File::Which qw"which";

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
    my $output_dir = qq"outputs/resfinder";
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
    my $output_dir = qq"outputs/phageterm";

    my $uncompress_string = qq"";
    my $input_string = qq"";
    my $delete_string = qq"";
    if ($options->{input} =~ /\:|\;|\,|\s+/) {
        my @in = split(/\:|\;|\,|\s+/, $options->{input});
        $uncompress_string = qq!
less $in[0] > r1.fastq && less $in[1] > r2.fastq
!;
        $input_string = qq! -f r1.fastq -p r2.fastq !;
        $delete_string = qq!rm r1.fastq && rm r2.fastq!;
    } else {
        $uncompress_string = qq!
less $options->{input} > r1.fastq
!;
        $input_string = qq! -f r1.fastq !;
        $delete_string = qq!rm r1.fastq!;
    }
    
    my $comment = qq!## This is a script to run phageterm.
!;
    my $jstring = qq!mkdir -p ${output_dir}
${uncompress_string}
PhageTerm.py ${input_string} \\
  -r $options->{library} \\
  -c $options->{cpus}
${delete_string}
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
    my $jstring = qq!mkdir -p ${output_dir}
long-orfs -n -t 1.15 $options->{input} ${output_dir}/first_run_longorfs.txt \\ 
  2> ${output_dir}/first_run_longorfs.out 1>&2
extract -t $options->{input} ${output_dir}/first_run_longorfs.txt > \\
  ${output_dir}/first_run_training.txt
build-icm -r ${output_dir}/first_run.icm < ${output_dir}/first_run_training.txt
glimmer3 -o50 -g110 -t30 \\
  $options->{input} \\
  ${output_dir}/first_run.icm \\
  ${output_dir}/first_run.out \\
  2>${output_dir}first_run_glimmer3.out 1>&2

## Use this first run output to go again

extract -t $options->{input} train.coords > ${output_dir}/second_run_training.txt
build-icm -r ${output_dir}/second_run.icm < ${output_dir}/second_run_training.txt
upstream-coords.awk 25 0 train.coords | extract $options->{input} - > \\
  ${output_dir}/second_run_upstream.txt
elph ${output_dir}/second_run_upstream.txt LEN=6 | get-motif-counts.awk > \\
  ${output_dir}/second_run_motif.txt
startuse='start-codon-distrib -3 $options->{input} train.coords'
glmmer3 -o50 -g110 -t30 -b ${output_dir}/second_run_motif.txt -P \${startuse} \\
  $options->{input} ${output_dir}/second_run.icm \\
  ${output_dir}/second_run.out

## Final run

long-orfs -n -t 1.15 $options->{input} ${output_dir}/third_run_longorfs.txt
extract -t $options->{input} ${output_dir}/third_run_longorfs.txt > \\
  ${output_dir}/third_run_training.txt
buid-icm -r ${output_dir}/third_run.icm < ${output_dir}/third_run_training.txt
glimmer3 -o50 -g110 -t30 $options->{input} \\
  ${output_dir}/third_run.icm \\
  ${output_dir}/third_run_intermediate.txt
tail +2 ${output_dir}/third_run_intermediate.txt > ${output_dir}/third_run.coords
upstream-coords.awk 25 0 ${output_dir}/third_run.coords | extract $options->{input} - \\
  ${output_dir}/third_run_upstream.txt
elph ${output_dir}/third_run_upstream.txt LEN=6 | get-motif-counts.awk > \\
  ${output_dir}/third_run_motif.txt
startuse='start-codon-distrib -3 $options->{input} ${output_dir}/third_run.coords'
glimmer3 -o50 -g110 -t30 -b ${output_dir}/third_run_motif.txt -P \${startuse} \\
  $options->{input} ${output_dir}/third_run.icm \\
  ${output_dir}/third_run_output.txt
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

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

=cut

1;
