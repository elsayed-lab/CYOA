package Bio::Adventure::Resistance;
## LICENSE: gplv2
## ABSTRACT:  Kitty!
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use Cwd qw(abs_path getcwd);
use File::Basename;
use File::Spec;
use File::Path qw"make_path";
use File::Which qw"which";
use File::ShareDir qw":ALL";

=head1 NAME

Bio::Adventure::Annotation - Do some searches to help annotate genes.

=head1 SYNOPSIS

=head2 C<Abricate>

Invoke abricate, ask it for the set of available databases, and try them all against the
full assembled contigs (not ORFs).

=cut
sub Abricate {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jprefix => '18',
        coverage => 80,
        identity => 80,
        modules => ['any2fasta', 'abricate', 'blast'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('abricate');
    die("Could not find abricate in your PATH.") unless($check);

    my $coverage = 70;
    $coverage = $options->{coverage} if (defined($options->{coverage}));
    my $identity = 70;
    $identity = $options->{identity} if (defined($options->{identity}));

    my $job_name = $class->Get_Job_Name();
    my $input_paths = $class->Get_Paths($options->{input});
    my $input_name;
    if (defined($input_paths->{dirname})) {
        $input_name = $input_paths->{dirname};
    } else {
        $input_name = $input_paths->{filebase_extension};
    }
    my $output_dir = qq"outputs/$options->{jprefix}abricate_${input_name}";
    my $species_string = qq"";
    my $comment = qq!## This is a script to run abricate.
!;
    my $jstring = qq!mkdir -p ${output_dir}
## First get the list of available databases:
dbs=\$(abricate --list | grep -v "^DATABASE" | awk '{print \$1}')
for db in \${dbs}; do
  abricate $options->{input} --db \${db} --nopath --noheader \\
  --minid ${identity} --mincov ${coverage} \\
  2>${output_dir}/abricate_\${db}.err \\
  1>${output_dir}/abricate_\${db}.tsv
  cat ${output_dir}/abricate_\${db}.tsv >> ${output_dir}/abricate_combined.tsv
done
abricate --summary ${output_dir}/*.tsv \\
  2>${output_dir}/abricate_summary.err \\
  1>${output_dir}/abricate_summary.txt
!;

    my $abricate = $class->Submit(
        cpus => 6,
        comment => $comment,
        jdepends => $options->{jdepends},
        jmem => 24,
        jname => "abricate_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        modules => $options->{modules},
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        output => qq"${output_dir}/abricate_combined.tsv",
        output_argannot => qq"${output_dir}/abricate_argannot.tsv",
        output_card => qq"${output_dir}/abricate_card.tsv",
        output_ecoh => qq"${output_dir}/abricate_ecoh.tsv",
        output_ecoli => qq"${output_dir}/abricate_ecoli_vf.tsv",
        output_megares => qq"${output_dir}/abricate_megares.tsv",
        output_ncbi => qq"${output_dir}/abricate_ncbi.tsv",
        output_plasmidfinder => qq"${output_dir}/abricate_plasmidfinder.tsv",
        output_resfinder => qq"${output_dir}/abricate_resfinder.tsv",
        output_summary => qq"${output_dir}/abricate_summary.txt",
        output_vfdb => qq"${output_dir}/abricate_vfdb.tsv",);
    return($abricate);
}

=head2 C<Resfinder>

Resfinder provides a database of resistance genes and search function
so that one may relatively quickly check a sequence database/assembly
for potentially problematic genes/point mutations.

=cut
sub Resfinder {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        species => undef,
        jprefix => '18',
        arbitrary => ' -l 0.6 -t 0.8 --acquired ',
        modules => ['resfinder'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('run_resfinder.py');
    die("Could not find resfinder in your PATH.") unless($check);

    my $resfinder_args = $options->{arbitrary};
    my $job_name = $class->Get_Job_Name();
    my $assembly_name = basename(dirname($options->{input}));
    my $output_dir = qq"outputs/$options->{jprefix}resfinder";
    my $species_string = qq"";
    if (defined($options->{species})) {
        $species_string = qq" --species $options->{species} ";
    }

    my $comment = qq!## This is a script to run resfinder.
!;
    my $jstring = qq!mkdir -p ${output_dir}
run_resfinder.py -ifa $options->{input} \\
  -o ${output_dir} \\
  ${resfinder_args} ${species_string}
!;
    my $resfinder = $class->Submit(
        cpus => 6,
        comment => $comment,
        jdepends => $options->{jdepends},
        jmem => 24,
        jname => "resfinder_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        modules => $options->{modules},
        output => qq"${output_dir}/resfinder.txt",
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($resfinder);
}

=head2 C<Rgi>

RGI is an alternative to Resfinder, I have not really explored it yet,
but it appears to provide a somewhat more in-depth database of
interesting genes than resfinder.  Its database structure is a bit unwieldy.

=cut
sub Rgi {
    my ($class, %args) = @_;
    my $check = which('rgi');
    die("Could not find rgi in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        arbitrary => '',
        jprefix => '15',
        modules => ['kma', 'jellyfish', 'bowtie', 'bwa', 'diamond', 'rgi'],);
    my $rgi_args = $options->{arbitrary};
    my $job_name = $class->Get_Job_Name();
    my $assembly_name = basename(dirname($options->{input}));
    my $output_dir = qq"outputs/$options->{jprefix}rgi";
    my $species_string = qq"";
    my $comment = qq!## This is a script to run rgi.
!;
    my $jstring = qq!mkdir -p ${output_dir}
rgi main --input_sequence $options->{input} \\
  --output_file ${output_dir}/rgi_result.txt --input_type protein \\
  --include_loose --clean
!;
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $rgi = $class->Submit(
        cpus => 6,
        comment => $comment,
        jdepends => $options->{jdepends},
        jmem => 24,
        jname => "rgi_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        modules => $options->{modules},
        output => qq"${output_dir}/rgi_result.txt",
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($rgi);
}

1;

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

=cut
