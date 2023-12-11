package Bio::Adventure::Resistance;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';
use Bio::Adventure::Config;
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

 Abricate provides a series of blast-based searches for resistance databases.
 https://github.com/tseemann/abricate

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
        identity => 80,);
    my $coverage = 70;
    $coverage = $options->{coverage} if (defined($options->{coverage}));
    my $identity = 70;
    $identity = $options->{identity} if (defined($options->{identity}));

    my $job_name = $class->Get_Job_Name();
    my $input_paths = $class->Get_Paths($options->{input});
    my $input_dir = $input_paths->[0]->{dirname};
    my $input_ext = $input_paths->[0]->{filebase_extension};
    my $input_name;
    if (defined($input_dir)) {
        $input_name = $input_dir;
    } else {
        $input_name = $input_ext;
    }
    my $output_dir = qq"outputs/$options->{jprefix}abricate_${input_name}";
    my $species_string = '';
    my $comment = '## This is a script to run abricate.';
    my $stderr = qq"${output_dir}/abricate.stderr";
    my $output_txt = qq"${output_dir}/abricate_summary.txt";
    my $jstring = qq!mkdir -p ${output_dir}
## First get the list of available databases:
dbs=\$(abricate --list | { grep -v "^DATABASE" || test \$? = 1; } | awk '{print \$1}')
for db in \${dbs}; do
  abricate $options->{input} --db \${db} --nopath --noheader \\
  --minid ${identity} --mincov ${coverage} \\
  2>>${stderr} \\
  1>${output_dir}/abricate_\${db}.tsv
  cat ${output_dir}/abricate_\${db}.tsv >> ${output_dir}/abricate_combined.tsv
done
abricate --summary ${output_dir}/*.tsv \\
  2>>${stderr} \\
  1>${output_txt}
!;

    my $abricate = $class->Submit(
        comment => $comment,
        jcpu => 4,
        jdepends => $options->{jdepends},
        jmem => 24,
        jname => qq"abricate_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        output => qq"${output_dir}/abricate_combined.tsv",
        output_argannot => qq"${output_dir}/abricate_argannot.tsv",
        output_card => qq"${output_dir}/abricate_card.tsv",
        output_ecoh => qq"${output_dir}/abricate_ecoh.tsv",
        output_ecoli => qq"${output_dir}/abricate_ecoli_vf.tsv",
        output_megares => qq"${output_dir}/abricate_megares.tsv",
        output_ncbi => qq"${output_dir}/abricate_ncbi.tsv",
        output_plasmidfinder => qq"${output_dir}/abricate_plasmidfinder.tsv",
        output_resfinder => qq"${output_dir}/abricate_resfinder.tsv",
        output_txt => ${output_txt},
        output_vfdb => qq"${output_dir}/abricate_vfdb.tsv",
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($abricate);
}

=head2 C<Resfinder>

 Perform a resfinder search.
 10.1093/jac/dkaa345

 Resfinder provides a database of resistance genes and search function
 so that one may relatively quickly check a sequence database/assembly
 for potentially problematic genes/point mutations.

=cut
sub Resfinder {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jprefix => '18',
        arbitrary => ' -l 0.6 -t 0.8 --acquired ',);
    my $resfinder_args = $options->{arbitrary};
    my $job_name = $class->Get_Job_Name();
    my $assembly_name = basename(dirname($options->{input}));
    my $output_dir = qq"outputs/$options->{jprefix}resfinder";
    my $species_string = qq"";
    if (defined($options->{species})) {
        $species_string = qq" --species $options->{species} ";
    }

    my $comment = '## This is a script to run resfinder.';
    my $jstring = qq!mkdir -p ${output_dir}
run_resfinder.py -ifa $options->{input} \\
  -o ${output_dir} \\
  ${resfinder_args} ${species_string}
!;
    my $resfinder = $class->Submit(
        comment => $comment,
        jcpu => 4,
        jdepends => $options->{jdepends},
        jmem => 24,
        jname => qq"resfinder_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        output => qq"${output_dir}/resfinder.txt",
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($resfinder);
}

=head2 C<Rgi>

 Perform a RGI search.
 10.1093/nar/gkz935

 RGI is an alternative to Resfinder, I have not really explored it yet,
 but it appears to provide a somewhat more in-depth database of
 interesting genes than resfinder.  Its database structure is a bit unwieldy.

=cut
sub Rgi {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        arbitrary => '',
        jprefix => '15',);
    my $rgi_args = $options->{arbitrary};
    my $job_name = $class->Get_Job_Name();
    my $assembly_name = basename(dirname($options->{input}));
    my $output_dir = qq"outputs/$options->{jprefix}rgi";
    my $species_string = qq"";
    my $comment = '## This is a script to run rgi.';
    my $jstring = qq!mkdir -p ${output_dir}
rgi main --input_sequence $options->{input} \\
  --output_file ${output_dir}/rgi_result.txt --input_type protein \\
  --include_loose --clean
!;
    my $rgi = $class->Submit(
        comment => $comment,
        jcpu => 4,
        jdepends => $options->{jdepends},
        jmem => 24,
        jname => "rgi_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        output => qq"${output_dir}/rgi_result.txt",
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($rgi);
}

1;

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

=cut
