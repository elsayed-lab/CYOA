package Bio::Adventure::Local;
## LICENSE: gplv2
## ABSTRACT:  Kitty!
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use Cwd;
use File::Basename qw "basename";
use File::Path qw"make_path remove_tree";
use File::Which qw"which";
use IO::Handle;

=head1 NAME

Bio::Adventure::Local - Write out jobs using scripts for running jobs on the local computer.

=head1 SYNOPSIS

use Bio::Adventure;
my $hpgl = new Bio::Adventure::Local;
$hpgl->Submit(%job_args);

This should invoke scripts with appropriate parameters for various jobs on the local host.

=head1 Methods

=head2 C<Submit>

This is responsible for the final writeup and submission of a local script.

=cut
sub Submit {
    my ($class, $parent, %args) = @_;
    my $options = $parent->Get_Vars(
        args => \%args,
        jprefix => '',
        jname => 'unknown');
    ## For arguments to bash, start with the defaults in the constructor in $class
    ## then overwrite with any application specific requests from %args
    my $script_file = qq"$options->{basedir}/scripts/$options->{jprefix}$options->{jname}.sh";
    my $bash_cmd_line = qq"bash ${script_file}";
    my $mycwd = getcwd();
    make_path("$options->{logdir}", {verbose => 0}) unless (-r qq"$options->{logdir}");
    make_path("$options->{basedir}/scripts", {verbose => 0}) unless (-r qq"$options->{basedir}/scripts");
    my $script_base = basename($script_file);
    my $bash_log = 'outputs/log.txt';

    ## Remove the need for two functions that do the same thing except one for perl and one for bash
    if ($options->{language} eq 'perl') {
        my $perl_file = qq"$options->{basedir}/scripts/$options->{jprefix}$options->{jname}.pl";
        my $perl_start = qq?#!/usr/bin/env perl
use strict;
use FileHandle;
use Bio::Adventure;
my \$out = FileHandle->new(">>${bash_log}");
my \$d = qx'date';
print \$out "###Started $script_file at \${d}";
chdir("$options->{basedir}");
my \$h = Bio::Adventure->new();
?;
##        if ($options->{option_file}) {
##            $perl_start .= qq!
##use Storable qw "freeze thaw store retrieve";
##my \$options = retrieve('$options->{option_file}');
##\$h->{options} = \$options;
##!;
##}
        my $perl_end = qq!
## The following lines give status codes and some logging
my \$jobid = "";
\$jobid = \$ENV{SLURM_JOBID} if (\$ENV{SLURM_JOBID});
my \$end_d = qx'date';
print \$out "####Finished \${jobid} ${script_base} at \${end_d}.";
close(\$out);
!;
        ## $perl_end .= qq"unlink($class->{option_file});\n" if ($options->{options_file});
        print "The job is:
$args{jstring}" if ($options->{verbose});
        my $total_perl_string = "";
        $total_perl_string .= "$perl_start\n";
        $total_perl_string .= "$args{comment}\n" if ($args{comment});
        $total_perl_string .= "$args{prescript}\n" if ($args{prescript});
        $total_perl_string .= "$args{jstring}\n" if ($args{jstring});
        $total_perl_string .= "$perl_end\n";

        my $perl_script = FileHandle->new(">$perl_file");
        print $perl_script $total_perl_string;
        $perl_script->close();
        chmod(0755, $perl_file);
        $args{jstring} = qq"${perl_file}\n";
    } ## End extra processing for submission of a perl script (perhaps not needed for slurm?

    my $module_string = '';
    if (scalar(@{$options->{modules}} > 0)) {
        $module_string = 'module add';
        for my $m (@{$options->{modules}}) {
            $module_string .= qq" ${m}";
        }
    }

    my $script_start = qq?#!/usr/bin/env bash
set -o errexit
set -o pipefail
echo "####Started ${script_file} at \$(date)" >> ${bash_log}
${module_string}
?;
    my $script_end = qq!## The following lines give status codes and some logging
echo "###Job status:\$?" >> ${bash_log}
echo "###Finished ${script_base} at \$(date), it took \$(( SECONDS / 60 )) minutes." >> ${bash_log}!;
    if ($options->{verbose}) {
        print qq"The job is:
$args{jstring}
";
    }

    my $total_script_string = '';
    $total_script_string .= qq"${script_start}\n";
    $total_script_string .= qq"$args{comment}\n" if ($args{comment});
    $total_script_string .= qq"$args{prescript}\n" if ($args{prescript});
    $total_script_string .= qq"$args{jstring}\n" if ($args{jstring});
    if ($args{postscript}) {
        $total_script_string .= qq!if [ \$? == "0" ]; then
   $args{postscript}
fi
!;
    }
    $total_script_string .= qq"${script_end}\n";

    my $script = FileHandle->new(">$script_file");
    print $script $total_script_string;
    $script->close();
    chmod(0755, $script_file);
    my $job_text = "";
    my $handle = IO::Handle->new;
    my $bash_pid = open($handle, "${bash_cmd_line} |");
    while (my $line = <$handle>) {
        $job_text = $job_text . $line;
    }
    my $job;
    sleep(1);
    my $job_id = $bash_pid;

    print "Starting a new job: ${job_id} $options->{jname}";
    if ($options->{jdepends}) {
        print ", depending on $options->{jdepends}.";
    }
    print "\n";
    $job = {
        job_id => $job_id,
        log => $bash_log,
        pid => $bash_pid,
        script_file => $script_file,
    };
    foreach my $k (keys %args) {
        next if ($k eq 'jstring');
        next if ($k eq 'comment');
        $job->{$k} = $args{$k};
    }
    my @wanted_vars = ('basedir', 'cpus', 'input', 'jname', 'jmem',
                       'jqueue', 'jwalltime', 'output');
    foreach my $w (@wanted_vars) {
        $job->{$w} = $options->{$w} if (!defined($job->{$w}));
    }

    if ($options->{verbose}) {
        use Data::Dumper;
        print Dumper $job;
    }
    ## Reset the environment in case we left any cruft behind
    my $reset = Bio::Adventure::Reset_Vars($parent);
    return($job);
}

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

    L<Bio::Adventure::Slurm> L<Bio::Adventure::Torque>

=cut

1;
