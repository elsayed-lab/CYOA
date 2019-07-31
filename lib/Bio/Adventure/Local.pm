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

has bash_args => (is => 'rw', default => '');
has bash_logdir => (is => 'rw', default => getcwd());
has basedir => (is => 'rw', default => getcwd());
has cpus => (is => 'rw', default => '4');
has depends_prefix => (is => 'rw', default => '--dependency=afterok');
has jname => (is => 'rw', default => 'unnamed');
has language => (is => 'rw', default => 'bash');
has loghost => (is => 'rw', default => 'localhost');
has mem => (is => 'rw', default => '6');
has options => (is => 'rw');
has option_file => (is => 'rw', default => '');
has partition => (is => 'rw', default => 'dpart');
has queue => (is => 'rw', default => 'workstation');
has queues => (is => 'rw', default => qq"throughput,workstation,long,large");
##has shell => (is => 'rw', default => '/usr/bin/bash');
has verbose => (is => 'rw', default => 0);
has walltime => (is => 'rw', default => '10:00:00');

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
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    ## For arguments to bash, start with the defaults in the constructor in $class
    ## then overwrite with any application specific requests from %args
    my $bash_log = qq"$options->{bash_logdir}/outputs/$options->{jname}.out";

    my $depends_string = "";
    if ($options->{depends}) {
        $depends_string = qq"$options->{depends_prefix}:$options->{depends}";
    }
    my $script_file = qq"$options->{basedir}/scripts/$options->{jprefix}$options->{jname}.sh";
    my $bash_cmd_line = qq"bash ${script_file}";
    my $mycwd = getcwd();
    make_path("$options->{basedir}/outputs/status", {verbose => 0}) unless (-r qq"$options->{basedir}/outputs/status");
    make_path("$options->{bash_logdir}", {verbose => 0}) unless (-r qq"$options->{bash_logdir}");
    make_path("$options->{basedir}/scripts", {verbose => 0}) unless (-r qq"$options->{basedir}/scripts");
    my $script_base = basename($script_file);

    ## Remove the need for two functions that do the same thing except one for perl and one for bash
    if ($options->{language} eq 'perl') {
        my $perl_file = qq"$options->{basedir}/scripts/$options->{jprefix}$options->{jname}.pl";
        my $perl_start = qq?#!/usr/bin/env perl
use strict;
use FileHandle;
use Bio::Adventure;
my \$out = FileHandle->new(">>outputs/log.txt");
my \$d = qx'date';
print \$out "###Started $script_file at \${d}";
chdir("$options->{basedir}");
my \$h = Bio::Adventure->new();
?;
        if ($options->{option_file}) {
            $perl_start .= qq!
use Storable qw "freeze thaw store retrieve";
my \$options = retrieve('$options->{option_file}');
\$h->{options} = \$options;
!;
}
        my $perl_end = qq!
## The following lines give status codes and some logging
my \$jobid = "";
\$jobid = \$ENV{SLURM_JOBID} if (\$ENV{SLURM_JOBID});
my \$end_d = qx'date';
print \$out "####Finished \${jobid} ${script_base} at \${end_d}.";
close(\$out);
!;
        $perl_end .= qq"unlink($class->{option_file});\n" if ($options->{options_file});
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

    my $script_start = qq?#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --mail-type=NONE
#SBATCH --workdir=$options->{basedir}
#SBATCH --partition=$options->{partition}
#SBATCH --qos=$options->{queue}
#SBATCH --nodes=1
#SBATCH --time=$options->{walltime}
#SBATCH --job-name=$options->{jname}
#SBATCH --mem=$options->{mem}G
#SBATCH --cpus-per-task=$options->{cpus}
#SBATCH --output=${bash_log}

echo "####Started ${script_file} at \$(date)" >> outputs/log.txt
cd $options->{basedir} || exit
?;
    my $script_end = qq!## The following lines give status codes and some logging
echo \$? > outputs/status/$options->{jname}.status
echo "###Finished \${SLURM_JOBID} ${script_base} at \$(date), it took \$(( SECONDS / 60 )) minutes." >> outputs/log.txt
!;
    ## It turns out that if a job was an array (-t) job, then the following does not work because
    ## It doesn't get filled in properly by qstat -f...

    ## The following lines used to be in the shell script postscript
    ## Copying the full job into the log is confusing, removing this at least temporarily.
    ##echo "####This job consisted of the following:" >> outputs/log.txt
    ##cat "\$0" >> outputs/log.txt

    $script_end .= qq!
##walltime=\$(qstat -f -t \"\${SLURM_JOBID}\" | grep 'resources_used.walltime' | awk -F ' = ' '{print \$2}')
##echo "#### walltime used by \${SLURM_JOBID} was: \${walltime:-null}" >> outputs/log.txt
##mem=\$(qstat -f -t | grep \"\${SLURM_JOBID}\" | grep 'resources_used.mem' | awk -F ' = ' '{print \$2}')
##echo "#### memory used by \${SLURM_JOBID} was: \${mem:-null}" >> outputs/log.txt
##vmmemory=\$(qstat -f -t \"\${SLURM_JOBID}\" | grep 'resources_used.vmem' | awk -F ' = ' '{print \$2}')
##echo "#### vmemory used by \${SLURM_JOBID} was: \${vmmemory:-null}" >> outputs/log.txt
##cputime=\$(qstat -f -t \"\${SLURM_JOBID}\" | grep 'resources_used.cput' | awk -F ' = ' '{print \$2}')
##echo "#### cputime used by \${SLURM_JOBID} was: \${cputime:-null}" >> outputs/log.txt
####qstat -f -t \${SLURM_JOBID} >> outputs/log.txt
!;

    if ($options->{verbose}) {
        print qq"The job is:
$args{jstring}
";
    }

    my $total_script_string = "";
    $total_script_string .= "${script_start}\n";
    $total_script_string .= "$args{comment}\n" if ($args{comment});
    $total_script_string .= "$args{prescript}\n" if ($args{prescript});
    $total_script_string .= "$args{jstring}\n" if ($args{jstring});
    if ($args{postscript}) {
        $total_script_string .= qq!if [ \$? == "0" ]; then
   $args{postscript}
fi
!;
    }
    $total_script_string .= "${script_end}\n";

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
    if ($args{depends}) {
        print ", depending on $args{depends}.";
    }
    print "\n";

    $job = {
        basedir => $options->{basedir},
        cpus => $options->{cpus},
        depends_string => $depends_string,
        job_args => \%args,
        job_id => $job_id,
        job_input => $options->{job_input},
        jname => $options->{jname},
        job_output => $options->{job_output},
        job_text => $job_text,
        log => $bash_log,
        mem => $options->{mem},
        queue => $options->{queue},
        pbs_id => $job_id,
        bash_args => $options->{bash_args},
        script_body => $args{jstring},
        script_file => $script_file,
        script_start => $script_start,
        submitter => $bash_cmd_line,
        walltime => $options->{walltime},
    };
    return($job);
}

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

    L<Bio::Adventure::Slurm> L<Bio::Adventure::Torque>

=cut

1;
