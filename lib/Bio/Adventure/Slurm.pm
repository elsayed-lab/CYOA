package Bio::Adventure::Slurm;
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

sub Check_Sbatch {
    my ($class, %args) = @_;
    my $path = which('sbatch');
    if (!defined($path)) {
        $path = which('bash');
    }
    return($path);
}

=head1 NAME

Bio::Adventure::Slurm - Submit jobs to the Slurm cluster.

=head1 SYNOPSIS

use Bio::Adventure;
my $hpgl = new Bio::Adventure::Slurm;

=head1 METHODS

=head2 C<Submit>

$hpgl->Submit(); invokes sbatch with (hopefully) appropriate
parameters for various jobs on our Slurm cluster.

=cut
sub Submit {
    my ($class, $parent, %args) = @_;
    my $options = $parent->Get_Vars(
        args => \%args);
    my $sbatch = $class->Check_Sbatch();
    my $depends_prefix = '--dependency=afterok';
    ## For arguments to sbatch, start with the defaults in the constructor in $class
    ## then overwrite with any application specific requests from %args
    my $sbatch_log = qq"$options->{logdir}/$options->{jname}.sbatchout";

    my $depends_string = "";
    if ($options->{jdepends}) {
        $depends_string = qq"${depends_prefix}:$options->{jdepends}";
    }
    if (!defined($options->{jname})) {
        $options->{jname} = $class->Get_Job_Name();
    }
    my $script_file = qq"$options->{basedir}/scripts/$options->{jprefix}$options->{jname}.sh";
    my $sbatch_cmd_line = qq"${sbatch} ${depends_string} ${script_file}";
    my $mycwd = getcwd();
    make_path("$options->{logdir}", {verbose => 0}) unless (-r qq"$options->{logdir}");
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
print \$out "## Started $script_file at \${d}";
chdir("$options->{basedir}");
my \$h = Bio::Adventure->new();
?;
        if ($options->{option_file}) {
            $perl_start .= qq!
use Storable qw "lock_retrieve";
local \$Storable::Eval = 1;
use FileHandle;
## Pull options from the option file: $parent->{option_file}
my \$options = lock_retrieve('$parent->{option_file}');
\$h->{options} = \$options;
my \$result;
!;
        }
        my $perl_end = qq!## The following lines give status codes and some logging
my \$jobid = "";
\$jobid = \$ENV{SLURM_JOBID} if (\$ENV{SLURM_JOBID});
my \$end_d = qx'date';
my \$host = qx'hostname';
print \$out "## \${host} finished \${jobid} ${script_base} at \${end_d}.\n";
close(\$out);
!;
        my $total_perl_string = "$perl_start\n";
        $total_perl_string .= "$options->{comment}\n" if ($options->{comment});
        $total_perl_string .= "$options->{prescript}\n" if ($options->{prescript});
        $total_perl_string .= "$options->{jstring}\n" if ($options->{jstring});
        $total_perl_string .= "$perl_end\n";

        my $perl_script = FileHandle->new(">$perl_file");
        print $perl_script $total_perl_string;
        $perl_script->close();
        chmod(0775, $perl_file);
        ## If I get this working properly, change this to:
        ## qq"${perl_file} && rm $options->{option_file}\n";
        $options->{jstring} = qq"set -o errexit
set -o pipefail

${perl_file}

## At least in theory, this rm should not happen if the script run successfully.
## For the moment, however, I will not run it so that I can continue testing my scripts.
## rm $parent->{option_file}
";
    } ## End extra processing for submission of a perl script (perhaps not needed for slurm?

    my $jname = qq"$options->{jprefix}$options->{jname}";

    my $nice_string = '';
    $nice_string = qq"--nice=$options->{jnice}" if (defined($options->{jnice}));

    my $script_start = qq?#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --mail-type=NONE
#SBATCH --chdir=$options->{basedir}
#SBATCH --partition=$options->{jpartition}
#SBATCH --qos=$options->{jqueue} ${nice_string}
#SBATCH --nodes=1 --requeue
#SBATCH --time=$options->{jwalltime}
#SBATCH --job-name=${jname}
#SBATCH --mem=$options->{jmem}G
#SBATCH --cpus-per-task=$options->{cpus}
#SBATCH --output=${sbatch_log}
?;
    if ($options->{array_string}) {
        $script_start .= qq"#SBATCH --array=$options->{array_string}
";
    }

    my $module_string = '';
    if (scalar(@{$options->{modules}} > 0)) {
        $module_string = qq"module add";
        for my $m (@{$options->{modules}}) {
            $module_string .= qq" $m";
        }
    }
    $script_start .= qq?
echo "## Started ${script_file} at \$(date) on \$(hostname) with id \${SLURM_JOBID}." >> outputs/log.txt

${module_string}
?;


    my $script_end = qq!
## The following lines give status codes and some logging
echo "## Job status: \$? " >> outputs/log.txt
echo "## \$(hostname) Finished \${SLURM_JOBID} ${script_base} at \$(date), it took \$(( SECONDS / 60 )) minutes." >> outputs/log.txt
!;
    ## It turns out that if a job was an array (-t) job, then the following does not work because
    ## It doesn't get filled in properly by qstat -f...

    ## The following lines used to be in the shell script postscript
    ## Copying the full job into the log is confusing, removing this at least temporarily.
    ##echo "####This job consisted of the following:" >> outputs/log.txt
    ##cat "\$0" >> outputs/log.txt

    $script_end .= qq!
walltime=\$(scontrol show job "\${SLURM_JOBID}" | grep RunTime | perl -F'/\\s+|=/' -lane '{print \$F[2]}')
echo "#### walltime used by \${SLURM_JOBID} was: \${walltime:-null}" >> outputs/log.txt
maxmem=\$(sstat --format=MaxVMSize -n "\${SLURM_JOBID}.batch")
echo "#### maximum memory used by \${SLURM_JOBID} was: \${maxmem:-null}" >> outputs/log.txt
avecpu=\$(sstat --format=AveCPU -n "\${SLURM_JOBID}.batch")
echo "#### average cpu used by \${SLURM_JOBID} was: \${avecpu:-null}" >> outputs/log.txt
!;

    my $total_script_string = "";
    $total_script_string .= "${script_start}\n";
    $total_script_string .= "$options->{comment}\n" if ($options->{comment});
    $total_script_string .= "$options->{prescript}\n" if ($options->{prescript});
    $total_script_string .= "$options->{jstring}\n" if ($options->{jstring});
    if ($options->{postscript}) {
        $total_script_string .= qq!if [ \$? == "0" ]; then
   $options->{postscript}
fi
!;
    }
    $total_script_string .= "${script_end}\n";

    my $script = FileHandle->new(">$script_file");
    if (!defined($script)) {
        die("Could not write the script: $script_file, check its permissions.")
    }
    print $script $total_script_string;
    $script->close();
    chmod(0755, $script_file);
    my $job_id = undef;
    my $handle = IO::Handle->new;
    my $sbatch_pid = open($handle, "${sbatch_cmd_line} |");
    while (my $line = <$handle>) {
        chomp($line);
        next unless($line =~ /Submitted/);
        $job_id = $line;
        $job_id =~ s/^.*Submitted batch job (\d+)/$1/g;
    }
    my $job;
    if (!defined($job_id)) {
        warn("The job id did not get defined.  sbatch likely failed.");
        return(undef);
    }
    sleep($options->{jsleep});
    my @jobid_list = split(/\./, $job_id);
    my $short_jobid = shift(@jobid_list);

    print "Starting a new job: ${short_jobid} $options->{jname}";
    if ($options->{jdepends}) {
        print ", depending on $options->{jdepends}.";
    }
    print "\n";
    $job = {
        job_id => $job_id,
        log => $sbatch_log,
        pid => $sbatch_pid,
        script_file => $script_file,
    };
    foreach my $k (keys %args) {
        next if ($k eq 'jstring');
        next if ($k eq 'comment');
        $job->{$k} = $args{$k};
    }
    my @wanted_vars = ('basedir', 'cpus', 'depends_string', 'input',
                       'jname', 'jmem', 'jqueue', 'jwalltime', 'output');
    foreach my $w (@wanted_vars) {
        $job->{$w} = $options->{$w} if (!defined($job->{$w}));
    }

    if ($options->{verbose}) {
        use Data::Dumper;
        print Dumper $job;
    }
    return($job);
}

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

L<Bio::Adventure::Torque> L<Bio::Adventure::Local> L<sbatch>

=cut

1;
