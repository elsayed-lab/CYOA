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

has basedir => (is => 'rw', default => getcwd());
has cpus => (is => 'rw', default => '4');
has depends_prefix => (is => 'rw', default => '--dependency=afterok');
has jname => (is => 'rw', default => 'unnamed');
has language => (is => 'rw', default => 'bash');
has loghost => (is => 'rw', default => 'localhost');
has mem => (is => 'rw', default => '6');
has sbatch => (is => 'rw', default => Check_Sbatch());
has sbatch_args => (is => 'rw', default => '--export=ALL --mail-type=NONE');
has sbatch_logdir => (is => 'rw', default => getcwd());
has options => (is => 'rw');
has options_file => (is => 'rw');
has partition => (is => 'rw', default => 'dpart');
has queue => (is => 'rw', default => 'workstation');
has queues => (is => 'rw', default => qq"throughput,workstation,long,large");
##has shell => (is => 'rw', default => '/usr/bin/bash');
has verbose => (is => 'rw', default => 0);
has walltime => (is => 'rw', default => '10:00:00');

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

=head2 Methods

=over 4

=item C<Submit>

    $hpgl->Submit(); invokes sbatch with (hopefully) appropriate
    parameters for various jobs on our Slurm cluster.

=cut
sub Submit {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    ## For arguments to sbatch, start with the defaults in the constructor in $class
    ## then overwrite with any application specific requests from %args
    my $sbatch_log = qq"$options->{sbatch_logdir}/outputs/$options->{jname}.sbatchout";

    my $depends_string = "";
    if ($options->{depends}) {
        $depends_string = qq"$options->{depends_prefix}:$options->{depends}";
    }
    my $script_file = qq"$options->{basedir}/scripts/$options->{jprefix}$options->{jname}.sh";
    my $sbatch_cmd_line = qq"$options->{sbatch} ${depends_string} ${script_file}";
    my $mycwd = getcwd();
    make_path("$options->{basedir}/outputs/status", {verbose => 0}) unless (-r qq"$options->{basedir}/outputs/status");
    make_path("$options->{sbatch_logdir}", {verbose => 0}) unless (-r qq"$options->{sbatch_logdir}");
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
        if ($class->{options_file}) {
            $perl_start .= qq!
use Storable qw "freeze thaw store retrieve";
local $Storable::Eval = 1;
my \$options = retrieve($class->{option_file});
\$h->{options} = \$options;
!;
        }
        my $perl_end = qq!## The following lines give status codes and some logging
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
        chmod(0775, $perl_file);
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
#SBATCH --output=${sbatch_log}
?;
    if ($options->{array_string}) {
        $script_start .= qq"#SBATCH --array=$options->{array_string}
";
    }
    $script_start .= qq?

echo "####Started ${script_file} at \$(date) on \$(hostname)" >> outputs/log.txt
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
walltime=\$(scontrol show job "\${SLURM_JOBID}" | grep RunTime | perl -F'/\\s+|=/' -lane '{print \$F[2]}')
echo "#### walltime used by \${SLURM_JOBID} was: \${walltime:-null}" >> outputs/log.txt
maxmem=\$(sstat --format=MaxVMSize -n "\${SLURM_JOBID}.batch")
echo "#### maximum memory used by \${SLURM_JOBID} was: \${maxmem:-null}" >> outputs/log.txt
avecpu=\$(sstat --format=AveCPU -n "\${SLURM_JOBID}.batch")
echo "#### average cpu used by \${SLURM_JOBID} was: \${avecpu:-null}" >> outputs/log.txt
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
    sleep(1);
    my @jobid_list = split(/\./, $job_id);
    my $short_jobid = shift(@jobid_list);

    print "Starting a new job: ${short_jobid} $options->{jname}";
    if ($options->{depends}) {
        print ", depending on $options->{depends}.";
    }
    print "\n";

    $job = {basedir => $options->{basedir},
            cpus => $options->{cpus},
            depends_string => $depends_string,
            job_args => \%args,
            job_id => $short_jobid,
            job_input => $options->{job_input},
            jname => $options->{jname},
            job_output => $options->{job_output},
            log => $sbatch_log,
            mem => $options->{mem},
            queue => $options->{queue},
            pbs_id => $job_id,
            sbatch_args => $options->{sbatch_args},
            script_body => $options->{jstring},
            script_file => $script_file,
            script_start => $script_start,
            submitter => $sbatch_cmd_line,
            walltime => $options->{walltime},
        };
    return($job);
}

=back

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

    L<Bio::Tools::Run::StandAloneBlast>

=cut

1;
