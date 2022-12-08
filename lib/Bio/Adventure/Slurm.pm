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
use File::Basename qw "basename dirname";
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

sub Get_QOS {
    my $avail_qos = {};
    my $qos = FileHandle->new("sacctmgr -p show qos |");
    my $count = 0;
  QOS: while (my $line = <$qos>) {
      $count++;
      next QOS if ($count == 1);

      my ($name, $priority, $gracetime, $preempt, $preempt_exempt, $preempt_mode, $flags,
          $usage_thresh, $usage_factor, $group_tres, $group_tres_min, $group_tres_run_min,
          $group_jobs, $group_submit, $group_wall, $max_resources_per_job, $max_tres_per_node,
          $max_tres_min, $max_wall, $max_resources_per_user, $max_jobs_pu, $max_submit_pu, $max_tres_pa,
          $max_jobs_pa, $max_submit_pa, $min_tres) = split(/\|/, $line);
      my $max_job_cpu = 0;
      my $max_job_gpu = 0;
      my $max_job_mem = 0;
      my $max_user_cpu = 0;
      my $max_user_gpu = 0;
      my $max_user_mem = 0;
      my $max_jobs = -1;
      my $max_hours = 0;

      if ($max_resources_per_job ne '') {
          $max_job_cpu = $max_resources_per_job;
          if ($max_job_cpu =~ /cpu=/) {
              $max_job_cpu =~ s/.*cpu=(\d+),.*$/$1/g;
          } else {
              $max_job_cpu = 0;
          }
          $max_job_gpu = $max_resources_per_job;
          if ($max_job_gpu =~ m/gpu=/) {
              $max_job_gpu =~ s/.*gpu=(\d+),.*$/$1/g;
          } else {
              $max_job_gpu = 0;
          }
          $max_job_mem = $max_resources_per_job;
          if ($max_job_mem =~ /mem=/) {
              $max_job_mem =~ s/.*mem=(\d+)G.*$/$1/g;
          } else {
              $max_job_mem = 0;
          }
      }

      if ($max_resources_per_user ne '') {
          $max_user_cpu = $max_resources_per_user;
          if ($max_user_cpu =~ /cpu=/) {
              $max_user_cpu =~ s/.*cpu=(\d+),.*$/$1/g;
          } else {
              $max_user_cpu = $max_job_cpu;
          }

          $max_user_gpu = $max_resources_per_user;
          if ($max_user_gpu =~ m/gpu=/) {
              $max_user_gpu =~ s/.*gpu=(\d+),.*$/$1/g;
          } else {
              $max_user_gpu = $max_job_cpu;
          }
          $max_user_mem = $max_resources_per_user;
          if ($max_user_mem =~ /mem=/) {
              $max_user_mem =~ s/.*mem=(\d+)G.*$/$1/g;
          } else {
              $max_user_mem = $max_job_mem;
          }
      }

      if ($max_wall ne '') {
          my ($days, $hms) = split(/\-/, $max_wall);
          my ($hours, $min, $sec) = split(/:/, $hms);
          $days = 0 if (!defined($days));
          $hours = 0 if (!defined($hours));
          $max_hours = ($days * 24) + $hours;
      }

      $avail_qos->{$name} = {
          max_job_cpu => $max_job_cpu,
          max_job_gpu => $max_job_gpu,
          max_job_mem => $max_job_mem,
          max_user_cpu => $max_user_cpu,
          max_user_gpu => $max_user_gpu,
          max_user_mem => $max_user_mem,
          max_jobs => $max_jobs,
          max_hours => $max_hours,
      };
      $avail_qos->{$name}->{used_cpu} = 0 if (!defined($avail_qos->{$name}->{used_cpu}));
      $avail_qos->{$name}->{used_gpu} = 0 if (!defined($avail_qos->{$name}->{used_gpu}));
      $avail_qos->{$name}->{used_mem} = 0 if (!defined($avail_qos->{$name}->{used_mem}));
      $avail_qos->{$name}->{used_hours} = 0 if (!defined($avail_qos->{$name}->{used_hours}));

  }
    $qos->close();

    ## Run through the qos keys and set undefined options to those from the default.
    for my $n (keys %{$avail_qos}) {
        if ($avail_qos->{$n}->{max_job_cpu} == 0) {
            $avail_qos->{$n}->{max_job_cpu} = $avail_qos->{default}->{max_job_cpu};
        }
        if ($avail_qos->{$n}->{max_job_gpu} == 0) {
            $avail_qos->{$n}->{max_job_gpu} = $avail_qos->{default}->{max_job_gpu};
        }
        if ($avail_qos->{$n}->{max_job_mem} == 0) {
            $avail_qos->{$n}->{max_job_mem} = $avail_qos->{default}->{max_job_mem};
        }
        if ($avail_qos->{$n}->{max_user_cpu} == 0) {
            $avail_qos->{$n}->{max_user_cpu} = $avail_qos->{default}->{max_user_cpu};
        }
        if ($avail_qos->{$n}->{max_user_gpu} == 0) {
            $avail_qos->{$n}->{max_user_gpu} = $avail_qos->{default}->{max_user_gpu};
        }
        if ($avail_qos->{$n}->{max_user_mem} == 0) {
            $avail_qos->{$n}->{max_user_mem} = $avail_qos->{default}->{max_user_mem};
        }
        if ($avail_qos->{$n}->{max_hours} == 0) {
            $avail_qos->{$n}->{max_hours} = $avail_qos->{default}->{max_hours};
        }
    }
}

sub Get_Associations {
    my $associations = {};
    my $assoc = FileHandle->new("sacctmgr -p show associations | grep $ENV{USER} |");
    while (my $line = <$assoc>) {
        my ($cluster, $account, $user, $partition, $share, $priority, $group_jobs, $group_tres,
            $group_wall, $group_tresmin, $max_jobs, $max_tres, $max_tres_per_node, $max_submit,
            $max_wall, $max_tres_min, $undef, $qos_lst, $def_qos, $group_tres_run_min) = split(/\|/, $line);
        my @possible_qos = split(/\,/, $qos_lst);
        my $inner_hash = {
            share => $share,
            qos => \@possible_qos,
            default_qos => $def_qos };
        $associations->{$cluster}->{$account} = $inner_hash;
    }
    $assoc->close();
}

sub Choose_Spec {
    my %args = @_;
    my $associations = $args{associations};
    my $qos_info = $args{qos_info};
    my $wanted_spec = $args{wanted_spec};

    my $chosen_account = '';
    my $chosen_qos = '';
  TOP: for my $cluster (keys %{$associations}) {
    ACCOUNT: for my $account (keys %{$associations->{$cluster}}) {
        my @qos = @{$associations->{$cluster}->{$account}->{qos}};
        my $found_qos = 0;
      QOS: for my $q (@qos) {
          my $qos_info = $avail_qos->{$q};
          my $stringent_mem = $qos_info->{max_job_mem} + $qos_info->{used_mem};
          my $stringent_cpu = $qos_info->{max_job_cpu} + $qos_info->{used_cpu};
          my $stringent_gpu = $qos_info->{max_job_gpu} + $qos_info->{used_gpu};
          my $stringent_hours = $qos_info->{max_hours} + $qos_info->{used_hours};
          ## If we pass this initial test, then the job should start immediately.
          if ($wanted_spec->{mem} <= $stringent_mem &&
              $wanted_spec->{cpu} <= $stringent_cpu &&
              $wanted_spec->{gpu} <= $stringent_gpu &&
              $wanted_spec->{walltime} <= $stringent_hours) {
              print "Found qos in first pass: $q wanted: $wanted_spec->{mem} $wanted_spec->{cpu} $wanted_spec->{walltime} vs. $stringent_mem $stringent_cpu $stringent_hours\n";
              $found_qos++;
              $qos_info->{used_mem} = $qos_info->{used_mem} + $wanted_spec->{mem};
              $qos_info->{used_cpu} = $qos_info->{used_cpu} + $wanted_spec->{cpu};
              $qos_info->{used_gpu} = $qos_info->{used_gpu} + $wanted_spec->{gpu};
              $qos_info->{used_hours} = $qos_info->{used_hours} + $wanted_spec->{walltime};
              $chosen_account = $account;
              $chosen_qos = $q;
              last TOP;
          } ## End found a suitable qos
      } ## End iterating over every qos

        ## If we get here, then there is no place to immediately queue the job
        ## because there are already jobs queued, so just pick a qos which is big enough.
        my $found_qos2 = 0;
      QOS2: for my $q (@qos) {
          my $qos_info = $avail_qos->{$q};
          if ($wanted_spec->{mem} <= $qos_info->{max_job_mem} &&
              $wanted_spec->{cpu} <= $qos_info->{max_job_cpu} &&
              $wanted_spec->{gpu} <= $qos_info->{max_job_gpu} &&
              $wanted_spec->{walltime} <= $qos_info->{max_mem}) {
              print "Found qos in second pass: $q wanted: $wanted_spec->{mem} $wanted_spec->{cpu} $wanted_spec->{walltime}\n";
              $found_qos2++;
              $qos_info->{used_mem} = $qos_info->{used_mem} + $wanted_spec->{mem};
              $qos_info->{used_cpu} = $qos_info->{used_cpu} + $wanted_spec->{cpu};
              $qos_info->{used_gpu} = $qos_info->{used_gpu} + $wanted_spec->{gpu};
              $qos_info->{used_hours} = $qos_info->{used_hours} + $wanted_spec->{walltime};
              $chosen_account = $account;
              $chosen_qos = $q;
              last TOP;
          } ## End found a suitable qos
      } ## End iterating over every qos
    }
  }
    print "Got $chosen_qos\n";
    my $ret = {
        qos_info => $qos_info,
        choice => $chosen_qos,
    };
    return($ret);
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
        args => \%args,
        jprefix => '',
        jname => 'unknown',);
    my $sbatch = $class->Check_Sbatch();
    my $depends_prefix = '--dependency=afterok';
    ## For arguments to sbatch, start with the defaults in the constructor in $class
    ## then overwrite with any application specific requests from %args
    my $sbatch_log = 'outputs/log.txt';

    my $depends_string = "";
    if ($options->{jdepends}) {
        $depends_string = qq"${depends_prefix}:$options->{jdepends}";
    }
    if (!defined($options->{jname})) {
        $options->{jname} = $class->Get_Job_Name();
    }
    my $jname = qq"$options->{jprefix}$options->{jname}";
    my $finished_file = qq"outputs/logs/${jname}.finished";
    my $job = {};
    foreach my $k (keys %args) {
        next if ($k eq 'jstring');
        next if ($k eq 'comment');
        $job->{$k} = $args{$k};
    }
    my @wanted_vars = ('basedir', 'depends_string', 'input',
                       'jcpus', 'jmem', 'jname', 'jqueue', 'jwalltime', 'output');
    foreach my $w (@wanted_vars) {
        $job->{$w} = $options->{$w} if (!defined($job->{$w}));
    }
    if ($options->{restart} && -e $finished_file) {
        print "The restart option is on, and this job appears to have finished.\n";
        return($job);
    }

    my $script_file = qq"$options->{basedir}/scripts/$options->{jprefix}$options->{jname}.sh";
    my $sbatch_cmd_line = qq"${sbatch} ${depends_string} ${script_file}";
    my $mycwd = getcwd();
    make_path("$options->{logdir}", {verbose => 0}) unless (-r qq"$options->{logdir}");
    make_path("$options->{basedir}/scripts", {verbose => 0}) unless (-r qq"$options->{basedir}/scripts");
    my $script_base = basename($script_file);

    ## Remove the need for two functions that do the same thing except one for perl and one for bash
    if ($options->{language} eq 'perl') {
        my $perl_base = qq"$options->{basedir}/scripts";
        my $perl_file = qq"${perl_base}/$options->{jprefix}$options->{jname}.pl";
        if (defined($options->{output_dir})) {
            $perl_base = $options->{output_dir};
        } elsif (defined($options->{output})) {
            $perl_base = dirname($options->{output});
        }
        my $perl_stderr = qq"${perl_base}/$options->{jprefix}$options->{jname}.stderr";
        my $perl_stdout = qq"${perl_base}/$options->{jprefix}$options->{jname}.stdout";
        my $perl_start = qq?#!/usr/bin/env perl
use strict;
use FileHandle;
use Bio::Adventure;
my \$out = FileHandle->new(">>${sbatch_log}");
my \$d = qx'date';
print \$out "## Started $script_file at \${d}";
chdir("$options->{basedir}");
my \$h = Bio::Adventure->new();
?;
        if (defined($parent->{option_file})) {
            $perl_start .= qq!
## Pull options from the option file: $parent->{option_file}
my \$loaded = \$h->Load_Vars(input => '$parent->{option_file}');
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
        my $total_perl_string = qq"${perl_start}\n";
        $total_perl_string .= qq"$options->{comment}\n" if ($options->{comment});
        $total_perl_string .= qq"$options->{prescript}\n" if ($options->{prescript});
        $total_perl_string .= qq"$options->{jstring}\n" if ($options->{jstring});
        $total_perl_string .= qq"${perl_end}\n";

        my $perl_script = FileHandle->new(">${perl_file}");
        print $perl_script $total_perl_string;
        $perl_script->close();
        chmod(0775, $perl_file);
        $options->{jstring} = qq"
${perl_file} \\
  2>${perl_stderr} \\
  1>${perl_stdout}
";
    } ## End extra processing for submission of a perl script (perhaps not needed for slurm?

    my $nice_string = '';
    $nice_string = qq"--nice=$options->{jnice}" if (defined($options->{jnice}));

    my $module_string = '';
    if (defined($options->{modules}) && scalar(@{$options->{modules}} > 0)) {
        $module_string = 'module add';
        for my $m (@{$options->{modules}}) {
            $module_string .= qq" ${m}";
        }
    }
    my $array_string = '';
    $array_string = qq"#SBATCH --array=$options->{array_string}" if ($options->{array_string});

    my $script_start = qq?#!$options->{shell}
#SBATCH --export=ALL --requeue --mail-type=NONE --open-mode=append
#SBATCH --chdir=$options->{basedir}
#SBATCH --partition=$options->{jpartition}
#SBATCH --qos=$options->{jqueue} ${nice_string}
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=$options->{jcpus}
#SBATCH --time=$options->{jwalltime}
#SBATCH --job-name=${jname}
#SBATCH --mem=$options->{jmem}G
#SBATCH --output=${sbatch_log}.sbatch
${array_string}
${module_string}
set -o errexit
set -o errtrace
set -o pipefail
export LESS='$ENV{LESS}'
echo "## Started ${script_file} at \$(date) on \$(hostname) with id \${SLURM_JOBID}." >> ${sbatch_log}

?;

    my $script_end = qq!
## The following lines give status codes and some logging
echo "Job status: \$? " >> ${sbatch_log}
minutes_used=\$(( SECONDS / 60 ))
echo "  \$(hostname) Finished \${SLURM_JOBID} ${script_base} at \$(date), it took \${minutes_used} minutes." >> ${sbatch_log}
if [[ -x "\$(command -v sstat)" && \! -z "\${SLURM_JOBID}" ]]; then
  echo "  walltime used by \${SLURM_JOBID} was: \${minutes_used:-null} minutes." >> ${sbatch_log}
  maxmem=\$(sstat -n -P --format=MaxVMSize -j "\${SLURM_JOBID}")
  echo "  maximum memory used by \${SLURM_JOBID} was: \${maxmem:-null}." >> ${sbatch_log}
  echo "" >> ${sbatch_log}
fi
touch ${finished_file}
!;

    my $total_script_string = '';
    $total_script_string .= qq"${script_start}\n";
    $total_script_string .= qq"$options->{comment}\n" if ($options->{comment});
    $total_script_string .= qq"$options->{prescript}\n" if ($options->{prescript});
    $total_script_string .= qq"$options->{jstring}\n" if ($options->{jstring});
    $total_script_string .= qq"$options->{postscript}\n" if ($options->{postscript});
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
    my $sbatch_pid = open($handle, qq"${sbatch_cmd_line} |");
    while (my $line = <$handle>) {
        chomp($line);
        next unless($line =~ /Submitted/);
        $job_id = $line;
        $job_id =~ s/^.*Submitted batch job (\d+)/$1/g;
    }
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
    $job->{log} = $sbatch_log;
    $job->{job_id} = $job_id;
    $job->{pid} = $sbatch_pid;
    $job->{script_file} = $script_file;

    if ($options->{verbose}) {
        use Data::Dumper;
        print Dumper $job;
    }
    ## Take a moment to reset the shell and language
    ##my $reset = Bio::Adventure::Reset_Vars($class);
    ##$reset = Bio::Adventure::Reset_Vars($parent);
    return($job);
}

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

L<Bio::Adventure::Torque> L<Bio::Adventure::Local> L<sbatch>

=cut

1;
