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

## List of accounts for this user
has accounts => (is => 'rw', default => undef);
## hash of the associations for this person
has association_data => (is => 'rw', default => undef);
## List of clusters available to this user
has cluster => (is => 'rw', default => undef);
## List of qos names visible to this user
has qos => (is => 'rw', default => undef);
## hash of the qos and their attributes.
has qos_data => (is => 'rw', default => undef);
has language => (is => 'rw', default => 'bash');
## Location of the sbatch executable.
has partitions => (is => 'rw', default => undef);
has sbatch => (is => 'rw', default => 'sbatch');
has slurm_test => (is => 'rw', default => 'testing_slurm_instance_variable_value');
## Current usage stats
has usage => (is => 'rw', default => undef);


sub BUILD {
    my ($class, $args) = @_;
    my $assoc = $class->Get_Associations();
    my $qos = $class->Get_QOS();
    $class->{qos_data} = $qos;
    my $partitions = $class->Get_Partitions();
    $class->{partitions} = $partitions;
    my @qos_names = sort keys %{$qos};
    $class->{qos} = \@qos_names;
    my $cluster = undef;
    my @accounts = ();
    my @clusters = keys %{$assoc};
    if (scalar(@clusters) > 0) {
        $cluster = $clusters[0];
        @accounts = keys %{$assoc->{$cluster}};
        ## Merge the qos information into the user's assocations
        ## In the hopes that this makes it easier to pick and choose queues
        ## So, when we wish to pull the current usage in a qos, we will
        ## ask for: $assoc->{$cluster}->{$account}->{$qos}->{used_mem} or
        ## whatever...  This does assume that the counters accross accounts
        ## are not shared, something which I have not yet tested.
        for my $iterate (keys @clusters) {
            for my $account (keys %{$assoc->{$iterate}}) {
                my @account_qos = @{$assoc->{$iterate}->{$account}->{qos}};
                for my $q (@account_qos) {
                    my %info = %{$qos->{$q}};
                    $assoc->{$iterate}->{$account}->{$q} = \%info;
                }
            }
        }
    } ## End checking for associations.

    $class->{cluster} = $cluster;
    $class->{accounts} = \@accounts;
    $class->{association_data} = $assoc;
}

sub Check_Sbatch {
    my ($class, %args) = @_;
    my $path = which('sbatch');
    if (!defined($path)) {
        $path = which('bash');
    }
    return($path);
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
    use Data::Dumper;
    print Dumper $chosen_qos;
    return($ret);
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
            partition => $partition,
            default_qos => $def_qos };
        $associations->{$cluster}->{$account} = $inner_hash;
    }
    $assoc->close();
    return($associations);
}

sub Get_Partitions {
    my $avail_part = {};
    my $part = FileHandle->new("sinfo -a --long -r |");
    my $count = 0;
  PART: while (my $line = <$part>) {
      $count++;
      next QOS if ($count < 2); ## First line is a timestamp, second is a header
      my ($part, $avail, $timelimit, $jobsize, $root,
          $over, $groups, $nodes, $state, $nodelist) = split(/\s+/, $line);
      my $default = 0;
      if ($part =~ /\*$/) {
          $default = 1;
      }
      $part =~ s/\*$//g;
      $avail_part->{$part} = {
          avail => $avail,
          timelimit => $timelimit,
          jobsize => $jobsize,
          root => $root,
          oversub => $over,
          groups => $groups,
          nodes => $nodes,
          state => $state,
          nodelist => $nodelist,
          default => $default,
      }
  }
    $part->close();
    return($avail_part);
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
    return($avail_qos);
}

sub Get_Usage {
    my $usage = FileHandle->new("squeue --me -o '%all' |");
    my $current = {};  ## Hash of the current usage by user.
    my $all_jobs = {};
    my $count = 0;
  USAGE: while (my $line = <$usage>) {
      chomp $line;
      next USAGE if ($line =~ /^ACCOUNT/);
      $count++;
      my ($account, $tres_per_node, $min_cpus, $min_tmp_disk, $end_time, $features, $group,
          $over_sub, $job_id, $name, $comment, $time_limit, $min_memory, $req_nodes, $command,
          $priority, $qos, $reason, $blank, $st, $user, $reservation, $wckey, $exclude_nodes,
          $nice, $sct, $job_id2, $exec_host, $cpus, $nodes, $dependency, $array_job_id, $group2,
          $sockets_per_node, $cores_per_socket, $threads_per_core, $array_task_id, $time_left,
          $time, $nodeslist, $contiguous, $partition, $priority2, $nodelist_reason, $start_time,
          $state, $uid, $submit_time, $licenses, $core_spec, $sched_nodes, $work_dir) = split(/\|/, $line);
      ## Some notes about what some of the above variables look like:
      ## tres_per_node: 'gres:gpu:3'
      ## I suspect we can do something fun with the priority information, priority2 looks
      ## like it would be easier to parse though.
      ## the group variable looks like it currently contains what I think of as user?
      $min_memory =~ s/G$//g;
      my $internal = {
          account => $account,
          array_job_id => $array_job_id,
          array_task_id => $array_task_id,
          blank => $blank,
          command => $command,
          comment => $comment,
          core_spec => $core_spec,
          contiguous => $contiguous,
          cores_per_socket => $cores_per_socket,
          cpus => $cpus,
          dependency => $dependency,
          end_time => $end_time,
          exclude_nodes => $exclude_nodes,
          exec_host => $exec_host,
          features => $features,
          group => $group,
          group2 => $group2,
          jobid2 => $job_id2,
          licenses => $licenses,
          min_cpus => $min_cpus,
          min_memory => $min_memory,
          min_tmp => $min_tmp_disk,
          name => $name,
          nice => $nice,
          nodelist_reason => $nodelist_reason,
          nodes => $nodes,
          nodes_list => $nodeslist,
          over_sub => $over_sub,
          partition => $partition,
          priority => $priority,
          priority2 => $priority2,
          qos => $qos,
          reason => $reason,
          req_nodes => $req_nodes,
          reservation => $reservation,
          sched_nodes => $sched_nodes,
          sct => $sct,
          sockets_per_node => $sockets_per_node,
          st => $st,
          start_time => $start_time,
          state => $state,
          submit_time => $submit_time,
          threads_per_core => $threads_per_core,
          time => $time,
          time_left => $time_left,
          time_limit => $time_limit,
          tres_per_node => $tres_per_node,
          uid => $uid,
          user => $user,
          wckey => $wckey,
          work_dir => $work_dir,
      };
      $all_jobs->{$job_id} = $internal;

      my $instance = {
          $partition => {
              $account => {
                  $qos => {
                      mem => $min_memory,
                      min_cpu => $min_cpus,
                      jobs => 1,
                  }, }, },
      };

      if ($state eq 'RUNNING') {
          $instance->{$partition}->{$account}->{$qos}->{running} = 1;
          $instance->{$partition}->{$account}->{$qos}->{queued} = 0;
          $instance->{$partition}->{$account}->{$qos}->{failed} = 0;
      } elsif ($state eq 'PENDING') {
          $instance->{$partition}->{$account}->{$qos}->{running} = 0;
          $instance->{$partition}->{$account}->{$qos}->{queued} = 1;
          $instance->{$partition}->{$account}->{$qos}->{failed} = 0;
      } elsif ($state eq 'COMPLETING') {
          $instance->{$partition}->{$account}->{$qos}->{running} = 1;
          $instance->{$partition}->{$account}->{$qos}->{queued} = 0;
          $instance->{$partition}->{$account}->{$qos}->{failed} = 0;
      } else {
          $instance->{$partition}->{$account}->{$qos}->{running} = 0;
          $instance->{$partition}->{$account}->{$qos}->{queued} = 0;
          $instance->{$partition}->{$account}->{$qos}->{failed} = 1;
      }
      if (!defined($current->{$user}->{$partition}->{$account}->{$qos})) {
          $current->{$user} = $instance;
      } else {
          $current->{$user}->{$partition}->{$account}->{$qos}->{running} += $instance->{$partition}->{$account}->{$qos}->{running};
          $current->{$user}->{$partition}->{$account}->{$qos}->{queued} += $instance->{$partition}->{$account}->{$qos}->{queued};
          $current->{$user}->{$partition}->{$account}->{$qos}->{failed} += $instance->{$partition}->{$account}->{$qos}->{failed};
          $current->{$user}->{$partition}->{$account}->{$qos}->{mem} += $instance->{$partition}->{$account}->{$qos}->{mem};
          $current->{$user}->{$partition}->{$account}->{$qos}->{min_cpu} += $instance->{$partition}->{$account}->{$qos}->{min_cpu};
          $current->{$user}->{$partition}->{$account}->{$qos}->{jobs} += 1;
      }
  }
    $usage->close();
    return($current);
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
        jname => 'unknown',);
    my $sbatch = $class->Check_Sbatch();
    my $depends_prefix = '--dependency=afterok';
    ## For arguments to sbatch, start with the defaults in the constructor in $class
    ## then overwrite with any application specific requests from %args
    my $sbatch_log = 'outputs/log.txt';

    my $wanted = {
        mem => $options->{jmem},
        cpu => $options->{jcpu},
        gpu => $options->{gpu},
        walltime => $options->{walltime},
    };
    my $chosen_qos = Choose_Spec(associations => $class->{association_data},
                                 qos_info => $class->{qos_data},
                                 wanted_spec => $wanted,);


    my $depends_string = '';
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
#SBATCH --job-name=${jname} ${nice_string}
#SBATCH --output=${sbatch_log}.sbatch
?;
    if (defined($chosen_qos->{partition})) {
        $script_start .= qq?
#SBATCH --partition=$chosen_qos->{partition}
?;
    }
    if (defined($chosen_qos->{qos})) {
        $script_start .= qq?
#SBATCH --partition=$chosen_qos->{qos}
?;
    }
    if (defined($chosen_qos->{cpu})) {
        $script_start .= qq?
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=$chosen_qos->{cpu}
?;
    }
    if (defined($chosen_qos->{time})) {
        $script_start .= qq?
#SBATCH --time=$chosen_qos->{time}
?;
    }
    if (defined($chosen_qos->{mem})) {
        $script_start .= qq?
#SBATCH --time=$chosen_qos->{mem}
?;
    }
    if (defined($chosen_qos->{account})) {
        $script_start .= qq?
#SBATCH --account=$chosen_qos->{account}
?;
    }
$script_start .= qq?
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

    ## Take a moment to reset the shell and language
    ##my $reset = Bio::Adventure::Reset_Vars($class);
    ##$reset = Bio::Adventure::Reset_Vars($parent);
    $parent->{language} = 'bash';
    return($job);
}

sub Test_Job {
    my ($class, %args) = @_;
    my $slurm = Bio::Adventure::Slurm->new();
    my $test = $slurm->{slurm_test};
    $slurm->Get_Usage();
    my $wanted = {
        mem => 40,
        cpu => 4,
        gpu => undef,
        walltime => '4:00:00',
    };
    use Data::Dumper;
    print "TESTME: association_data\n";
    print Dumper $slurm->{association_data};
    print "TESTME: qos_data\n";
    print Dumper $slurm->{qos_data};
    ##my $chosen_qos = Choose_Spec(associations => $slurm->{association_data},
    ##                             qos_info = $slurm->{qos_data},
    ##                             wanted_spec => $wanted,);

}

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

L<Bio::Adventure::Torque> L<Bio::Adventure::Local> L<sbatch>

=cut

1;
