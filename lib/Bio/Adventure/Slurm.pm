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
use Data::Dumper;
use File::Basename qw"basename dirname";
use File::Path qw"make_path remove_tree";
use File::Which qw"which";
use IO::Handle;
use POSIX qw"floor";

## List of accounts, partitions, qos for this user
has accounts => (is => 'rw', default => undef);
has partitions => (is => 'rw', default => undef);
has partition_names => (is => 'rw', default => '');
has qos_names => (is => 'rw', default => '');
## hash of the associations for this person
has association_data => (is => 'rw', default => undef);
has allowed_qos => (is => 'rw', default => undef);

## The name of the chosen account, cluster, partition and qos for this job.
has chosen_account => (is => 'rw', default => '');
has chosen_cluster => (is => 'rw', default => '');
has chosen_partition => (is => 'rw', default => '');
has chosen_qos => (is => 'rw', default => '');

## hash of the qos and their attributes.
has qos_data => (is => 'rw', default => undef);
## Any attributes here which are also in Adventure.pm get set to the values
## from Adventure.pm which is confusing to me.
##has jprefix => (is => 'rw', default => '00');
##has jname => (is => 'rw', default => 'unknown');
##has language => (is => 'rw', default => 'bash');
## Location of the sbatch executable.

##has sbatch => (is => 'rw', default => 'sbatch');
has slurm_test => (is => 'rw', default => 'testing_slurm_instance_variable_value');
## Current usage stats
has usage => (is => 'rw', default => undef);

sub BUILD {
    my ($class, $args) = @_;
    if (!defined($class->{qos_data})) {
        my $qos_data = Get_QOS();
        my @qos_names = sort keys %{$qos_data};
        $class->{qos_names} = \@qos_names;
        $class->{qos_data} = $qos_data;
    }

    if (!defined($class->{usage})) {
        $class->{usage} = Get_Usage();
    }

    ## Give me the set of partitions/clusters/qos available to my current user.
    if (!defined($class->{assciation_data})) {
        my $qos_names = $class->{qos_names};
        my $qos = $class->{qos_data};
        my $assoc = $class->Get_Associations();
        $class->{association_data} = $assoc;
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

        if (!defined($class->{partitions})) {
            my $partition_info = $class->Get_Partitions();
            $class->{partitions} = $partition_info;
            my $partition_names_string = '';
            my @partition_names = keys %{$partition_info};
            for my $k (@partition_names) {
                $partition_names_string .= qq"${k},";
            }
            $class->{partitions} = $partition_info;
            $class->{partition_names} = $partition_names_string;
        }

        if (!defined($class->{allowed_qos})) {
            print Dumper $class->{association_data};
            $class->{allowed_qos} = $class->Get_My_QOS(
                partitions => $class->{partitions},
                qos => $class->{qos_data},
                assoc => $class->{association_data});
        }

        $class->{cluster} = $cluster;
        $class->{accounts} = \@accounts;
        $class->{association_data} = $assoc;
    }
    return($class);
}

sub Get_My_QOS {
    my ($class, %args) = @_;
    ## I want this datastructure to be
    ## my_qos->{partition}->{cluster}->{account} = list_of_qoses
    my $my_qos = {};
    my $all_partitions = $args{partitions};
    my $my_assoc = $args{assoc};
    my $all_qos = $args{qos};
    ## The partition hash is keyed off the partition names, the
    ## association hash is keyed off slurm accounts.
    ## Thus, lets iterate over the partitions in the set of all partitions.
  PART: for my $part (keys %{$all_partitions}) {
      my $part_info = $all_partitions->{$part};
      my $allowed_accounts = $part_info->{allowaccounts};
      my @allowed_qos = split(/,/, $part_info->{allowqos});
      my @partition_allowed_accounts = split(/,/, $allowed_accounts);
      ## Then over the set of clusters to which I am associated.
    CL: for my $cluster (keys %{$my_assoc}) {
        my $account_info = $my_assoc->{$cluster};
        ## Followed by the set of accounts to which I am associated.
      ACC: for my $account (keys %{$account_info}) {
          ## Add a check to see if this account is allowed on this partition.
          ## This is separate from checking the set of QOSes allowed to the account.
          my $allowed_partition = 0;
          for my $part_allowed_account (@partition_allowed_accounts) {
              if ($account eq $part_allowed_account) {
                  $allowed_partition = 1;
              }
          }
          next ACC unless ($allowed_partition);
          ## We now have two sets of partitions/accounts->allowedqos (@allowed_qos)
          ## and accounts->partitions (@personal_qos), so we can use a
          ## key lookup to get the shared qos between them.
          my @personal_qos = @{$my_assoc->{$cluster}->{$account}->{qos}};

          my %shared_qos = map { $_ => 1 } @allowed_qos;
          my @allowed_shared = grep { $shared_qos{$_} } @personal_qos;
          $my_qos->{$part}->{$cluster}->{$account} = \@allowed_shared;
      }
    }
  }
    print "We should now have the set of qos to which I am permitted access by account/partition
        partition        cluster        account    qos\n";
    print Dumper $my_qos;
    return($my_qos);
}


sub Check_Sbatch {
    my ($class, %args) = @_;
    my $path = which('sbatch');
    if (!defined($path)) {
        $path = which('bash');
    }
    return($path);
}

sub Choose_QOS {
    my ($class, %args) = @_;
    my $wanted_spec = $args{wanted_spec};
    my $current_usage = $args{current_usage};
    my $associations = $class->{association_data};
    my $allowed = $class->{allowed_qos};

    my $qos_info = $class->{qos_data};
    $wanted_spec->{walltime_hours} = Convert_to_Hours($wanted_spec->{walltime});
    my $chosen_cluster = '';
    my $chosen_account = '';
    my $chosen_partition = '';
    my $chosen_qos = '';
    my $found_qos = 0;
    my @all_qos = @{$class->{qos_names}};
    my @qos;
    my $skip_default = 1;
  PART: for my $partition (keys %{$allowed}) {
      my %inner_part = %{$allowed->{$partition}};
    CLUST: for my $cluster (keys %inner_part) {
        my %inner_clust = %{$allowed->{$partition}->{$cluster}};
      ACCOUNT: for my $account (keys %inner_clust) {
          my @allowed_qos = @{$inner_clust{$account}};
          ## print "TESTME: On cluster: $cluster account: $account we have qos: @qos\n";
          my $potential_qos = {};
          ## We need a little logic which says: if the only qos for a person is 'default'
          ## then that person is likely associated with all qoses.
          my @qos;
          if (defined($associations->{$cluster}->{$account}->{qos})) {
              @qos = @{$associations->{$cluster}->{$account}->{qos}};
          } else {
              next ACCOUNT;
          }
          ## It appears that the default qos is a weird special case in
          ## a bunch of ways: 1.  It fills in the values for all other queues,
          ## 2.  at least at umiacs, it may not be used by all of one's associated accounts
          ## e.g. if I use it on nexus I get a failed job allocation with:
          ## 'Job's_account_not_permitted_to_use_this_partition_(tron_allows_nexus_not_cbcb)'
          ## So, for the moment, let us explicitly not use default.
          ## 3.  It appears that if a user is associated with only 'default', then that means
          ## the user is actually associated with every qos...
          ## Thus the check on the length of the qos array here, if it is comprised of only
          ## ('default'), then we will try and choose from every qos provided by
          ## sacctmgr show qos.
          if (scalar(@qos) == 1 && $qos[0] eq 'default') {
              @qos = @all_qos;
          }
        QOS: for my $q (@qos) {
            next QOS if ($q eq 'default' && $skip_default);
            my $info = $qos_info->{$q};
            ## As currently written, this info->{used_mem} is incorrect because it _should_ be
            ## getting that information from %current_usage, once we add those up and compare
            ## to the maximum, then we should be able to match up to a correct qos
            ## Ideally, I would like to have some knowledge about which qos/nodes are idle
            ## and choose those, but I think I don't sufficiently understand the cluster to do that yet...

            ## We need to make separate calls on each of these criteria because some clusters do not
            ## provide a maximum cpus/mem/gpus on a per-qos basis.

            my $stringent_mem = 0;
            my $stringent_cpu = 0;
            my $stringent_gpu = 0;
            my $stringent_hours = 0;
            if ($info->{max_job_mem}) {
                $stringent_mem = $info->{max_job_mem} + $info->{used_mem};
                # print "Checking $q stringent $stringent_mem vs wanted $wanted_spec->{mem} memory.\n";
                if ($wanted_spec->{mem} > $stringent_mem) {
                    # print "Stringent: This job wants $wanted_spec->{mem} which is more than ${stringent_mem}, not using qos ${q}\n";
                    next QOS;
                }
            }
            if ($info->{max_job_cpu}) {
                $stringent_cpu = $info->{max_job_cpu} + $info->{used_cpu};
                if ($wanted_spec->{cpu} > $stringent_cpu) {
                    # print "Stringent: This job wants $wanted_spec->{cpu} which is more than ${stringent_cpu}, not using qos ${q}\n";
                    next QOS;
                }
            }
            if ($info->{max_job_gpu}) {
                $stringent_gpu = $info->{max_job_gpu} + $info->{used_gpu};
                if ($wanted_spec->{gpu} > $stringent_gpu) {
                    ## print "Stringent: This job wants $wanted_spec->{cpu} which is more than ${stringent_cpu}, not using qos ${q}\n";
                    next QOS;
                }
            }
            if ($info->{max_hours}) {
                $stringent_hours = $info->{max_hours} + $info->{used_hours};
                # print "Checking $q stringent $stringent_hours vs wanted $wanted_spec->{walltime_hours} hours.\n";
                if ($wanted_spec->{walltime_hours} > $stringent_hours) {
                    print "Stringent: This job wants $wanted_spec->{walltime_hours} which is more than ${stringent_hours}, not using qos ${q}\n";
                    next QOS;
                }
            }
            my $potential_metrics = {
                delta_cpu => $info->{max_job_cpu} - $info->{used_cpu},
                delta_gpu => $info->{max_job_gpu} - $info->{used_gpu},
                delta_mem => $info->{max_job_mem} - $info->{used_mem},
                delta_hours => $info->{max_hours} - $info->{used_hours},
            };
            $potential_qos->{$q} = $potential_metrics;
          } ## End iterating over stringent qos.
          my $num_potential = scalar(keys(%{$potential_qos}));
          $found_qos = $num_potential;
          if ($num_potential) {
              $chosen_qos = Choose_Among_Potential_QOS($potential_qos);
          }

          ## print "TESTME: QOS $chosen_qos has the closest memory requirements.\n";
          $qos_info->{$chosen_qos}->{used_mem} = $qos_info->{$chosen_qos}->{used_mem} + $wanted_spec->{mem};
          $qos_info->{$chosen_qos}->{used_cpu} = $qos_info->{$chosen_qos}->{used_cpu} + $wanted_spec->{cpu};
          $qos_info->{$chosen_qos}->{used_gpu} = $qos_info->{$chosen_qos}->{used_gpu} + $wanted_spec->{gpu};
          $qos_info->{$chosen_qos}->{used_hours} = $qos_info->{$chosen_qos}->{used_hours} + $wanted_spec->{walltime_hours};
          $chosen_account = $account;
          $chosen_cluster = $cluster;
          $chosen_partition = $partition;
          last PART;
      } ## End iterating over accounts

        unless ($found_qos) {
          ACCOUNT2: for my $account (keys %{$associations->{$cluster}}) {
              ## print "TESTME2: On cluster: $cluster account: $account we have qos: @qos\n";
              my $potential_qos = {};
              ## If we get here, then there is no place to immediately queue the job
              ## because there are already jobs queued, so just pick a qos which is big enough.
              my $found_qos2 = 0;
            QOS2: for my $q (@qos) {
                next QOS2 if ($q eq 'default' && $skip_default);
                print "mem: $wanted_spec->{mem} vs $qos_info->{$q}->{max_job_mem}
cpu: $wanted_spec->{  cpu} vs $qos_info->{$q}->{max_job_cpu}
gpu: $wanted_spec->{gp     u} vs $qos_info->{$q}->{max_job_gpu}
time: $wanted_spec->{wallti   me_hours} vs $qos_info->{$q}->{max_hours}\n";
                my $info = $qos_info->{$q};
                if ($info->{max_job_mem}) {
                    if ($wanted_spec->{mem} > $info->{max_job_mem}) {
                        # print "This job wants $wanted_spec->{mem} which is more than $info->{max_job_mem}, not using qos ${q}\n";
                        next QOS2;
                    }
                }
                if ($info->{max_job_cpu}) {
                    if ($wanted_spec->{cpu} > $info->{max_job_cpu}) {
                        # print "This job wants $wanted_spec->{cpu} which is more than $info->{max_job_cpu}, not using qos ${q}\n";
                        next QOS2;
                    }
                }
                if ($info->{max_job_gpu}) {
                    if ($wanted_spec->{gpu} > $info->{max_job_gpu}) {
                        # print "This job wants $wanted_spec->{gpu} which is more than $info->{max_job_gpu}, not using qos ${q}\n";
                        next QOS2;
                    }
                }
                if ($info->{max_hours}) {
                    if ($wanted_spec->{walltime_hours} > $info->{max_hours}) {
                        # print "This job wants $wanted_spec->{walltime_hours} which is more than $info->{max_hours}, not using qos ${q}\n";
                        next QOS2;
                    }
                }
                my $potential_metrics = {
                    delta_cpu => $info->{max_job_cpu} - $info->{used_cpu},
                    delta_gpu => $info->{max_job_gpu} - $info->{used_gpu},
                    delta_mem => $info->{max_job_mem} - $info->{used_mem},
                    delta_hours => $info->{max_hours} - $info->{used_hours},
                };
                $potential_qos->{$q} = $potential_metrics;
            } ## End iterating over the non-stringent qos.
              my $num_potential = scalar(keys(%{$potential_qos}));
              if ($num_potential) {
                  $chosen_qos = Choose_Among_Potential_QOS($potential_qos);
                  $chosen_account = $account;
                  $chosen_cluster = $cluster;
              }

              print "Found qos in second pass: $chosen_qos\n";
              $found_qos2++;
              $qos_info->{used_mem} = $qos_info->{$chosen_qos}->{used_mem} + $wanted_spec->{mem};
              $qos_info->{used_cpu} = $qos_info->{$chosen_qos}->{used_cpu} + $wanted_spec->{cpu};
              $qos_info->{used_gpu} = $qos_info->{$chosen_qos}->{used_gpu} + $wanted_spec->{gpu};
              $qos_info->{used_hours} = $qos_info->{$chosen_qos}->{used_hours} + $wanted_spec->{walltime_hours};
              last TOP;
          }
      } ## End iterating over a second attempt of accounts.
    } ## End a second pass if we didn't find anything the first time.
  } ## End iterating over every association

    if ($chosen_qos eq '') {
        print "Something went wrong, no qos was chosen.
The job wanted $wanted_spec->{mem} mem, $wanted_spec->{cpu} cpu, and $wanted_spec->{walltime_hours} hours.
Setting it to scavenger.\n";
        $chosen_qos = 'scavenger';
        $chosen_account = 'scavenger';
  }

    print "Choose_QOS: Got $chosen_qos\n";
    my $ret = {
        qos_info => $qos_info,
        chosen_qos => $chosen_qos,
        chosen_partition => $chosen_partition,
        chosen_account => $chosen_account,
        chosen_cluster => $chosen_cluster,
    };
    return($ret);
}

sub Convert_to_Hours {
    my $string = shift;
    my $days = 0;
    my $hms = '00:00:00';
    if ($string =~ m/\-/) {
        ($days, $hms) = split(/\-/, $string);
    } else {
        $hms = $string;
    }

    my $hours = 1;
    my $min = 0;
    my $sec = 0;
    if ($hms =~ m/:/) {
        ($hours, $min, $sec) = split(/:/, $hms);
    } else {
        $hours = $hms;
    }
    my $final_hours = ($days * 24) + $hours;
    if ($sec > 0) {
        $min++;
    }
    if ($min > 0) {
        $final_hours++;
    }
    ## E.g. round up, effectively ignoring minutes and seconds.
    return($final_hours);
}

sub Convert_to_Walltime {
    my $hours = shift;
    my $walltime = qq"${hours}:00:00";
    if ($hours > 24) {
        my $days = floor($hours / 24);
        my $hours = $hours % 24;
        $walltime = qq"${days}-${hours}:00:00";
    }
    return($walltime);
}

sub Get_Associations {
    my $associations = {};
    my $assoc = FileHandle->new("sacctmgr -p show user $ENV{USER} --associations |");
  ASSOC: while (my $line = <$assoc>) {
      chomp $line;
      next ASSOC if ($line =~ /^User/);
      ## Interesting, on the cbcb cluster, every user has only 1 associated QOS (default).
      ## This breaks my assumption that show associations would have every qos with every user.
      my ($user, $def_account, $admin, $cluster, $account, $partition, $share, $priority,
          $max_jobs, $max_nodes, $max_cpus, $max_submit, $max_wall, $max_cpu_min, $qos_lst,
          $def_qos) = split(/\|/, $line);
      my @avail_qos = split(/,/, $qos_lst);
      my $inner_hash = {
          default_account => $def_account,
          admin_access => $admin,
          cluster => $cluster,
          account => $account,
          partition => $partition,
          share => $share,
          priority => $priority,
          max_jobs => $max_jobs,
          max_nodes => $max_nodes,
          max_cpus => $max_cpus,
          max_submit => $max_submit,
          max_wall => $max_wall,
          max_cpu_min => $max_cpu_min,
          qos => \@avail_qos,
          default_qos => $def_qos,
      };
      $associations->{$cluster}->{$account} = $inner_hash;
    }
    $assoc->close();
    return($associations);
}

sub Get_Partitions {
    my ($class, %args) = @_;
    my $partitions = {};
    my $count = 0;
    my $partition_info = {};
    my $current_partition = '';
    my $part = FileHandle->new("scontrol show partitions |");
  PART: while (my $line = <$part>) {
      chomp $line;
      $count++;
      ## Each new partition entry begins with partitionname.
      if ($line =~ /^PartitionName/) {
          my @name_info = split(/=/, $line);
          $current_partition = $name_info[1];
      } elsif ($line =~ /^$/) {
          ## A blank line ends a partition.
          my %current_partition_info = %{$partition_info};
          $partitions->{$current_partition} = \%current_partition_info;
          $partition_info = {};
      } else {
          $line =~ s/^\s+//g;
          ## Everything else is a set of name:value pairs separated by = and space.
          my @pairs = split(/\s+/, $line);
          for my $pair (@pairs) {
              my @name_value = split(/=/, $pair);
              my $num_elements = scalar(@name_value);
              if ($num_elements == 2) {
                  my $name = $name_value[0];
                  my $value = $name_value[1];
                  $name = lc($name);
                  $partition_info->{$name} = $value;
              } elsif ($num_elements > 2) {
                  my $name = shift @name_value;
                  my $tmp_line = $line;
                  my $remove_string = qq"${name}=";
                  $tmp_line =~ s/.*${remove_string}//g;
                  my @new_array = split(/,/, $tmp_line);
                  for my $n (@new_array) {
                      my ($inner_key, $inner_value) = split(/=/, $n);
                      my $new_key = qq"${name}_${inner_key}";
                      $partition_info->{$new_key} = $inner_value;
                  }
              } else {
                  print "I do not think we should have other element numbers: ${num_elements}\n";
                  sleep(10);
              }
          } ## End looking at the equal sign separated pairs.
      } ## End the else{}
  } ## End reading scontrol show partitions
    ## The names I expect to find in each inner hash are:
    ## allowgroups,allowaccounts, allowqos,allocnodes,default,qos,
    ## defaulttime,disablerootjobs,exclusiveuser,gracetime,hidden,
    ## maxnodes,maxtime,minnodes,lln,maxcpuspernode,nodes,legacy,
    ## priorityjobfactor,prioritytier,rootonly,reqresv,oversubscribe,
    ## overtimelimit,preemptmode,state,totalcpus,totalnodes,
    ## selecttypeparameters,jobdefaults,defmempercpu,maxmempernode,tres,
    ## hmm, tres and tresbillingweights are weird.
    $part->close();
    return($partitions);
}

sub Get_QOS {
    my $avail_qos = {};
    my $qos = FileHandle->new("sacctmgr -p show qos |");
    my $count = 0;
  QOS: while (my $line = <$qos>) {
      $count++;
      next QOS if ($count == 1);
      ## Note that the cbcb cluster uses MaxJobsPU to limit the number of concurrent jobs allowed to run per QOS
      ## and that it appears to not set a maximum number of cpus
      my ($name, $priority, $gracetime, $preempt, $preempt_exempt, $preempt_mode, $flags,
          $usage_thresh, $usage_factor, $group_tres, $group_tres_min, $group_tres_run_min,
          $group_jobs, $group_submit, $group_wall, $max_resources_per_job, $max_tres_per_node,
          $max_tres_min, $max_wall, $max_resources_per_user, $max_jobs_pu, $max_submit_pu, $max_tres_pa,
          $max_jobs_pa, $max_submit_pa, $min_tres) = split(/\|/, $line);
      ## print "TESTME: name: $name pr: $priority gr: $gracetime pr $preempt wa $max_wall\n";
      my $max_job_cpu = 0;
      my $max_job_gpu = 0;
      my $max_job_mem = 0;
      my $max_user_cpu = 0;
      my $max_user_gpu = 0;
      my $max_user_mem = 0;
      my $max_jobs = -1;
      my $max_hours = 0;

      if ($max_resources_per_job) {
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

      if ($max_resources_per_user) {
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

      if ($max_wall) {
          my $days = 0;
          my $hms = '';
          if ($max_wall =~ /\-/) {
              ($days, $hms) = split(/\-/, $max_wall);
          } else {
              $hms = $max_wall;
          }
          my $hours = 0;
          my $min = 0;
          my $sec = 0;
          if ($hms) {
              ($hours, $min, $sec) = split(/:/, $hms);
          }
          $days = 0 if (!defined($days));
          $hours = 1 if (!defined($hours));
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
      ## This is perhaps a little silly/redundant, given that I am immediately
      ## filling in this information from the default queue.
      $avail_qos->{$name}->{used_cpu} = 0 if (!defined($avail_qos->{$name}->{used_cpu}));
      $avail_qos->{$name}->{used_gpu} = 0 if (!defined($avail_qos->{$name}->{used_gpu}));
      $avail_qos->{$name}->{used_mem} = 0 if (!defined($avail_qos->{$name}->{used_mem}));
      $avail_qos->{$name}->{used_hours} = 0 if (!defined($avail_qos->{$name}->{used_hours}));
      $avail_qos->{$name}->{max_job_cpu} = 0 if (!defined($avail_qos->{$name}->{max_job_cpu}));
      $avail_qos->{$name}->{max_job_gpu} = 0 if (!defined($avail_qos->{$name}->{max_job_gpu}));
      $avail_qos->{$name}->{max_job_mem} = 0 if (!defined($avail_qos->{$name}->{max_job_mem}));
      $avail_qos->{$name}->{max_user_cpu} = 0 if (!defined($avail_qos->{$name}->{max_user_cpu}));
      $avail_qos->{$name}->{max_user_gpu} = 0 if (!defined($avail_qos->{$name}->{max_user_gpu}));
      $avail_qos->{$name}->{max_user_mem} = 0 if (!defined($avail_qos->{$name}->{max_user_mem}));
      $avail_qos->{$name}->{max_jobs} = 0 if (!defined($avail_qos->{$name}->{max_jobs}));
      $avail_qos->{$name}->{max_hours} = 0 if (!defined($avail_qos->{$name}->{max_hours}));

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

sub Get_Spec {
    my ($class, %args) = @_;
    my $options = $args{options};
    my $wanted = {
        mem => 1,
        cpu => 1,
        gpu => 0,
        walltime => 1
    };

    if (defined($options->{mem}) && defined($options->{jmem})) {
        print "Both mem and jmem are defined, that is confusing, using jmem.\n";
        $wanted->{mem} = $options->{jmem};
    } elsif (defined($options->{mem})) {
        print "Mem is defined, ideally this should be jmem.\n";
        $wanted->{mem} = $options->{mem};
    } elsif (defined($options->{jmem})) {
        $wanted->{mem} = $options->{jmem};
    } else {
        print "Neither mem nor jmem is defined, defaulting to 10G.\n";
        $wanted->{mem} = 10;
    }

    my $walltime_string = '00:40:00';
    if (defined($options->{walltime}) && defined($options->{jwalltime})) {
        print "Both walltime and jwalltime are defined, that is confusing, using jwalltime.\n";
        $walltime_string = $options->{jwalltime};
    } elsif (defined($options->{walltime})) {
        print "Walltime is defined, ideally this should be jwalltime.\n";
        $walltime_string = $options->{walltime};
    } elsif (defined($options->{jwalltime})) {
        $walltime_string = $options->{jwalltime};
    } else {
        print "Neither walltime nor jwalltime is defined, defaulting to 40 minutes.\n";
    }
    # print "About to convert to hours: $walltime_string\n";
    my $walltime_hours = Convert_to_Hours($walltime_string);
    $wanted->{walltime_hours} = $walltime_hours;
    $wanted->{walltime} = Convert_to_Walltime($walltime_hours);

    if (defined($options->{cpu}) && defined($options->{jcpu})) {
        print "Both cpu and jcpu are defined, that is confusing, using jcpu.\n";
        $wanted->{cpu} = $options->{jcpu};
    } elsif (defined($options->{cpu})) {
        print "Cpu is defined, ideally this should be jcpu.\n";
        $wanted->{cpu} = $options->{cpu};
    } elsif (defined($options->{jcpu})) {
        $wanted->{cpu} = $options->{jcpu};
    } else {
        print "Neither cpu nor jcpu is defined, defaulting to 1.\n";
        $wanted->{cpu} = 1;
    }

    if (defined($options->{gpu}) && defined($options->{jgpu})) {
        print "Both gpu and jgpu are defined, that is confusing, using jgpu.\n";
        $wanted->{gpu} = $options->{jgpu};
    } elsif (defined($options->{gpu})) {
        print "Gpu is defined, ideally this should be jgpu.\n";
        $wanted->{gpu} = $options->{gpu};
    } elsif (defined($options->{jgpu})) {
        $wanted->{gpu} = $options->{jgpu};
    } else {
        $wanted->{gpu} = 0;
    }

    return($wanted);
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
      ## FIXME: If we are getting a M, then we need to divide the memory by
      ## 1000 if we want to keep counting by gigs.
      $min_memory =~ s/M$//g;
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
          ## I think I would like this to print some information about failed jobs perhaps here?
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
    my $class_jprefix = $class->{jprefix};
    my $options = $parent->Get_Vars(
        args => \%args,
        jname => 'unknown',);
    my $sbatch = $class->Check_Sbatch();
    my $depends_prefix = '--dependency=afterok';
    ## For arguments to sbatch, start with the defaults in the constructor in $class
    ## then overwrite with any application specific requests from %args
    my $sbatch_log = 'outputs/log.txt';

    ## Get my current usage, which fills in a hash of:
    ## usage->{user}->{partition}->{account}->{qos}->{mem, cpu, jobs, running, queued, failed}
    ## This at least in theory could also hold the various 50+ pieces of information acquired from
    ## $(squeue --me -o '%all')
    my $usage = $class->Get_Usage();
    ## Get_Spec takes the putative requirements of the job, which reside in $options->{jmem, jcpu, jwalltime}
    ## and return a fully numeric specification of time, cpus, gpus, memory.
    ## I state 'fully numeric' because the walltime is written as a string and I convert it to hours
    ## with a minimum of 1.
    my $wanted = $class->Get_Spec(options => $options);
    ## Choose_QOS should take the current usage and spec for this job and find a combination of
    ## cluster/partition/account/qos which will work and return that information.
    ## with the caveat that scavenger is a weirdo special case.
    ## In order to do this, it pulls the user's associations from the constructor.
    ## Reminder, this is a hash of assoc->{cluster}->{account}->{qos} -- which may be a bad organizational
    ## choice because the partition is inside this innermost hash.  I don't yet fully understand
    ## slurm queueing and organization, so keep this in mind until I do...
    my $chosen_qos = $class->Choose_QOS(wanted_spec => $wanted,
                                        current_usage => $usage);

    $class->{chosen_cluster} = $chosen_qos->{chosen_cluster};
    $class->{chosen_partition} = $chosen_qos->{chosen_partition};
    $class->{chosen_account} = $chosen_qos->{chosen_account};
    $class->{chosen_qos} = $chosen_qos->{chosen_qos};

    my $depends_string = '';
    if ($options->{jdepends}) {
        $depends_string = qq"${depends_prefix}:$options->{jdepends}";
    }
    if (!defined($options->{jname})) {
        $options->{jname} = $class->Get_Job_Name();
    }
    if (!defined($options->{jprefix})) {
        $options->{jprefix} = '';
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
                       'jcpu', 'jmem', 'jname', 'jqueue', 'jwalltime', 'output');
    foreach my $w (@wanted_vars) {
        $job->{$w} = $options->{$w} if (!defined($job->{$w}));
    }
    ## Now fill in the job's qos etc from the slurm instance
    my $qos_string = '';
    if ($class->{chosen_qos}) {
        $qos_string = $class->{chosen_qos};
    } else {
        print "QOS is not defined, setting it to the empty string\n";
    }
    my $account_string = '';
    if ($class->{chosen_account}) {
        $account_string = $class->{chosen_account};
    } else {
        print "account is not defined, setting it to the empty string\n";
    }
    my $cluster_string = '';
    if ($class->{chosen_cluster}) {
        $cluster_string = $class->{chosen_cluster};
    } else {
        print "cluster is not defined, setting it to the empty string.\n";
    }
    my $partition_string = '';
    if ($class->{chosen_partition}) {
        $partition_string = $class->{chosen_partition};
    } else {
        ## Partition is often empty, I don't quite know why yet.
        print "partition is not defined, setting it to the empty string\n";
        $partition_string = '';
    }
    ##  Need to catch the special case of scavenger
    if ($qos_string eq 'scavenger') {
        $account_string = 'scavenger';
        $partition_string = 'scavenger';
    }
    if ($account_string eq 'scavenger') {
        $qos_string = 'scavenger';
        $partition_string = 'scavenger';
    }
    if ($partition_string eq 'scavenger') {
        $qos_string = 'scavenger';
        $account_string = 'scavenger';
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
    my $walltime_string = qq"$wanted->{walltime_hours}:00:00";
    my $mem_string = qq"$wanted->{mem}G";
    my $script_start = qq?#!$options->{shell}
#SBATCH --export=ALL --requeue --mail-type=NONE --open-mode=append
#SBATCH --chdir=$options->{basedir}
#SBATCH --job-name=${jname} ${nice_string}
#SBATCH --output=${sbatch_log}.sbatch
?;
    $script_start .= qq?#SBATCH --account=${account_string}\n? if ($account_string);
    $script_start .= qq?#SBATCH --partition=${partition_string}\n? if ($partition_string);
    $script_start .= qq?#SBATCH --qos=${qos_string}\n? if ($qos_string);
    ## FIXME: This should get smarter and be able to request multiple tasks and nodes.
    $script_start .= qq?#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=$wanted->{cpu}\n? if (defined($wanted->{cpu}));
    $script_start .= qq?#SBATCH --time=${walltime_string}\n? if (defined($wanted->{walltime}));
    $script_start .= qq?#SBATCH --mem=${mem_string}\n? if (defined(${mem_string}));
    $script_start .= qq"${array_string}\n" if ($array_string);
    $script_start .= qq?set -o errexit
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
  ## I am not sure why, but when I run a script in an interactive session, the maxmem variable
  ## gets set correctly everytime, but when it is run by another node, sometimes it does not.
  ## Lets try and figure that out...
  echo "TESTME: \${maxmem}"
  if [[ \! -z "\${maxmem}" ]]; then
    echo "  maximum memory used by \${SLURM_JOBID} was: \${maxmem}." >> ${sbatch_log}
  else
    echo "  The maximum memory did not get set for this job: \${SLURM_JOBID}." >> ${sbatch_log}
    sstat -P -j "\${SLURM_JOBID}" >> ${sbatch_log}
  fi
  echo "" >> ${sbatch_log}
fi
touch ${finished_file}
!;

    my $total_script_string = '';
    $total_script_string .= qq"${script_start}\n";
    $total_script_string .= qq"$options->{comment}\n" if ($options->{comment});
    $total_script_string .= qq"$options->{prescript}\n" if ($options->{prescript});
    ## The prescript contains the module() definition when needed.
    ## So put that after it.
    $total_script_string .= qq"${module_string}\n" if ($module_string);

    $total_script_string .= qq"$options->{jstring}\n" if ($options->{jstring});
    $total_script_string .= qq"$options->{postscript}\n" if ($options->{postscript});
    $total_script_string .= "${script_end}\n";

    my $script = FileHandle->new(">$script_file");
    if (!defined($script)) {
        die("Could not write the script: $script_file
  $!");
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
        warn("The job id did not get defined, submission likely failed.");
        return(undef);
    }
    sleep($options->{jsleep});
    my @jobid_list = split(/\./, $job_id);
    my $short_jobid = shift(@jobid_list);

    print "Starting a new job: ${short_jobid} $options->{jname} with
cluster: ${cluster_string} partition: ${partition_string} account: ${account_string} qos: ${qos_string}
";
    if ($options->{jdepends}) {
        print "This job depends on $options->{jdepends}.\n\n";
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
    my $cyoa = Bio::Adventure->new();
    my $options = $cyoa->Get_Vars(
        args => \%args,
        jmem => 10,
        jwalltime => 4,);

    my $jstring = qq?echo "Starting test job."
sleep 300
echo "Ending test job."
?;
    my $job = $cyoa->Submit(
        input => 'test',
        jstring => $jstring,
        jmem => $options->{jmem},
        jcpu => $options->{jcpu},
        jwalltime => $options->{jwalltime});
}

sub Choose_Among_Potential_QOS {
    my $choices = shift;
    my $favorites = {
        hours => '',
        delta_hours => 1e9,
        mem => '',
        delta_mem => 1e9,
        cpu => '',
        delta_cpu => 1e9,
        gpu => '',
        delta_gpu => 1e9,
    };
    for my $qos (keys %{$choices}) {
        my $info = $choices->{$qos};
        if ($info->{delta_hours} < $favorites->{delta_hours}) {
            $favorites->{hours} = $qos;
            $favorites->{delta_hours} = $info->{delta_hours};
        }
        if ($info->{delta_mem} < $favorites->{delta_mem}) {
            $favorites->{mem} = $qos;
            $favorites->{delta_mem} = $info->{delta_mem};
        }
        if ($info->{delta_cpu} < $favorites->{delta_cpu}) {
            $favorites->{cpu} = $qos;
            $favorites->{delta_cpu} = $info->{delta_cpu};
        }
        if ($info->{delta_gpu} < $favorites->{delta_gpu}) {
            $favorites->{gpu} = $qos;
            $favorites->{delta_gpu} = $info->{delta_gpu};
        }
    }
    ## I think what I want to do at some point is take the ratio of deltamem/totalmem and
    ## deltahours/totalhours and choose the one where that is smallest.
    ## but for now I will just choose the one with the tightest memory.
    my $choice = $favorites->{mem};
    return($choice);
}

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

L<Bio::Adventure::Torque> L<Bio::Adventure::Local> L<sbatch>

=cut

1;
