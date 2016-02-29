package CYOA;

sub Test_Job {
    my $me = shift;
    my %args = @_;
    my $job_string = qq"/bin/true";
    my $job = $me->Qsub(job_name => 'test',
                        job_string => $job_string,
                        comment => '## hi!',
        );
}

=head1 NAME
    CYOA::Qsub - Submit jobs to the torque cluster.

=head1 SYNOPSIS

    use CYOA;
    my $hpgl = new CYOA;
    $hpgl->Cutadapt();

=head2 Methods

=over 4

=item C<Qsub>

    $hpgl->Qsub(); invokes qsub with (hopefully) appropriate
    parameters for various jobs on our Torque cluster.

=cut
sub Qsub {
    my $me = shift;
    my %args = @_;
    my $queue = $args{queue};

    ## Figure out if this job depends on anything
    ## The source of a dependency is usually via a parameter
    ## but might also be passed as a command line argument
    ## In the former case, it will be via $args{}, while
    ## in the latter case, it is via $me->{depends}
    ## Remember that dependency lists are colon separated.
    my $depends = $args{depends} if ($args{depends});
    $depends = $me->{depends} if (!$depends and $me->{depends});
    ## Set up the job name.  This is only allowed to come from %args.
    my $job_name;
    my $name_suffix = substr($me->{basename}, 0, 8);
    $name_suffix =~ s/_|\-//g;
    if ($args{job_name}) {
        my $name_prefix = $args{job_name};
        $job_name = qq"${name_prefix}-${name_suffix}";
    } else {
        $job_name = $name_suffix;
    }

    ## For arguments to qsub, start with the defaults in the constructor in $me
    ## then overwrite with any application specific requests from %args
    my $qsub_args = $me->{qsub_args};
    $qsub_args = $args{qsub_args} if ($args{qsub_args});
    my $qsub_queue = $me->{qsub_queue};
    $qsub_queue = $args{qsub_queue} if ($args{qsub_queue});
    my $qsub_shell = $me->{qsub_shell};
    $qsub_shell = $args{qsub_shell} if ($args{qsub_shell});
    my $qsub_mem = $me->{qsub_mem};
    $qsub_mem = $args{qsub_mem} if ($args{qsub_mem});
    my $qsub_wall = $me->{qsub_wall};
    $qsub_wall = $args{qsub_wall} if ($args{qsub_wall});
    my $qsub_cpus = $me->{qsub_cpus};
    $qsub_cpus = $args{qsub_cpus} if ($args{qsub_cpus});
    my $qsub_logdir = $me->{qsub_logdir};
    $qsub_logdir = $args{qsub_logdir} if ($args{qsub_logdir});
    my $qsub_log = qq"${qsub_logdir}/${job_name}.qsubout";
    my $job_output = "";
    $job_output = $args{output} if (defined($args{output}));
    my $job_input = "";
    $job_input = $args{input} if (defined($args{input}));
    my $depends_string = $me->{qsub_depends};
    if (defined($args{depends_type})) {
        if ($args{depends_type} eq 'array') {
            $depends_string = $me->{qsub_dependsarray};
        }
    }
    $depends_string .= $depends if (defined($depends));

    my $jobid_name = qq"$me->{basename}-$args{job_name}";
    my $script_file = qq"$me->{basedir}/scripts/${job_name}.sh";
    my $mycwd = cwd();
    make_path("$me->{basedir}/outputs/status", {verbose => 0}) unless (-r qq"$me->{basedir}/outputs/status");
    make_path("$qsub_logdir", {verbose => 0}) unless (-r qq"$qsub_logdir");
    make_path("$me->{basedir}/sequences", {verbose => 0}) unless (-r qq"$me->{basedir}/sequences");
    make_path("$me->{basedir}/scripts", {verbose => 0}) unless (-r qq"$me->{basedir}/scripts");
    my $script_base = basename($script_file);
    my $script_start = qq?#PBS -V -S ${qsub_shell} -q ${qsub_queue}
#PBS -d $me->{basedir}
#PBS -N ${job_name} -l mem=${qsub_mem}gb -l walltime=${qsub_wall} -l ncpus=${qsub_cpus}
#PBS -o ${qsub_log} ${qsub_args}
echo "####Started ${script_file} at \$(date)" >> outputs/log.txt
cd $me->{basedir} || exit
?;
    my $script_end = qq!## The following lines give status codes and some logging
echo \$? > outputs/status/${job_name}.status
echo "###Finished \${PBS_JODID} $script_base at \$(date), it took \$(( \$SECONDS / 60 )) minutes." >> outputs/log.txt
!;
    ## It turns out that if a job was an array (-t) job, then the following does not work because
    ## It doesn't get filled in properly by qstat -f...

    ## The following lines used to be in the shell script postscript
    ## Copying the full job into the log is confusing, removing this at least temporarily.
    ##echo "####This job consisted of the following:" >> outputs/log.txt
    ##cat "\$0" >> outputs/log.txt

    if ($me->{pbs}) {
        $script_end .= qq!
walltime=\$(qstat -f -t \"\${PBS_JOBID}\" | grep 'resources_used.walltime' | awk -F ' = ' '{print \$2}')
echo "####PBS walltime used by \${PBS_JOBID} was: \${walltime:-null}" >> outputs/log.txt
mem=\$(qstat -f -t | grep \"\${PBS_JOBID}\" | grep 'resources_used.mem' | awk -F ' = ' '{print \$2}')
echo "####PBS memory used by \${PBS_JOBID} was: \${mem:-null}" >> outputs/log.txt
vmmemory=\$(qstat -f -t \"\${PBS_JOBID}\" | grep 'resources_used.vmem' | awk -F ' = ' '{print \$2}')
echo "####PBS vmemory used by \${PBS_JOBID} was: \${vmmemory:-null}" >> outputs/log.txt
cputime=\$(qstat -f -t \"\${PBS_JOBID}\" | grep 'resources_used.cput' | awk -F ' = ' '{print \$2}')
echo "####PBS cputime used by \${PBS_JOBID} was: \${cputime:-null}" >> outputs/log.txt
##qstat -f -t \${PBS_JOBID} >> outputs/log.txt
!;
    }

    if ($me->{verbose}) {
        print qq"The job is:
$args{job_string}
";
    }

    my $total_script_string = "";
    $total_script_string .= "$script_start\n";
    $total_script_string .= "$args{comment}\n" if ($args{comment});
    $total_script_string .= "$args{prescript}\n" if ($args{prescript});
    $total_script_string .= "$args{job_string}\n" if ($args{job_string});
    if ($args{postscript}) {
        $total_script_string .= qq!if [ \$? == "0" ]; then
   $args{postscript}
fi
!;
    }
    $total_script_string .= "$script_end\n";

    my $script = new FileHandle;
    $script->open(">$script_file");
    print $script $total_script_string;
    $script->close();

    my $qsub = qq"${qsub_shell} $script_file";
    if ($me->{pbs}) {
        $qsub = qq"qsub -W $depends_string $script_file";
    }
    my $job_id = undef;
    my $qsub_pid = open(my $fh, "-|", $qsub);
    while(my $line = <$fh>) {
      chomp($line);
      $job_id = $line;
    }
    close($fh);
    if (!defined($job_id)) {
        warn("The job id did not get defined.  qsub likely failed.");
        return(undef);
    }
    my @jobid_list = split(/\./, $job_id);
    my $short_jobid = shift(@jobid_list);

    print "Starting a new job: ${short_jobid} ${jobid_name}";
    if ($depends) {
        my @short_dep = split(/\./, $depends);
        my $shortened_dep = shift(@short_dep);
        print ", depending on ${shortened_dep}.";
    }
    print "\n";

    my $job = { id => $jobid_name,
                submitter => $qsub,
                mem => $qsub_mem,
                walltime => $qsub_wall,
                cpus => $qsub_cpus,
                jobname => $job_name,
                log => $qsub_log,
                depends_string => $depends_string,
                queue => $qsub_queue,
                qsub_args => $qsub_args,
                basedir => $me->{basedir},
                pbs_id => $job_id,
                script_file => $script_file,
                script_start => $script_start,
                script_body => $args{job_string},
                output => $job_output,
                input => $job_input,
    };
    return($job);
}

sub Qsub_Arbitrary {
    my $me = shift;
    my %args = @_;
    my $queue = $args{queue};
    my $depends = $me->{depends};
    if (!$depends) {
        $depends = $args{depends};
    }
    $me->Check_Options(['job_string']);
    my $job_name = $me->{job_name};
    $job_name = 'arbitrary' if (!$job_name);
    my $job = $me->Qsub(job_name => $job_name,
                        depends => $depends,
                        job_string => $me->{job_string},
        );
    return($job);
}


sub Qsub_Perl {
    my $me = shift;
    my %args = @_;
    my $queue = $args{queue};
    if ($args{depends}) {
        print "This job depends on $args{depends}\n";
    }

    my $job_name;
    my $name_suffix = substr($me->{basename}, 0, 8);
    if ($args{job_name}) {
        my $name_prefix = $args{job_name};
        $job_name = qq"${name_prefix}-${name_suffix}";
    } else {
        $job_name = $name_suffix;
    }
    my $qsub_args = $me->{qsub_args};
    $qsub_args = $args{qsub_args} if ($args{qsub_args});
    my $qsub_queue = $me->{qsub_queue};
    $qsub_queue = $args{qsub_queue} if ($args{qsub_queue});
    my $qsub_shell = '/usr/bin/env perl';
    my $qsub_mem = $me->{qsub_mem};
    $qsub_mem = $args{qsub_mem} if ($args{qsub_mem});
    my $qsub_wall = $me->{qsub_wall};
    $qsub_wall = $args{qsub_wall} if ($args{qsub_wall});
    my $qsub_cpus = $me->{qsub_cpus};
    $qsub_cpus = $args{qsub_cpus} if ($args{qsub_cpus});
    my $qsub_logdir = $me->{qsub_logdir};
    $qsub_logdir = $args{qsub_logdir} if ($args{qsub_logdir});
    my $qsub_log = qq"${qsub_logdir}/${job_name}.qsubout";
    my $depends_string = $me->{qsub_depends};
    $depends_string .= $args{depends} if (defined($args{depends}));
    my $job_output = "";
    $job_output = $args{output} if (defined($args{output}));
    my $job_input = "";
    $job_input = $args{input} if (defined($args{input}));
    my $script_file = qq"$me->{basedir}/scripts/${job_name}.pl";
    my $mycwd = cwd();
    make_path("$me->{basedir}/outputs", {verbose => 0}) unless (-r qq"$me->{basedir}/outputs");
    make_path("$me->{basedir}/sequences", {verbose => 0}) unless (-r qq"$me->{basedir}/sequences");
    make_path("$me->{basedir}/scripts", {verbose => 0}) unless (-r qq"$me->{basedir}/scripts");
    my $script_start = qq?#!${qsub_shell}
#PBS -V -S ${qsub_shell} -q ${qsub_queue} -d $me->{basedir}
#PBS -N $job_name -l mem=${qsub_mem}gb -l walltime=$qsub_wall -l ncpus=$qsub_cpus
#PBS -o $qsub_log $qsub_args
open(OUT, ">>outputs/log.txt");
my \$d = qx'date';
print OUT "###Started $script_file at \${d}";
chdir("$me->{basedir}");
?;
    my $script_end = qq!## The following lines give status codes and some logging
print OUT "####This job consisted of the following:";
!;
    print "The job is:
$args{job_string}" if ($me->{debug});

    my $total_script_string = "";
    $total_script_string .= "$script_start\n";
    $total_script_string .= "$args{comment}\n" if ($args{comment});
    $total_script_string .= "$args{prescript}\n" if ($args{prescript});
    $total_script_string .= "$args{job_string}\n" if ($args{job_string});
    $total_script_string .= "$script_end\n";

    my $script = new FileHandle;
    $script->open(">$script_file");
    print $script $total_script_string;
    $script->close();

    my $qsub = qq"$qsub_shell $script_file";
    if ($me->{pbs}) {
        $qsub = qq"qsub -W $depends_string $script_file";
    }
    my $job_id;
    ##sleep(1);
    my $qsub_pid = open(my $fh, "-|", $qsub);
    while(my $line = <$fh>) {
      chomp($line);
      $job_id = $line;
    }
    close($fh);
    my $jobid_name = qq"$me->{basename}-$args{job_name}";
    print "Starting a new job: $jobid_name\nEOJ\n\n";
    my $job = { id => $jobid_name,
                submitter => $qsub,
                mem => $qsub_mem,
                walltime => $qsub_wall,
                cpus => $qsub_cpus,
                jobname => $job_name,
                log => $qsub_log,
                depends_string => $depends_string,
                queue => $qsub_queue,
                qsub_args => $qsub_args,
                basedir => $me->{basedir},
                pbs_id => $job_id,
                script_file => $script_file,
                script_start => $script_start,
                script_body => $args{job_string},
                output => $job_output,
                input => $job_input,
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
