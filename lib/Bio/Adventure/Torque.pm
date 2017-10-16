package Bio::Adventure::Torque;
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

has args => (is => 'rw', default => '-j eo -V -m n');
has basedir => (is => 'rw', default => getcwd());
has cpus => (is => 'rw', default => '4');
has depends_prefix => (is => 'rw', default => 'depend=afterok:');
has depends => (is => 'rw', default => undef);
has dependsarray_prefix => (is => 'rw', default => 'depend=afterokarray:');
has job_name => (is => 'rw', default => 'unnamed');
has job_type => (is => 'rw', default => 'unknown');
has logdir => (is => 'rw', default => getcwd());
has loghost => (is => 'rw', default => 'localhost');
has mem => (is => 'rw', default => '6');
has qsub => (is => 'rw', default => Check_Qsub());
has queue => (is => 'rw', default => 'workstation');
has queues => (is => 'rw', default => qq"throughput,workstation,long,large");
has shell => (is => 'rw', default => '/usr/bin/bash');
has verbose => (is => 'rw', default => 0);
has walltime => (is => 'rw', default => '10:00:00');

sub Check_Qsub {
  my ($class, %args) = @_;
  my $path = which('qsub');
  if (!defined($path)) {
    $path = which('bash');
  }
  return($path);
}

=head1 NAME

    Bio::Adventure::Qsub - Submit jobs to the torque cluster.

=head1 SYNOPSIS

    use Bio::Adventure;
    my $hpgl = new Bio::Adventure;
    $hpgl->Cutadapt();

=head2 Methods

=over 4

=item C<Qsub>

    $hpgl->Qsub(); invokes qsub with (hopefully) appropriate
    parameters for various jobs on our Torque cluster.

=cut
sub Submit {
  my ($class, %args) = @_;
  my $options = $class->Get_Vars(args => \%args);

  my $log = qq"$options->{logdir}/outputs/$options->{job_name}.qsubout";

  my $depends_string = "";
  if ($options->{depends}) {
      if ($options->{depends_type} eq 'array') {
          $depends_string = qq"$options->{dependsarray_prefix}$options->{depends}";
      } else {
          $depends_string = qq"$options->{depends_prefix}:$options->{depends}";
      }
  }
  my $script_file = qq"$options->{basedir}/scripts/$options->{job_prefix}$options->{job_name}.sh";
  my $qsub_cmd_line = qq"$options->{qsub} ${script_file}";
  my $has_pbs = 0;
  if ($options->{qsub} =~ /qsub$/) {
    $qsub_cmd_line = qq"$options->{qsub} -W ${depends_string} ${script_file}";
    $has_pbs = 1;
  }

  my $mycwd = getcwd();
  make_path("$options->{basedir}/outputs/status", {verbose => 0}) unless (-r qq"$options->{basedir}/outputs/status");
  make_path("$options->{logdir}", {verbose => 0}) unless (-r qq"$options->{logdir}");
  make_path("$options->{basedir}/scripts", {verbose => 0}) unless (-r qq"$options->{basedir}/scripts");
  my $script_base = basename($script_file);
  my $script_start = qq?#!/usr/bin/env bash
#PBS -V -S $options->{shell} -q $options->{queue}
#PBS -d $options->{basedir}
#PBS -N $options->{job_name} -l mem=$options->{mem}gb -l walltime=$options->{walltime} -l ncpus=$options->{cpus}
#PBS -o ${log} $options->{args}
echo "####Started ${script_file} at \$(date)" >> outputs/log.txt
cd $options->{basedir} || exit
?;
  my $script_end = qq!## The following lines give status codes and some logging
echo \$? > outputs/status/$options->{job_name}.status
echo "###Finished \${PBS_JOBID} $script_base at \$(date), it took \$(( SECONDS / 60 )) minutes." >> outputs/log.txt
!;
  ## It turns out that if a job was an array (-t) job, then the following does not work because
  ## It doesn't get filled in properly by qstat -f...

  ## The following lines used to be in the shell script postscript
  ## Copying the full job into the log is confusing, removing this at least temporarily.
  ##echo "####This job consisted of the following:" >> outputs/log.txt
  ##cat "\$0" >> outputs/log.txt

  if ($has_pbs) {
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

  if ($class->{verbose}) {
    print qq"The job is:
$args{job_string}
";
  }

  my $total_script_string = "";
  $total_script_string .= "${script_start}\n";
  $total_script_string .= "$options->{comment}\n" if ($options->{comment});
  $total_script_string .= "$options->{prescript}\n" if ($options->{prescript});
  $total_script_string .= "$options->{job_string}\n" if ($options->{job_string});
  if ($options->{postscript}) {
    $total_script_string .= qq!if [ \$? == "0" ]; then
   $options->{postscript}
fi
!;
  }
  $total_script_string .= "${script_end}\n";

  my $script = FileHandle->new(">$script_file");
  print $script $total_script_string;
  $script->close();
  chmod(0755, $script_file);

  my $job_id = undef;
  my $qsub_string = qq"${qsub_cmd_line} |";
  my $handle = IO::Handle->new;
  my $qsub_pid = open($handle, "${qsub_cmd_line} |");
  while(my $line = <$handle>) {
    chomp($line);
    $job_id = $line;
  }

  my $job;
  my $short_jobid = "";
  if ($has_pbs) {
    if (!defined($job_id)) {
      warn("The job id did not get defined.  qsub likely failed.");
      return(undef);
    }
    my @jobid_list = split(/\./, $job_id);
    my $short_jobid = shift(@jobid_list);

    print "Starting a new job: ${short_jobid} $options->{job_name}";
    if ($options->{depends}) {
      my @short_dep = split(/\./, $options->{depends});
      my $shortened_dep = shift(@short_dep);
      print ", depending on ${shortened_dep}.";
    }
    print "\n";
  }

  $job = {
      basedir => $options->{basedir},
      cpus => $options->{cpus},
      depends_string => $depends_string,
      job_id => $short_jobid,
      job_input => $options->{job_input},
      job_name => $options->{job_name},
      job_output => $options->{job_output},
      log => $log,
      mem => $options->{mem},
      queue => $options->{queue},
      pbs_id => $job_id,
      qsub_args => $options->{args},
      script_body => $args{job_string},
      script_file => $script_file,
      script_start => $script_start,
      submitter => $qsub_cmd_line,
      walltime => $options->{wall},
  };
  return($job);
}

sub Submit_Perl {
  my ($class, %args) = @_;

  my $job_type;
  if ($args{job_type}) {
    $job_type = $args{job_type};
  } else {
    $job_type = $class->{job_type};
  }
  my $job_name;
  if ($args{job_name}) {
    $job_name = $args{job_name};
  } else {
    $job_name = $class->{job_name};
  }
  $job_name = qq"${job_type}-${job_name}";
  my $job_prefix = "";
  $job_prefix = $args{job_prefix} if ($args{job_prefix});
  ## For arguments to qsub, start with the defaults in the constructor in $class
  ## then overwrite with any application specific requests from %args
  my $log = qq"$class->{logdir}/${job_name}.qsubout";
  my $job_output = "";
  $job_output = $args{output} if (defined($args{output}));
  my $job_input = "";
  $job_input = $args{input} if (defined($args{input}));

  my $depends_string = $class->{depends_prefix};
  if (defined($args{depends_type})) {
    if ($args{depends_type} eq 'array') {
      $depends_string = $class->{dependsarray_prefix};
    }
  }
  $depends_string .= $class->{job_depends} if (defined($class->{job_depends}));
  my $script_file = qq"$class->{basedir}/scripts/${job_prefix}${job_name}.sh";

  my $qsub = qq"$class->{qsub} ${script_file}";
  my $has_pbs = 0;
  if ($class->{qsub} =~ /qsub$/) {
    $qsub = qq"$class->{qsub} -W ${depends_string} ${script_file}";
    $has_pbs = 1;
  }

  my $mycwd = getcwd();
  make_path("$class->{basedir}/outputs/status", {verbose => 0}) unless (-r qq"$class->{basedir}/outputs/status");
  make_path("$class->{logdir}", {verbose => 0}) unless (-r qq"$class->{logdir}");
  make_path("$class->{basedir}/scripts", {verbose => 0}) unless (-r qq"$class->{basedir}/scripts");
  my $script_base = basename($script_file);

  $script_file = qq"$class->{basedir}/scripts/${job_prefix}${job_name}.pl";
  my $script_start = qq?#!/usr/bin/env perl
use strict;
use FileHandle;
my \$out = FileHandle->new(">>outputs/log.txt");
my \$d = qx'date';
print \$out "###Started $script_file at \${d}";
chdir("$class->{basedir}");
?;
  my $script_end = qq!## The following lines give status codes and some logging
my \$jobid = "";
\$jobid = \$ENV{PBS_JOBID} if (\$ENV{PBS_JOBID});
my \$end_d = qx'date';
print \$out "####Finished \${jobid} ${script_base} at \${end_d}.";
close(\$out);
!;
  print "The job is:
$args{job_string}" if ($class->{verbose});

  my $total_script_string = "";
  $total_script_string .= "$script_start\n";
  $total_script_string .= "$args{comment}\n" if ($args{comment});
  $total_script_string .= "$args{prescript}\n" if ($args{prescript});
  $total_script_string .= "$args{job_string}\n" if ($args{job_string});
  $total_script_string .= "$script_end\n";

  my $script = FileHandle->new(">$script_file");
  print $script $total_script_string;
  $script->close();
  chmod(0755, $script_file);

  my $bash_script = $script_file;
  $bash_script =~ s/\.pl/\.sh/g;
  my %new_args = %args;
  $new_args{job_string} = qq"${script_file}\n";
  my $shell_job = $class->Submit(%new_args);

  my $perl_job = {id => $shell_job->{id},
                  submitter => $qsub,
                  mem => $class->{mem},
                  walltime => $class->{walltime},
                  cpus => $class->{cpus},
                  jobname => $shell_job->{jobname},
                  log => $log,
                  depends_string => $depends_string,
                  queue => $class->{queue},
                  qsub_args => $class->{args},
                  basedir => $class->{basedir},
                  pbs_id => $shell_job->{pbs_id},
                  script_file => $shell_job->{script_file},
                  script_start => $shell_job->{script_start},
                  script_body => $args{job_string},
                  output => $shell_job->{output},
                  input => $shell_job->{input}};

  return($shell_job);
  ##return([$shell_job, $perl_job]);
}

=back

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

    L<Bio::Tools::Run::StandAloneBlast>

=cut

1;
