#!/usr/bin/env perl
# -*-Perl-*-
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";

use File::Basename;
use Getopt::Long qw"GetOptionsFromArray";
use Switch;

use Bio::Adventure;

=head1 NAME

    cyoa - Choose Your Own Adventure!  (In Sequence Processing)

=head1 SYNOPSIS

    ## Starts the menu-based function chooser
    cyoa
    ## Perform an rnaseq specific job to use bwa with the lmajor reference and 2 paired sequence files.
    cyoa --task rnaseq --method bwa --species lmajor --input fwd_reads.fastq.gz:rev_reads.fastq.gz

=head2 Methods & Globals

=over 4

=item C<menus>

    The global variable $menus keeps a set of silly choose your own
    adventure quotes along with the possible task->method listings.
    The choices slot in each key of the hash is a list of function
    names found in Adventure::something.  The assumption is that all the
    functions use Check_Options() to ensure that they receive all the
    relevant parameters.

=cut

## Going to pull an initial instance of cyoa in order to handle command line arguments
our $cyoa = new Bio::Adventure;

## The general idea is to have a toplevel 'task' to perform
## Something like TNSeq, RNASeq, etc
##  Then a second-level method to perform therein, which will provide a match to the appropriate function
if (!defined($cyoa->{options}->{task})) {  ## If task is not defined, then we need to run Main()
    Main($cyoa);
} elsif (!defined($cyoa->{options}->{method})) { ## If task is defined, but method is not, use the menu system to get that information
    Iterate($cyoa, task => $cyoa->{options}->{task});
} else {
    Run_Method($cyoa, task => $cyoa->{options}->{task}, method => $cyoa->{options}->{method});
}


sub Main {
  my ($class, %args) = @_;
  my $term = $class->{options}->{term};
  if (!defined($term)) {
      $term = Bio::Adventure::Get_Term();
  }
  my $menus = $class->{options}->{menus};

  my $finished = 0;
  while ($finished != 1) {
    my @choices = sort keys(%{$menus});
    my $top_level = $term->get_reply(prompt => 'Choose your adventure!',
                                     choices => \@choices,
                                     default => 'TNSeq',);
    my $tasks = Iterate($class, task => $top_level, interactive => 1,);
    my $bool = $term->ask_yn(prompt => 'Are you done?',
                             default => 'n',);
    $finished = 1 if ($bool eq '1');
    #my $string = q[some_command -option --no-foo --quux='this thing'];
    #my ($options,$munged_input) = $term->parse_options($string);
    ### don't have Term::UI issue warnings -- default is '1'
    ### Retrieve the entire session as a printable string:
    #$hist = Term::UI::History->history_as_string;
    #$hist = $term->history_as_string;
  } ## End while not finished.
}

sub Iterate {
  my ($class, %args) = @_;
  my $options = $class->{options};
  my $menus = $options->{menus};
  my $term = $options->{term};
  if (!$term) {
      $term = Bio::Adventure::Get_Term();
  }
  my $type = $args{task};
  my $finished = 0;
  my $tasks_done = 0;
  while ($finished != 1) {
    my $choices = $menus->{$type}->{choices};
    my @choice_list = sort keys(%{$choices});
    my $task_name = $term->get_reply(prompt => $menus->{$type}->{message},
                                     choices => \@choice_list,);
    my $task = $choices->{$task_name};
    if ($task eq 'Cutadapt' or $task =~ /Trim/) {
      if ($cyoa->{input}) {
        my $new_input = basename($options->{input}, @{$options->{suffixes}}) . "-trimmed.fastq";
        $cyoa->{options}->{input} = $new_input;
        print "I perceive that thou dost perform a trimming operation, resetting the input to $new_input.\n";
      }
    }
    my $cyoa_result;
    {
      $tasks_done = $tasks_done++;
      no strict 'refs';
      $cyoa_result = &{$task}($cyoa, type => $type, interactive => $args{interactive});
    }
    my $bool = $term->ask_yn(prompt => 'Go back to the toplevel?',
                             default => 'n',);
    if ($bool eq '1') {
      $finished = 1;
      my $reset = $term->get_reply(prompt => 'How many of your parameters do you want to reset?',
                                   choices => ['Everything', 'Just input', 'Nothing'],
                                   default => 'Just input');
      if ($reset eq 'Just input') {
        $cyoa->{input} = undef;
      } elsif ($reset eq 'Everything') {
        $cyoa->{input} = undef;
        $cyoa->{species} = undef;
        ## Some other stuff which I can't remember now
      } else {
        ## Leave $cyoa alone
      }
    }
  } ## End the task while loop
  return($tasks_done);
}

sub Run_Method {
  my ($class, %args) = @_;
  my $options = $class->{options};
  my $type = $args{task};
  my $method = $args{method};

  my @choices = @{$options->{methods_to_run}};
  my $job_count = 0;
  for my $job (@choices) {
    $job_count = $job_count++;
    print "About to run: $job\n";
    my $process;
    {
      no strict 'refs';
      $process = &${job}($cyoa, type => $type, interactive => $args{interactive});
    }
    ##my $process = $cyoa->$job(type => $type, interactive => $args{interactive},);

    ##
    ## I have decided that this is a dumb thing to do, I may change my mind and so am not deleting it yet.
    ##
    ## my $print_job = new FileHandle(">last_job.txt");
    ## my $process_class = ref($process);
    ## if ($process_class eq 'ARRAY') {
    ##   ## This happens in the case of bowtie mappings because we make a list of jobs.
    ##   for my $p (@{$process}) {
    ##     print $print_job "$p->{job_id}\n" if (defined($p->{job_id}));
    ##   }
    ## } else {
    ##   print $print_job "$process->{job_id}\n" if (defined($process->{job_id}));
    ## }
    ## $print_job->close();

  }
  return($job_count);
}

=back

=head1 AUTHOR - atb

Email abelew@gmail.com

=head1 SEE ALSO

    L<Adventure> L<Adventure::Convert> L<Adventure::RNASeq_Trim>

=cut

## EOF
