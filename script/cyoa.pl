#!/usr/bin/env perl
use Cwd qw"abs_path cwd";
use lib abs_path('lib');
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

=head2 C<menus>

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
## Then a second-level method to perform therein,
## which will provide a match to the appropriate function
my $overrides = $cyoa->{variable_getopt_overrides};
if (defined($overrides->{help})) {
    my $helped = $cyoa->Help();
}
if (!defined($overrides->{task}) && !defined($overrides->{method})) {
    ## If task is not defined, then we need to run Main()
    $cyoa = Main($cyoa);
} elsif (!defined($overrides->{method})) {
    ## If task is defined, but method is not, use the menu system to get that information
    $cyoa = Iterate($cyoa, task => $overrides->{task});
} else {
    my $task = $cyoa->{task};
    $task = $overrides->{task} if (defined($overrides->{task}));
    ## If both task and method are defined, run whatever is requested.
    $cyoa = Run_Method($cyoa, task => $task,
                       method => $overrides->{method});
}

=head2 C<Main>

The Main() function's job is to actually perform the various processes.

=cut
sub Main {
    my ($class, %args) = @_;
    my $term = $class->{term};
    if (!defined($term)) {
        $term = Bio::Adventure::Get_Term();
    }
    my $menus = Bio::Adventure::Get_Menus();
    my $finished = 0;
    while ($finished != 1) {
        my @choices = sort keys(%{$menus});
        my $top_level = $term->get_reply(prompt => 'Choose your adventure!',
                                         choices => \@choices, default => 'TNSeq',);
        my $tasks = Iterate($cyoa, term => $term, task => $top_level, interactive => 1,);
        my $bool = $term->ask_yn(prompt => 'Are you done?', default => 'n',);
        $finished = 1 if ($bool eq '1');
        ## my $string = q[some_command -option --no-foo --quux='this thing'];
        ## my ($options,$munged_input) = $term->parse_options($string);
        ## don't have Term::UI issue warnings -- default is '1'
        ## Retrieve the entire session as a printable string:
        ## $hist = Term::UI::History->history_as_string;
        ## $hist = $term->history_as_string;
    } ## End while not finished.
    return($class)
}

=head2 C<Iterate>

The Iterate runs multiple methods and comes back to the menu interface after each.

=cut
sub Iterate {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars();
    my $term = $args{term};
    my $menus = Bio::Adventure::Get_Menus();
    my $task = $args{task};
    my $finished = 0;
    my $process;
    while ($finished != 1) {
        my $choices = $menus->{$task}->{choices};
        my @choice_list = sort keys(%{$choices});
        my $task_name = $term->get_reply(prompt => $menus->{$task}->{message},
                                         choices => \@choice_list, default => 1);
        ## This now returns a function reference.
        my $task = $choices->{$task_name};
        my $cyoa_result;
        my $start_dir = cwd();
        my @dirs = ($start_dir);
        if (defined($options->{directories})) {
            my @dirs = split(/:|\s+|\,/, $options->{directories});
        }
        foreach my $d (@dirs) {
            chdir($d);
            $process = $task->($class, task => $task, interactive => $args{interactive});
            chdir($start_dir);
        }
        my $bool = $term->ask_yn(prompt => 'Go back to the toplevel?', default => 'n',);
        if ($bool eq '1') {
            $finished = 1;
            my $reset = $term->get_reply(
                prompt => 'How many of your parameters do you want to reset?',
                choices => ['Everything', 'Just input', 'Nothing'], default => 'Just input');
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
    return($class);
}

=head2 C<Run_Method>

The Run_Method() runs an individual method, kind of like it says on the tin.

=cut
sub Run_Method {
    my ($class, %args) = @_;
    my $task = $args{task};
    my $method = $args{method};

    my @choices = @{$class->{methods_to_run}};
    my $job_count = 0;
    my $process;
    for my $job (@choices) {
        $job_count = $job_count++;
        my $cyoa_result;
        my $start_dir = cwd();
        my @dirs = ($start_dir);
        if (defined($class->{directories})) {
            @dirs = split(/:|\s+|\,/, $class->{directories});
        }
        my $class_copy = $class;
        foreach my $d (@dirs) {
            sleep(1);
            chdir($d);
            $process = $job->($class_copy, task => $task, interactive => $args{interactive}, %args);
            chdir($start_dir);
        }
    }
    return($class);
}

=head1 AUTHOR - atb

Email abelew@gmail.com

=head1 SEE ALSO

    L<Adventure> L<Adventure::Convert> L<Adventure::Trim>

=cut

## EOF
