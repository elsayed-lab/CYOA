package HPGL::Job;
our %jobs;

=head 1 NAME
    HPGL::Job
    A super-simple state-keeper for torque jobs.
=cut
sub new {
    my ($class, %args) = @_;
    my $me = bless {}, $class;
    foreach my $key (keys %args) {
        $me->{$key} = $args{$key} if ($args{$key});
    }
    my $id = $me->{id};
    $jobs{$id} = $me;
    return($me);
}
