=head2
    Cleanup_Seq()
=cut
sub Cleanup_Seq {
    my $me = shift;
    my %args = @_;
    my $depends = "";
    if ($args{depends}) {
        $depends = $args{depends};
    }
    my $job_string = qq!rm -rf status outputs scripts sequences && find . -name '*.fastq*' -exec rm {} ';'
!;
    my $clean = $me->Qsub(job_name => "th_clean",
                          qsub_queue => 'throughput',
                          depends => $depends,
                          job_string => $job_string,
			  prescript => $args{prescript},
			  postscript => $args{postscript},
        );
    return($clean);
}

=head2
    Cleanup_Tophat()
=cut
sub Cleanup_Tophat {
    my $me = shift;
    my %args = @_;
    my $depends = "";
    if ($args{depends}) {
        $depends = $args{depends};
    }
    my $job_string = qq!rm -rf tophat_out
!;
    my $clean = $me->Qsub(job_name => "th_clean",
                          qsub_queue => 'throughput',
                          depends => $depends,
                          job_string => $job_string,
			  prescript => $args{prescript},
			  postscript => $args{postscript},
        );
    return($clean);
}
