=head2
    Recompress()
=cut
sub Recompress {
    my $me = shift;
    my %args = @_;
    my $input;
    if ($args{input}) {
        $input = $args{input};
    } else {
        $input = $me->{input};
    }
    my ($input1, $input2);
    if ($input =~ /\:|\,/) {
	($input1, $input2) = split(/\:|\,/, $input);
    } else {
	$input1 = $input;
    }
    my $job_name = "xz";
    $job_name = $args{job_name} if ($args{job_name});
    my $depends = "";
    $depends = $args{depends} if ($args{depends});
    my $comment = "";
    $comment = $args{comment} if ($args{comment});
    my $basename = $me->{basename};
    my $job_string = qq!pxz -f -9e $input1 $input2
!;
    if ($args{output}) {
        $job_string .= qq!mv ${input1}.xz $args{output}
!;
	if ($input2) {
	    $job_string .= qq!mv ${input2}.xz $args{output2}
!;
	}
    }

    my $trim_jobid = qq"${basename}_xz";
    my $compression = $me->Qsub(job_name => $job_name,
                                depends => $depends,
                                qsub_queue => "throughput",
                                qsub_cpus => 1,
                                qsub_mem => 4,
                                qsub_wall => "18:00:00",
                                job_string => $job_string,
                                input => $input,
                                comment => $comment,
				prescript => $args{prescript},
				postscript => $args{postscript},
        );
    return($compression);
}

=head2
    Uncompress()
=cut
sub Uncompress {
    my $me = shift;
    my %args = @_;
    my $input;
    if ($args{input}) {
        $input = $args{input};
    } else {
        $input = $me->{input};
    }
    my ($input1, $input2);
    if ($input =~ /\:|\,/) {
	($input1, $input2) = split(/\:|\,/, $input);
    } else {
	$input1 = $input;
    }
    my $job_name = "unxz";
    $job_name = $args{job_name} if ($args{job_name});
    my $depends = "";
    $depends = $args{depends} if ($args{depends});
    my $comment = "";
    $comment = $args{comment} if ($args{comment});
    my $basename = $me->{basename};
    my $job_string = qq!pxz -f -d $input1 $input2
!;

    my $trim_jobid = qq"${basename}_unxz";
    my $compression = $me->Qsub(job_name => $job_name,
                                depends => $depends,
                                qsub_queue => "throughput",
                                qsub_cpus => 1,
                                qsub_mem => 4,
                                qsub_wall => "18:00:00",
                                job_string => $job_string,
                                input => $input,
                                comment => $comment,
				prescript => $args{prescript},
				postscript => $args{postscript},
        );
    return($compression);
}
