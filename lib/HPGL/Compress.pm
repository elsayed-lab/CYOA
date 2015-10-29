package HPGL;
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
        $input2 = "";
    }
    my $job_name = "xz";
    $job_name = $args{job_name} if ($args{job_name});
    my $depends = "";
    $depends = $args{depends} if ($args{depends});
    my $comment = "";
    $comment = $args{comment} if ($args{comment});
    my $basename = $me->{basename};
    my $job_string = "";
    my $indir = dirname($input1);
    my ($in1, $in2);
    $in1 = dirname($input1) . '/' . basename($input1, ('.gz','.bz2'));
    $in2 = dirname($input2) . '/' . basename($input2, ('.gz','.bz2')) if ($input2);

    if ($input1 =~ /\.gz/) {
        $job_string = qq!gunzip -f $input1 $input2 && xz -f -9e $in1 $in2\n!;
    } elsif ($input1 =~ /\.bz2/) {
        $job_string = qq!bunzip2 -f $input1 $input2 && xz -f -9e $in1 $in2\n!;
    } elsif ($input1 =~ /\.fast[a|q]$/) {
        $job_string = qq!pxz -f -9e $in1 $in2\n!;
    } else {
        $job_string = qq!pxz -f -9e $in1 $in2\n!;
    }
    if ($args{output}) {
        $job_string .= qq!mv ${in1}.xz $args{output}\n!;
	if ($input2) {
	    $job_string .= qq!mv ${in2}.xz $args{output}\n!;
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
    my $job_string;
    if ($input =~ /\.gz$/) {
        $job_string = qq!gunzip -f $input1 $input2\n!;
    } elsif ($input =~ /\.bz2$/) {
        $job_string = qq!bunzip2 -f $input1 $input2\n!;
    } elsif ($input =~ /\.xz$/) {
        $job_string = qq!xz -f -d $input1 $input2\n!;
    } else {
        $job_string = qq!xz -f -d $input1 $input2\n!;
    }

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

1;
