package CYOA;

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
    my ($in1, $in2) = "";
    $in1 = dirname($input1) . '/' . basename($input1, ('.gz','.bz2'));
    $in2 = dirname($input2) . '/' . basename($input2, ('.gz','.bz2')) if ($input2);
    $in2 = "" if (!defined($in2));
    if ($input1 =~ /\.gz/) {
        $job_string = qq!gunzip -f ${input1} ${input2} && nice xz -f -9e ${in1} ${in2} !;
    } elsif ($input1 =~ /\.bz2/) {
        $job_string = qq!bunzip2 -f ${input1} ${input2} && nice xz -f -9e ${in1} ${in2} !;
    } elsif ($input1 =~ /\.fast[a|q]$/) {
        $job_string = qq!nice -n 20 xz -f -9e ${in1} ${in2} !;
    } else {
        $job_string = qq!nice -n 20 xz -f -9e ${in1} ${in2} !;
    }
    my $input_dir = dirname(${in1});
    $job_string .= qq" 2>${input_dir}/${job_name}.log 1>&2";
    if ($args{output}) {
        $job_string .= qq! && mv ${in1}.xz $args{output}\n!;
	if ($input2) {
	    $job_string .= qq! && mv ${in2}.xz $args{output2}\n!;
	}
    }

    my $trim_jobid = qq"${basename}_xz";
    my $compression = $me->Qsub(
        comment => $comment,
        depends => $depends,
        input => $input,
        job_name => $job_name,
        job_prefix => $args{job_prefix},
        job_string => $job_string,
        prescript => $args{prescript},
        postscript => $args{postscript},
        qsub_cpus => 1,
        qsub_mem => 4,
        qsub_queue => "workstation",
        qsub_wall => "18:00:00",
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
    my $input1 = "";
    my $input2 = "";
    if ($input =~ /\:|\;|\,|\s+/) {
        ($input1, $input2) = split(/\:|\;|\,|\s+/, $input);
        ## Switch various separators with a single space.
        $input =~ s/\:|\;|\,|\s+/ /g;
    } else {
        $input1 = $input;
    }
    my $job_name = "unzip";
    $job_name = $args{job_name} if ($args{job_name});
    my $depends = "";
    $depends = $args{depends} if ($args{depends});
    my $comment = "";
    $comment = $args{comment} if ($args{comment});
    my $basename = $me->{basename};
    my $job_string;
    if ($input =~ /\.gz$/) {
        $job_string = qq!gunzip -f ${input}\n!;
    } elsif ($input =~ /\.bz2$/) {
        $job_string = qq!bunzip2 -f ${input}\n!;
    } elsif ($input =~ /\.xz$/) {
        $job_string = qq!xz -f -d ${input}\n!;
    } else {
        $job_string = qq!xz -f -d ${input}\n!;
    }
    my $input_dir = dirname(${input1});
    $job_string .= qq" 2>${input_dir}/${job_name}.log 1>&2";

    my $trim_jobid = qq"${basename}_unxz";
    my $compression = $me->Qsub(
        comment => $comment,
        depends => $depends,
        input => $input,
        job_name => $job_name,
        job_prefix => $args{job_prefix},
        job_string => $job_string,
        prescript => $args{prescript},
        postscript => $args{postscript},
        qsub_cpus => 1,
        qsub_mem => 4,
        qsub_queue => "workstation",
        qsub_wall => "18:00:00",
        );
    return($compression);
}

1;
