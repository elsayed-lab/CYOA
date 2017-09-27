package Bio::Adventure::Compress;
## LICENSE: gplv2
## ABSTRACT:  Kitty!
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use File::Basename;

=head2

    Recompress()

=cut
sub Recompress {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        xz_input => $class->{options}->{input},
        basename => "recomp",
        job_name => "xz",
        job_depends => "",
        comment => "",
    );
    my $input = $options->{input};
    my ($input1, $input2) = "";
    if ($input =~ /\:|\,/) {
        ($input1, $input2) = split(/\:|\,/, $input);
    } else {
        $input1 = $input;
        $input2 = "";
    }
    my $job_name = $options->{job_name};
    my $job_depends = $options->{job_depends};
    my $comment = $options->{comment};
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

    my $trim_jobid = qq"$options->{basename}_xz";
    my $compression = $class->Submit(
        comment => $comment,
        cpus => 1,
        input => $input,
        job_depends => $job_depends,
        job_name => $job_name,
        job_prefix => $args{job_prefix},
        job_string => $job_string,
        mem => 4,
        prescript => $args{prescript},
        postscript => $args{postscript},
        queue => "workstation",
        wall => "18:00:00",
    );
    return($compression);
}

=head2

    Uncompress()

=cut
sub Uncompress {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        xz_input => $class->{options}->{input},
        basename => "uncomp",
        job_name => "unzip",
        job_depends => "",
        comment => "",
    );
    my $input = $options->{input};
    my ($input1, $input2) = "";
    if ($input =~ /\:|\;|\,|\s+/) {
        ($input1, $input2) = split(/\:|\;|\,|\s+/, $input);
        ## Switch various separators with a single space.
        $input =~ s/\:|\;|\,|\s+/ /g;
    } else {
        $input1 = $input;
        $input2 = "";
    }
    my $job_string = "";
    my $indir = dirname($input1);
    my ($in1, $in2) = "";
    $in1 = dirname($input1) . '/' . basename($input1, ('.gz','.bz2'));
    $in2 = dirname($input2) . '/' . basename($input2, ('.gz','.bz2')) if ($input2);
    $in2 = "" if (!defined($in2));
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
    $job_string .= qq" 2>${input_dir}/$options->{job_name}.log 1>&2";

    my $trim_jobid = qq"$options->{basename}_unxz";
    my $compression = $class->Submit(
        comment => $options->{comment},
        cpus => 1,
        input => $input,
        job_depends => $options->{job_depends},
        job_name => $options->{job_name},
        job_prefix => $options->{job_prefix},
        job_string => $job_string,
        mem => 4,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => "workstation",
        wall => "18:00:00",
    );
    return($compression);
}

1;
