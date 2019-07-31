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

=head1 NAME

Bio::Adventure::Compress - Handle the (de/re)compression of data files.

=head1 SYNOPSIS

Some tools can handle compressed input, some cannot.  Bio::Adventure makes heavy use
of bash subshells <(less filename) to automagically handle any format, but at times
one must still decompress/recompress some data.

=head1 METHODS

=head2 C<Recompress>

Invoke xz to recompress a given input file.

=cut
sub Recompress {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['xz_input'],
        basename => 'recomp',
        comment => '',
        depends => '',
        jname => 'xz',
        queue => 'long'
    );
    my $input = $options->{xz_input};
    my ($input1, $input2) = "";
    if ($input =~ /\:|\,/) {
        ($input1, $input2) = split(/\:|\,/, $input);
    } else {
        $input1 = $input;
        $input2 = "";
    }
    my $jname = $options->{jname};
    my $depends = $options->{depends};
    my $comment = $options->{comment};
    my $jstring = "";
    my $indir = dirname($input1);
    my ($in1, $in2) = "";
    $in1 = dirname($input1) . '/' . basename($input1, ('.gz','.bz2'));
    $in2 = dirname($input2) . '/' . basename($input2, ('.gz','.bz2')) if ($input2);
    $in2 = "" if (!defined($in2));
    if ($input1 =~ /\.gz/) {
        $jstring = qq!gunzip -f ${input1} ${input2} && nice xz -f -9e ${in1} ${in2} !;
    } elsif ($input1 =~ /\.bz2/) {
        $jstring = qq!bunzip2 -f ${input1} ${input2} && nice xz -f -9e ${in1} ${in2} !;
    } elsif ($input1 =~ /\.fast[a|q]$/) {
        $jstring = qq!nice -n 20 xz -f -9e ${in1} ${in2} !;
    } else {
        $jstring = qq!nice -n 20 xz -f -9e ${in1} ${in2} !;
    }
    my $input_dir = dirname(${in1});
    $jstring .= qq" 2>${input_dir}/${jname}.err \\
   1>${input_dir}/${jname}.out";
    if ($args{output}) {
        $jstring .= qq! && mv ${in1}.xz $args{output}\n!;
        if ($input2) {
            $jstring .= qq! && mv ${in2}.xz $args{output2}\n!;
        }
    }

    my $trim_jobid = qq"$options->{basename}_xz";
    my $compression = $class->Submit(
        comment => $comment,
        cpus => 1,
        input => $input,
        depends => $depends,
        jname => $jname,
        jprefix => $args{jprefix},
        jstring => $jstring,
        mem => 4,
        prescript => $args{prescript},
        postscript => $args{postscript},
        queue => 'long',
        walltime => '18:00:00',
    );
    return($compression);
}

=head2 C<Uncompress>

This really should not be necessary for anything anymore given my use of bash's <().

=cut
sub Uncompress {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['xz_input'],
        basename => 'uncomp',
        comment => '',
        depends => '',
        jname => 'unzip',
        queue => 'long',
    );
    my $input = $options->{xz_input};
    my ($input1, $input2) = "";
    if ($input =~ /\:|\;|\,|\s+/) {
        ($input1, $input2) = split(/\:|\;|\,|\s+/, $input);
        ## Switch various separators with a single space.
        $input =~ s/\:|\;|\,|\s+/ /g;
    } else {
        $input1 = $input;
        $input2 = "";
    }
    my $jstring = "";
    my $indir = dirname($input1);
    my ($in1, $in2) = "";
    $in1 = dirname($input1) . '/' . basename($input1, ('.gz','.bz2'));
    $in2 = dirname($input2) . '/' . basename($input2, ('.gz','.bz2')) if ($input2);
    $in2 = "" if (!defined($in2));
    if ($input =~ /\.gz$/) {
        $jstring = qq!nice gunzip -f ${input}\n!;
    } elsif ($input =~ /\.bz2$/) {
        $jstring = qq!nice bunzip2 -f ${input}\n!;
    } elsif ($input =~ /\.xz$/) {
        $jstring = qq!nice xz -f -d ${input}\n!;
    } else {
        $jstring = qq!nice xz -f -d ${input}\n!;
    }
    my $input_dir = dirname(${input1});
    $jstring .= qq" 2>${input_dir}/$options->{jname}.err \\
  1>${input_dir}/$options->{jname}.out";

    my $trim_jobid = qq"$options->{basename}_unxz";
    my $compression = $class->Submit(
        comment => $options->{comment},
        cpus => 1,
        input => $input,
        depends => $options->{depends},
        jname => $options->{jname},
        jprefix => $options->{jprefix},
        jstring => $jstring,
        mem => 4,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        queue => 'long',
        walltime => '18:00:00',
    );
    return($compression);
}

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

L<xz> L<gzip> L<bash>

=cut

1;
