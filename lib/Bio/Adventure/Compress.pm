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
sub Compress {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        comment => '## Compressing files.',
        jmem => 8,
        jname => 'xz',
        jprefix => '',
        jqueue => 'long',
        jwalltime => '24:00:00',
        modules => undef,);
    my $input_paths = $class->Get_Paths($options->{input});
    my $check = which('xz');
    die('Could not find xz in your PATH.') unless($check);
    my $jstring = "";
    my $output_string = '';
    for my $in (@{$input_paths}) {
        my $in_dir = $in->{directory};
        my $in_base = $in->{filebase_compress};
        my $in_full = $in->{fullpath};
        my $output_file = qq"${in_dir}/${in_base}.xz";
        $output_string .= qq"${output_file}:";
        $jstring .= qq!
## Compressing ${in_full}
echo "Compressing ${in_full}"
if [ -f "${in_full}" ]; then
  xz -9e -f ${in_full}
  if [ "\$?" -ne "0" ]; then
    echo "The compression of ${in_full} failed."
  fi
else
  echo "The input: ${in_full} does not exist."
fi

!;
    }
    $output_string =~ s/:$//g;

    my $compression = $class->Submit(
        comment => $options->{comment},
        jdepends => $options->{jdepends},
        input => $options->{input},
        jcpu => 1,
        jmem => $options->{jmem},
        jname => $options->{jname},
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jwalltime => $options->{jwalltime},
        output => $output_string,);
    return($compression);
}

sub Recompress {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        comment => '## Recompressing files.',
        jname => 'xz',
        jmem => 8,
        jprefix => '99',
        jqueue => 'long',
        jwalltime => '12:00:00',);
    my $input_paths = $class->Get_Paths($options->{input});

    my $jstring = "";
    my $output_string = '';
    for my $in (@{$input_paths}) {
        my $in_dir = $in->{directory};
        my $in_base = $in->{filebase_compress};
        my $in_full = $in->{fullpath};
        my $output_file = qq"${in_dir}/${in_base}.xz";
        $output_string .= qq"${output_file}:";
        if ($in_full =~ /\.gz|\.bz2|\.zip/) {
            $jstring .= qq!
less ${in_full} 2>${in_dir}/${in_base}_decompress.stderr 1>${in_dir}/${in_base}
xz -9e -f ${in_dir}/${in_base}
rm ${in_full}
rm ${in_dir}/${in_base}_decompress.stderr!;
        } else {
            $jstring .= qq!
xz -9e -f ${in_full}
!;
        }
    }
    $output_string =~ s/:$//g;

    my $compression = $class->Submit(
        comment => $options->{comment},
        jdepends => $options->{jdepends},
        input => $options->{input},
        jcpu => 1,
        jmem => $options->{jmem},
        jname => $options->{jname},
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jwalltime => $options->{jwalltime},
        output => $output_string,);
    return($compression);
}

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

L<xz> L<gzip> L<bash>

=cut

1;
