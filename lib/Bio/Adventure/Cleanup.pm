package Bio::Adventure::Cleanup;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

=head1 NAME

Bio::Adventure::Cleanup - Delete the various directories/files created by the
tools invoked by Bio::Adventure.

=head1 SYNOPSIS

use Bio::Adventure;
my $hpgl = new Bio::Adventure;
$hpgl->Cleanup();

=head1 METHODS

=head2 C<Cleanup>

Delete the various directories/files from the RNASeq tools invoked by Bio::Adventure.pm.

=cut
sub Cleanup {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars();
    my $depends = $options->{depends};
    my $input_files = $class->{input};
    $input_files =~ s/:/ /g;
    my $unzipped = $input_files;
    $unzipped =~ s/\.gz//g;
    my $trimmed_input = $unzipped;
    $trimmed_input =~ s/\.fastq/-trimmed\.fastq/g;
    my $trimmed_paired = $unzipped;
    $trimmed_paired =~ s/\.fastq/-trimmed_paired\.fastq/g;
    my $trimmed_unpaired = $unzipped;
    $trimmed_unpaired =~ s/\.fastq/-trimmed_unpaired\.fastq/g;
    my $jstring = qq!rm -rf outputs scripts sequences ${input_files} ${unzipped} \\
  ${trimmed_input} ${trimmed_paired} ${trimmed_unpaired}
!;
    print "Execute: $jstring\n";
    return(1);
}

sub Cleanup_Phage_Assembly {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        jprefix => '99',
        jname => 'cleanup_phage');

    my $job_name = $class->Get_Job_Name();
    my $input_paths = $class->Get_Paths($options->{input});

    my $jstring = qq!## Rando fastq files
stuff=\$(find . -type f -name '*.fastq')
for file in \${stuff}; do
  if [[ \! -f "\${file}" ]]; then
    echo "Removing \${file}"
    rm -f "\${file}"
  else
   echo "\${file} does not appear to exist."
  fi
done

## Core dumps
stuff=\$(find . -type f -name core)
for file in \${stuff}; do
  if [[ \! -f "\${file}" ]]; then
    echo "Removing \${file}"
    rm -f "\${file}"
  else
    echo "\${file} does not appear to exist."
  fi
done

## Hisat indexes
stuff=\$(find . -type f -name '*.ht2')
for file in \${stuff}; do
  if [[ \! -f "\${file}" ]]; then
    rm -f "\${file}"
  else
    echo "\${file} does not appear to exist."
  fi
done

## tmp files from the various prediction tools
stuff=\$(find . -type f -name '*.tmp.*')
for file in \${stuff}; do
  if [[ \! -f "\${file}" ]]; then
    rm -f "\${file}"
  else
    echo "\${file} does not appear to exist."
  fi
done

## Trinotate junk
stuff=\$(find . -type d -name '*_dir*')
for file in \${stuff}; do
  if [[ \! -f "\${file}" ]]; then
    echo "Removing \${file}"
    rm -rf "\${file}"
  else
    echo "\${file} does not appear to exist."
  fi
done

## Recompress random fastq files, but not symlinks
stuff=\$(find . -type f -size +0 -name '*.fastq')
for file in \${stuff}; do
  if [[ \! -f "\${file}" ]]; then
    echo "Recompressig \${file}"
    xz -9e -f "\${file}"
  else
    echo "\${file} does not appear to exist."
  fi
done
!;
    my $comment = '## Cleanup some of the mess.';
    my $clean = $class->Submit(
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => $options->{jname},
        jprefix => $options->{jprefix},
        jstring => $jstring,);
    return($clean);
}

=head1 AUTHOR - atb

Email <abelew@gmail.com>

=head1 SEE ALSO

    L<Bio::Adventure>

=cut

1;
