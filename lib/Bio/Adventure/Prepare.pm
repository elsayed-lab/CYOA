package Bio::Adventure::Prepare;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';
use Bio::Adventure::Torque;

use Text::CSV;

=head1 NAME

Bio::Adventure::Prepare - Get raw data ready for processing.

=head1 SYNOPSIS

use Bio::Adventure;
my $hpgl = new Bio::Adventure;
$hpgl->Prepare(csv => 'all_samples.csv');

=head1 Methods

=head2 C<Copy_Raw>

Copy from the raw_data store.  This is probably no longer useful.

=cut
sub Copy_Raw {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args,
                                   required => ['raw_dir', 'sampleid'],);
    if ($args{interactive}) {
        print "Run with: cyoa --task prepare --method copy --csv all_samples.csv --raw_dir /some/directory
      or: cyoa --task prepare --method copy --raw_dir /some/directory --sampleid HPGL0xxx\n";
    }
    my $data = {};
    if ($options->{sampleid}) {
        my $id = $options->{sampleid};
        $data->{$id} = "raw";
    } else {
        $data = Bio::Adventure::Prepare::Read_Samples($class, %args);
    }

    my $files_read = 0;
    my $forward_files = 0;
    my $reverse_files = 0;
    foreach my $id (keys %{$data}) {
        my $lc_id = lc($id);
        my $uc_id = uc($id);
        my $input_dir = qq"$options->{raw_dir}/${uc_id}/unprocessed";
        opendir(my $dh, $input_dir) or die;
      READ: while (readdir $dh) {
            next READ unless ($_ =~ /\.gz$/);
            $files_read++;
            my $cat = qx"mkdir -p ${lc_id}";
            my $cmd = "";
            if ($_ =~ /_R1_/) {
                $cmd = qq"cat ${input_dir}/$_ >> ${lc_id}/${lc_id}_forward.fastq.gz";
                $forward_files++;
            } elsif ($_ =~ /_R2_/) {
                $cmd = qq"cat ${input_dir}/$_ >> ${lc_id}/${lc_id}_reverse.fastq.gz";
                $reverse_files++;
            } else {
                print "Didn't copy $_\n";
            }
            print "$cmd\n";
            my $ret = qx($cmd);
        }                       ## End reading the directory
        closedir $dh;
        print "Sample ${id} finished.\n";
    }                     ## End each sample in the csv or commandline arguments
    print "Read ${files_read} files including ${forward_files} forward and ${reverse_files} reverse sequences.\n";
    return($files_read);
}

=head2 C<Read_Samples>

Read the sample csv file.

=cut
sub Read_Samples {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args, required => ['csv_file']);
    if (!-r $options->{csv_file}) {
        print "Unable to find $options->{csv_file}.\n";
        die();
    }
    my $fh = new FileHandle("<$options->{csv_file}");
    my $csv = Text::CSV->new({binary => 1});
    my $row_count = 0;
    my @headers = ();
    my $data = {};
    while (my $row = $csv->getline($fh)) {
        $row_count++;
        if ($row_count == 1) {
            @headers = @{$row};
        } else {
            my $id = $row->[0];
            foreach my $c (1 .. $#headers) {
                my $key = $headers[$c];
                $data->{$id}->{$key} = $row->[$c];
            }
        }
    }
    $csv->eof or $csv->error_diag();
    $fh->close();
    return($data);
}

1;
