package Bio::Adventure::Ceph;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

=head1 NAME

    Bio::Adventure::Ceph - upload/download  sequences from the umiacs amazon S3 service

=head1 SYNOPSIS

    use CYOA;
    my $hpgl = new CYOA;
    $hpgl->Dump_Reads();

=over 4

=Methods

=cut

=item C<Ceph_Upload>

    Upload data

=cut

=item C<Dump_Reads>

    $hpgl->Dump_Reads()
    dumps reads from the given hpgl directory into the current working
    directory into concatenated files.  This is in opposition to
    Download_Reads which keeps them as separate files.

=cut
sub Dump_Reads {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args, required => ['hpgl',]);
    my $casava_passed = 0;
    my $line_count = 0;
    my $connection = Bio::Adventure::Ceph::Ceph_Connect($class, %args);
    my $result = Bio::Adventure::Ceph::Download_Files($class, %args, connection => $connection);
    my $seq = $line_count / 4;
    print "$seq sequences were downloaded, of which $casava_passed passed casava's filtering.\n";
    return($seq)
}

sub Ceph_Connect {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $connection = Net::Amazon::S3->new({aws_access_key_id => $options->{access_key},
                                           aws_secret_access_key => $options->{secret_key},
                                           host => $options->{host},
                                           secure => 1,
                                           retry => 0,});
    die("The Amazon::S3 connection failed.") if (!defined($connection));
    return($connection);
}

sub Download_Files {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $connection = $options->{connection};
    my $fwd = "$options->{bucket}_forward.fastq";
    my $rev = "$options->{bucket}_reverse.fastq";
    my ($fwd_result, $rev_result);
    if ($options->{orientation} eq 'both') {
        $fwd_result = Bio::Adventure::Ceph::Perform_Download($class, %args,
                                                            orientation => 'forward',
                                                            dest => $fwd);
        $rev_result = Bio::Adventure::Ceph::Perform_Download($class, %args,
                                                             orientation => 'reverse',
                                                             dest => $rev);
    } elsif ($options->{orientation} eq 'forward') {
        $fwd_result = Bio::Adventure::Ceph::Perform_Download($class, %args,
                                                             orientation => 'forward',
                                                             dest => $fwd);
    } elsif ($options->{orientation} eq 'reverse') {
        $rev_result = Bio::Adventure::Ceph::Perform_Download($class, %args,
                                                             orientation => 'reverse',
                                                             dest => $rev);
    } else {
        die("I can only download forward reads, reverse reads, or both.");
    }
    return(0);
}

sub Perform_Download {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $connection = $options->{connection};
    my $dest = $options->{dest};
    my $orientation = $options->{orientation};
    if (-r $dest) {
        print "$dest already exists, deleting it.\n" if ($options->{verbose});
        unlink($dest);
    }

    my $response = $connection->buckets;
    my @buckets = @{$response->{buckets}};
    my $c = 0;
    my $found_bucket = 0;
  BUCKETS: foreach my $bucket (@buckets) {
        my $bucket_name = $bucket->bucket;
        if ($bucket_name eq $options->{bucket}) {
            $found_bucket++;
            my $res = Bio::Adventure::Ceph::Download_Bucket($class, %args,
                                                            bucket => $bucket,
                                                            dest => $dest,
                                                            orientation => $orientation);

            if ($options->{multi}) {
              INNER: foreach my $buck (@buckets) {
                  LETTERS: foreach my $letter ("a".."z") {
                        my $test = join("", $options->{bucket}, $letter);
                        my $new_name = $buck->bucket();
                        print "Now testing for $new_name and $test\n";
                        if ($new_name eq $test) {
                            Bio::Adventure::Ceph::Download_Bucket($class, %args,
                                                                  bucket => $buck,
                                                                  dest => $dest,
                                                                  orientation => $orientation);
                            next LETTERS;
                        } else {
                            next INNER;
                        }
                    }     ## End rechecking the list of buckets for xxxa -> xxxz
                }         ## End making the names xxxa -> xxxz
            }             ## End if multi-buckets are on
        } else {          ## End if we find hpgl0xxx
            next BUCKETS;
        }
    }
    if ($found_bucket == 0) {
        die("Did not find the bucket: $options->{bucket}");
    }
    return(@buckets);
}

sub Download_Bucket {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => %args);
    my $bucket = $options->{bucket};
    my $dest = $options->{dest};
    my $orientation = $options->{orientation};
    my $out = FileHandle->new(">>${dest}");
    my $connection = $options->{connection};
    my $test = $bucket->list_all or die(qq"$connection->err : $connection->errstr");
    my @keys = @{$bucket->list_all->{keys} || []};

    ## The following for loop creates a list of files to download from ceph.
    ## The logic therein is very specific to Najib's data and should be changed for other
    ## applications.
    my @files_to_download = ();
    foreach my $key (@keys) {
        my $keyname = $key->{key};
        next if ($keyname =~ /\.md5$/);
        if ($orientation eq 'forward') {
            ## The difference between forward/reverse is R1 vs R2
            if ($keyname =~ /.*_[A|T|G|C]+_L\d+_R1_\d+\.fastq/) {
                push(@files_to_download, $keyname);
            }
        } else {
            if ($keyname =~ /.*_[A|T|G|C]+_L\d+_R2_\d+\.fastq/) {
                push(@files_to_download, $keyname);
            }
        }
    }

    ## Now that we have a list of files to download, do it!
    foreach my $file (@files_to_download) {
        my $tmp_dest = basename($file);
        if (-r $tmp_dest) {
            print "$tmp_dest already exists, deleting it.\n" if ($options->{verbose});
            unlink($tmp_dest);
        }
        print "Downloading $file.\n";
        $bucket->get_key_filename($file, undef, $tmp_dest);
        my $md5_file = "$file.md5";
        my $md5_dest = "$tmp_dest.md5";
        if (-r $md5_dest) {
            print "$md5_dest already exists, deleting it.\n" if ($options->{verbose});
            unlink($md5_dest);
        }
        $bucket->get_key_filename($md5_file, undef, $md5_dest);
        print "Checking md5 Sum\n";
        my $m = Bio::Adventure::Ceph::Check_MD5($class, %args,
                                                dest => $tmp_dest,
                                                md5 => $md5_dest);
        my @suffixes = (".gz",".bz2",".xz");
        my $inter_dest = basename($file, @suffixes);
        if ($file =~ /\.fastq$/) { ## Then it wasn't compressed
            print "The file was not compressed.\n";
        } else {
            if (-r $inter_dest) {
                print "$inter_dest already exists, deleting it.\n" if ($options->{verbose});
                unlink($inter_dest);
            }
            print "Extracting $tmp_dest.\n";
            my $ae = Archive::Extract->new(archive => $tmp_dest);
            my $ok = $ae->extract(to => $inter_dest);
        }
        print "Concatenating $inter_dest onto $dest.\n" if ($options->{verbose});
        my $in = FileHandle->new("<$inter_dest");
        my $line_count = 0;
        my $casava_passed = undef;
        while (my $line = <$in>) {
            $line_count++;
            $casava_passed++ if ($line =~ m/^@.* [^:]*:N:[^:]*:/);
            print $out;
        }
        $in->close();
        print "Cleaning up $inter_dest and $tmp_dest.\n" if ($options->{verbose});
        unlink($inter_dest) if (-r $inter_dest);
        unlink($tmp_dest) if (-r $tmp_dest);
    }
    $out->close();
    return($bucket);
}

sub Ceph_List {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $connection = $options->{connection};
    my @buckets = ();
    if (!$options->{bucket}) {
        my $buck = Bio::Adventure::Ceph::List_Buckets($class, %args);
    } else {
        my $response = $connection->buckets;
        @buckets = @{$response->{buckets}};
        my $found_bucket = 0;
        foreach my $bucket (@buckets) {
            my $name = $bucket->bucket();
            if ($name eq $options->{bucket}) {
                $found_bucket++;
            } else {
                next;
            }
            my @keys = @{$bucket->list_all->{keys} || []};
            my @files_to_download = ();
            foreach my $key (@keys) {
                print "$key->{key}\t$key->{size}\t$key->{last_modified}\n";
            }
        }
        if ($found_bucket == 0) {
            print "Did not find $options->{bucket}\n";
            my $list = Bio::Adventure::Ceph::List_Buckets($class, %args);
        }
    }
    return(@buckets);
}

sub List_Buckets {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $connection = $options->{connection};
    print "Available buckets are:\n";
    my $response = $connection->buckets;
    my @buckets = @{$response->{buckets}};
    foreach my $bucket (@buckets) {
        my $name = $bucket->bucket();
        print "$name\n";
    }
    return(@buckets);
}

sub Check_MD5 {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(args => \%args);
    my $md5_file = $options->{md5};
    my $downloaded_file = $options->{dest};
    my $md5_in = FileHandle->new("<$md5_file");
    my $original_md5 = "";
    while (my $line = <$md5_in>) {
        chomp($line);
        my ($md5, $filename) = split(/\s+/, $line);
        $original_md5 = $md5;
    }
    $md5_in->close();
    unlink($md5_file);
    my $checksum = Digest::MD5->new();
    my $md5_check = FileHandle->new("<$downloaded_file");
    my $chk = \$md5_check;
    ## $chk = \*MD5_CHECK;  ## I don't understand how this works...
    $checksum->addfile($chk);
    my $new_md5_checksum = $checksum->hexdigest;
    my ($new_md5, $new_filename) = split(/\s+/, $new_md5_checksum);
    if ($new_md5 eq $original_md5) {
        return(1);
    } else {
        print "The checksums no longer agree.  This file has changed since it was uploaded,
or the download failed and should be reperformed.\n";
        return(undef);
    }
    return($chk);
}

1;
