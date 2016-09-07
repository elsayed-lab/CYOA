package CYOA;
use common::sense;
use autodie qw":all";

=head1 NAME
    CYOA::Ceph - upload/download  sequences from the umiacs amazon S3 service

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
sub Ceph_Upload {
    ## ?
}



=item C<Dump_Reads>

    $hpgl->Dump_Reads()
    dumps reads from the given hpgl directory into the current working
    directory into concatenated files.  This is in opposition to
    Download_Reads which keeps them as separate files.

=cut
sub Dump_Reads {
    my $me = shift;
    my %args = @_;
    $me->Check_Options(args => \%args, needed => ['hpgl',]);
    my $casava_passed = 0;
    my $line_count = 0;
    my $connection;

    sub main {
        %options = Parse_Options();
        $connection = Connect();
        Download_Files();
        my $seq = $line_count / 4;
        print "$seq sequences were downloaded, of which $casava_passed passed casava's filtering.\n";
    }

    sub Connect {
        my $connection = new Net::Amazon::S3({aws_access_key_id => $options{access_key},
                                              aws_secret_access_key => $options{secret_key},
                                              host => $options{host},
                                              secure => 1,
                                              retry => 0,});
        die("The Amazon::S3 connection failed.") if (!defined($connection));
        return($connection);
    }


    sub Download_Files {
        my $fwd = "$options{bucket}_forward.fastq";
        my $rev = "$options{bucket}_reverse.fastq";
        if ($options{orientation} eq 'both') {
            Perform_Download('forward', $fwd);
            Perform_Download('reverse', $rev);
        } elsif ($options{orientation} eq 'forward') {
            Perform_Download('forward', $fwd);
        } elsif ($options{orientation} eq 'reverse') {
            Perform_Download('reverse', $rev);
        } else {
            die("I can only download forward reads, reverse reads, or both.");
        }
    }

    sub Perform_Download {
        my $orientation = shift;
        my $dest = shift;
        if (-r $dest) {
            print "$dest already exists, deleting it.\n" if ($options{verbose});
            unlink($dest);
        }

        my $response = $connection->buckets;
        my @buckets = @{$response->{buckets}};
        my $c = 0;
        my $found_bucket = 0;
      BUCKETS: foreach my $bucket (@buckets) {
          my $bucket_name = $bucket->bucket;
          if ($bucket_name eq $options{bucket}) {
              $found_bucket++;
              Download_Bucket($bucket, $dest, $orientation);

              if ($options{multi}) {
                INNER: foreach my $buck (@buckets) {
                  LETTERS: foreach my $letter ("a".."z") {
                      my $test = join("", $options{bucket}, $letter);
                      my $new_name = $buck->bucket();
                      print "Now testing for $new_name and $test\n";
                      if ($new_name eq $test) {
                          Download_Bucket($buck, $dest, $orientation);
                          next LETTERS;
                      } else {
                          next INNER;
                      }
                  } ## End rechecking the list of buckets for xxxa -> xxxz
                } ## End making the names xxxa -> xxxz
              } ## End if multi-buckets are on
          } else { ## End if we find hpgl0xxx
              next BUCKETS;
          }
      }
        if ($found_bucket == 0) {
            die("Did not find the bucket: $options{bucket}");
        }
    }

    sub Download_Bucket {
        my $bucket = shift;
        my $dest = shift;
        my $orientation = shift;
        open(OUT, ">>$dest");
        my $test = $bucket->list_all or die $connection->err . ": " . $connection->errstr;
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
                print "$tmp_dest already exists, deleting it.\n" if ($options{verbose});
                unlink($tmp_dest);
            }
            print "Downloading $file.\n";
            $bucket->get_key_filename($file, undef, $tmp_dest);
            my $md5_file = "$file.md5";
            my $md5_dest = "$tmp_dest.md5";
            if (-r $md5_dest) {
                print "$md5_dest already exists, deleting it.\n" if ($options{verbose});
                unlink($md5_dest);
            }
            $bucket->get_key_filename($md5_file, undef, $md5_dest);
            print "Checking md5 Sum\n";
            my $m = Check_MD5($tmp_dest, $md5_dest);
            my @suffixes = (".gz",".bz2",".xz");
            my $inter_dest = basename($file, @suffixes);
            if ($file =~ /\.fastq$/) {  ## Then it wasn't compressed
                print "The file was not compressed.\n";
            } else {
                if (-r $inter_dest) {
                    print "$inter_dest already exists, deleting it.\n" if ($options{verbose});
                    unlink($inter_dest);
                }
                print "Extracting $tmp_dest.\n";
                my $ae = new Archive::Extract(archive => $tmp_dest);
                my $ok = $ae->extract(to => $inter_dest);
            }
            print "Concatenating $inter_dest onto $dest.\n" if ($options{verbose});
            open(IN, "<:mmap", $inter_dest);
            while(<IN>) {
                $line_count++;
                $casava_passed++ if ($_ =~ m/^@.* [^:]*:N:[^:]*:/);
                print OUT;
            }
            close(IN);
            print "Cleaning up $inter_dest and $tmp_dest.\n" if ($options{verbose});
            unlink($inter_dest) if (-r $inter_dest);
            unlink($tmp_dest) if (-r $tmp_dest);
        }
        close(OUT);
    }

    sub List {
        if (!$options{bucket}) {
            List_Buckets();
        } else {
            my $response = $connection->buckets;
            my @buckets = @{$response->{buckets}};
            my $found_bucket = 0;
            foreach my $bucket (@buckets) {
                my $name = $bucket->bucket();
                if ($name eq $options{bucket}) {
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
                print "Did not find $options{bucket}\n";
                List_Buckets();
            }
        }
    }

    sub List_Buckets {
        print "Available buckets are:\n";
        my $response = $connection->buckets;
        my @buckets = @{$response->{buckets}};
        foreach my $bucket (@buckets) {
            my $name = $bucket->bucket();
            print "$name\n";
        }
    }

    sub Check_MD5 {
        my $downloaded_file = shift;
        my $md5_file = shift;
        open(MD5_IN, "<$md5_file");
        my $original_md5 = "";
        while(my $line = <MD5_IN>) {
            chomp($line);
            my ($md5, $filename) = split(/\s+/, $line);
            $original_md5 = $md5;
        }
        close(MD5_IN);
        unlink($md5_file);
        my $checksum = new Digest::MD5;
        open(MD5_CHECK, "<$downloaded_file");
        my $chk = \*MD5_CHECK;
        $chk = \*MD5_CHECK;  ## I don't understand how this works...
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
    }

    sub Parse_Options {
        my %options = (
            bucket => undef,
            orientation => 'both',
            debug => 0,
            md5 => 0,
            sample => 0,
            help => 0,
            list => 0,
            file => undef,
            verbose => 1,
            multi => 0,
            host => $ENV{CEPH_HOST},
            access_key => $ENV{CEPH_RO_ID},
            secret_key => $ENV{CEPH_RO_KEY},
            );
        my $opt = GetOptions(
            "bucket|b:s" => \$options{bucket},
            "orientation|o:s" => \$options{orientation},
            "debug|d" => \$options{debug},
            "md5|m" => \$options{md5},
            "sample|s" => \$options{sample},
            "help|h" => \$options{help},
            "list|l" => \$options{list},
            "file|f:s" => \$options{file},
            "verbose|v:s" => \$options{verbose},
            "multi|m:s" => \$options{multi},
            );

        if ($options{help}) {
            print "This script is intended to make downloading and viewing reads from ceph easier.
An example invocation might be:
> dump_reads.pl -b hpgl0223 -o r
That will dump the reverse (-o r) reads from the 'bucket' hpgl0223 (forward is default)
> dump_reads.pl -l
This will list the available buckets.
";
            exit(0);
        }
        if ($options{orientation} =~ /^r|^R/) {
            $options{orientation} = 'reverse';
        } elsif ($options{orientation} =~ /^f|^F/) {
            $options{orientation} = 'forward';
        } else {
            $options{orientation} = 'both';
        }
        if ($options{list}) {
            List();
            exit(0);
        }

        ## First ensure that we have all the information required to download sequence files.
        my $term = new Term::ReadLine('>');
        my $attribs = $term->Attribs;
        $attribs->{completion_suppress_append} = 1;
        my $OUT = $term->OUT || \*STDOUT;
        if (!defined($options{bucket})) {
            $options{bucket} = $term->readline("Please provide a hpgl to download (lowercase):");
            $options{bucket} =~ s/\s+$//g;
        }
        return(%options);
    }
}

1;
