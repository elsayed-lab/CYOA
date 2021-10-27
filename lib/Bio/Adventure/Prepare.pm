package Bio::Adventure::Prepare;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use Text::CSV;

=head1 NAME

Bio::Adventure::Prepare - Get raw data ready for processing.

=head1 SYNOPSIS

use Bio::Adventure;
my $hpgl = new Bio::Adventure;
$hpgl->Prepare(csv => 'all_samples.csv');

=head1 Methods

=head2 C<Download_NCBI_Accession>

Given an accession, get the fasta/genbank/etc file.

=cut
sub Download_NCBI_Accession {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jprefix => '11',);
    my $job_name = $class->Get_Job_Name();
    ## Make an array of the accession(s)
    my $unique = [$options->{input}, ];

    my $eutil = Bio::DB::EUtilities->new(-eutil => 'esummary',
                                         -email => 'abelew@gmail.com',
                                         -db => 'nucleotide',
                                         -id => $unique,);

    while (my $docsum = $eutil->next_DocSum) {
        my $acc_version = '';
        my $accession = '';
      ITEMS: while (my $item = $docsum->next_Item) {
          my $item_name = $item->get_name;
          if ($item_name eq 'AccessionVersion') {
              $acc_version = $item->get_content();
          } elsif ($item_name eq 'Caption') {
              $accession = $item->get_content();
          } else {
              next ITEMS;
          }
      } ## End checking the document summary

        ## Now check if we already have this file
        if (-r qq"${acc_version}.gb") {
            print "Already have: ${acc_version}.gb\n";
      } else {
          print "Downloading ${accession}\n";
          my $download = Bio::DB::EUtilities->new(-eutil => 'efetch',
                                                  -db => 'nucleotide',
                                                  -rettype => 'gb',
                                                  -email => 'abelew@umd.edu',
                                                  -id => $accession,);
          my $output_file = qq"${acc_version}.gb";
          $download->get_Response(-file => $output_file);
          sleep(1);

          my @current_files = glob(qq"${acc_version}*");
          my ($first_acc, $first_ver, $first_ext);
          my ($second_acc, $second_ver, $second_ext);
          if (scalar(@current_files) > 1) {
              ($first_acc, $first_ver, $first_ext) = split(/\./, $current_files[0]);
              ($second_acc, $second_ver, $second_ext) = split(/\./, $current_files[1]);
              if ($first_ver > $second_ver) {
                  unlink($current_files[1]);
              } else {
                  unlink($current_files[0]);
              }
          }
      } ## Finished checking if we already have this accession
    } ## Finished iterating over every phage ID

}

=head2 C<Fastq_Dump>

Invoke fastq_dump to download some data from sra.

=cut
sub Fastq_Dump {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        output => undef);

    my $fastq_comment = qq"## This script should download an sra accession to local fastq.gz files.
";
    my $job_basename = $class->Get_Job_Name();
    my @inputs = split(/\,/, $options->{input});
    my @outputs = ();
    my $first_output = undef;
    if (defined($options->{output})) {
        @outputs = split(/\,/, $options->{output});
        $first_output = $outputs[0];
    }

    my %fastq_jobs = ();
    my $count = 0;
    my $fastq_job;
    for my $i (0 .. $#inputs) {
        $count++;
        my $in = $inputs[$i];
        print "Invoking fastq-dump for ${in}.\n";
        my $jstring = "";
        if (defined($outputs[$i]) || defined($first_output)) {
            $outputs[$i] = $first_output if (!defined($outputs[$i]));
            $jstring = qq"mkdir -p $outputs[$i] && \\
  fastq-dump --outdir $outputs[$i] \\
    --gzip --skip-technical --readids \\
    --read-filter pass --dumpbase \\
    --split-3 --clip ${in} && \\
";
        } else {
            $jstring = qq"mkdir -p ${in} && \\
  fastq-dump --outdir ${in} --gzip --skip-technical --readids \\
    --read-filter pass --dumpbase \\
    --split-3 --clip ${in}
";
        }

        my $current_fastq_job = $class->Submit(
            comment => $fastq_comment,
            input => $in,
            jdepends => $options->{jdepends},
            jname => qq"fqd_${in}",
            jstring => $jstring,
            jprefix => "01",
            jmem => 12,
            jwalltime => '6:00:00',);

        if (defined($fastq_job)) {
            $fastq_job->{$count} = $current_fastq_job;
        } else {
            $fastq_job = $current_fastq_job;
        }
    } ## Foreach my $input
    return($fastq_job);
}

=head2 C<Read_Samples>

Read the sample csv file.

=cut
sub Read_Samples {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['csv_file']);
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

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

L<fastq-dump>

=cut

1;
