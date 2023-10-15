package Bio::Adventure::Structure;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';
use feature 'try';
no warnings 'experimental::try';

use Bio::SeqFeature::Generic;
use Bio::SeqIO;
use Capture::Tiny qw":all";
use Cwd qw"abs_path getcwd cwd";
use File::Basename;
use File::Spec;
use File::Path qw"make_path rmtree";
use File::Which qw"which";
use File::ShareDir qw":ALL";
use IPC::Open2;
use Spreadsheet::Read;
use Symbol qw"gensym";

=head2 C<RNAFold_Windows>

 Run RNAFold on a rolling window across a sequence.
 10.1093/nar/gkn188

=cut
sub RNAFold_Windows {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        jcpu => 1,
        jprefix => 80,
        jname => 'vienna',
        length => 201,
        required => ['input'],
        step => 3,);
    my $output_name = basename($options->{input}, ('.gbk', '.fsa', '.fasta',));
    my $output_dir = qq"outputs/$options->{jprefix}rnafold";
    my $output = qq"${output_dir}/${output_name}.tsv.xz";
    my $comment = '## Iterate over a sequence with RNAfold.';
    my $jstring = qq?
use Bio::Adventure::Structure;
\$h->Bio::Adventure::Structure::RNAFold_Windows_Worker(
  input => '$options->{input}',
  length => '$options->{length}',
  step => '$options->{step}',
  output => '${output}',
  jprefix => '$options->{jprefix}',);
?;
    my $folder = $class->Submit(
        input => $options->{input},
        output => ${output},
        jname => 'vienna',
        jprefix => $options->{jprefix},
        length => $options->{length},
        step => $options->{step},
        jstring => $jstring,
        comment => $comment,
        language => 'perl',);
    return($folder);
}

=head2 C<RNAFold_Windows_Worker>

 Does the actual work of rolling across a sequence and running rnafold.

=cut
sub RNAFold_Windows_Worker {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        length => 201,
        step => 3,
        jname => 'vienna',
        jprefix => 80,
        required => ['input', 'output'],
        modules => ['vienna']);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $input_paths = $class->Get_Paths($options->{output});
    ## Put the data here!  First key is location, second is structure/mfe.
    my $output = $options->{output};
    my $out_dir = dirname($output);
    my $log = FileHandle->new(">${out_dir}/rnafold.log");
    print $log "Starting rolling window fold of $options->{input} with length $options->{length} and step $options->{step}.\n";
    my $out_txt = basename($output, ('.xz'));
    my $out_txt_path = qq"${out_dir}/${out_txt}";
    make_path($out_dir);
    my $txt_writer = FileHandle->new(">${out_txt_path}");
    print $txt_writer "contig\tstart\tend\tA\tU\tG\tC\tGC\tpaired\tbp_percent\tmfe\tmfe_bp\tmfe_gc\tstructure\n";

    my $length = $options->{length};
    my $step = $options->{step};
    my $in = FileHandle->new("less $options->{input} |");
    my $seqio;
    if ($options->{input} =~ /\.fasta|\.fsa/) {
        $seqio = Bio::SeqIO->new(-format => 'fasta', -fh => $in);
    } elsif ($options->{input} =~ /\.gb|\.gbff|\.gbk|\.gbf/) {
        $seqio = Bio::SeqIO->new(-format => 'genbank', -fh => $in);
    } else {
        die("I do not understand the format for $options->{input}");
    }
    my $results = {};
    my $reader = gensym();
    my $writer = gensym();
    my $pid = open2($reader, $writer, 'RNAfold -t 4 --noPS');
    my $output_line;
  CONTIGS: while (my $seq = $seqio->next_seq) {
      ## Set up the first sequence and the stepper
      my $seq_length = $seq->length;
      my $circular = $seq;
      my $id = $seq->id;
      my $make_circular = $circular->is_circular(1);
      ## Make the sequence length divisible by 3 to simplify logic later
      if (($seq_length % 3) == 1) {
          $circular = $circular->trunc(1, ($seq_length - 1));
      } elsif (($seq_length % 3) == 2) {
          $circular = $circular->trunc(1, ($seq_length - 2));
      }
      my $post_length = $circular->length;
      ## Start the iterator a little bit before the beginning of the sequence.
      my $start = (1 - $length) + $step;
      my $end = $start + $length;
      my $continue = 1;
      my %nt_counts = (A => 0, U => 0, G => 0, C => 0);
      my $gc_content = 0;
      my $bp_count = 0;
      my $bp_percent = 0;
      my $key = qq"${id}_${start}_${end}";
    STEP: while ($continue) {
        if ($end > $post_length) {
            $start = $start - $post_length;
            $end = $end - $post_length;
        }

        my $sub_sequence = '';
        if ($start <= 0 && $end <= 0) {
            ## This should not happen, but we can handle it easily enough
            print "Somehow got a zero on start or end, this should not happen\n";
        } elsif ($start <= 0) {
            my $pre_start = $post_length + $start;
            my $pre_end = $post_length;
            my $remaining = ($length + $start);
            my $post_start = 1;
            my $post_end = $post_start + $remaining;
            my $pre = $circular->subseq($pre_start, $pre_end);
            my $post = $circular->subseq($post_start, $post_end);
            $sub_sequence = $pre . $post;
        }   elsif ($end <= 0) {
            ## This really shouldn't happen.
            print "end is less than zero, that is weird.\n";
        } else  {
            $sub_sequence = $circular->subseq($start, $end);
        }

        ## Make sure to send a return character to RNAfold
        $sub_sequence .= "\n";
        ## Send the subsequence of interest to RNAfold.
        print $writer $sub_sequence;
        ## And read its first line of output.
        $output_line = <$reader>;
        ## If the line is comprised entirely of word characters, then it is printing the sequence.
        if ($output_line =~ /^\w+$/) {
            my $seq_line = $output_line;
            ## Count up the nucleotides.
            $nt_counts{A} = $seq_line =~ tr/A//;
            $nt_counts{U} = $seq_line =~ tr/U//;
            $nt_counts{G} = $seq_line =~ tr/G//;
            $nt_counts{C} = $seq_line =~ tr/C//;
            $gc_content = sprintf("%0.4f", ($nt_counts{G} + $nt_counts{C}) / $options->{length});
            ## Now skip down to the next line of output from RNAfold.
            $output_line = <$reader>;
        }

        ## Pick out the basepairs and calculated MFE value.
        my ($structure, $mfe) = split(/\s+/, $output_line);
        ## Count up the base pairs
        my $structure_string = $structure;
        $bp_count = $structure_string =~ tr/\.//;
        $bp_percent = sprintf("%0.4f", $bp_count / $options->{length});

        $mfe =~ s/\(|\)//g;
        my $normalized_mfe_bp = sprintf("%0.4f", $mfe / $bp_count);
        my $normalized_mfe_gc = sprintf("%0.4f", $mfe / ($nt_counts{G} + $nt_counts{C}));
        my $txt_string = qq"${id}\t${start}\t${end}\t$nt_counts{A}\t$nt_counts{U}\t$nt_counts{G}\t$nt_counts{C}\t$gc_content\t$bp_count\t$bp_percent\t${mfe}\t${normalized_mfe_bp}\t${normalized_mfe_gc}\t${structure}\n";
        print $txt_writer $txt_string;
        ## print $txt_string;
        ## Make a record of the contig/start/end so we can tell when we are finished.
        $results->{$key} = 1;

        ## Done with the step, set new start/end/key and check to see if we are done.
        $start = $start + $step;
        $end = $end + $step;
        $key = qq"${id}_${start}_${end}";
        if (defined($results->{$key})) {
            last STEP;
        }
    }
  }
    $writer->close();
    my $compressed = qx"xz -9e -f ${out_txt_path}";
    return($results);
}

=head2 C<Suppa>

  Set up and invoke Suppa, a differential transcript and splicing
  event calculator.  10.1186/s13059-018-1417-1

  I am putting this function in this file because it is relevant to
  the gene structure, which is admittedly a bit of a stretch.

  This function should be able to invoke the tool which creates the
  catalog of splicing events (suppa.py generateEvents), create the
  table of TPM values by transcript.  I am thinking to use the sample
  sheet as the input variable with a parameter containing the
  salmon/etc filename and have this generate the associated table.
  Given that information, this should be able to run the psiPerIsoform
  function and/or psiPerEvent.  Note, it turns out that suppa provides
  a helper function for creating the expression table 'joinFiles'; so
  that will simplify things quite a bit.

=cut
sub Suppa {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        type => 'SE',
        condition_column => 'drug',
        file_column => 'hg38100salmonfile',
        jprefix => '90',);
    my $suppa_dir = qq"outputs/$options->{jprefix}suppa_$options->{species}";
    my $gff_file = qq"$options->{libpath}/$options->{libtype}/$options->{species}.gtf";
    my $reader = Spreadsheet::Read->new($options->{input});
    my $sheet = $reader->sheet(1);
    ## Strangely, Spreadsheet::Read only uses numeric values for columns, so grab the first row
    ## and get the number of my column from it...
    my @column_names = $sheet->cellrow(1);
    print "TESTME: @column_names\n";
    my $file_number = undef;
    my $condition_number = undef;
    my $count = 0;
  COLUMNS: for my $name (@column_names) {
      $count++;
      if ($name eq $options->{file_column}) {
          print "Found $options->{file_column} as number: $count\n";
          $file_number = $count;
      }
      if ($name eq $options->{condition_column}) {
          print "Found $options->{condition_column} as number: $count\n";
          $condition_number = $count;
      }
  }
    my @filenames = $sheet->cellcolumn($file_number);
    my @conditions = $sheet->cellcolumn($condition_number);
    my $files_by_condition = {};
    my $cond_count = -1;
    for my $cond (@conditions) {
        print "Working on ${cond}\n";
        $cond_count++;
        my @cond_filenames;
        if (defined($files_by_condition->{$cond})) {
            @cond_filenames = @{$files_by_condition->{$cond}};
        } else {
            @cond_filenames = ();
        }
        my $this_file = $filenames[$cond_count];
        push(@cond_filenames, $this_file);
        print "TESTME: @cond_filenames\n";
        $files_by_condition->{$cond} = \@cond_filenames;
    }
    my $join_string = '';
    my $psi_string = '';
    my $psi_local_string = '';
    print "Creating joinfiles, psiperisoform, and psiperevent strings.\n";
    for my $cond (sort keys %{$files_by_condition}) {
        my $file_string = '';
        for my $f (@{$files_by_condition->{$cond}}) {
            $file_string .= qq" ${f} ";
        }
        $join_string .= qq"suppa.py joinFiles -f tpm -i ${file_string} -o ${cond}.tpm
";
        print "TESTME: $join_string\n";
        $psi_string .= qq"suppa.py psiPerIsoform -g ${gff_file} -e ${cond}.tpm -o ${cond}.psi
";
        $psi_local_string .= qq"suppa.py psiPerEvent --ioe-file local_as_events.txt --expression-file ${cond}.tpm -o ${cond}_local.psi
";
    }

    my $file_string = '';
    for my $f (@filenames) {
        $file_string .= " ${f} ";
    }
    my $jstring = qq!mkdir -p ${suppa_dir}
suppa.py generateEvents -i ${gff_file} -f ioe -e SE -o local_as_events_SE.txt
suppa.py generateEvents -i ${gff_file} -f ioe -e SS -o local_as_events_SS.txt
suppa.py generateEvents -i ${gff_file} -f ioe -e MX -o local_as_events_MX.txt
suppa.py generateEvents -i ${gff_file} -f ioe -e RI -o local_as_events_RI.txt
suppa.py generateEvents -i ${gff_file} -f ioe -e FL -o local_as_events_FL.txt
suppa.py generateEvents -i ${gff_file} -f ioi -o transcript_events.txt
${join_string}
${psi_string}
${psi_local_string}
!;
    my $suppa_job = $class->Submit(
        jstring => $jstring,
        );
    return($suppa_job);
}

1;
