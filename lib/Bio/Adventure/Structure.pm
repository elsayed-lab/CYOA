package Bio::Adventure::Structure;
## LICENSE: gplv2
## ABSTRACT:  Kitty!
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
use Symbol qw"gensym";

=head2 C<RNAFold_Windows>

 Run RNAFold on a rolling window across a sequence.
 10.1093/nar/gkn188

=cut
sub RNAFold_Windows {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        length => 201,
        step => 3,
        jprefix => 80,
        jname => 'vienna',
        required => ['input'],
        modules => ['vienna']);
    my $output_name = basename($options->{input}, ('.gbk'));
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
        output => $options->{output},
        jname => 'vienna',
        jprefix => $options->{jprefix},
        length => $options->{length},
        step => $options->{step},
        jstring => $jstring,
        comment => $comment,
        language => 'perl',
        shell => '/usr/bin/env perl',);
    $class->{language} = 'bash';
    $class->{shell} = '/usr/bin/env bash';
    return($folder);
}

=head2 C<RNAFold_Windows_Worker

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
    my $check = which('RNAfold');
    die("Could not find RNAfold in your PATH.") unless($check);
    my $input_paths = $class->Get_Paths($options->{output});
    ## Put the data here!  First key is location, second is structure/mfe.
    my $output = $options->{output};
    my $out_dir = dirname($output);
    my $out_txt = basename($output, ('.xz'));
    my $out_txt_path = qq"${out_dir}/${out_txt}";
    make_path($out_dir);
    my $txt_writer = FileHandle->new(">${out_txt_path}");
    print $txt_writer "contig\tstart\tend\tA\tU\tG\tC\tGC\tpaired\tbp_percent\tmfe\tmfe_bp\tmfe_gc\tstructure\n";

    my $length = $options->{length};
    my $step = $options->{step};
    my $in = FileHandle->new("less $options->{input} |");
    my $seqio = Bio::SeqIO->new(-format => 'genbank', -fh => $in);
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

1;
