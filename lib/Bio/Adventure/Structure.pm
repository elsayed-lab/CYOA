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

sub RNAFold_Windows {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        length => 201,
        overlap => 198,
        required => ['input'],
        modules => ['vienna']);

    ## Put the data here!  First key is location, second is structure/mfe.
    my $results = {};

    my $length = $options->{length};
    my $overlap = $options->{overlap};
    my $step = $options->{length} - $options->{overlap};
    my $in = FileHandle->new("less $options->{input} |");
    my $seqio = Bio::SeqIO->new(-format => 'genbank', -fh => $in);

    my $reader = gensym();
    my $writer = gensym();
    my $pid = open2($reader, $writer, 'RNAfold');
    my $output_line;
  CONTIGS: while (my $seq = $seqio->next_seq) {
      ## Set up the first sequence and the stepper
      ## Make the sequence length divisible by 3 to simplify logic later
      my $seq_length = $seq->length;
      my $circular = $seq;
      my $id = $seq->id;
      my $make_circular = $circular->is_circular(1);
      if (($seq_length % 3) == 1) {
          $circular = $circular->trunc(1, ($seq_length - 1));
      } elsif (($seq_length % 3) == 2) {
          $circular = $circular->trunc(1, ($seq_length - 2));
      }
      my $post_length = $circular->length;
      my $start = 1 - $overlap;
      my $end = $start + $length;
      my $continue = 1;
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

        $sub_sequence .= "\n";
        ## Now start reading/writing to RNAfold.
        print $writer $sub_sequence;
        $output_line = <$reader>;
        unless ($output_line =~ /^\w+$/) {
            my ($structure, $mfe) = split(/\s+/, $output_line);
            $mfe =~ s/\(|\)//g;
            $results->{$key}->{structure} = $structure;
            $results->{$key}->{mfe} = $mfe;
            print "Set structure and mfe to $structure and $mfe\n";
        }

        ## Done with the step, set new start/end/key and check to see if we are done.
        $start = $start + $step;
        $end = $end + $step;
        $key = qq"${id}_${start}_${end}";
        if (defined($results->{$key})) {
            last STEP;
        }
    }
  }
    return($results);
}

1;
