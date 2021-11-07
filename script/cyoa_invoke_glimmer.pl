#!/usr/bin/env perl
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Cwd;

use File::Basename qw"dirname basename";
use File::Path qw"make_path remove_tree";
use File::Which qw"which";
use Getopt::Long qw"GetOptionsFromArray";

use Bio::FeatureIO;
use Bio::DB::Fasta;
use Bio::Tools::GFF;
use Bio::Adventure;

## This script will attempt to invoke glimmer3 in a sane fashion,
## putting the various outputs into a consistent format and returning
## a final multi-fasta of the detected ORFs

## The commands run follow the suggestions at:
## https://ccb.jhu.edu/software/glimmer/glim302notes.pdf

## Here is a shell script which basically performs everything:

=head1 EXAMPLE

mkdir -p outputs/glimmer
long-orfs -n -t 1.15 outputs/shovill_r1_trimmed-corrected/contigs.fa outputs/glimmer/first_run_longorfs.txt
extract outputs/shovill_r1_trimmed-corrected/contigs.fa outputs/glimmer/first_run_longorfs.txt \
  1>outputs/glimmer/first_run_training.txt 2>outputs/glimmer/first_run_training.err
build-icm -r outputs/glimmer/first_run.icm < outputs/glimmer/first_run_training.txt \
  2>outputs/glimmer/first_run_icm.out 1>&2
glimmer3 -o50 -g110 -t30 \
  outputs/shovill_r1_trimmed-corrected/contigs.fa \
  outputs/glimmer/first_run.icm \
  outputs/glimmer/first_run.out

## Use this first run output to go again

extract -t outputs/shovill_r1_trimmed-corrected/contigs.fa \
        outputs/glimmer/first_run.out.predict train.coords \
        1>outputs/glimmer/second_run_training.txt \
        2>outputs/glimmer/second_run_training.err
build-icm -r outputs/glimmer/second_run.icm < outputs/glimmer/second_run_training.txt
upstream-coords.awk 25 0 outputs/glimmer/first_run.out.predict | extract outputs/shovill_r1_trimmed-corrected/contigs.fa - > \
  outputs/glimmer/second_run_upstream.txt
elph outputs/glimmer/second_run_upstream.txt LEN=6 2>outputs/glimmer/elph.err | \
    get-motif-counts.awk 2>outputs/glimmer/second_run_motif.txt 1>&2

startuse=$(start-codon-distrib -3 outputs/shovill_r1_trimmed-corrected/contigs.fa outputs/glimmer/first_run.out.predict)
glimmer3 -o50 -g110 -t30 -b outputs/glimmer/second_run_motif.txt -P ${startuse} \
  outputs/shovill_r1_trimmed-corrected/contigs.fa outputs/glimmer/second_run.icm \
  outputs/glimmer/second_run.out \
  2>outputs/glimmer/second_run.err 1>&2

=cut

## So yeah, that is a bit of a PITA.  It is difficult to follow, the
## intermediate steps are poorly documented, and the final output is
## not trivially usable by anything.  So let us try and fix it up a
## little.

our $cyoa = new Bio::Adventure;
my $loaded = $cyoa->Module_Loader(modules => 'glimmer');
my $done = Run_Glimmer($cyoa);
$loaded = $cyoa->Module_Loader(modules => 'glimmer',
                               action => 'unload');

sub Run_Glimmer {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        long_t => '1.15',
        glimmer_o => '50',
        glimmer_g => '110',
        glimmer_t => '30',
        upstream_args => '25 0',
        elph_len => '6',
        jprefix => '16',
        );
    my $input_base = basename(dirname($options->{input}));
    my $outdir = qq"$options->{basedir}/outputs/$options->{jprefix}glimmer_${input_base}";
    make_path($outdir) unless(-d $outdir);
    my $first_run = First_Run(
        outdir => $outdir,
        input => $options->{input},
        long_t => $options->{long_t},
        glimmer_o => $options->{glimmer_o},
        glimmer_g => $options->{glimmer_g},
        glimmer_t => $options->{glimmer_t},
        );
    my $second_run = Second_Run(
        outdir => $outdir,
        input => $options->{input},
        predict => $first_run->{glimmer},
        upstream_args => $options->{upstream_args},
        elph_len => $options->{elph_len},
        glimmer_o => $options->{glimmer_o},
        glimmer_g => $options->{glimmer_g},
        glimmer_t => $options->{glimmer_t},
        );

    ## The filename returned by Second_Run() is almost, but not quite useful.
    my $gff = Make_Predict_Useful(
        input => $options->{input},
        glimmer_output => $second_run->{glimmer},
        outdir => $outdir);
}

sub First_Run {
    my %args = @_;
    my $run = 'first';
    my $long_orfs = Invoke_Longorfs(
        t => $args{long_t},
        run => $run,
        outdir => $args{outdir},
        input => $args{input},
        );
    print "The ${run} run long orfs file is: ${long_orfs}\n";
    my $extracted = Invoke_Extract(
        run => $run,
        outdir => $args{outdir},
        input => $args{input},
        long_orfs => $long_orfs,
        name => 'training',
        );
    print "The ${run} run extracted file is: ${extracted}\n";
    my $icm = Invoke_Buildicm(
        run => $run,
        outdir => $args{outdir},
        input => $extracted,
        );
    print "The ${run} run icm file is: ${icm}\n";
    my $glimmer = Invoke_Glimmer3(
        run => $run,
        outdir => $args{outdir},
        input => $args{input},
        icm => $icm,
        t => $args{glimmer_t},
        g => $args{glimmer_g},
        o => $args{glimmer_o},
        );
    print "The ${run} run prediction file is: ${glimmer}\n";
    my $ret = {
        run => $run,
        longorfs => $long_orfs,
        extracted => $extracted,
        icm => $icm,
        glimmer => $glimmer,
    };
    return($ret);
}

sub Second_Run {
    my %args = @_;
    my $run = 'second';
    my $extracted = Invoke_Extract(
        run => $run,
        outdir => $args{outdir},
        input => $args{input},
        long_orfs => $args{predict},
        extract_args => '-t',
        name => 'training',
        );
    print "The ${run} run extracted file is: ${extracted}\n";
    my $icm = Invoke_Buildicm(
        run => $run,
        outdir => $args{outdir},
        input => $extracted,
        );
    print "The ${run} run icm file is: ${icm}\n";
    my $upstream = Invoke_Upstream(
        run => $run,
        input => $args{input},
        outdir => $args{outdir},
        predict => $args{predict},
        upstream_args => $args{upstream_args},
        );
    my $elph = Invoke_Elph(
        run => $run,
        outdir => $args{outdir},
        upstream => $upstream,
        elph_len => $args{elph_len},
        );
    my $codon_distribution = Invoke_Startdistrib(
        outdir => $args{outdir},
        input => $args{input},
        predict => $args{predict}
        );
    my $glimmer = Invoke_Glimmer3(
        run => $run,
        outdir => $args{outdir},
        input => $args{input},
        icm => $icm,
        t => $args{glimmer_t},
        g => $args{glimmer_g},
        o => $args{glimmer_o},
        P => $codon_distribution,
        b => $elph,
        );
    print "The ${run} run prediction file is: ${glimmer}\n";
    my $ret = {
        extracted => $extracted,
        icm => $icm,
        glimmer => $glimmer,
    };
    return($ret);
}


sub Invoke_Longorfs {
    my %args = @_;
    my $program = 'long-orfs';
    my $check = which($program);
    die("Could not find ${program} in your PATH.") unless($check);
    my $t = '1.15';
    if ($args{t}) {
        $t = $args{t};
    }

    my $longorfs_outfile = qq"$args{outdir}/$args{run}_longorfs.txt";
    my $argstring = qq" -n -t ${t} $args{input} ${longorfs_outfile} ";
    my $stdout_err = $longorfs_outfile;
    $stdout_err =~ s/\.txt$/.out/g;

    my $shell_fh = FileHandle->new;
    print "Longorfs: About to run: ${program} ${argstring}\n";
    my $longorfs_stdout = FileHandle->new(">${stdout_err}");
    my $longorfs_pid = open($shell_fh, "${program} ${argstring} >$stdout_err 2>&1 |");
    while (my $line = <$shell_fh>) {
        print $longorfs_stdout $line;
    }
    $longorfs_stdout->close();
    $shell_fh->close();

    ## This is a good place to put some of the information
    ## printed by longorfs and returning it.
    return($longorfs_outfile);
}

sub Invoke_Extract {
    my %args = @_;
    my $program = 'extract';
    my $check = which($program);
    die("Could not find ${program} in your PATH.") unless($check);

    my $extract_args = qq"";
    if ($args{extract_args}) {
        $extract_args = $args{extract_args};
    }

    my $stdout_file = qq"$args{outdir}/$args{run}_$args{name}.txt";
    my $stderr_file = qq"$args{outdir}/$args{run}_$args{name}.stderr";
    my $outstring = qq" 1>${stdout_file} 2>${stderr_file}";
    my $argstring = qq" ${extract_args} $args{input} $args{long_orfs} ${outstring}";

    my $shell_fh = FileHandle->new;
    print "Extract: About to run: ${program} ${argstring}\n";
    my $shell_pid = open($shell_fh, "${program} ${argstring} |");
    while (my $line = <$shell_fh>) {
        chomp($line);
        print $line;
    }
    $shell_fh->close();

    ## This is a good place to put some of the information
    ## printed by longorfs and returning it.
    return($stdout_file);
}

sub Invoke_Buildicm {
    my %args = @_;
    my $program = 'build-icm';
    my $check = which($program);
    die("Could not find ${program} in your PATH.") unless($check);

    my $output_file = qq"$args{outdir}/$args{run}.icm";
    my $argstring = qq" -r ${output_file} < $args{input}";

    my $shell_fh = FileHandle->new;
    print "Build-icm: About to run: ${program} ${argstring}\n";

    my $stdout_file = qq"$args{outdir}/$args{run}_build-icm.stdout";
    my $stderr_file = qq"$args{outdir}/$args{run}_build-icm.stderr";
    my $outstring = qq" 1>${stdout_file} 2>${stderr_file}";

    my $shell_pid = open($shell_fh, "${program} ${argstring} ${outstring}|");
    while (my $line = <$shell_fh>) {
        chomp($line);
        print $line;
    }
    $shell_fh->close();

    ## This is a good place to put some of the information
    ## printed by longorfs and returning it.
    return($output_file);
}


sub Invoke_Glimmer3 {
    my %args = @_;
    my $program = 'glimmer3';
    my $check = which($program);
    die("Could not find ${program} in your PATH.") unless($check);

    my $P_arg = "";
    if ($args{P}) {
        $P_arg = qq" -P $args{P} ";
    }

    my $b_arg = "";
    if ($args{b}) {
        $b_arg = qq" -b $args{b} ";
    }

    my $output_file = qq"$args{outdir}/$args{run}_glimmer.out";
    my $argstring = qq" -o$args{o} -g$args{g} -t$args{t} ${b_arg} ${P_arg} $args{input} $args{icm} ${output_file} ";

    my $shell_fh = FileHandle->new;
    print "Glimmer3: About to run: ${program} ${argstring}\n";

    my $stdout_file = qq"$args{outdir}/$args{run}_glimmer.stdout";
    my $stderr_file = qq"$args{outdir}/$args{run}_glimmer.stderr";
    my $outstring = qq" 1>${stdout_file} 2>${stderr_file}";

    my $shell_pid = open($shell_fh, "${program} ${argstring} ${outstring}|");
    while (my $line = <$shell_fh>) {
        chomp($line);
        print $line;
    }
    $shell_fh->close();

    my $predict_file = qq"${output_file}.predict";
    return($predict_file);
}

sub Invoke_Upstream {
    my %args = @_;
    my $program = 'upstream-coords.awk';
    my $check = which($program);
    die("Could not find ${program} in your PATH.") unless($check);

    my $orfs_file = qq"$args{outdir}/$args{run}_longorfs.txt";
    my $output_file = qq"$args{outdir}/$args{run}_upstream.txt";
    my $stderr_file = qq"$args{outdir}/$args{run}_upstream.stderr";
    my $argstring = qq" $args{upstream_args} $args{predict} 2>$stderr_file 1>${orfs_file}";

    my $shell_fh = FileHandle->new;
    print "Upstream: About to run: ${program} ${argstring}\n";
    my $shell_pid = open($shell_fh, "${program} ${argstring}|");
    while (my $line = <$shell_fh>) {
        chomp($line);
        print $line;
    }
    $shell_fh->close();

    my $upstream_extracted = Invoke_Extract(
        input => $args{input},
        outdir => $args{outdir},
        run => $args{run},
        long_orfs => $orfs_file,
        name => 'upstream',
        );

    ## This is a good place to put some of the information
    ## printed by longorfs and returning it.
    return($upstream_extracted);
}

sub Invoke_Elph {
    my %args = @_;
    my $program = 'elph';
    my $check = which($program);
    die("Could not find ${program} in your PATH.") unless($check);

    my $orfs_file = qq"$args{outdir}/$args{run}_elph.txt";
    my $output_file = qq"$args{outdir}/$args{run}_motif.txt";
    my $stderr_file = qq"$args{outdir}/$args{run}_motif.stderr";
    my $argstring = qq" $args{upstream} LEN=$args{elph_len} 2>$args{outdir}/$args{run}_elph.stderr | get-motif-counts.awk 2>${stderr_file} 1>${output_file}";

    my $shell_fh = FileHandle->new;
    print "Elph: About to run: ${program} ${argstring}\n";
    my $shell_pid = open($shell_fh, "${program} ${argstring}|");
    while (my $line = <$shell_fh>) {
        chomp($line);
        print $line;
    }
    $shell_fh->close();

    ## This is a good place to put some of the information
    ## printed by longorfs and returning it.
    return($output_file);
}


sub Invoke_Startdistrib {
    my %args = @_;
    my $program = 'start-codon-distrib';
    my $check = which($program);
    die("Could not find ${program} in your PATH.") unless($check);

    my $argstring = qq" -3 $args{input} $args{predict} 2>$args{outdir}/startdistrib.err";

    my $shell_fh = FileHandle->new;
    print "Startdistrib: About to run: ${program} ${argstring}\n";
    my $shell_pid = open($shell_fh, "${program} ${argstring}|");
    my $out_string = qq"";
    while (my $line = <$shell_fh>) {
        chomp($line);
        $out_string .= $line;
    }
    $shell_fh->close();

    ## This is a good place to put some of the information
    ## printed by longorfs and returning it.
    return($out_string);
}


sub Make_Predict_Useful {
    my %args = @_;
    my $outdir = $args{outdir};
    my $in = FileHandle->new("<$args{glimmer_output}");
    my $in_fasta = Bio::DB::Fasta->new($args{input});
    my @ids = $in_fasta->get_all_primary_ids;
    my $out_fasta = Bio::SeqIO->new(-file => qq">${outdir}/predicted_cds.fasta", -format => 'fasta');
    my $out_gff = Bio::Tools::GFF->new(-file => qq">${outdir}/predicted_cds.gff", -gff_version => 3);

    my $current_contig = '';
    my $contig_id = '';
  ENTRIES: while (my $line = <$in>) {
      chomp $line;
      if ($line =~ /^>/) {
          ## Then this tells us the contig
          $current_contig = $line;
          $contig_id = $line;
          $contig_id =~ s/^>(\w+)?\s+(.*)$/$1/g;
          if ($contig_id =~ /^>/) {
              $contig_id =~ s/^>//g;
          }
          next ENTRIES;
      }
      ## Example prediction
      ## orf00001   109351       25  +2    12.06
      my ($orf_id, $start, $end, $frame, $score) = split(/\s+/, $line);
      my $seqobj = $in_fasta->get_Seq_by_id($contig_id);
      my $max_length = $seqobj->length;
      my $seq;
      my $strand = 1;
      if ($frame =~ /^\-/) {
          $strand = -1;
      }
      my $suffix = '';
      my $gff_start = 0;
      my $gff_end = 0;
      if ($strand == 1) {
          if ($end < $start) {
              $end = $max_length;
              $seq = $seqobj->subseq($start, $max_length);
              $suffix = 'passed end of contig';
          } else {
              $seq = $seqobj->subseq($start, $end);
          }
          $gff_start = $start;
          $gff_end = $end;
      } else {
          if ($end > $start) {
              $end = 0;
              $seq = $seqobj->subseq($end, $start);
              $seq =~ tr/ATGC/TACG/;
              $seq = reverse($seq);
              $suffix = 'passed end of contig';
          } else {
              $seq = $seqobj->subseq($end, $start);
              $seq =~ tr/ATGC/TACG/;
              $seq = reverse($seq);
          }
          $gff_start = $end;
          $gff_end = $start;
      }
      my $fasta_id = qq"${contig_id}_${orf_id} ";
      my $fasta_desc = qq"${start} ${end} ${suffix}";
      $fasta_id =~ s/\s+//g;

      my $seq_object = Bio::PrimarySeq->new(
          -id => $fasta_id, -seq => $seq,
          -description => $fasta_desc,
          -is_circular => 0,
          -alphabet => 'dna');
      $out_fasta->write_seq($seq_object);

      my $start_phase = $frame;
      $start_phase =~ s/^.{1}(\d)$/$1/g;
      my $phase = $start_phase - 1;
      my $gff_feature = Bio::SeqFeature::Lite->new(
          -seq_id => $contig_id,
          -name => $orf_id,
          -start => $gff_start,
          -end => $gff_end,
          -score => $score,
          -phase => qq"$phase",
          -strand => $strand,
          -type => 'glimmer',
          -attributes => { comment => $suffix, glimmer_phase => $start_phase },
          );
      $out_gff->write_feature($gff_feature);
  } ## End iterating over every line of the glimmer3 output.

    $in->close();
    ## close($in_fasta);
    ## close($out_fasta);
    ## close($out_gff);
}
