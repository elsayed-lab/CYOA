package Bio::Adventure::Index;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use File::Basename;
use File::Copy qw"cp";
use File::Path qw"make_path";
use File::Spec;
use File::Which qw"which";

=head1 NAME

 Bio::Adventure::Index - Invoke the various indexing utilities.

=head1 SYNOPSIS

=head2 C<BT1_Index>

 Create a bowtie1 index.

 This requires an 'input' argument.  This is a fasta file which will be used
 to write indexes in the libdir/libtype/indexes directory.  It assumes that the
 input fasta file is named in a way which is human-readable as the species.

 Thus, if the input fasta is 'lmajor_v46.fasta' and libtype is 'genome' it will
 be copied to $libdir/genome/lmajor_v46.fasta and indexes will have the basename
 'lmajor_v46'.

 This also uses less(1) to send an uncompressed copy of the input file to
 libdir/libtype/species.fasta.

=cut
sub BT1_Index {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        modules => ['bowtie1'],);
    my $species = basename($options->{input}, ('.gz', '.bz2', '.xz'));
    $species = basename($species, ('.fasta', '.fa'));
    my $copied_location = qq"$options->{libpath}/$options->{libtype}/${species}.fasta";
    if (!-f $copied_location) {
        my $copied = qx"less $options->{input} > ${copied_location}";
    }
    my $output_dir = qq"$options->{basedir}/outputs/$options->{jprefix}bt1index";
    my $stdout = qq"${output_dir}/bt1_index.stdout";
    my $stderr = qq"${output_dir}/bt1_index.stderr";

    my $jstring = qq!mkdir -p ${output_dir}
bowtie-build $options->{input} \\
  $options->{libdir}/$options->{libtype}/indexes/${species} \\
  2>${stderr} 1>${stdout}
!;
    my $comment = qq!## Generating bowtie1 indexes for species: ${species}
## in $options->{libdir}/$options->{libtype}/indexes!;
    my $bt1_index = $class->Submit(
        comment => $comment,
        jname => qq"bt1idx_${species}",
        jdepends => $options->{jdepends},
        jstring => $jstring,
        jprefix => '10',
        output => $output_dir,
        stderr => $stderr,
        stdout => $stdout,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($bt1_index);
}

=head2 C<BT2_Index>

 Create a bowtie2 index.

 This requires an 'input' argument.  This is a fasta file which will be used
 to write indexes in the libdir/libtype/indexes directory.  It assumes that the
 input fasta file is named in a way which is human-readable as the species.

 Thus, if the input fasta is 'lmajor_v46.fasta' and libtype is 'genome' it will
 be copied to $libdir/genome/lmajor_v46.fasta and indexes will have the basename
 'lmajor_v46'.

 This also uses less(1) to send an uncompressed copy of the input file to
 libdir/libtype/species.fasta.

=cut
sub BT2_Index {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jprefix => '',
        modules => ['bowtie2'],);
    my $libtype = $options->{libtype};
    my $libdir = File::Spec->rel2abs($options->{libpath});
    my $species = basename($options->{input}, ('.gz', '.bz2', '.xz'));
    $species = basename($species, ('.fasta', '.fa'));
    my $copied_location = qq"$options->{libpath}/$options->{libtype}/${species}.fasta";
    if (!-f $copied_location) {
        my $copied = qx"less $options->{input} > ${copied_location}";
    }
    my $output_dir = qq"$options->{basedir}/outputs/$options->{jprefix}bt2_index";
    my $stdout = qq"${output_dir}/index.stdout";
    my $stderr = qq"${output_dir}/index.stderr";
    my $jstring = qq!mkdir -p ${output_dir}
bowtie2-build $options->{input} \\
  $options->{libdir}/${libtype}/indexes/${species} \\
  2>${stderr} 1>${stdout}
!;
    my $comment = qq!## Generating bowtie2 indexes for species: ${species}
## in $options->{libdir}/${libtype}/indexes!;
    my $indexer = $class->Submit(
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => qq"bt2idx_${species}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        modules => $options->{modules},
        output => $output_dir,
        stderr => $stderr,
        stdout => $stdout,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($indexer);
}

=head2 C<BWA_Index>

 Create bwa indexes.

=cut
sub BWA_Index {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jprefix => '',
        modules => ['bwa'],);
    my $species = basename($options->{input}, ('.fasta', '.fa'));
    my $copied_location = qq"$options->{libpath}/$options->{libtype}/${species}.fa";
    if (!-f $copied_location) {
        cp($options->{input}, $copied_location);
    }
    my $output_dir = qq"$options->{basedir}/outputs/$options->{jprefix}bwa_index";
    my $stdout = qq"${output_dir}/bwa_index.stdout";
    my $stderr = qq"${output_dir}/bwa_index.stderr";
    my $jstring = qq!mkdir -p ${output_dir}
start=\$(pwd)
cd $options->{libdir}/$options->{libtype}/indexes
ln -sf ../${species}.fa .
bwa index ${species}.fa \\
  2>${stderr} \\
  1>${stdout}
cd \$start
!;
    my $basedir = qq"$options->{libpath}/$options->{libtype}/indexes";
    my $index_sa = qq"${basedir}/${species}.fa.sa";
    my $index_pac = qq"${basedir}/${species}.fa.pac";
    my $index_bwt = qq"${basedir}/${species}.fa.bwt";
    my $index_ann = qq"${basedir}/${species}.fa.ann";
    my $index_amb = qq"${basedir}/${species}.fa.amb";
    my $index_fa = qq"${basedir}/${species}.fa";

    my $comment = qq!## Generating bwa indexes for species: ${species}
## in $options->{libdir}/$options->{libtype}/indexes!;
    my $bwa_index = $class->Submit(
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => qq"$options->{jprefix}bwaidx",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        modules => $options->{modules},
        output => $output_dir,
        output_sa => $index_sa,
        output_pac => $index_pac,
        output_bwt => $index_bwt,
        output_ann => $index_ann,
        output_amb => $index_amb,
        output_fa => $index_fa,
        stderr => $stderr,
        stdout => $stdout,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($bwa_index);
}

=head2 C<Check_Blastdb>

 Check_Blastdb makes sure that there is an appropriately formatted blastdb for
 the library and type of blast performed.

=cut
sub Check_Blastdb {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        blast_tool => 'blastn',
        type => 'prot',
        required => ['input'],
        modules => ['blast'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $libname = basename($options->{input}, ('.fasta', '.fa', '.faa', '.fsa', '.fna'));
    if (defined($options->{library})) {
        $libname = $options->{library};
    }
    ## First check for the relevant library in $ENV{BLASTDB}
    ## If it isn't there, make one in basedir/blastdb/
    my $foundlib = 0;
    if ($options->{type} ne 'prot' && $options->{type} ne 'nucl') {
        die(qw"makeblastdb requires either a type of 'prot' or 'nucl'.");
    }

    my $mismatches = 0;
    my $matches = 0;
    my $test_in = Bio::SeqIO->new(-file => $options->{input}, -format => 'Fasta');
    while (my $seq = $test_in->next_seq) {
        my $guess = $seq->alphabet;
        if ($guess eq 'protein') {
            $guess = 'prot';
        } else {
            $guess = 'nucl';
        }
        if ($options->{type} eq $guess) {
            $matches++;
        } else {
            $mismatches++;
        }
    }
    my $sum = $matches + $mismatches;
    print "Out of ${sum} sequences, ${matches} were guessed to be $options->{type} and ${mismatches} were not.\n";
    my $checklib = qq"${libname}.psq";
    my $checklib_zero = qq"${libname}.00.psq";
    if ($options->{type} eq 'nucl') {
        $checklib_zero = qq"${libname}.00.nsq";
        $checklib = qq"${libname}.nsq";
    }
    my $db_directory = $ENV{BLASTDB};
    my $lib = '';
    my $relative_directory = 'blastdb';
    if (!defined($ENV{BLASTDB})) {
        $ENV{BLASTDB} = "$options->{basedir}/blastdb";
        $db_directory = "$options->{basedir}/blastdb";
    } else {
        $relative_directory = qq"$ENV{BLASTDB}";
    }

    print "Looking for ${checklib} / ${checklib_zero} in either $ENV{BLASTDB} or $options->{basedir}/blastdb.\n";
    ## Start with BLASTDB
    if (-f "$ENV{BLASTDB}/${checklib}" or
        -f "$ENV{BLASTDB}/${checklib_zero}") {
        $foundlib++;
        $lib = qq"$ENV{BLASTDB}/${libname}";
        print "Found an existing blast database at ${lib}.\n";
    } else {
        print "Did not find an existing blast database.\n";
    }

    ## If we do not find the blast database, create it in the basedir.
    if (!$foundlib) {
        if (!-d qq"$options->{basedir}/blastdb") {
            make_path(qq"$options->{basedir}/blastdb");
        }
        my $formatdb_command = qq"makeblastdb \\
  -in $options->{input} \\
  -dbtype $options->{type} \\
  -out ${db_directory}/${libname} \\
  2>${db_directory}/makeblastdb.stderr \\
  1>${db_directory}/makeblastdb.stdout";
        print "The makeblastdb command run is: ${formatdb_command}\n";
        my $formatdb_ret = qx"${formatdb_command}";
    }
    my $final_directory = qq"${relative_directory}/${libname}";
    return($final_directory);
}

=head2 C<Extend_Kraken_DB>

 Add some more sequences to an existing kraken2 database.

=cut
sub Extend_Kraken_DB {
    my ($class, %args) = @_;
    my $check = which('kraken2');
    die('Could not find kraken2 in your PATH.') unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        library => 'viral',
        modules => ['kraken'],);
    ## kraken2 --db ${DBNAME} --paired --classified-out cseqs#.fq seqs_1.fq seqs_2.fq
    my $job_name = $class->Get_Job_Name();
    my $output_dir = qw"outputs/extend_kraken";

    my $stderr = qq"${output_dir}/kraken2-build.stderr";
    my $stdout = qq"${output_dir}/kraken2-build.stdout";
    my $comment = '## This is a script to extend an existing kraken2 library with some new sequences.';
    my $jstring = qq!mkdir -p ${output_dir}
kraken2-build --download-taxonomy --db \${KRAKEN_DB_PATH}/$options->{library} \\
              2>${stderr} \\
              1>${stdout}
kraken2-build --add-to-library $options->{input} --db \${KRAKEN_DB_PATH}/$options->{library} \\
              2>>${stderr} \\
              1>>${stdout}
kraken2-build --download-library $options->{library} \\
              --db \${KRAKEN_DB_PATH}/$options->{library} \\
              2>>${stderr} \\
              1>>${stdout}
kraken2-build --build --db \${KRAKEN_DB_PATH}/$options->{library} \\
              2>>${stderr} \\
              1>>${stdout}
!;
    my $kraken = $class->Submit(
        comment => $comment,
        jcpus => 6,
        jdepends => $options->{jdepends},
        jname => qq"kraken_${job_name}",
        jprefix => '99',
        jstring => $jstring,
        jmem => 96,
        modules => $options->{modules},
        output => qq"${output_dir}/kraken2-build.out",
        stderr => $stderr,
        stdout => $stdout,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => 'large',);
    return($kraken);
}

=head2 C<Hisat2_Index>

 Create a hisat2 index using ${species}.fasta and leaves it in the indexes/
 directory.

=cut
sub Hisat2_Index {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        output_dir => undef,
        required => ['input'],
        modules => ['hisat2'],
        jprefix => '21',);
    my $libtype = $options->{libtype};
    my $libdir = File::Spec->rel2abs($options->{libpath});
    my $species = basename($options->{input}, ('.fasta', '.fa', '.fsa'));
    my $copied_location = qq"$options->{libpath}/$options->{libtype}/${species}.fasta";
    my $stdout = qq"hisat2_index_${species}.stdout";
    my $stderr = qq"hisat2_index_${species}.stderr";
    my $output_dir;
    if (defined($options->{output_dir})) {
        $output_dir = $options->{output_dir};
    } else {
        $output_dir = qq"$options->{basedir}/$options->{jprefix}hisat2_index";
    }
    $stdout = qq"${output_dir}/${stdout}";
    $stderr = qq"${output_dir}/${stderr}";
    make_path($output_dir, {verbose => 0}) unless (-r $output_dir);

    my $copied = undef;
    if (-r $copied_location) {
        print "The indexes appear to exist at: ${copied_location}.\n";
    } else {
        print "Copying $options->{input} to ${copied_location}\n";
        $copied = cp($options->{input}, $copied_location);
    }
    my $jstring = qq!mkdir -p ${output_dir}
hisat2-build $options->{input} \\
  $options->{libdir}/${libtype}/indexes/${species} \\
  2>${stderr} \\
  1>${stdout}
!;
    my $comment = qq!## Generating hisat2 indexes for species: ${species}
## in $options->{libdir}/${libtype}/indexes!;
    my $indexer = $class->Submit(
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => qq"ht2idx_${species}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        output => $stdout,
        stderr => $stderr,
        stdout => $stdout,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($indexer);
}

=head2 C<Kallisto_Index

 Use kallisto and an annotated_CDS fasta sequence library to create an index.

=cut
sub Kallisto_Index {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        modules => ['kallisto'],
        jprefix => '45',
        required => ['input'],);
    my $cds = basename($options->{input}, ('.fasta', '.fa'));
    my $species = $cds;
    $species =~ s/_cds//g;
    $species =~ s/_nt//g;
    my $copied_location = qq"$options->{libpath}/$options->{libtype}/${cds}.fasta";
    if (!-f $copied_location) {
        cp($options->{input}, $copied_location);
    }
    my $libtype = $options->{libtype};
    my $output_dir = qq"$options->{basedir}/outputs/$options->{jprefix}kallisto_index";
    my $stdout = qq"${output_dir}/index.stdout";
    my $stderr = qq"${output_dir}/index.stderr";
    my $input = File::Spec->rel2abs($options->{input});
    my $jstring = qq!mkdir -i ${output_dir}
kallisto index -i $options->{libdir}/${libtype}/indexes/${species}.idx \\
  ${input} \\
  2>${stderr} 1>${stdout}
!;
    my $comment = qq!## Generating kallisto indexes for species: ${species}
## in $options->{libdir}/${libtype}/indexes!;
    my $ka_index = $class->Submit(
        comment => $comment,
        jdepends => $options->{jdepends},
        jstring => $jstring,
        jname => 'kalidx',
        jprefix => $options->{jprefix},
        modules => $options->{modules},
        stderr => $stderr,
        output => $output_dir,
        stdout => $stdout,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($ka_index);
}

=head2 C<Make_Codon_Table>

 Assist caical by creating a codon table for an arbitrary species in
 the expected format.

 Given a reference assembly, print out a codon table for use with a
 codon adapatation calculator.  I wrote a little piece of code which
 works from a fasta/gff pair, but currently it is dumb.  I ought to
 improve its logging and have it return a useful/interesting data
 structure.

=cut
sub Make_Codon_Table {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species'],
        jprefix => '80',);
    my $out_table = qq"$options->{libpath}/codon_tables/$options->{species}.txt";
    my $out_dir = dirname($out_table);
    unless (-d $out_dir) {
        my $made = make_path($out_dir);
    }

    ## I have a few suffixes for writing genbank files.
    my @potential_suffixes = ('gbff', 'gbk', 'gbf', 'gb', 'genbank');
    my @compressions = ('gz', 'xz', 'bz2');
    my $in_gbff = '';
  POTENTIAL: for my $potential (@potential_suffixes) {
      my $start = qq"$options->{libpath}/$options->{libtype}/$options->{species}.${potential}";
      if (-r $start) {
          $in_gbff = $start;
          last POTENTIAL;
      }
      for my $c (@compressions) {
          my $comp_test = qq"${start}.${c}";
          if (-r $comp_test) {
              $in_gbff = $comp_test;
              last POTENTIAL;
          }
      }
  }
    ## The table already exists, move on -- or maybe I should have it overwrite?
    if (-r $out_table) {
        return(1);
  }

    ## This is a little dumb, but an easy way to think through writing out the table.
    ## E.g. I will do a for(for(for())) over these to get the 64 codons.
    my @first = ('T', 'C', 'A', 'G');
    my @second = ('T', 'C', 'A', 'G');
    my @third = ('T', 'C', 'A', 'G');

    ## Set up the pieces which will hold the data of interest.
    my $total_codons = 0;
    my %codon_counts = ();
    my $in = FileHandle->new("less ${in_gbff} |");
    my $seqio = Bio::SeqIO->new(-format => 'genbank', -fh => $in);
    my $seq_count = 0;
  SEQ: while (my $seq = $seqio->next_seq) {
      $seq_count++;
      ## print "Starting sequence $seq_count\n";
      my @feature_list = $seq->get_SeqFeatures();
      my $f_count = 0;
    FEAT: for my $f (@feature_list) {
        next FEAT unless ($f->primary_tag eq 'CDS');
        $f_count++;
        ## print "Feature: $f_count\n";
        my $sequence_string = $f->seq->seq;
        my $trans = $f->seq->translate->seq;
        my @seq_str = split(//, $sequence_string);
        my @tr_str = split(//, $trans);
      CHOMP: while (scalar(@tr_str) > 0) {
          my $amino = shift @tr_str;
          my $nt1 = shift @seq_str;
          my $nt2 = shift @seq_str;
          my $nt3 = shift @seq_str;
          my $codon = qq"${nt1}${nt2}${nt3}";
          if (defined($codon_counts{$codon})) {
              $codon_counts{$codon}++;
          } else {
              $codon_counts{$codon} = 1;
          }
          $total_codons++;
      } ## Done pulling apart the sequence arrays.
    } ## Iterating over the features in this sequence.
  } ## End going through the sequences of the assembly.
    $in->close();

    if ($total_codons == 0) {
        print "This failed: $options->{species}\n";
        return(undef);
    }
    my $divisor = 1000.0 / $total_codons;
    my $table = FileHandle->new(">${out_table}");
    for my $f (@first) {
        for my $s (@second) {
            my $string = '';
            for my $t (@third) {
                my $codon = qq"${f}${s}${t}";
                my $per_thousand = sprintf("%.1f", $codon_counts{$codon} * $divisor);
                $string .= qq"${codon} ${per_thousand}($codon_counts{$codon})   ";

            }
            print $table qq"${string}\n";
        }
        print $table "\n";
    }
    $table->close();
}

=item C<RSEM_Index

 Use RSEM and an annotated_CDS fasta sequence library to create a transcript index.

=cut
sub RSEM_Index {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        modules => ['rsem', 'bowtie2']);

    my $species = basename($options->{input}, ('.fasta', '.fa'));
    $species =~ s/_cds//g;
    my $copied_location = qq"$options->{libpath}/$options->{libtype}/${species}.fasta";
    if (!-f $copied_location) {
        cp($options->{input}, $copied_location);
    }
    my $output_dir = qq"$options->{basedir}/outputs/$options->{jprefix}rsem_index";
    my $stdout = qq"${output_dir}/index.stdout";
    my $stderr = qq"${output_dir}/index.stderr";
    my $comment = '## RSEM Index creation.';
    my $jstring = qq!mkdir -p ${output_dir}
rsem-prepare-reference --bowtie2 $options->{input} ${species} \\
  2>${stderr} 1>${stdout}
!;
    my $jobid = $class->Submit(
        comment => $comment,
        jdepends => $options->{jdepends},
        jstring => $jstring,
        jname => 'rsemidx',
        jprefix => $options->{jprefix},
        modules => $options->{modules},
        stderr => $stderr,
        stdout => $stdout,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($jobid);
}

=head2 C<Salmon_Index>

 Invoke salmon with an annotated_CDS fasta sequence library to create a
 transcript index.  Note that newer version of salmon would like a set
 of decoys, which may be performed in a couple of ways, the second of
 which I am copy/pasting from the documentation.

 The second is to use the entire genome of the organism as the decoy
 sequence. This can be done by concatenating the genome to the end of
 the transcriptome you want to index and populating the decoys.txt
 file with the chromosome names. Detailed instructions on how to
 prepare this type of decoy sequence is available here. This scheme
 provides a more comprehensive set of decoys, but, obviously, requires
 considerably more memory to build the index

=cut
sub Salmon_Index {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        modules => ['salmon'],);
    my $libtype = $options->{libtype};
    my $genome = File::Spec->rel2abs($options->{input});

    my $cds = basename($options->{input}, ('.fasta', '.fa'));
    my $cds_dir = dirname($options->{input});
    my $species = $cds;
    ## Drop the suffixes which might be annoying.
    $species =~ s/_cds//g;
    $species =~ s/_nt//g;
    my $species_file = qq"${cds_dir}/${species}.fasta";
    my $copied_location = qq"$options->{libpath}/$options->{libtype}/${cds}.fasta";
    my $species_location = qq"$options->{libpath}/$options->{libtype}/${species}.fasta";

    if (!-r $copied_location) {
        cp($options->{input}, $copied_location);
    }
    my $decoy_copy_string = qq'';
    my $jstring = qq'';

    my $output_dir = qq"$options->{basedir}/outputs/$options->{jprefix}salmon_index";
    my $stdout = qq"${output_dir}/index.stdout";
    my $stderr = qq"${output_dir}/index.stderr";

    my $index_input = $options->{input};
    my $index_string = qq!
salmon index -t ${index_input} -i $options->{libdir}/${libtype}/indexes/${species}_salmon_index!;
    if (-f $species_location or -f $species_file) {
        if (!-f $species_location) {
            cp($species_file, $species_location);
        }
        my $decoy_location = qq"$options->{libdir}/${libtype}/${species}_decoys.fasta";
        $decoy_copy_string = qq!less $options->{input} > ${decoy_location} && less ${species_file} >> ${decoy_location}
less ${species_file} | { grep '^>' || test \$? = 1; } | sed 's/^>//g' >> ${decoy_location}.txt
!;
        $index_input = $decoy_location;
        $index_string = qq!mkdir -p ${output_dir}
salmon index \\
  -t ${index_input} \\
  -i $options->{libdir}/${libtype}/indexes/${species}_salmon_index \\
!;
        $jstring = qq!${decoy_copy_string}
${index_string} \\
  --decoys ${decoy_location}.txt \\
  2>${stderr} 1>${stdout}
!;
    } else {
        warn("This function would prefer to make a decoy aware index set which requires the full genome.");
        say("Waiting 10 seconds to see if you want to quit and gather that genome,
otherwise a decoy-less index will be generated.");
        sleep(10);
        $jstring = qq!${index_string} \\
  2>${stderr} 1>${stdout}
!;
    }

    my $comment = qq!## Generating salmon indexes for species: ${species}
## in $options->{libdir}/${libtype}/indexes!;
    my $jobid = $class->Submit(
        comment => $comment,
        jdepends => $options->{jdepends},
        jstring => $jstring,
        jname => qq"salidx_${species}",
        jmem => 24,
        jprefix => '15',
        modules => $options->{modules},
        output => $output_dir,
        stderr => $stderr,
        stdout => $stdout,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($jobid);
}

=head2 C<STAR_Index>

 Create indexes appropriate for STAR.

=cut
sub STAR_Index {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input',],
        modules => ['star'],);
    my $comment = '## STAR Index creation.';
    my $libtype = 'genome';
    $libtype = $options->{libtype} if ($options->{libtype});
    my $species = basename($options->{input}, ('.fasta', '.fa'));
    my $copied_location = qq"$options->{libpath}/$options->{libtype}/${species}.fasta";
    if (!-f $copied_location) {
        cp($options->{input}, $copied_location);
    }
    my $star_refdir = "$options->{libdir}/${libtype}/indexes/${species}_star_index";
    my $output_dir = qq"$options->{basedir}/outputs/$options->{jprefix}star_index";
    my $stdout = qq"${output_dir}/index.stdout";
    my $stderr = qq"${output_dir}/index.stderr";
    my $jstring = qq!mkdir -p ${output_dir}
STAR \\
  --runMode genomeGenerate \\
  --runThreadN 12 \\
  --genomeDir ${star_refdir} \\
  --genomeFastaFiles $options->{libdir}/${libtype}/${species}.fasta \\
  --sjdbGTFfile $options->{libdir}/${libtype}/${species}.gtf \\
  --limitGenomeGenerateRAM 160000000000 \\
  2>${stderr} 1>${stdout}

!;
    my $jobid = $class->Submit(
        comment => $comment,
        jdepends => $options->{depends},
        jstring => $jstring,
        jname => 'staridx',
        jprefix => $options->{jprefix},
        jmem => 180,
        modules => $options->{modules},
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        stderr => $stderr,
        stdout => $stdout,
        jqueue => 'xlarge',);
    return($jobid);
}

1;

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

=cut
