package Bio::Adventure::Index;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use File::Basename;
use File::Copy qw"cp";
use File::Spec;
use File::Which qw"which";

=head1 NAME

Bio::Adventure::Index - Invoke the various indexing utilities.

=head1 SYNOPSIS

=cut

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
    my $copied_location = qq"$options->{libdir}/$options->{libtype}/${species}.fasta";
    if (!-f $copied_location) {
        my $copied = qx"less $options->{input} > ${copied_location}";
    }

    my $jstring = qq!bowtie-build $options->{input} \\
  $options->{libdir}/$options->{libtype}/indexes/${species}
!;
    my $comment = qq!## Generating bowtie1 indexes for species: ${species} in $options->{libdir}/$options->{libtype}/indexes!;
    my $bt1_index = $class->Submit(
        comment => $comment,
        jname => qq"bt1idx_${species}",
        jdepends => $options->{jdepends},
        jstring => $jstring,
        jprefix => '10',
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
        modules => ['bowtie2'],);
    my $libtype = $options->{libtype};
    my $libdir = File::Spec->rel2abs($options->{libdir});
    my $species = basename($options->{input}, ('.gz', '.bz2', '.xz'));
    $species = basename($species, ('.fasta', '.fa'));
    my $copied_location = qq"$options->{libdir}/$options->{libtype}/${species}.fasta";
    if (!-f $copied_location) {
        my $copied = qx"less $options->{input} > ${copied_location}";
    }
    my $jstring = qq!bowtie2-build $options->{input} \\
  $options->{libdir}/${libtype}/indexes/${species}
!;
    my $comment = qq!## Generating bowtie2 indexes for species: ${species} in $options->{libdir}/${libtype}/indexes!;
    my $indexer = $class->Submit(
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => "bt2idx_${species}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        modules => $options->{modules},
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
        modules => ['bwa'],);
    my $species = basename($options->{input}, ('.fasta', '.fa'));
    my $copied_location = qq"$options->{libdir}/$options->{libtype}/${species}.fa";
    if (!-f $copied_location) {
        cp($options->{input}, $copied_location);
    }
    my $jstring = qq!
start=\$(pwd)
cd $options->{libdir}/$options->{libtype}/indexes
ln -sf ../${species}.fa .
bwa index ${species}.fa \\
  2>bwa_index.err \\
  1>bwa_index.out
cd \$start
!;
    my $comment = qq!## Generating bwa indexes for species: ${species} in $options->{libdir}/$options->{libtype}/indexes!;
    my $bwa_index = $class->Submit(
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => "bwaidx",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        modules => $options->{modules},
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
        required => ['input'],);
    my $libname = basename($options->{input}, '.fasta');
    ## First check for the relevant library in $ENV{BLASTDB}
    ## If it isn't there, make one in basedir/blastdb/
    my $foundlib = 0;
    if ($options->{type} ne 'prot' && $options->{type} ne 'nucl') {
        die("makeblastdb requires either a type of 'prot' or 'nucl'.");
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
    my $relative_directory = qq"blastdb";
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
  -out ${db_directory}/${libname}";
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
    die("Could not find kraken2 in your PATH.") unless($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        library => 'viral',
        modules => ['kraken'],);
    ## kraken2 --db ${DBNAME} --paired --classified-out cseqs#.fq seqs_1.fq seqs_2.fq
    my $job_name = $class->Get_Job_Name();
    my $output_dir = qq"outputs/extend_kraken";

    my $comment = qq!## This is a script to extend an existing kraken2 library with some new sequences.
!;
    my $jstring = qq!mkdir -p ${output_dir}
kraken2-build --download-taxonomy --db \${KRAKEN_DB_PATH}/$options->{library} \\
              2>${output_dir}/kraken2-build.out 1>&2
kraken2-build --add-to-library $options->{input} --db \${KRAKEN_DB_PATH}/$options->{library} \\
              2>>${output_dir}/kraken2-build.out 1>&2
kraken2-build --download-library $options->{library} \\
              --db \${KRAKEN_DB_PATH}/$options->{library} \\
              2>>${output_dir}/kraken2-build.out 1>&2
kraken2-build --build --db \${KRAKEN_DB_PATH}/$options->{library} \\
              2>>${output_dir}/kraken2-build.out 1>&2
!;
    my $kraken = $class->Submit(
        cpus => 6,
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => "kraken_${job_name}",
        jprefix => "99",
        jstring => $jstring,
        jmem => 96,
        modules => $options->{modules},
        output => qq"${output_dir}/kraken2-build.out",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => "large",);
    return($kraken);
}

=head2 C<HT2_Index>

Create a hisat2 index using ${species}.fasta and leaves it in the indexes/
directory.

=cut
sub HT2_Index {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        modules => ['hisat2'],
        jprefix => '21',);
    my $libtype = $options->{libtype};
    my $libdir = File::Spec->rel2abs($options->{libdir});
    my $species = basename($options->{input}, ('.fasta', '.fa'));
    my $copied_location = qq"$options->{libdir}/$options->{libtype}/${species}.fasta";
    if (!-f $copied_location) {
        cp($options->{input}, $copied_location);
    }
    my $jstring = qq!
hisat2-build $options->{input} \\
  $options->{libdir}/${libtype}/indexes/${species}
!;
    my $comment = qq!## Generating hisat2 indexes for species: ${species} in $options->{libdir}/${libtype}/indexes!;
    my $indexer = $class->Submit(
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => qq"ht2idx_${species}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
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
    my $copied_location = qq"$options->{libdir}/$options->{libtype}/${cds}.fasta";
    if (!-f $copied_location) {
        cp($options->{input}, $copied_location);
    }
    my $libtype = $options->{libtype};
    my $input = File::Spec->rel2abs($options->{input});
    my $jstring = qq!
kallisto index -i $options->{libdir}/${libtype}/indexes/${species}.idx \\
  ${input}
!;
    my $comment = qq!## Generating kallisto indexes for species: ${species} in $options->{libdir}/${libtype}/indexes!;
    my $ka_index = $class->Submit(
        comment => $comment,
        jdepends => $options->{jdepends},
        jstring => $jstring,
        jname => "kalidx",
        jprefix => $options->{jprefix},
        modules => $options->{modules},
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    return($ka_index);
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
    my $copied_location = qq"$options->{libdir}/$options->{libtype}/${species}.fasta";
    if (!-f $copied_location) {
        cp($options->{input}, $copied_location);
    }

    my $comment = qq"## RSEM Index creation.";
    my $jstring = qq!
rsem-prepare-reference --bowtie2 $options->{input} ${species}
!;
    my $jobid = $class->Submit(
        comment => $comment,
        jdepends => $options->{jdepends},
        jstring => $jstring,
        jname => "rsemidx",
        jprefix => $options->{jprefix},
        modules => $options->{modules},
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
    my $copied_location = qq"$options->{libdir}/$options->{libtype}/${cds}.fasta";
    my $species_location = qq"$options->{libdir}/$options->{libtype}/${species}.fasta";

    if (!-f $copied_location) {
        cp($options->{input}, $copied_location);
    }
    my $decoy_copy_string = qq'';
    my $jstring = qq'';
    my $index_input = $options->{input};
    my $index_string = qq!
salmon index -t ${index_input} -i $options->{libdir}/${libtype}/indexes/${species}_salmon_index!;
    if (-f $species_location or -f $species_file) {
        if (!-f $species_location) {
            cp($species_file, $species_location);
        }
        my $decoy_location = qq"$options->{libdir}/${libtype}/${species}_decoys.fasta";
        $decoy_copy_string = qq!less $options->{input} > ${decoy_location} && less ${species_file} >> ${decoy_location}
less ${species_file} | grep '^>' | sed 's/^>//g' > ${decoy_location}.txt
!;
        $index_input = $decoy_location;
        $index_string = qq!
salmon index -t ${index_input} -i $options->{libdir}/${libtype}/indexes/${species}_salmon_index!;
        $jstring = qq!${decoy_copy_string}
${index_string} --decoys ${decoy_location}.txt
!;
    } else {
        warn("This function would prefer to make a decoy aware index set which requires the full genome.");
        say("Waiting 10 seconds to see if you want to quit and gather that genome,
otherwise a decoy-less index will be generated.");
        sleep(10);
    }

    my $comment = qq!## Generating salmon indexes for species: ${species} in $options->{libdir}/${libtype}/indexes!;
    my $jobid = $class->Submit(
        comment => $comment,
        jdepends => $options->{jdepends},
        jstring => $jstring,
        jname => "salidx_${species}",
        jmem => 24,
        jprefix => "15",
        modules => $options->{modules},
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
    my $comment = qq"## STAR Index creation.";
    my $libtype = 'genome';
    $libtype = $options->{libtype} if ($options->{libtype});
    my $species = basename($options->{input}, ('.fasta', '.fa'));
    my $copied_location = qq"$options->{libdir}/$options->{libtype}/${species}.fasta";
    if (!-f $copied_location) {
        cp($options->{input}, $copied_location);
    }
    my $star_refdir = "$options->{libdir}/${libtype}/indexes/${species}_star_index";
    my $jstring = qq!
STAR \\
  --runMode genomeGenerate \\
  --runThreadN 12 \\
  --genomeDir ${star_refdir} \\
  --genomeFastaFiles $options->{libdir}/${libtype}/${species}.fasta \\
  --sjdbGTFfile $options->{libdir}/${libtype}/${species}.gtf \\
  --limitGenomeGenerateRAM 160000000000
!;
    my $jobid = $class->Submit(
        comment => $comment,
        jdepends => $options->{depends},
        jstring => $jstring,
        jname => "staridx",
        jprefix => $options->{jprefix},
        jmem => 180,
        modules => $options->{modules},
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => 'xlarge',);
    return($jobid);
}

1;

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

=cut
