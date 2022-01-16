package Bio::Adventure::Annotation;
## LICENSE: gplv2
## ABSTRACT:  Kitty!
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";
use Moo;
extends 'Bio::Adventure';

use Bio::SeqIO;
use Bio::Seq;
use Bio::SeqFeature::Generic;
use Bio::Tools::GuessSeqFormat;
use Cwd qw"abs_path getcwd cwd";
use Data::Table;
use Data::Table::Excel qw"tables2xlsx";
use File::Basename;
use File::Copy qw"cp";
use File::Spec;
use File::Path qw"make_path";
use File::Which qw"which";
use File::ShareDir qw":ALL";
use List::MoreUtils qw"any";
use Template;
use Text::CSV_XS::TSV;

=head2 C<Interproscan>

 Invoke interproscan on a set of sequences.

 Interproscan is (I think) the gold-standard of similarity-based search
 tools.  It provides a single interface for searching against a
 comprehensive array of databases.

=item C<Arguments>

 input(required): Fasta file containing amino acid sequences.
 jprefix(21): Prefix of the job name/output directory.
 jmem(8): Expected memory consumption.
 cpus(4): Limit the number of cpus per job with this.
 modules('interproscan'): Environment module to load.

=item C<Invocation>

> cyoa --task annot --method interpro --input amino_acid_data.fasta

=cut
sub Interproscan {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jprefix => '21',
        jmem => 8,
        cpus => 4,
        modules => ['interproscan'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('interproscan.sh');
    die('Could not find interproscan in your PATH.') unless($check);

    my $job_name = $class->Get_Job_Name();
    my $inputs = $class->Get_Paths($options->{input});
    my $cwd_name = basename(cwd());

    my $interproscan_exe_dir = dirname($check);
    ## Hey, don't forget abs_path requires a file which already exists.
    my $input_filename = basename($options->{input});
    my $output_filename = qq"${input_filename}.tsv";
    my $input_dir = dirname($options->{input});
    my $input_dirname = basename($input_dir);
    my $input_path = abs_path($input_dir);
    $input_path = qq"${input_path}/${input_filename}";
    my $output_dir = qq"outputs/$options->{jprefix}interproscan_${input_dirname}";
    my $comment = qq!## This is a interproscan submission script
!;
    my $stdout = qq"${output_dir}/interproscan.stdout";
    my $stderr = qq"${output_dir}/interproscan.stderr";
    my $jstring = qq!mkdir -p ${output_dir}
start=\$(pwd)
cd ${output_dir}
perl -pe 's/\\*//g' ${input_path} > ${input_filename}
interproscan.sh -i ${input_filename} \\
  2>interproscan.stderr \\
  1>interproscan.stdout
ln -sf ${output_filename} interproscan.tsv
cd \${start}
!;
    my $interproscan = $class->Submit(
        cpus => 8,
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => "interproscan_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => 16,
        modules => $options->{modules},
        output => qq"${output_dir}/interproscan.tsv",
        output_gff => qq"${output_dir}/${input_filename}.gff3",
        output_tsv => qq"${output_dir}/interproscan.tsv",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => 'large',
        stdout => $stdout,
        stderr => $stderr,
        walltime => '144:00:00',);

    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($interproscan);
}

=head2 C<Kraken>

 Use kraken2 to taxonomically classify reads.

 Kraken2 is a kmer-based read classifier.  It is quite fast once its
 databases are built.

=item C<Arguments>

 input(required): Set of fastq reads.
 library('viral'): Kraken2 database to classify against.
 jprefix(11): Prefix for the job name and output directory.
 modules('kraken'): Environment module to load.

=item C<Invocation>

> cyoa --task annot --method kraken --input read1.fastq.xz:read2.fastq.xz --library standard

=cut
sub Kraken {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        library => 'viral',
        jprefix => '11',
        modules => ['kraken'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('kraken2');
    die('Could not find kraken2 in your PATH.') unless($check);
    ## kraken2 --db ${DBNAME} --paired --classified-out cseqs#.fq seqs_1.fq seqs_2.fq
    my $job_name = $class->Get_Job_Name();
    my $input_directory = basename(cwd());
    my $output_dir = qq"outputs/$options->{jprefix}kraken_$options->{library}";
    make_path($output_dir);
    my $input_string = "";
    if ($options->{input} =~ /\:|\;|\,|\s+/) {
        my @in = split(/\:|\;|\,|\s+/, $options->{input});
        $input_string = qq" --paired <(less $in[0]) <(less $in[1]) ";
    } else {
        $input_string = qq"<(less $options->{input}) ";
    }
    my $comment = qq!## This is a kraken2 submission script
!;
    my $stdout = qq"${output_dir}/kraken.stdout";
    my $stderr = qq"${output_dir}/kraken.stderr";
    my $jstring = qq!kraken2 --db $ENV{KRAKEN2_DB_PATH}/$options->{library} \\
  --report ${output_dir}/kraken_report.txt --use-mpa-style \\
  --use-names ${input_string} \\
  --classified-out ${output_dir}/classified#.fastq.gz \\
  --unclassified-out ${output_dir}/unclassified#.fastq.gz \\
  2>${stderr} \\
  1>${stdout}
!;
    my $kraken = $class->Submit(
        cpus => 6,
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => "kraken_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => 96,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => 'large',
        modules => $options->{modules},
        output => qq"${output_dir}/kraken_report.txt",
        stdout => $stdout,
        stderr => $stderr);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($kraken);
}

=head2 C<Prokka>

 Perform automated assembly annotation, intended for bacteria.

 Prokka is an automagic annotation tool for bacterial assemblies.  It
 seems useful for other relatively small genomes.  It is one of the
 suite of tools written by Torsten Seeman, and probably my
 favorite. It takes an assembly, runs prodigal, some tRNA searches,
 some BLAST searches, and generates fasta/gbk/etc files from the results.

=item C<Arguments>

 input(required): Input assembly.
 arbitrary(''): Add arbitrary arguments here.
 coverage(30): Intended for a coverage filter, not currently used.
 evalue('1e-05'): Intended as a confidence filter, not currently used.
 gcode('11'): Use this codon code (bacterial is 11).
 genus('phage'): Set the genus tag in the output files.
 kingdom('bacteria'): Set the kingdom tag in the output files.
 locus_tag('unknownphage'): Set the locus tag in the output files.
 species('virus'): Set the species tag in the output files.
 jprefix('19'): Use this prefix for the job/output.
 modules('prokka'): Use this environment module.

=item C<Invocation>

> cyoa --task annot --method prokka --input assembled_genomic_sequence.fasta

=cut
sub Prokka {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        arbitrary => '',
        coverage => 30,
        evalue => '1e-05',
        gcode => '11',
        genus => 'phage',
        kingdom => 'bacteria',
        locus_tag => 'unknownphage',
        species => 'virus',
        jmem => 12,
        jprefix => '19',
        modules => ['prokka'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('prokka');
    die('Could not find prokka in your PATH.') unless($check);

    my $job_name = $class->Get_Job_Name();
    my $input_paths = $class->Get_Paths($options->{input});
    my $kingdom_string = '';
    if ($options->{kingdom} ne '') {
        $kingdom_string = qq" --kingdom $options->{kingdom} --gcode $options->{gcode} ";
    }

    my $training_file = qq"$options->{libdir}/hmm/$options->{species}_gc$options->{gcode}.training";
    my $cwd_name = basename(cwd());
    ## If this is being run a part of a pipeline, the basename(dirname(input)) will be something
    ## useful/interesting.  If not, it will not likely be of use but instead be '.'
    ## Let us check for that and set $input_name to something useful in that case.
    my $input_name;
    if (defined($input_paths->{dirname})) {
        $input_name = $input_paths->{dirname};
    } else {
        $input_name = $input_paths->{filebase_extension};
    }
    my $locus_tag;
    if ($options->{locus_tag}) {
        $locus_tag = $options->{locus_tag};
    } else {
        $locus_tag = basename($input_name, ('.fasta'));
    }
    my $output_dir = qq"outputs/$options->{jprefix}prokka";
    my $comment = qq!## This is a script to run prokka.
!;

##    prokka --addgenes -rfam --force --locustag EAb03 --genus Phage --species virus --cdsrnaolap --usegenus --kingdom Bacteria --gcode 11 --prodigaltf /home/trey/libraries/hmm/abaumannii_phage_samples_gc11.training --outdir outputs/19prokka_15termreorder_13unicycler --prefix EAb03 outputs/15termreorder_13unicycler/final_assembly_reordered.fasta --evalue '1e-05' --coverage 30 2>&1 | less

    my $jstring = qq!mkdir -p ${output_dir}
prokka --addgenes --rfam --force ${kingdom_string} \\
  --locustag ${locus_tag} --genus $options->{genus} \\
  --compliant --cdsrnaolap --usegenus \\
  --prodigaltf ${training_file} \\
  --outdir ${output_dir} \\
  --prefix ${cwd_name} \\
  $options->{input} \\
  2>${output_dir}/prokka.stderr \\
  1>${output_dir}/prokka.stdout
!;

    my $error_file =  qq"${output_dir}/${cwd_name}.err";
    my $peptide_file = qq"${output_dir}/${cwd_name}.faa";
    my $cds_file = qq"${output_dir}/${cwd_name}.ffn";
    my $assembly_copy = qq"${output_dir}/${cwd_name}.fna";
    my $assembly_renamed_contigs = qq"${output_dir}/${cwd_name}.fsa";
    my $genbank_file = qq"${output_dir}/${cwd_name}.gbk";
    my $gff_file = qq"${output_dir}/${cwd_name}.gff";
    my $sqn_file = qq"${output_dir}/${cwd_name}.sqn";
    my $tbl_file = qq"${output_dir}/${cwd_name}.tbl";
    my $tsv_file = qq"${output_dir}/${cwd_name}.tsv";

    my $prokka = $class->Submit(
        cpus => 6,
        comment => $comment,
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
        jname => "prokka_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => 24,
        modules => $options->{modules},
        output => $cds_file,
        output_error => $error_file,
        output_peptide => $peptide_file,
        output_cds => $cds_file,
        output_assembly => $assembly_copy,
        output_assemblyrenamed => $assembly_renamed_contigs,
        output_fsa => $assembly_renamed_contigs,
        output_genbank => $genbank_file,
        output_gff => $gff_file,
        output_sqn => $sqn_file,
        output_tbl => $tbl_file,
        output_tsv => $tsv_file,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => 'workstation',);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($prokka);
}

sub Transposonpsi {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        input_faa => '',
        jprefix => '21',
        jmem => 8,
        cpus => 4,
        modules => ['transposonpsi'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('transposonPSI.pl');
    die('Could not find transposonPSI in your PATH.') unless($check);

    my $job_name = $class->Get_Job_Name();
    my $inputs = $class->Get_Paths($options->{input});
    my $cwd_name = basename(cwd());

    my $interproscan_exe_dir = dirname($check);
    ## Hey, don't forget abs_path requires a file which already exists.
    my $input_filename = basename($options->{input});
    my $output_filename = qq"${input_filename}.tsv";
    my $input_dir = dirname($options->{input});
    my $input_dirname = basename($input_dir);
    my $input_path = abs_path($input_dir);
    $input_path = qq"${input_path}/${input_filename}";
    my $output_dir = qq"outputs/$options->{jprefix}interproscan_${input_dirname}";
    my $comment = qq!## This is a interproscan submission script
!;
    my $stdout = qq"${output_dir}/interproscan.stdout";
    my $stderr = qq"${output_dir}/interproscan.stderr";
    my $jstring = qq!mkdir -p ${output_dir}
start=\$(pwd)
cd ${output_dir}
perl -pe 's/\\*//g' ${input_path} > ${input_filename}
interproscan.sh -i ${input_filename} \\
  2>interproscan.stderr \\
  1>interproscan.stdout
ln -sf ${output_filename} interproscan.tsv
cd \${start}
!;
    my $interproscan = $class->Submit(
        cpus => 8,
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => "interproscan_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => 16,
        modules => $options->{modules},
        output => qq"${output_dir}/interproscan.tsv",
        output_gff => qq"${output_dir}/${input_filename}.gff3",
        output_tsv => qq"${output_dir}/interproscan.tsv",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => 'large',
        stdout => $stdout,
        stderr => $stderr,
        walltime => '144:00:00',);

    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($interproscan);

}

=head2 C<Extract_Annotations>

 Pull apart the encoded trinotate annotations into a more readable format.

 Trinotate uses its own format to encode annotations, it wraps
 multiple hits into one cell of a tsv output file with a combination
 of backticks(`) and caret(^).  This function reads that and splits
 them up into separate cells.  It is called by Extract_Trinotate().

=item C<Arguments>

 identity: Cutoff for percent identity between a query and hit.
 evalue: Cutoff by evalue between query and a hit.
 ids: Set of IDs to explicitly extract from the trinotate output.
 fh: filehandle to read.

=cut
sub Trinotate_Extract_Annotations {
    my ($class, $datum, %args) = @_;
    my $min_identity = $args{identity};
    my $eval_max = $args{evalue};
    my $ids = $args{ids};
    my $fh = $args{fh};
    my $options = $class->Get_Vars();

    my $id = $datum->{prot_id};
    $id = $datum->{transcript_id} if ($id eq '.');
    if (defined($options->{species})) {
        my $species = $options->{species};
        $species =~ s/\s+/_/g;
        $id =~ s/TRINITY/${species}/g;
    }
    $id =~ s/\:\:/; /g;

    if (defined($ids->{$id})) {
        ## Then this record has been processed, just increment it and move on.
        $ids->{$id}++;
    } elsif (!defined($datum->{swissprot} && !defined($datum->{rnammer}))) {
        $ids->{$id} = 1;
    } elsif (defined($datum->{rnammer})) {
        my $seq = $datum->{sequence};
        $seq = join("\n", ($seq =~ m/.{1,80}/g));
        $ids->{$id} = 1;
        my $header_string = qq"${id}; ${id}; undef; undef; $datum->{rnammer}->[0]->{name}, $datum->{rnammer}->[0]->{region}; undef; undef";
        print $fh qq">${header_string}
${seq}
";
    } elsif ($datum->{swissprot}->[0]->{identity} >= $min_identity &&
                 $datum->{swissprot}->[0]->{e_value} <= $eval_max) {
        $ids->{$id} = 1;
        my $seq = $datum->{sequence};
        $seq = join("\n", ($seq =~ m/.{1,80}/g));
        my $header_string = qq"${id}; $datum->{swissprot}->[0]->{fullname}; $datum->{swissprot}->[0]->{species}; $datum->{swissprot}->[0]->{e_value}; $datum->{swissprot}->[0]->{identity}";
        print $fh qq">${header_string}
${seq}
";
    }
    return($ids);
}

=head2 C<Extract_Trinotate>

 Parse the trinotate encoded blast results into a simpler table.

 The trinotate output format is a bit... unwieldy.  This seeks to
 parse out the useful information from it.

=item C<Arguments>

 input(required): Input csv file from trinotate.
 output(interesting.fasta): Output for the parsed csv.
 evalue(1e-10): Evalue cutoff for trinotate output.
 identity(70): Minimum percent identity cutoff for trinotate output.

=item C<Invocation>

> cyoa --task assembly --method extract --input trinotate_output.csv

=cut
sub Extract_Trinotate {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jname => 'trin_rsem',
        output => 'interesting.fasta',
        evalue => 1e-10,
        identity => 70,);
    my $job_name = $class->Get_Job_Name();
    my $trinity_out_dir = qq"outputs/trinity_${job_name}";
    my $input = FileHandle->new("<$options->{input}");
    my $parser = Parse::CSV->new(
        handle => $input,
        sep_char => "\t",
        names => 1,);

    my $count = 0;
    my $all_annotations = [];
    my $out = FileHandle->new(">$options->{output}");
    while (my $object = $parser->fetch) {
        my $filled_in = Read_Write_Annotation($class, $object,
                                              db => $all_annotations,
                                              count => $count,
                                              evalue => $options->{evalue},
                                              identity => $options->{identity},
                                              fh => $out);
        $count = $count + 1;
    }
    $out->close();
    $input->close();
    return($count);
}

=head2 C<Read_Write_Annotation>

 Read trinotate output and write out extracted annotations.

 Called by Extract_Trinotate() to help write out the trinotate csv information
 into an easier-to-read format.

=cut
sub Read_Write_Annotation {
    my ($class, $object, %args) = @_;
    my $annotations = $args{db};
    my $element_number = $args{count};
    my $fh = $args{fh};
    ## $object is a hash reference containing the various csv fields.
    ## Extract them by name and dump them to the array of material
    my $element = {};
    my $filled_in = 0;
    my $ids = {};

    foreach my $key (keys %{$object}) {
        my @result_list = ();
        my @field_elements = ();
        my $field_text = $object->{$key};
        if ($field_text =~ /\`/) {
            @field_elements = split(/\`/, $field_text);
        } else {
            $field_elements[0] = $field_text;
        }

        for my $t (0..$#field_elements) {
            my $ret = {};
            if ($key eq 'sprot_Top_BLASTX_hit') {
                if ($field_text eq '.') {
                    ## Null case: if a . then just drop out.
                    $element->{swissprot} = undef;
                } else {
                    ## The swissprot elements look something like:
                    ## Name        ## ARIA_ARATH^
                    ## Name again? ## ARIA_ARATH^
                    ## Query range ## Q:7-573,H:513-701^
                    ## %identity   ## 75.66%ID^
                    ## E-value     ## E:3e-102^
                    ## Full name   ## RecName: Full=ARM REPEAT PROTEIN INTERACTING WITH ABF2;^
                    ## Ontology    ## Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliophyta; eudicotyledons; Gunneridae; Pentapetalae; rosids; malvids; Brassicales; Brassicaceae; Camelineae; Arabidopsis',
                    ## Start by pulling apart the data using the chosen (somewhat odd) separator '^'.
                    my ($name, $name_again, $query_coords, $identity, $e_value, $fullname, $ontology) = split(/\^/, $field_text);
                    ## Clean up the fields a little
                    $identity =~ s/\%ID//g;
                    $e_value =~ s/E://g;
                    $fullname =~ s/RecName: //g;
                    $fullname =~ s/Full=//g;
                    $fullname =~ s/\;//g;
                    my ($query, $h) = split(/\,/, $query_coords);
                    $query =~ s/Q\://g;
                    $h =~ s/H\://g;
                    $name =~ s/Full=//g;
                    my @ontology_list = split(/; /, $ontology);
                    my $species = qq"$ontology_list[$#ontology_list - 1] $ontology_list[$#ontology_list]";
                    ## Put them back into the $ret
                    $ret->{name} = $name;
                    $ret->{coords} = $query_coords;
                    $ret->{identity} = $identity;
                    $ret->{e_value} = $e_value;
                    $ret->{fullname} = $fullname;
                    $ret->{species} = $species;
                    $filled_in = $filled_in + 6;
                    ## And refill field_elements with this information.
                    $field_elements[$t] = $ret;
                    ## Finally, fill in the swissprot entry with this information.
                    $element->{swissprot} = \@field_elements;
                }
            } elsif ($key eq 'gene_ontology_pfam') {
                if ($field_text eq '.') {
                    $element->{pfam_go} = undef;
                } else {
                    my ($id, $ontology, $name) = split(/\^/, $field_text);
                    $ret->{id} = $id;
                    $ret->{ontology} = $ontology;
                    $ret->{name} = $name;
                    $filled_in = $filled_in + 3;
                    $field_elements[$t] = $ret;
                    $element->{pfam_go} = \@field_elements;
                }
            } elsif ($key eq 'gene_ontology_blast') {
                if ($field_text eq '.') {
                    $element->{blast_go} = undef;
                } else {
                    my ($id, $ontology, $name) = split(/\^/, $field_text);
                    $ret->{id} = $id;
                    $ret->{ontology} = $ontology;
                    $ret->{name} = $name;
                    $filled_in = $filled_in + 3;
                    $field_elements[$t] = $ret;
                    $element->{blast_go} = \@field_elements;
                }
            } elsif ($key eq 'RNAMMER') {
                if ($field_text eq '.') {
                    $element->{rnammer} = undef;
                } else {
                    my ($name, $region) = split(/\^/, $field_text);
                    $ret->{name} = $name;
                    $ret->{region} = $region;
                    $filled_in = $filled_in + 2;
                    $field_elements[$t] = $ret;
                    $element->{rnammer} = \@field_elements;
                }
            } elsif ($key eq 'eggnog') {
                if ($field_text eq '.') {
                    $element->{eggnog} = undef;
                } else {
                    my ($id, $name) = split(/\^/, $field_text);
                    $ret->{id} = $id;
                    $ret->{name} = $name;
                    $filled_in = $filled_in + 2;
                    $field_elements[$t] = $ret;
                    $element->{eggnog} = \@field_elements;
                }
            } elsif ($key eq 'transcript_id') {
                $filled_in = $filled_in + 1;
                $element->{transcript_id} = $field_text;
            } elsif ($key eq 'Kegg') {
                if ($field_text eq '.') {
                    $element->{kegg} = undef;
                } else {
                    my ($id, $name) = split(/\^/, $field_text);
                    $ret->{id} = $id;
                    $ret->{name} = $name;
                    $filled_in = $filled_in + 2;
                    $field_elements[$t] = $ret;
                    $element->{kegg} = \@field_elements;
                }
            } elsif ($key eq 'prot_coords') {
                if ($field_text eq '.') {
                    $element->{prot_coords} = undef;
                } else {
                    $filled_in = $filled_in + 1;
                    $element->{prot_coords} = $field_text;
                }
            } elsif ($key eq 'prot_id') {
                $element->{prot_id} = $field_text;
            } elsif ($key eq 'transcript') {
                $element->{sequence} = $field_text;
            } elsif ($key eq 'TmHMM') {
                ## An example field
                ## ExpAA=124.31^PredHel=6^Topology=i59-78o83-105i126-148o152-174i187-209o229-251i
                if ($field_text eq '.') {
                    $element->{tmhmm} = undef;
                } else {
                    my ($expaa, $pred_helixes, $topology) = split(/\^/, $field_text);
                    $expaa =~ s/ExpAA=//g;
                    $pred_helixes =~ s/PredHel=//g;
                    $topology =~ s/Topology=//g;
                    $ret->{expaa} = $expaa;
                    $ret->{pred_helixes} = $pred_helixes;
                    $ret->{topology} = $topology;
                    $filled_in = $filled_in + 3;
                    $field_elements[$t] = $ret;
                    $element->{tmhmm} = \@field_elements;
                }
            } elsif ($key eq 'SignalP') {
                ## sigP:1^18^0.613^YES
                if ($field_text eq '.') {
                    $element->{signalp} = undef;
                } else {
                    my ($first, $second, $third, $boolean) = split(/\^/, $field_text);
                    $first =~ s/sigP\://g;
                    $ret->{first} = $first;
                    $ret->{second} = $second;
                    $ret->{third} = $third;
                    $ret->{boolean} = $boolean;
                    $filled_in = $filled_in + 4;
                    $field_elements[$t] = $ret;
                    $element->{tmhmm} = \@field_elements;
                }
            }
            ## Set the nth annotation to its annotation.
            $annotations->[$element_number] = $element;
        } ## End iterating over every element in a multi-element field
    } ## End the foreach every key in the object
    $ids = Trinotate_Extract_Annotations($class, $element,
                                         ids => $ids,
                                         fh => $fh,
                                         evalue => $args{evalue},
                                         identity => $args{identity},);
    return($filled_in);
}

=head2 C<Transdecoder>

 Run transdecoder on a transcriptome.

 Submit a trinity denovo sequence assembly to transdecoder.

=item C<Arguments>

 input(required): Output from trinity for post processing.

=item C<Invocation>

> cyoa --task assembly --method transdecoder --input trinity.fasta

=cut
sub Transdecoder {
    my ($class, %args) = @_;
    my $check = which('TransDecoder.LongOrfs');
    die('Could not find transdecoder in your PATH.') unless($check);
    my $transdecoder_exe_dir = dirname($check);
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        modules => ['transdecoder'],);
    my $transdecoder_input = File::Spec->rel2abs($options->{input});
    my $job_name = $class->Get_Job_Name();
    my $output_dir = qq"outputs/trinity_${job_name}";
    my $comment = qq!## This is a transdecoder submission script
!;
    my $jstring = qq!mkdir -p ${output_dir} && cd ${output_dir} && \\
  TransDecoder.LongOrfs -t ${transdecoder_input} \\
    2>${output_dir}/transdecoder_longorfs_${job_name}.err \\
    1>${output_dir}/transdecoder_longorfs_${job_name}.out
TransDecoder.Predict -t ${transdecoder_input} \\
  2>${output_dir}/transdecoder_predict_${job_name}.err \\
  1>${output_dir}/transdecoder_predict_${job_name}.out
${transdecoder_exe_dir}/util/cdna_alignment_orf_to_genome_orf.pl \\
  ${output_dir}/transcripts.fasta.transdecoder.gff3 \\
  ${output_dir}/transcripts.gff3 \\
  ${transdecoder_input} \\
  2>${output_dir}/cdna_alignment_${job_name}.err \\
  1>${output_dir}/transcripts.fasta.transdecoder.genome.gff3
!;
    my $transdecoder = $class->Submit(
        comment => $comment,
        jmem => 4,
        jname => "transdecoder_${job_name}",
        jprefix => '47',
        jstring => $jstring,
        modules => $options->{modules},
        output => qq"${output_dir}/transcripts.fasta.transdecoder.genome.gff3",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => 'workstation',);
    return($transdecoder);
}

=head2 C<Trinotate>

 Submit a trinity denovo sequence assembly to trinotate.

 In the time since writing this, I added a cheesy hack to allow it to
 run on genomic assemblies as well.

=item C<Arguments>

 input(required): Input fasta from trinity.

=item C<Invocation>

> cyoa --task assembly --method trinotate --input trinity.fasta

=cut
sub Trinotate {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jprefix => '20',
        modules => ['divsufsort', 'transdecoder', 'blast', 'blastdb', 'signalp', 'hmmer',
                    'tmhmm', 'rnammer', 'trinotate', ],
        trinotate => 'autoTrinotate.pl',
        config => 'conf.txt',
        );

    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('Trinotate');
    die("Could not find Trinotate in your PATH:
 $ENV{PATH}.") unless($check);

    my $job_name = $class->Get_Job_Name();
    my $cwd_name = basename(cwd());
    my $trinotate_exe_dir = dirname($check);
    ## Once again, abs_path only works on stuff which already exists.
    ## So create the output directory, and use that.
    my $input_paths = $class->Get_Paths($options->{input});
    my $input_full = $input_paths->{fullpath};
    my $output_name = basename($input_full, ('.fasta', '.fa', '.fna', '.fsa', '.ffn'));
    $output_name = qq"${output_name}.tsv";
    my $output_dir = qq"outputs/$options->{jprefix}trinotate";
    my $stdout = qq"${output_dir}/trinotate_${job_name}.stdout";
    my $stderr = qq"${output_dir}/trinotate_${job_name}.stderr";
    $output_dir .= qq"$input_paths->{dirname}" if (defined($input_paths->{dirname}));
    my $comment = qq!## This is a trinotate submission script
!;
    my $jstring = qq!mkdir -p ${output_dir}
start=\$(pwd)
cd ${output_dir}
ln -sf $input_paths->{fullpath} .
if [ -f $input_paths->{filename}.gene_trans_map ]; then
  ln -sf $input_paths->{filename}.gene_trans_map .
else
  ids=\$({ grep "^>" $input_paths->{fullpath} || test \$? = 1; } | sed 's/>//g' | awk '{print \$1}')
  rm -f $input_paths->{filename}.gene_trans_map
  for i in \${ids}; do
    echo "\${i}	\${i}" >> $input_paths->{filename}.gene_trans_map
  done
fi

${trinotate_exe_dir}/auto/$options->{trinotate} \\
  --conf ${trinotate_exe_dir}/auto/$options->{config} \\
  --Trinotate_sqlite ${trinotate_exe_dir}/sample_data/Trinotate.boilerplate.sqlite \\
  --transcripts $input_paths->{filename} \\
  --gene_to_trans_map $input_paths->{filename}.gene_trans_map \\
  --CPU 6 \\
  2>trinotate_${job_name}.stderr \\
  1>trinotate_${job_name}.stdout
mv Trinotate.tsv ${output_name}
cd \${start}
!;
    my $trinotate = $class->Submit(
        cpus => 6,
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => qq"trinotate_$input_paths->{filename}_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => 12,
        modules => $options->{modules},
        output => qq"${output_dir}/${output_name}",
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        jqueue => 'large',
        jwalltime => '144:00:00',
        stdout => $stdout,
        stderr => $stderr);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($trinotate);
}

=head2 C<Rosalind_Plus>

 Attempt to ensure that the plus strand has the most putative ORFs.

 Use prodigal to count up putative ORFs on the plus and minus strands.
 If the number of minus strand ORFs is larger than plus, reverse
 complement the sequence.  This function just puts the actual function
 which does the work onto the dependency chain.

 Just a little note from a friend: "This is my formal proposal to start
 calling the 'Watson' and 'Crick' strands 'Rosalind' and 'Franklin'."

=item C<Arguments>

 input(required): Input assembly.
 This should handle all the prodigal options too, FIXME!

=cut
sub Rosalind_Plus {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jprefix => '44',);

    my $input_seq = $options->{input};
    my $job_name = 'rosalindplus';
    my $output_dir = qq"outputs/$options->{jprefix}${job_name}";
    make_path($output_dir);
    my $rosalind_prodigal = Bio::Adventure::Feature_Prediction::Prodigal($class,
        gcode => '11',
        input => $options->{input},
        jdepends => $options->{jdepends},
        jname => $job_name,
        jprefix => $options->{jprefix},
        modules => ['prodigal'],
        output_dir => $output_dir,
        species => 'phages',);
    $options->{jdepends} = $rosalind_prodigal->{job_id};
    my $input_gff = $rosalind_prodigal->{output_gff};
    my $output_file = basename($options->{input});
    $output_file = qq"${output_dir}/${output_file}";
    my $comment_string = qq"## This takes the prodigal output and checks to see which strand has
## more ORFs, if it is the Franklin strand, then the chromsome is reverse complemented.
";
    my $jstring = qq!
use Bio::Adventure::Annotation;
my \$result = \$h->Bio::Adventure::Annotation::Rosalind_Plus_Worker(
  gff => '${input_gff}',
  input => '$options->{input}',
  jdepends => '$options->{jdepends}',
  jname => '$options->{jname}',
  jprefix => '$options->{jprefix}',
  output => '${output_file}',
  output_dir => '${output_dir}',);
!;
    my $log = qq"${output_dir}/final_assembly.log";
    my $rewrite = $class->Submit(
        comment => $comment_string,
        gff => $options->{gff},
        input => $options->{input},
        jdepends => $options->{jdepends},
        jmem => 8,
        jname => $options->{jname},
        jprefix => $options->{jprefix},
        jstring => $jstring,
        language => 'perl',
        log => $log,
        output => $output_file,
        output_dir => $output_dir,
        shell => '/usr/bin/env perl',);
    return($rewrite);
}

=head2 C<Rosalind_Plus_Worker>

 Read results from prodigal, count the +/- strands, and flip the
 genome if the number of - is greater than +.  Does the actual work of
 rewriting a genome to put the majority ORFs on the plus strand.
 Rosalind_Plus() calls this.

=cut
sub Rosalind_Plus_Worker {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],);

    my $input_gff = $options->{gff};
    my $input_seq = $options->{input};
    my $job_name = 'rosalindplus';
    my $output_dir = qq"outputs/$options->{jprefix}${job_name}";
    make_path($output_dir);

    my $orf_count = {};
    my $annotation_in = FileHandle->new("<${input_gff}");
    my $write_file = basename($options->{input});
    $write_file = qq"${output_dir}/${write_file}";
    my $log_file = basename($write_file, ('.fasta'));
    $log_file = qq"${output_dir}/${log_file}.log";
    my $read_fasta = Bio::SeqIO->new(-file => qq"<$options->{input}", -format => 'Fasta');
    my $write_fasta = Bio::SeqIO->new(-file => qq">${write_file}", -format => 'Fasta');
    my $write_log = FileHandle->new(">${log_file}");
    print $write_log "Counting ORFs on the current Rosalind and Franklin strand.
If the Franklin strand is larger, flipping them.
";
    ## A prodigal line looks like:
    ## EAb06   Prodigal_v2.6.3 CDS     3205    3414    8.7     -       0       ID=1_10;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.248;conf=88.17;score=8.74;cscore=-0.38;sscore=9.12;rscore=8.34;uscore=-0.93;tscore=2.86;
    while (my $line = <$annotation_in>) {
        next if ($line =~ /^#/);
        my ($name, $prod, $type, $start, $end, $score, $strand, $phase, $annot) = split(/\t/, $line);
        my @annot_list = split(/;/, $annot);
        my $id_annot = $annot_list[0];
        my ($contig, $orf);
        if ($id_annot =~ /ID=(.*)?_(\d+)$/) {
            $contig = $1;
            $orf = $2;
        } else {
            die('Could not get contig and orf.');
        }
        if (!defined($orf_count->{$contig})) {
            $orf_count->{$contig} = {
                plus => 0,
                minus => 0,
            };
        }
        if ($strand eq '+' || $strand eq '1') {
            $orf_count->{$contig}->{plus}++;
        } else {
            $orf_count->{$contig}->{minus}++;
        }
    } ## End iterating over every detected ORF
    $annotation_in->close();
    ## Now we should have some idea of which strand has the most ORFs for each contig.

  TESTLOOP: while (my $seq = $read_fasta->next_seq) {
      my $id = $seq->id;
      my $test = $orf_count->{$id};
      $test->{plus} = 0 if (!defined($test->{plus}));
      $test->{minus} = 0 if (!defined($test->{minus}));
      if ($test->{plus} >= $test->{minus}) {
          ## Then leave this contig alone.
          $write_fasta->write_seq($seq);
          print $write_log "${id} was unchanged and has $test->{plus} plus and $test->{minus} minus ORFs.\n";
      } else {
          my $tmp = $seq->revcom();
          $write_fasta->write_seq($tmp);
          print $write_log "${id} was reverse-complemented and now has $test->{minus} plus and $test->{plus} minus ORFs.\n";
      }
  } ## Done flipping sequence files.
    $write_log->close();
    return($orf_count);
}

1;

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

=cut
