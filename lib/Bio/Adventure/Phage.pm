package Bio::Adventure::Phage;
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

use Bio::DB::EUtilities;
use Bio::Restriction::Analysis;
use Bio::Restriction::Enzyme;
use Bio::Restriction::EnzymeCollection;
use Bio::SeqFeature::Generic;
use Bio::Tools::Run::Alignment::StandAloneFasta;
use Bio::Tools::Run::StandAloneBlastPlus;
use Capture::Tiny qw":all";
use Cwd qw"abs_path getcwd cwd";
use File::Basename;
use File::Spec;
use File::Path qw"make_path rmtree";
use File::Which qw"which";
use File::ShareDir qw":ALL";
use Text::CSV qw"csv";
use WWW::Mechanize;

=head2 C<Bacphlip>

 Predict lytic vs. lysogenic phages via gene-based ML.

 Bacphlip seeks to use a set of training-generated gene and
 COG-defined categories in order to initialize a classifier which,
 given a single-contig assembly, searches all reading frames, looks
 for hits associated with lytic/lysogenic phages, and defines the most
 likely category thereof.
 It is therefore a cousin to the question I have been posed by
 the people at WRAIR.  Much more encouraging to my eyes, the author
 did a beautiful job of documenting his code and methodologies.  I am
 reasonably certain I can use this work to learn how to actually
 perform similar tasks in other contexts; as a result, I think any
 real ML tasks I succeed in will, by definition, require a reference
 to this.

=cut
sub Bacphlip {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        jmem => 12,
        jname => 'bacphlip',
        jprefix => '80',
        modules => ['bacphlip'],
        executables => ['bacphlip'],);
    my $loaded = $class->Module_Loader(
        modules => $options->{modules},
        executables => $options->{executables},);
    my @suffixes = split(/,/, $options->{suffixes});
    my $out_base = basename($options->{input}, @suffixes);
    my $in_base = basename($options->{input});
    my $output_dir = qq"outputs/$options->{jprefix}$options->{jname}";
    my $output = qq"${output_dir}/${in_base}.bacphlip";
    my $output_hmm = qq"${output_dir}/${in_base}.hmmsearch";
    my $output_frames = qq"${output_dir}/${in_base}.6frame";
    my $output_tsv = qq"${output_dir}/${in_base}.hmmsearch.tsv";
    my $stderr = qq"${output_dir}/${out_base}.stderr";
    my $stdout = qq"${output_dir}/${out_base}.stdout";
    my $comment = qq'## Running bacphlip using $options->{input}.';
    my $jstring = qq?mkdir -p ${output_dir}
cp $options->{input} ${output_dir}
start=\$(pwd)
contig_number=\$(grep "^>" $options->{input} | wc -l)
multi_arg=""
if (( "\${contig_number}" > 1 )); then
  multi_arg=" --multi_fasta "
fi
cd ${output_dir}
bacphlip -f \${multi_arg} -i ${in_base} \\
  2>${out_base}.stderr \\
  1>${out_base}.stdout
cd \${start}
?;
    my $job = $class->Submit(
        comment => $comment,
        output => $output,
        output_hmm => $output_hmm,
        output_frames => $output_frames,
        output_tsv => $output_tsv,
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
        jname => $options->{jname},
        jprefix => $options->{jprefix},
        jstring => $jstring,);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload',);
    return($job);
}

=head2 C<Caical>

 Use caical against a reference species.
 10.1186/1745-6150-3-38

 It takes a bit of hunting on their web server to find the actual
 script they use, it is somewhere at the bottom.  Once found, with a
 little work it can be made to run successfully.  This script does
 just that, given an input set of ORFs and host species, it will
 create a codon usage table of the host and run caical of the input
 against it.

 The caical perl script comes with a fair number of arguments, this
 just blindly uses the defaults.

=item C<Arguments>

 input(required): Input set of CDS sequences.
 species(required): Reference species to query.
 jmem(4): Expected memory.
 jprefix(80): Jobname/directory prefix.
 modules(caical): Environment module to load.

=cut

sub Caical {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        jmem => 4,
        jprefix => '80',
        modules => ['caical']);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('caical');
    die('Could not find caical in your PATH.') unless($check);
    my $test = ref($options->{suffixes});
    my @suffixes = split(/,/, $options->{suffixes});
    my $species = $options->{species};
    my $out_base = basename($options->{input}, @suffixes);
    my $jname = qq"caical_${out_base}_vs_${species}";
    my $output_dir = qq"outputs/$options->{jprefix}${jname}";
    make_path($output_dir);
    my $output_cai = qq"${output_dir}/${out_base}_cai.txt";
    my $random_sequences = qq"${output_dir}/${out_base}_random_sequences.txt";
    my $expected_cai = qq"${output_dir}/${out_base}_expected.txt";
    my $stderr = qq"${output_dir}/${out_base}.stderr";
    my $stdout = qq"${output_dir}/${out_base}.stdout";
    my $comment = qq'## Running caical against ${species}.';

    my $jstring = qq!use Bio::Adventure::Phage;
my \$result = \$h->Bio::Adventure::Phage::Caical_Worker(
  input => '$options->{input}',
  random_sequences => '$random_sequences',
  expected_cai => '$expected_cai',
  species => '$options->{species}',
  jdepends => '$options->{jdepends}',
  jname => '${jname}',
  jprefix => '$options->{jprefix}',
  output => '${output_cai}',
  stdout => '${stdout}',
  stderr => '${stderr}',
  output_dir => '${output_dir}',);
!;
    my $cai = $class->Submit(
        comment => $comment,
        input => $options->{input},
        random_sequences => $random_sequences,
        expected_cai => $expected_cai,
        species => $species,
        jdepends => $options->{jdepends},
        jname => $jname,
        jprefix => $options->{jprefix},
        jstring => $jstring,
        output => $output_cai,
        stdout => $stdout,
        stderr => $stderr,
        language => 'perl',
        output_dir => $output_dir);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload',);
    return($cai);
}

sub Caical_Worker {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        jmem => 4,
        random_sequences => '',
        expected_cai => '',
        output => '',
        jprefix => '80',
        modules => ['caical']);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('caical');
    die('Could not find caical in your PATH.') unless($check);
    my $species = $options->{species};
    ## Then the species argument is a filename, so pull the species from it.
    if (-r $species) {
        my $in = FileHandle->new("<$species");
        while (my $line = <$in>) {
            chomp $line;
            next if ($line =~ /^\s*$/);
            $species = $line;
        }
        $in->close();
    }

    my $index = qq"$options->{libpath}/codon_tables/${species}.txt";
    if (!-r $index) {
        my $written = $class->Bio::Adventure::Index::Make_Codon_Table(
            species => $species);
    }

    my $jstring = qq?caical -g 11 \\
  -f $options->{input} \\
  -h ${index} \\
  -o1 $options->{output} \\
  -o2 $options->{random_sequences} \\
  -o3 $options->{expected_cai} \\
  2>$options->{stderr} 1>$options->{stdout}
?;
    my $ran_caical = qx"$jstring";
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($ran_caical);
}

## We may need to rewrite the input fasta file because caical
## cries if a sequence is not %3 == 0.  If so, use this.
sub Caical_CDS {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'species'],
        jmem => 4,
        jprefix => '80',
        modules => ['caical']);
    my $seqio_in = Bio::SeqIO->new(-format => 'fasta', -file => $options->{input});
    my $output_dir = qq"";
    my $out_base = '';
    my $cai_fasta = qq">${output_dir}/${out_base}_input.fasta";
    my $seqio_out = Bio::SeqIO->new(-format => 'fasta', -file => $cai_fasta);
  SEQS: while (my $seq = $seqio_in->next_seq) {
      my $len = $seq->length;
      if (($len % 3) == 0) {
          $seqio_out->write($seq);
      }
  }
}

=head2 C<Classify_Phage>

 Use the ICTV data to attempt to classify a viral assembly.
 https://talk.ictvonline.org/taxonomy/

 Use tblastx against a database extracted from reference accessions
 provided by the ICTV.  This should hopefully provide a reasonably
 complete set of taxonomic classifications for viral assemblies.  This
 function pretty much just calls Blast_Classify().

=item C<Arguments>

 input(required): Input assembly.
 evalue(0.01): Cutoff evalue.
 blast_tool(tblastx): Choose a blast method when searching.
 library(ictv): Blast database basename.
 topn(5): Keep this number of hits.
 jmem(12): Expected memory usage.
 jprefix('18'): Job/output prefix.
 modules(blastdb, blast): Environment modules to load.

=cut
sub Classify_Phage {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        evalue => 0.01,
        blast_tool => 'tblastx',
        library => 'ictv',
        topn => 5,
        jmem => 12,
        jprefix => '18',
        modules => ['blastdb', 'blast'],);

    my $input_dir = basename(dirname($options->{input}));
    my $output_dir = qq"outputs/$options->{jprefix}classify_${input_dir}";
    if (-d $output_dir) {
        my $removed = rmtree($output_dir);
    }
    my $output_tsv = qq"${output_dir}/$options->{library}_filtered.tsv";
    my $output_blast = qq"${output_dir}/$options->{library}_hits.txt";
    my $output_log = qq"${output_dir}/classify.log";
    my $comment = qq"## This will perform a tblastx search against a manually curated set of
## ICTV taxonomies and attempt to provide the most likely classification for this assembly.\n";
    my $jstring = qq?
use Bio::Adventure;
use Bio::Adventure::Phage;
my \$result = Bio::Adventure::Phage::Classify_Phage_Worker(\$h,
  blast_tool => '$options->{blast_tool}',
  evalue => '$options->{evalue}',
  input => '$options->{input}',
  library => '$options->{library}',
  output => '${output_tsv}',
  output_blast => '${output_blast}',
  output_dir => '${output_dir}',
  output_log => '${output_log}',
  topn => '$options->{topn}',);
?;
    my $cjob = $class->Submit(
        blast_tool => $options->{blast_tool},
        evalue => $options->{evalue},
        input => $options->{input},
        comment => $comment,
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
        jname => 'classify_ictv',
        jprefix => $options->{jprefix},
        jstring => $jstring,
        language => 'perl',
        library => $options->{library},
        modules => $options->{modules},
        output => $output_tsv,
        output_blast => $output_blast,
        output_dir => $output_dir,
        output_log => $output_log,
        topn => $options->{topn},);
    return($cjob);
}

=head2 C<Classify_Phage_Worker>

 Does the actual work of attempting to classify a viral sequence
 against the ICTV database.

 In order for this to actually give a hit genus/species/accession/etc,
 it must read a CSV copy of the ICTV data, which was copied to the
 blast database directory as 'ictv.csv'.  It extracts from that csv
 file the various taxonomy data and will cross reference the results
 of the blast search against that data structure.

 It then invokes blast, parses the output, and pulls out the hits
 which pass the various cutoffs.  Upon completion, it collects the
 results, sorts them, and writes them out as a tab separated file.

=item C<Arguments>

 input(required): Input assembly to search against.
 evalue(0.01): Cutoff evalue.
 blast_tool(tblastx): I experimented with a few aligners and blast tools.
 library(ictv): Database name to search against.
 output_log(classify.log): Log filename to write.
 output_blast(ictv_hits.txt): Filename to write the blast hits.
 output_dir('.'): Directory to write the results.
 output('ictv_filtered.tsv'): TSV file to write the final results.
 score(1000): Cutoff score -- this is the stringent portion.
 topn(5): Keep this number of hits passing the score/evalue cutoffs.
 jprefix('18'): Job/directory prefix.
 jcpus(6): Use this number of cpus.
 modules('blast', 'blastdb'): Environment modules to load.

=cut
sub Classify_Phage_Worker {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input',],
        evalue => 0.01,
        blast_tool => 'tblastx',
        library => 'ictv',
        output_log => 'classify.log',
        output_blast => 'ictv_hits.txt',
        output_dir => '.',
        output => 'ictv_filtered.tsv',
        score => 1000,
        topn => 5,
        jprefix => '18',
        jcpus => 6,
        modules => ['blast', 'blastdb'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which($options->{blast_tool});
    die("Could not find $options->{blast_tool} in your PATH.") unless($check);

    ## Read the xref file, located in the blast database directory as ${library}.csv
    my $xref_file = qq"$ENV{BLASTDB}/$options->{library}.csv";
    ## Use Text::CSV to create an array of hashes with keynames:
    ## 'taxon', 'virusnames', 'virusgenbankaccession', 'virusrefseqaccession'
    my $xref_aoh = csv(in => $xref_file, headers => 'auto');
    my $log = FileHandle->new(">$options->{output_log}");
    ## First check for the test file, if it exists, then this should
    ## just symlink the input to the output until I think of a smarter
    ## way of thinking about this.
    my $blast_output = Bio::SearchIO->new(-format => 'blast', );
    my $number_hits = 0;
    my $blast_outfile = qq"$options->{output_blast}";
    my $final_fh = FileHandle->new(">$options->{output}");
    ## Print the tsv header: contig, description, taxon,length
    print $final_fh "contig\tquery_description\ttaxon\tname\tquery_length\thit_length\thit_accession\thit_description\thit_bit\thit_sig\thit_score\thit_family\thit_genus\n";
    my @params = (
        -e => $options->{evalue},
        -db_name => $options->{library},
        -outfile => $blast_outfile,
        -num_threads => $options->{jcpus},
        -program => $options->{blast_tool},);
    my $seq_count = 0;
    my $e = '';

    my $log_message = qq"Starting blast search of $options->{input}
against $options->{library} using tool: $options->{blast_tool} with $options->{jcpus} cpus.
Writing blast results to $options->{output_blast}.
Writing filtered results to $options->{output}.
";
    print $log $log_message;
    my $search = Bio::Tools::Run::StandAloneBlastPlus->new(@params);
    my @parameters = $search->get_parameters;
    print $log qq"Blast parameters: @parameters.\n";

    my $blast_report = $search->tblastx(
        -query => $options->{input},
        -outfile => $blast_outfile,
        -method_args => [ '-num_alignments' => $options->{topn},
                          '-num_threads' => $options->{jcpus}, ]);

    my $result_data = {};
    my $search_output = Bio::SearchIO->new(-file => $blast_outfile, -format => 'blast');
    ## Set some filters for the output.
    $search_output->min_score($options->{score});
    $search_output->max_significance($options->{evalue});
    my $element_count = 0;
    my $result_count = 0;

    my @hit_lst = ();
  RESULTLOOP: while (my $result = $search_output->next_result) {
      $result_count++;
      my $query_name = $result->query_name();
      my $query_length = $result->query_length();
      my $query_descr = $result->query_description();
      my $stats = $result->available_statistics();
      my $hits = $result->num_hits();
      print $log "The query being considered is: ${query_name}.\n";
      my $datum_ref = {
          description => $query_descr,
          name => $query_name,
          stats => $stats,
          num_hits => $hits,
          length => $query_length,
          hit_data => {},
      };
      $result_data->{$query_name} = $datum_ref;
      my $hit_count = 0;
    HITLOOP: while (my $hits = $result->next_hit()) {
        $number_hits = $number_hits++;
        my $hit_name = $hits->name();
        print $log "This result has a hit on ${hit_name}\n";
        ## The hit_name should cross reference to one of the two accession columns in xref_aoh
        ## from the beginning of this function.
        my $longname = '';
        my ($taxon, $family, $genus);
        my $xref_found = 0;
      XREFLOOP: for my $hashref (@{$xref_aoh}) {
          my $first_name = $hashref->{virusgenbankaccession};
          my $second_name = $hashref->{virusrefseqaccession};
          my $hash_taxon = $hashref->{taxon};
          my $hash_longname = $hashref->{virusnames};
          my $hash_family = $hashref->{family};
          my $hash_genus = $hashref->{genus};
          if ($hit_name eq $first_name or $hit_name eq $second_name) {
              $xref_found++;
              print $log "A cross reference was found to the ICTV taxon: ${hash_taxon}.\n";
              $longname = $hash_longname;
              $taxon = $hash_taxon;
              $family = $hash_family;
              $genus = $hash_genus;
              last XREFLOOP;
          }
      }
        my $hit_length = 0;
        if (defined($hits->length())) {
            $hit_length = $hits->length();
        }
        my $hit_acc = $hits->accession();
        my $hit_descr = $hits->description();
        my $hit_score = $hits->score();
        my $hit_sig = $hits->significance();
        my $hit_bits = $hits->bits();
        my %hit_datum = (
            query_name => $query_name,
            query_length => $query_length,
            query_descr => $query_descr,
            stats => $stats,
            length => $hit_length,
            acc => $hit_acc,
            description => $hit_descr,
            score => $hit_score,
            sig => $hit_sig,
            bit => $hit_bits,
            longname => $longname,
            taxon => $taxon,
            family => $family,
            genus => $genus);
        push (@hit_lst, \%hit_datum);
        $result_data->{$query_name}->{hit_data}->{$hit_name} = \%hit_datum;
        $hit_count++;
    } ## End the hitloop.
  } ## End of the blast search

    print $log "The number of results recorded: ${result_count}.\n";
    ## Sort largest to smallest score.
    my @sorted = sort { $b->{score} <=> $a->{score} } @hit_lst;
    for my $d (@sorted) {
        my $hit_string = qq"$d->{query_name}\t$d->{query_descr}\t$d->{taxon}\t$d->{longname}\t$d->{length}\t$d->{query_length}\t$d->{acc}\t$d->{description}\t$d->{bit}\t$d->{sig}\t$d->{score}\t$d->{family}\t$d->{genus}\n";
        print $final_fh $hit_string;
    }
    $final_fh->close();
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    $log->close();
    return($result_data);
}

=head2 C<Download_NCBI_Assembly_UID>

 Perform a web crawl to download a reference genome from NCBI.

 In theory, the various entrez modules should make downloading
 genomes/etc from NCBI relatively easy.  In practice, they are a right
 pain in the arse.  I fought for hours with them to download reference
 genomes and only ever got it working about 1 time in 4.  So, I
 decided to say 'screw it' and write a simple web crawler to do it
 instead.

=item C<Arguments>

 uid: a UID to download.

=cut
sub Download_NCBI_Assembly_UID {
    my %args = @_;
    my $out_fh = $args{out};

    ## If this assembly has not already been downloaded, get it.
    my $mech = WWW::Mechanize->new(autocheck => 1);
    my $url = qq"https://www.ncbi.nlm.nih.gov/assembly?LinkName=genome_assembly&from_uid=$args{uid}";
    my $data = $mech->get($url);
    my @links = $mech->find_all_links(
        tag => 'a', text_regex => qr/ASM/i);
  FIRST: for my $k (@links) {
      my ($url2, $title2) = @{$k};
  }
    my $first = $links[0];
    my $title = '';
    ($url, $title) = @{$first};
    my $accession = basename($url);
    return($accession);
}

=head2 C<Download_NCBI_Assembly_Accession>

 Perform a web crawl to download a reference genome from NCBI.  This
 is similar in idea to the download_uid function above, but the logic
 is slightly different to get an accession.

=cut
sub Download_NCBI_Assembly_Accession {
    my %args = @_;
    my $out_fh = $args{out};
    my $accession = $args{accession};
    ## If this assembly has not already been downloaded, get it.
    my $mech = WWW::Mechanize->new(autocheck => 1);
    my $url = qq"https://www.ncbi.nlm.nih.gov/assembly/${accession}";
    print $out_fh "Searching ${url} for appropriate download links.\n";

    my $data = $mech->get($url);
    my @links = $mech->find_all_links(
        tag => 'a', text_regex => qr/FTP/i);
    my $title;
  LINKS: foreach my $l (@links) {
      ($url, $title) = @{$l};
      if ($title eq 'FTP directory for GenBank assembly') {
          my $followed = $mech->follow_link(url => $url);
          my @final = $mech->find_all_links(
              tag => 'a', text_regex => qr/gbff.gz$/i);
          print $out_fh "Final download link: $final[0][0].\n";
          my $final_followed = $mech->follow_link(url => $final[0][0]);
          my $output = FileHandle->new(">$args{file}");
          my $printed = $mech->response->content();
          print $output $printed;
          $output->close();
          last LINKS;
      }
  } ## End looking at links hopefully containing an assembly.
    return($title);
}

=head2 C<Filter_Host_Kraken>

 Set the host species to the most represented by kraken, and filter.

 This function reads the results of a standard kraken run and pulls out
 the most represented species.  It then goes to NCBI and crawls to
 find a reference genome for that species, downloads it, and indexes
 it.  Finally, it runs hisat2 using the reads and this downloaded
 genome as a reference.  Then the unaligned reads are hopefully only
 from the phage and not the host bacterium.

=item C<Arguments>

 input(required): The kraken report.
 input_fastq(required): The trimmed/corrected reads.
 jdepends(''): Job this depends upon.
 jmem(8): Expected memory usage.
 jprefix('06'): Job/directory prefix.

=cut
sub Filter_Host_Kraken {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'input_fastq'],
        jdepends => '',
        jmem => 8,
        jprefix => '06',);
    my $comment = '## Use kraken results to choose a host species to filter against.';
    my $output_dir = qq"outputs/$options->{jprefix}filter_kraken_host";
    make_path($output_dir);
    my $out_r1_name = 'r1_host_filtered.fastq.xz';
    my $out_r2_name = 'r2_host_filtered.fastq.xz';
    my $output_files = qq"${output_dir}/${out_r1_name}:${output_dir}/${out_r2_name}";
    my $log = qq"${output_dir}/kraken_filter.log";
    my $jstring = qq!
use Bio::Adventure;
use Bio::Adventure::Phage;
my \$result = Bio::Adventure::Phage::Filter_Kraken_Worker(\$h,
  output => '${output_files}',
  output_dir => '${output_dir}',
  input => '$options->{input}',
  input_fastq => '$options->{input_fastq}',
  jdepends => '$options->{jdepends}',
  jname => 'kraken_host',
  jprefix => '$options->{jprefix}',
  job_log => '${log}');
!;
    my $host = $class->Submit(
        input => $options->{input},
        input_fastq => $options->{input_fastq},
        job_log => $log,
        output => $output_files,
        output_dir => $output_dir,
        comment => $comment,
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
        jname => 'hostfilter',
        jprefix => $options->{jprefix},
        jstring => $jstring,
        language => 'perl',);
    return($host);
}

=head2 C<Filter_Kraken_Worker>

 Do the work for Filter_Kraken().

 This does the actual work for the previous function.

=item C<Arguments>

 input(required): The kraken report file.
 output(required): The unaligned reads produced by hisat2.
 jname(krakenfilter): Job name (maybe not needed?)
 type(species): Which element from the kraken report to extract?

=cut
sub Filter_Kraken_Worker {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'output'],
        job_log => 'filter_kraken.log',
        jname => 'krakenfilter',
        jprefix => '06',
        type => 'species',);
    my $separator;
    if ($options->{type} eq 'domain') {
        $separator = 'd__';
    } elsif ($options->{type} eq 'phylum') {
        $separator = 'p__';
    } elsif ($options->{type} eq 'class') {
        $separator = 'c__';
    } elsif ($options->{type} eq 'order') {
        $separator = 'o__';
    } elsif ($options->{type} eq 'family') {
        $separator = 'f__';
    } elsif ($options->{type} eq 'genus') {
        $separator = 'g__';
    } elsif ($options->{type} eq 'species') {
        $separator = 's__';
    } else {
        $separator = 's__';
    }

    my $output_dir = $options->{output_dir};
    my $out = FileHandle->new(">$options->{job_log}");
    my $in = FileHandle->new("<$options->{input}");
    print $out "Starting search for best kraken host strain.\n";
    print $out "Reading kraken report: $options->{input}.\n";
    my %species_observed = ();
    my $most_species = '';
    my $most_observations = 0;
    my $line_count = 0;
    while (my $line = <$in>) {
        chomp $line;
        next unless ($line =~ m/$separator/);
        $line_count++;
        my ($entry, $number) = split(/\t/, $line);
        my ($stuff, $species) = split(/$separator/, $entry);
        $species_observed{$species} = $number;
        if ($number > $most_observations) {
            $most_species = $species;
            $most_observations = $number
        }
    }
    foreach my $s (keys %species_observed) {
        print $out "$options->{type} ${s} was observed $species_observed{$s} times.\n";
    }
    print $out "\n";
    $species_observed{most} = {
        species => $most_species,
        observations => $most_observations, };
    print $out "The most observed species was: ${most_species}.\n";
    print STDOUT "The most observed species was: ${most_species}.\n";
    my $escaped_species = $most_species;
    $escaped_species =~ s/\s/\+/g;

    my $host_species_accession = undef;
    if (-f 'host_species.txt') {
        my $host_file = FileHandle->new("<host_species.txt");
        while (my $line = <$host_file>) {
            chomp $line;
            $host_species_accession = $line;
        }
        $host_file->close();
        print $out "The host_species.txt file already exists and suggests using ${host_species_accession}.\n";
    }

    my $need_download = 0;
    if (defined($host_species_accession)) {
        my $accession_file = qq"$options->{libpath}/$options->{libtype}/${host_species_accession}.fasta";
        if (-r $accession_file) {
            print $out "Found fasta file for the putative host: ${host_species_accession}.\n";
            print "Found fasta file for the putative host: ${host_species_accession}.\n";
        } else {
            $need_download = 1;
        }
    } else {
        print $out "Downloading new assembly for ${most_species}.\n";
        print "Downloading new assembly for ${most_species}.\n";
        $need_download = 1;
    }

    if ($need_download) {
        my $search_url = qq"https://www.ncbi.nlm.nih.gov/assembly/?term=${escaped_species}";
        my $mech = WWW::Mechanize->new(autocheck => 1);
        print $out "Searching ${search_url} for appropriate download links.\n";
        my $search_data = $mech->get($search_url);
        my @search_links = $mech->find_all_links(
            tag => 'a', text_regex => qr/ASM/i);
        my $first_hit = $search_links[0];
        my ($assembly_link, $assembly_title) = @{$first_hit};
        my $accession = basename($assembly_link);
        if ($accession =~ /\,/) {
            my @multiple = split(/\,/, $accession);
            $accession = $multiple[0];
        }
        my $downloaded_file = qq"$options->{libpath}/$options->{libtype}/${accession}.gbff.gz";
        if (-r $downloaded_file) {
            print STDOUT "The file: ${downloaded_file} already exists.\n";
            print $out "The file: ${downloaded_file} already exists.\n";
        } else {
            my $assembly_url = qq"https://www.ncbi.nlm.nih.gov/assembly/${accession}";
            print $out "Searching ${assembly_url} for appropriate download links.\n";
            my $assembly_data = $mech->get($assembly_url);
            my @download_links = $mech->find_all_links(
                tag => 'a', text_regex => qr/FTP/i);
          LINKS: foreach my $l (@download_links) {
              my ($download_url, $download_title) = @{$l};
              if ($download_title eq 'FTP directory for GenBank assembly') {
                  my $ftp_followed = $mech->follow_link(url => $download_url);
                  my @download_links = $mech->find_all_links(
                      tag => 'a', text_regex => qr/gbff.gz$/i);
                  print $out "Final download link: $download_links[0][0].\n";
                  my $download_followed = $mech->follow_link(url => $download_links[0][0]);
                  my $output = FileHandle->new(">${downloaded_file}");
                  my $printed = $mech->response->content();
                  print $output $printed;
                  $output->close();
                  last LINKS;
              }
          } ## End looking at links hopefully containing an assembly.
        } ## End of when the download file does not exist.

        print $out "Writing a new host_species.txt file with: ${accession}.\n";
        my $host = FileHandle->new(">host_species.txt");
        print $host qq"${accession}\n";
        $host->close();
        $host_species_accession = $accession;

        ## Convert the assembly to fasta/gff/etc.
        print $out "Converting ${downloaded_file} assembly to fasta/gff.\n";
        my $test_name = basename($downloaded_file, ('.gz', '.gbff'));
        my $test_input = qq"$options->{libpath}/$options->{libtype}/${test_name}.fasta";
        my $converted;
        if (-r $test_input) {
            $converted = $test_input;
        } else {
            $converted = $class->Bio::Adventure::Convert::Gb2Gff(input => $downloaded_file);
        }
    } ## End checking if the host_species was defined.

    my $cyoa_shell = Bio::Adventure->new(cluster => 0);
    my $index_input = qq"$options->{libpath}/$options->{libtype}/${host_species_accession}.fasta";
    my $index_location = qq"$options->{libpath}/$options->{libtype}/indexes/${host_species_accession}.1.ht2";
    if (-r $index_location) {
        print $out "Found indexes at: ${index_location}\n";
        print "Found indexes at: ${index_location}\n";
    } else {
        print $out "Did not find indexes, running Hisat2_Index() now.\n";
        print "Did not find indexes, running Hisat2_Index() now.\n";
        my $indexed = $cyoa_shell->Bio::Adventure::Index::Hisat2_Index(
            input => $index_input,
            output_dir => $options->{output_dir});
        ## Check to see if there is a text file containing the putative host species
        ## If it does not exist, write it with the species name.
    }
    ## Now perform the filter using our temp_cyoa.
    print $out "Filtering out reads which map to ${host_species_accession}.\n";
    my $filter = $cyoa_shell->Bio::Adventure::Map::Hisat2(
        input => $options->{input_fastq},
        output_dir => $options->{output_dir},
        jprefix => $options->{jprefix},
        do_htseq => 0,
        species => $host_species_accession,);
    my $filtered_reads = $filter->{unaligned_comp};
    my ($in_r1, $in_r2) = split(/\:|\;|\,|\s+/, $filtered_reads);
    $in_r1 = File::Spec->rel2abs($in_r1);
    $in_r2 = File::Spec->rel2abs($in_r2);
    my ($out_r1, $out_r2) = split(/\:|\;|\,|\s+/, $options->{output});
    print $out "Symlinking final output files to $options->{output_dir}\n";
    unlink $out_r1 if (-l $out_r1);
    unlink $out_r2 if (-l $out_r2);
    if (-e $out_r1) {
        print "The file: $out_r1 already exists.\n";
    } else {
        my $s1 = symlink($in_r1, $out_r1);
    }
    if (-e $out_r2) {
        print "The file: $out_r2 already exists.\n";
    } else {
        my $s2 = symlink($in_r2, $out_r2);
    }
    return(%species_observed);
}

=head2 C<Get_DTR>

 Extract the direct-terminal-repeats from a phageterm run.

 Given the lengths we go in order to make phageterm actually run
 successfully, it makes some sense that we should also have to work
 to get the resulting direct terminal repeats from it.  This function
 seeks to extract them.  It reads the resulting fasta file and creates
 a SeqFeature from it.

 This is slightly complicated by the fact that we cannot trust
 phageterm to write a single copy, two, three, or four on the ends of
 the assembly; but instead it adds a variable number depending on the
 type of DTR.  As a result, this function includes a regex-based
 search of the output assembly to find the n copies of the DTR and set
 the DTR features accordingly.

=item C<Arguments>

 input_fsa: Fasta file of the reorganized genome produced by
  phageterm.
 input_dtr: Smaller fasta file containing the dtr sequence.

=cut
sub Get_DTR {
    my ($class, %args) = @_;
    my $input_fsa = $args{input_fsa};
    my $input_dtr = $args{input_dtr};
    my $log_fh = $args{log_fh};
    my $dtr_type_file = $input_dtr;
    $dtr_type_file =~ s/_dtr\.fasta/_nrt\.txt/g;
    return(undef) unless (-r $args{input_fsa});
    return(undef) unless (-r $args{input_dtr});
    return(undef) unless (-f $dtr_type_file);

    my $dtr_type_read = FileHandle->new("<${dtr_type_file}");
    my $dtr_type = '';
    while (my $line = <$dtr_type_read>) {
        chomp $line;
        $dtr_type = $line;
    }
    $dtr_type_read->close();
    print $log_fh "Got DTR type: ${dtr_type}.\n";

    my $dtr_read = Bio::SeqIO->new(-file => $input_dtr, -format => 'Fasta');
    my $dtr_sequence = '';
    my $dtr_length = 0;
    my $dtr_id = '';
  DTR: while (my $dtr_seq = $dtr_read->next_seq()) {
      next DTR unless(defined($dtr_seq->id));
      $dtr_id = $dtr_seq->id;
      $dtr_sequence = $dtr_seq->seq;
      $dtr_length = $dtr_seq->length;
  }

    my @dtr_features = ();
    my $fsa_read = Bio::SeqIO->new(-file => $input_fsa, -format => 'Fasta');
  FSA: while (my $genome_seq = $fsa_read->next_seq()) {
      my $contig_sequence = $genome_seq->seq;
      my $contig_id = $genome_seq->id;
    DTR_SEARCH: while ($contig_sequence =~ m/$dtr_sequence/g) {
        my $dtr_end = pos($contig_sequence);
        my $dtr_start = $dtr_end - ($dtr_length - 1);
        my $dtr_feature = Bio::SeqFeature::Generic->new(
            -primary => 'misc_feature',
            -seq_id => $contig_id,
            -source => 'PhageTerm',
            -start => $dtr_start,
            -end => $dtr_end,
            -strand => +1,
            -score => undef,
            -frame => 0,
            -tag => {
                'product' => 'Direct Terminal Repeat',
                'inference' => 'COORDINATES:profile:PhageTerm',
                'note' => qq"DTR type: ${dtr_type}",
            },);
        push(@dtr_features, $dtr_feature);
    } ## End matching on this contig
  } ## End iterating over teh contigs
    return(\@dtr_features);
}

=head2 C<Phageterm>

 Attempt to run phageterm and detect direct terminal repeats.
 10.1038/s41598-017-07910-5

 Use phageterm on a viral assembly to look for terminal repeat
 regions.  If they are found, rearrange the assembly to move them to
 the ends.

 This function goes to some trouble to try to work around some of the
 fragilities in the phageterm implementation.  Thus it splits
 multi-contig assemblies into separate files and runs on them.  In
 addition, it works around the nasty file type detection in phageterm
 by simply decompressing the input files.

=item C<Arguments>

 input(required): Input fastq reads.
 library(required): Assembly created from the input reads.
 cpus(8): Use this number of cpus.
 jmem(12): and this amount of memory.
 jprefix('14'): Output/jobname prefix.
 modules('phageterm'): load the phageterm environment module.

=cut
sub Phageterm {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'library'],
        cpus => 8,
        jmem => 12,
        jprefix => '14',
        modules => ['phageterm'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('PhageTerm.py');
    die('Could not find phageterm in your PATH.') unless($check);

    my $job_name = $class->Get_Job_Name();
    my $inputs = $class->Get_Paths($options->{input});
    my $cwd_name = basename(cwd());
    my $assembly_relative = $options->{library};
    my $assembly_full = abs_path($options->{library});
    my $assembly_name = basename(dirname($options->{library}));
    my $output_dir = qq"outputs/$options->{jprefix}phageterm";
    if ($assembly_name ne '.' || !defined($assembly_name)) {
        $output_dir .= qq"_${assembly_name}";
    }
    if (-d $output_dir) {
        my $removed = rmtree($output_dir);
    }

    my $comment = qq!## This is a script to run phageterm.
## Phageterm has some peculiarities which require one to be
## extra careful when running it.
!;
    my $test_file = qq"${output_dir}/${cwd_name}_direct-term-repeats.fasta";
    my $cwd_test_file = basename($test_file);
    my $output_file = qq"${output_dir}/phageterm_final_assembly.fasta";
    my $dtr_type_file = qq"${output_dir}/phageterm_final_nrt.txt";
    my $dtr_sequence_file = qq"${output_dir}/phageterm_final_dtr.fasta";
    my $dtr_test_file = qq"${output_dir}/phageterm_final_nodtr.txt";
    my $jstring = qq?
use Bio::Adventure::Phage;
my \$result = \$h->Bio::Adventure::Phage::Phageterm_Worker(
  input => '$options->{input}',
  library => '$options->{library}',
  jprefix => '$options->{jprefix}',
  jname => 'phageterm_${job_name}',
  output_dir => '${output_dir}',
  output => '${output_file}',);
?;
    my $phageterm_run = $class->Submit(
        input => $options->{input},
        output => $output_file,
        output_type => $dtr_type_file,
        output_dtr => $dtr_sequence_file,
        test_file => $dtr_test_file,
        language => 'perl',
        comment => $comment,
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
        jname => 'phageterm',
        jprefix => $options->{jprefix},
        jstring => $jstring,);
    return($phageterm_run);
}

=head2 C<Phageterm_Worker>

 Does the actual work of the phageterm function.

 This handles the logic of splitting the input assembly, decompressing
 the input files, and looking over the phageterm results to try to
 make sense of them.

=cut
sub Phageterm_Worker {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'library'],
        output_dir => '.',
        output => 'final_sequence.fasta',
        jprefix => '19',
        modules => ['phageterm']);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    ## I am going to use this function to make it easier for phageterm
    ## to run without shenanigans.  In order to accomplish this, I will
    ## make a separate run for each contig of the assembly, run it on them
    ## individually, interpret the results, and make a final assembly.
    my $workdir = $options->{output_dir};
    if (-d $workdir) {
        qx"rm -f ${workdir}/*";
    } else {
        my $made = make_path($workdir);
    }
    my $phage_log = FileHandle->new(">${workdir}/phageterm_worker.log");
    my $read_string = '';
    ## Step1: Decompress the reads so phageterm doesn't cry.
    if ($options->{input} =~ /\:|\;|\,|\s+/) {
        my @in = split(/\:|\;|\,|\s+/, $options->{input});
        ## I think abs_path only works on things which already exist.
        ## Which is a problem if we are running this in the middle of a pipeline
        my $in_dir = dirname($in[0]);
        my $prefix = abs_path($in_dir);
        my $r1_dirname = dirname($in[0]);
        my $r2_dirname = dirname($in[1]);
        my $r1_filename = basename($in[0]);
        my $r2_filename = basename($in[1]);
        my $r1_filebase = basename($r1_filename, ('.xz', '.gz', '.bz2'));
        my $r2_filebase = basename($r2_filename, ('.xz', '.gz', '.bz2'));
        my $decompressed_r1 = qx"less $in[0] > ${workdir}/r1.fastq";
        my $decompressed_r2 = qx"less $in[1] > ${workdir}/r2.fastq";
        print $phage_log "Copied reads to ${workdir} as r1.fastq and r2.fastq.\n";
        $read_string = 'r1.fastq r2.fastq';
    } else {
        my $r1_filename = basename($options->{input});
        my $r1_dirname = dirname($options->{input});
        my $r1 = abs_path($options->{input});
        my $r1_filebase = basename($r1_filename, ('.xz', '.gz', '.bz2'));
        my $decompressed = qx"less $options->{input} > ${workdir}/r1.fastq";
        print $phage_log "Copied reads to ${workdir} as r1.fastq.\n";
        $read_string = 'r1.fastq';
    }
    ## At this point we have the reads, so grab the assembly and split it up into individual sequences/file.
    my $in_assembly = Bio::SeqIO->new(-file => qq"<$options->{library}",
                                      -format => 'Fasta');
    print $phage_log "Separating contigs from assembly.\n";
    my $contig_number = 0;
  CONTIGS: while (my $contig = $in_assembly->next_seq()) {
      $contig_number++;
      my $id = $contig->id;
      my $contig_name = qq"contig${contig_number}";
      my $out_cwd_file = qq"${contig_name}.fasta";
      my $out_file = qq"${workdir}/${out_cwd_file}";
      my $out_contig = Bio::SeqIO->new(-file => ">${out_file}",
                                       -format => 'Fasta');
      print $phage_log "Extracted contig ${contig_number} as ${out_file}.\n";
      $out_contig->write_seq($contig);
      ## Now that we have the individual contigs separated, run phageterm on them.
      print $phage_log "Running phageterm on the individual contig.\n";
      my $job_string = qq"cd ${workdir} && \\
PhageTerm.py -f ${read_string} \\
  -r ${out_cwd_file} --nrt --report_title ${contig_name} \\
  2>${contig_name}.stderr 1>${contig_name}.stdout && \\
  mv nrt.txt ${contig_name}_nrt.txt";
      my $running_job = qx"${job_string}";
  }
    ## Now that we have run phageterm on each contig
    ## Go back over them and find which (if any) have
    ## a DTR and put the first one on the front of the assembly.
    print $phage_log "Looking through phageterm results for the first DTR.\n";
    my $dtr_contig = '0';
  DTR_HUNT: for my $num (1 .. $contig_number) {
      ## There are two conditions required for a good DTR:
      ## 1.  a direct-term-repeats.fasta file
      ## 2.  a rewritten genome file, which on occasion is unfortunately empty.
      my $test_dtr = qq"${workdir}/contig${num}_direct-term-repeats.fasta";
      my $test_genome = qq"${workdir}/contig${num}_sequence.fasta";
      my $test_lines = qx"wc -l ${test_genome} | awk '{print \$1}'";
      ## If we find the dtr file, take note of the contig and break out.
      print $phage_log "Looking for ${test_dtr} and counting lines of ${test_genome} which is ${test_lines}\n";
      if (-r $test_genome and $test_lines > 2) {
          print $phage_log "Using contig $num as the DTR containing contig.\n";
          $dtr_contig = $num;
          last DTR_HUNT;
      }
  } ## End the DTR hunt

    ## Last step, if we found a dtr, rewrite the assembly with it as the first contig
    ## But first, handle the case where there is no DTR.
    if ($dtr_contig eq '0') {
        print $phage_log "No DTR was observed, copying the input to the phageterm output file.\n";
        my $copied = qx"cp $options->{library} ${workdir}/phageterm_final_assembly.fasta && \\
  echo 'no dtr' > ${workdir}/phageterm_final_nodtr.txt";
  } else {
      ## Otherwise, put the DTR as the first contig and copy the rest
      print $phage_log "A DTR was found, writing to the phageterm output file.\n";
      my $out_assembly = Bio::SeqIO->new(-file => qq">${workdir}/phageterm_final_assembly.fasta",
                                         -format => 'Fasta');
      print $phage_log "Writing the first output assembly entry: ${workdir}/contig${dtr_contig}_sequence.fasta.\n";
      my $in_assembly = Bio::SeqIO->new(-file => qq"<${workdir}/contig${dtr_contig}_sequence.fasta",
                                        -format => 'Fasta');
      while (my $sequence = $in_assembly->next_seq()) {
          $out_assembly->write_seq($sequence);
      }
      ## For the rest of the contigs, copy the unmodified sequences to the final assembly.
      print $phage_log "Writing all other contigs to the final assembly.\n";
    CONTIG_LOOP: for my $num (1 .. $contig_number) {
        next CONTIG_LOOP if ($num eq $dtr_contig);
        print $phage_log "Writing contig ${num} to the final assembly.\n";
        my $last_contig = Bio::SeqIO->new(-file => qq"<${workdir}/contig${num}.fasta",
                                          -format => 'Fasta');
        while (my $seq = $last_contig->next_seq()) {
            $out_assembly->write_seq($seq);
        }
    } ## End the contig loop
      ## Now lets copy the DTR and nrt files
      print $phage_log "Copying the dtr from contig ${dtr_contig} to phageterm_final_dtr.fasta.\n";
      my $dtr_copy = qx"cp ${workdir}/contig${dtr_contig}_direct-term-repeats.fasta ${workdir}/phageterm_final_dtr.fasta";
      print $phage_log "Copying the NRT file from ${dtr_contig} to phageterm_final_nrt.txt.\n";
      my $ntr_copy = qx"cp ${workdir}/contig${dtr_contig}_nrt.txt ${workdir}/phageterm_final_nrt.txt";
  } ## End when we do find a DTR
    print $phage_log "Deleting the uncompressed input files.\n";
    my $delete_crap = qx"rm -f ${workdir}/r1.fastq ${workdir}/r2.fastq ${workdir}/contig*.fasta 2>/dev/null 1>/dev/null";
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
}

=head2 C<Phastaf>

 Invoke Torsten Seeman's phastaf tool to seek out phage-derived sequences.
 https://github.com/tseemann/phastaf

 Phastaf is (I assume) a joke on phaster, and a much more straight
 forward way to hunt for phage-derived sequences in an input assembly.
 It basically just runs tblastx on the input assembly and returns some
 files which describe the locations of high-quality hits.

=item C<Arguments>

 input(required): Input assembly.
 cpus(8): Use this number of cpus.
 jmem(12): And this amount of memory.
 jprefix('14'): with this prefix for the jobname/directory.
 modules('phastaf'): Load the phastaf environment module.

=cut
sub Phastaf {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        cpus => 8,
        jmem => 12,
        jprefix => '14',
        modules => ['phastaf'],);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('phastaf');
    die('Could not find phastaf in your PATH.') unless($check);

    my $job_name = $class->Get_Job_Name();
    my $input_paths = $class->Get_Paths($options->{input});
    my $input_full = $input_paths->{fullpath};
    my $output_dir = qq"outputs/$options->{jprefix}phastaf";
    $output_dir .= "_$input_paths->{dirname}" if (defined($input_paths->{dirname}));
    if (-d $output_dir) {
        my $removed = rmtree($output_dir);
    }
    my $comment = '## This is a script to run phastaf.';
    my $coords = qq"${output_dir}/diamond.coords";
    my $jstring = qq?
mkdir -p ${output_dir}
phastaf --force --outdir ${output_dir} \\
  --cpus $options->{cpus} \\
  $options->{input} \\
  2>${output_dir}/phastaf.stderr \\
  1>${output_dir}/phastaf.stdout
?;
    my $output_file = qq"${output_dir}/something.txt";
    my $phastaf = $class->Submit(
        cpus => $options->{cpus},
        coordinates => $coords,
        output => $output_file,
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => qq"phastaf_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => $options->{jmem},
        modules => $options->{modules},
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($phastaf);
}

=head2 C<Restriction_Catalog>

 Make a catalog of restriction enzyme hits from an input sequence.

 My initial goal for this was to also have it pull the restriction
 enzymes from a host genome and cross reference them.  Sadly, the
 references I checked against had only ~ 1 annotated restriction
 endonuclease -- which makes no damn sense to me.  I guess they are on
 plasmids or something?  In any event, I figured it would still be
 nice to have an idea of the restriction sites.

=item C<Arguments>

 input(required): Input assembly.
 library(host_species.txt): Either a species name or filename. (Currently not used)
 jmem(8): Memory expected.
 jprefix('29'): Job/directory prefix.

=cut
sub Restriction_Catalog {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        library => 'host_species.txt',
        jmem => 8,
        jprefix => '29',);
    my $output_dir = qq"outputs/$options->{jprefix}re_catalog";
    my $re_output = qq"${output_dir}/re_catalog.tsv";
    my $host_output = qq"${output_dir}/re_host_catalog.tsv";
    if (-d $output_dir) {
        my $removed = rmtree($output_dir);
    }
    my $paths = $class->Get_Paths($re_output);
    my $comment = '## Go on a restriction enzyme hunt!';
    my $jstring = qq!
use Bio::Adventure;
use Bio::Adventure::Phage;
my \$result = Bio::Adventure::Phage::Restriction_Catalog_Worker(\$h,
  input => '$options->{input}',
  jname => 're_catalog',
  jprefix => '$options->{jprefix}',
  library => '$options->{library}',
  output => '${re_output}',);
!;
    my $re_job = $class->Submit(
        input => $options->{input},
        comment => $comment,
        jmem => $options->{jmem},
        jname => 're_catalog',
        jprefix => $options->{jprefix},
        jstring => $jstring,
        language => 'perl',
        library => $options->{library},
        modules => $options->{modules},
        output => $re_output,);
    return($re_job);
}

=head2 C<Restriction_Catalog_Worker>

 This does the actual work for Restriction_Catalog().

 Bio::Restriction makes this really easy.

=cut
sub Restriction_Catalog_Worker {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'library'],
        jprefix => '29',
        output => 're_catalog.tsv',);
    my $output_dir = dirname($options->{output});
    my $log = FileHandle->new(">${output_dir}/re_catalog.log");

    my $collection = Bio::Restriction::EnzymeCollection->new();

    my $query = Bio::SeqIO->new(-file => $options->{input}, -format => 'Fasta');
    my %data = ();
  SEARCHES: while (my $contig = $query->next_seq()) {
      my $sequence = $contig->seq;
      my $analysis = Bio::Restriction::Analysis->new(-seq => $contig);
      my $internal = {};
      for my $en ($collection->each_enzyme()) {
          my $name = $en->name;
          my $site = $en->site;
          if (!defined($site)) {
              $site = '';
          }
          my $overhang = $en->overhang_seq;
          if (!defined($overhang)) {
              $overhang = '';
          }
          my $cuts = $analysis->cuts_by_enzyme($name);
          if (!defined($cuts)) {
              $cuts = 0;
          }
          my @fragments = $analysis->fragments($en);
          $internal = {
              site => $site,
              overhang => $overhang,
              cuts => $cuts,
              fragments => \@fragments,
          };
          if ($data{$name}) {
              $internal->{cuts} = $internal->{cuts} + $data{$name}->{cuts};
              push(@fragments, $internal->{fragments});
              $internal->{fragments} = \@fragments;
          }
          $data{$name} = $internal;
      }
  }

    my $re_tsv = FileHandle->new(">$options->{output}");
    print $re_tsv qq"RE\tSite\tOverhang\tCuts\n";
    for my $name (sort keys %data) {
        print $re_tsv qq"${name}\t$data{$name}->{site}\t$data{$name}->{overhang}\t$data{$name}->{cuts}\n";
    }
    $re_tsv->close();
}

=head2 C<Terminase_Reorder>

 Reorder an assembly so that a terminase gene is at the beginning.
 Terminase library: 10.1186/s12864-021-08029-8

 This will take an arbitrary assembly file, invoke prodigal on it to get ORFs,
 pass them to Reorder_Search (above), and copy the assembly to a new file
 with the most likely terminase sequence at the beginning.  It is
 worth noting that I took the terminase library directly out of MAS.

 We should potentially consider switching this from prodigal to phanotate.

=item C<Arguments>

 input(required): Input assembly.
 evalue(0.01): Cutoff evalue.
 fasta_tool(fastx36): Use this aligner.
 gcode(11): Assume the bacterial genetic code.
 library(terminase): Name of the library to search against.
 species('phage'): Used for prodigal training.
 test_file('direct-term-repeats.fasta'): File to test whether to run
  against (I think this is no longer needed).
 jmem(8): Expected memory usage.
 jprefix('15'): Job/directory prefix.
 modules('fasta', 'blast', 'blastdb'): Environment modules to load.

=cut
sub Terminase_ORF_Reorder {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        evalue => 0.01,
        fasta_tool => 'fastx36',
        gcode => '11',
        library => 'terminase',
        species => 'phages',
        test_file => 'direct-term-repeats.fasta',
        jmem => 8,
        jprefix => '15',
        modules => ['fasta', 'blast', 'blastdb'],);
    my $input_dir = basename(dirname($options->{input}));
    my $output_dir = qq"outputs/$options->{jprefix}termreorder_${input_dir}";
    my $final_output = qq"${output_dir}/final_assembly.fasta";
    if (-d $output_dir) {
        my $removed = rmtree($output_dir);
    }
    my $paths = $class->Get_Paths($final_output);

    my $prodigal_outname = 'prodigal';
    my $prodigal_cds = qq"${output_dir}/${prodigal_outname}_cds.fasta";

    my $term_prodigal = $class->Bio::Adventure::Feature_Prediction::Prodigal(
        gcode => $options->{gcode},
        input => $options->{input},
        jdepends => $options->{jdepends},
        jprefix => $options->{jprefix},
        output_dir => $output_dir,
        prodigal_outname => $prodigal_outname,
        species => $options->{species},);
    ## Once that is finished running, we should have files:
    ## 'translated.fasta' and 'cds.fasta' which we can use to go terminase hunting.
    ## I am duplicating a bunch of options here, that might be dumb.
    $options->{jdepends} = $term_prodigal->{job_id};
    my $comment = qq"## This should use the predicted prodigal ORFs to search against a
## local terminase sequence database.  Then for each contig, move the best
## terminase hit to the front of the sequence.\n";
    my $output_tsv = qq"${output_dir}/$options->{library}_summary.tsv";
    my $jstring = qq!
use Bio::Adventure;
use Bio::Adventure::Phage;
my \$result = Bio::Adventure::Phage::Terminase_ORF_Reorder_Worker(\$h,
  evalue => '$options->{evalue}',
  fasta_tool => '$options->{fasta_tool}',
  input => '$options->{input}',
  jdepends => '$options->{jdepends}',
  jname => 'terminase_reorder',
  jprefix => '$options->{jprefix}',
  library => '$options->{library}',
  query => '${prodigal_cds}',
  output => '${final_output}',
  output_dir => '${output_dir}',
  test_file => '$options->{test_file}',);
!;
    my $tjob = $class->Submit(
        evalue => $options->{evalue},
        fasta_tool => 'fastx36',
        input => $options->{input},
        library => $options->{library},
        query => $prodigal_cds,
        output => $final_output,
        output_dir => $output_dir,
        output_tsv => $output_tsv,
        test_file => $options->{test_file},
        language => 'perl',
        comment => $comment,
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
        jname => 'terminase_reorder',
        jprefix => $options->{jprefix},
        jstring => $jstring,
        modules => $options->{modules},);
    $tjob->{prodigal_job} = $term_prodigal;
    return($tjob);
}

=head2 C<Terminase_ORF_Reorder_Worker>

 Reorder an assembly from the result of a blast(or presumbly other) search.

 This does the work for the Terminase_ORF_Reorder().  It contains the
 logic for performing a blast/fasta search, parsing the results,
 filtering for high quality hits, and reordering the assembly accordingly.

=cut
sub Terminase_ORF_Reorder_Worker {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'query'],
        evalue => 0.01,
        fasta_tool => 'fastx36',
        jprefix => '15',
        library => 'terminase',
        modules => ['fasta', 'blast', 'blastdb'],
        output_file => 'final_assembly.fasta',
        test_file => '',);

    my $output_dir = dirname($options->{output});
    my $log = FileHandle->new(">${output_dir}/reorder.log");
    ## First check for the test file, if it exists, then this should
    ## just symlink the input to the output until I think of a smarter
    ## way of thinking about this.
    my $new_filename;
    if ($options->{output_file}) {
        ## Just in case a path is provided
        my $outfile = basename($options->{output_file});
        $new_filename = qq"${output_dir}/${outfile}";
        print "options->output_file exists: $new_filename\n";
    } else {
        $new_filename = basename($options->{input}, ('.fasta'));
        $new_filename = qq"${output_dir}/${new_filename}_reordered.fasta";
    }

    ## Load up fasta36 and make sure we are able to run a search against the
    ## translated nucleotide sequences.
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('fastx36');
    die("Could not find fasta36 in your PATH.") unless($check);

    ## Now perform the search for potential terminases.
    my $library_file = qq"$ENV{BLASTDB}/$options->{library}.fasta";
    my $fasta_output = Bio::SearchIO->new(-format => 'fasta', );
    my $number_hits = 0;
    print $log "Starting fasta search of $options->{input}
  against ${library_file} using tool: $options->{fasta_tool}.\n";
    my $query = Bio::SeqIO->new(-file => $options->{input}, -format => 'Fasta');
    my $fasta_outfile = qq"${output_dir}/$options->{library}_hits.txt";
    my @params = (
        b => 1,
        O => $fasta_outfile,
        program => $options->{fasta_tool},);
    my $seq_count = 0;
    my $e = '';
    my $search = Bio::Tools::Run::Alignment::StandAloneFasta->new(@params);
    $search->library($library_file);
    my @fasta_output;
    my ($stdout, $stderr, @returns)  = capture {
        try {
            @fasta_output = $search->run($options->{query});
        } catch ($e) {
            warn "An error occurred: $e";
        }
    };

    ## Collect the results from the fasta search and write them out.
    my $hit_fh = FileHandle->new(">${output_dir}/$options->{library}_summary.tsv");
    print $hit_fh "Name\tLength\tAccession\tDescription\tScore\tSignificance\tBit\tHitStrand\tQueryStrand\n";
    my $result_data = {};
    my $search_output = Bio::SearchIO->new(-file => $fasta_outfile, -format => 'fasta');
    my $element_count = 0;
    my $result_count = 0;
  RESULTS: while (my $result = $search_output->next_result) {
      $result_count++;
      my $query_name = $result->query_name();
      my $query_length = $result->query_length();
      my $query_descr = $result->query_description();
      my $stats = $result->available_statistics();
      my $hits = $result->num_hits();

      $result_data->{$query_name} = {
          description => $query_descr,
          name => $query_name,
          stats => $stats,
          num_hits => $hits,
          length => $query_length,
          hit_data => [],
      };
      my $hit_count = 0;
    HITLOOP: while (my $hits = $result->next_hit()) {
        $number_hits = $number_hits++;
        my $hit_name = $hits->name();
        my $hit_length = $hits->length();
        my $hit_strand = $hits->strand('hit');
        my $query_strand = $hits->strand('query');
        my $hit_acc = $hits->accession();
        my $hit_descr = $hits->description();
        my $hit_score = $hits->score();
        my $hit_sig = $hits->significance();
        my $hit_bits = $hits->bits();
        my @hit_data = @{$result_data->{$query_name}->{hit_data}};
        my $hit_datum = {
            name => $hit_name,
            length => $hit_length,
            acc => $hit_acc,
            description => $hit_descr,
            score => $hit_score,
            sig => $hit_sig,
            bit => $hit_bits,
            hit_strand => $hit_strand,
            query_strand => $query_strand,
        };
        print $hit_fh "${hit_name}\t${hit_length}\t${hit_acc}\t${hit_descr}\t${hit_score}\t${hit_sig}\t${hit_bits}\t${hit_strand}\t${query_strand}\n";
        push(@hit_data, $hit_datum);
        $result_data->{$query_name}->{hit_data} = \@hit_data;
        $hit_count++;
    } ## End the hitloop.
  } ## End of the individual fasta search  (maybe not needed?)
    $hit_fh->close();

    ## We wrote the hits, now check to see if we actually need to reorder the assembly.
    ## If so, continue, if not, exit out here.
    print $log "Testing for $options->{test_file}\n";
    if (-r $options->{test_file}) {
        print $log "The test file exists, reorder the genome to the terminase.\n";
    } else {
        my $full_input = abs_path($options->{input});
        my $pwd = getcwd();
        my $full_new = qq"${pwd}/${new_filename}";
        print $log "The test file does not exist, do not reorder the genome to the terminase.
Symlinking ${full_input} to ${full_new} and stopping.\n";
        if (-f $full_new or -l $full_new) {
            my $unlink_lst = ($full_new);
            my $removed = unlink($unlink_lst);
        }
        my $linked = symlink($full_input, $full_new);
        if (-f $full_new) {
            print $log "The symlink succeeded, the file ${full_new} should exist.\n";
        } else {
            print $log "Something went wrong with ln -s ${full_input} ${full_new}\n";
            my $linked = qx"ln -s ${full_input} ${full_new}";
        }
        return($full_new);
    } ## End checking for the test file.

    ## At this point we should have a data structure
    ## $result_data, with primary keys as the ORF IDs
    ## a few key-value pairs for the individual putative ORFs, then
    ## the key 'hit_data' which is an array of hits containing the scores.
    ## So we should now hunt down the best score
    ## The following loop will therefore go over the set of hits and look for the lowest E-value.
    my $best_score = 1;
    my $best_query = '';
    my $best_description = '';
    my $total_hits = 0;
    my $best_query_strand;
    my $best_hit_strand;
  SCANNER: foreach my $query_id (keys %{$result_data}) {
      my $query_description = $result_data->{$query_id}->{description};
      my @hit_arr = @{$result_data->{$query_id}->{hit_data}};
      my $hits = scalar(@hit_arr);
      $total_hits = $total_hits + $hits;
      if ($hits < 1) {
          next SCANNER;
      }
      for my $hit (@hit_arr) {
          if ($hit->{sig} < $best_score) {
              $best_query = $query_id;
              $best_query_strand = $hit->{query_strand};
              $best_hit_strand = $hit->{hit_strand};
              $best_description = $query_description;
              $best_score = $hit->{sig};
              my $hit_strand = $hit->{strand};
              print $log "Found a new best hit: ${query_id} with evalue $hit->{sig}\n";
              print $log "  This is on query strand: ${best_query_strand} and hit strand: ${best_hit_strand}\n";

          } elsif ($hit->{sig} == $best_score) {
              print $log "Found an equivalent hit: ${query_id} with evalue $hit->{sig}.\n";
          }
      }
  } ## Finished iterating over every potential result.
    print $log "Out of ${total_hits} hits, ${best_query} was chosen with an e-value of: ${best_score}.\n";
    print $log "Best terminase description: ${best_description}\n";
    ## print "Best terminase description: ${best_description}\n";
    ## Now things get a bit confusing: Keep in mind that prodigal's ORFs are printed as:
    ## ${contig}_${orf} # start # end # strand # stuff, so the first thing to do is pull out the best contig.
    my $best_contig = 'failed';
    my $best_ORF = 'failed';
    if ($best_query =~ /^(.*)?_(\d+)$/) {
        $best_contig = $1;
        $best_ORF = $2;
    }

    ## Now lets take a moment to set up the input and output sequence files.
    my $in_assembly = Bio::SeqIO->new(-file => $options->{input} ,
                                      -format => 'Fasta');
    print $log "Writing a new assembly to ${new_filename}\n";
    my $out_assembly = Bio::SeqIO->new(-file => ">${new_filename}" ,
                                       -format => 'Fasta');

    ## When reordering the sequence, if it is on the + strand we will take:
    ## S        T                    E
    ## Start    Terminase            End
    ## >>>>>>>>>|||||||||>>>>>>>>>>>>>>>
    ## From the T of terminase to end as the first subsequence
    ## Then from the S of start to T-1 as the second.

    ## Conversely, if it is on the - strand:
    ## S        T       3            E D
    ## Start    Terminas3            End
    ## <<<<<<<<<|||||||||<<<<<<<<<<<<<<<
    ## From the 3' of terminase to the S of start (revcomp'd) as the first subsequence
    ## Then from the D of end to 3+1 as the second.
    my @other_objects = ();
    my %assembly_lengths = ();
    my $first_object;
    while (my $seqobj = $in_assembly->next_seq()) {
        my $current_id = $seqobj->id();
        my $current_length = $seqobj->length();
        if ($current_id eq $best_contig) {
            $first_object = $seqobj;
        } else {
            push(@other_objects, $seqobj);
        }
        $assembly_lengths{$current_id} = $current_length;
    } ## End iterating over the sequences.

    ## The last thing to remember is that prodigal encodes its position information
    ## as a set of # thing # thing # things, so lets grab out the information of interest.
    my @start_end_array = split(/\s*#\s*/, $best_description);
    my $best_seq_start = $start_end_array[1];
    $best_seq_start =~ s/\s*//g;
    my $best_seq_end = $start_end_array[2];
    $best_seq_end =~ s/\s*//g;
    my $best_seq_strand = $start_end_array[3];
    $best_seq_strand =~ s/\s*//g;

    ## No matter what, we will need the sequence length:
    my $sequence_end = $first_object->length;
    ## Set the coordinates for the two pieces
    my ($first_start, $first_end, $second_start, $second_end);
    ## And the resulting sequences.
    my ($first_seq, $second_seq);
    if ($best_seq_strand == 1) {
        $first_start = $best_seq_start;
        $first_end = $sequence_end;
        $second_start = 1;
        $second_end = $best_seq_start - 1;
        $first_seq = $first_object->subseq($first_start, $first_end);
        $second_seq = $first_object->subseq($second_start, $second_end);
        print $log "Merging two + strand pieces:
  ${first_start}:${first_end} and ${second_start}:${second_end}.\n";
    } else {
        $first_start = 1;
        $first_end = $best_seq_end;
        $second_start = $best_seq_end + 1;
        ## $second_start = $best_seq_end;
        $second_end = $sequence_end;
        ## yeah yeah I know one is supposed to use revcomp, but damn
        ## it annoys me when it yells at me for not using a sequence object.
        $first_seq = $first_object->subseq($first_start, $first_end);
        print $log "Merging two - strand pieces:
  ${first_start}:${first_end}";

        $first_seq = reverse($first_seq);
        $first_seq =~ tr/ATGCUatgcu/TACGAtacga/;
        if ($second_start < $second_end) {
            $second_seq = $first_object->subseq($second_start, $second_end);
            $second_seq = reverse($second_seq);
            $second_seq =~ tr/ATGCUatgcu/TACGAtacga/;
            print $log " and ${second_start}:${second_end}";
        } else {
            print $log " and nothing, something is weird with the second sequence.";
        }
        print $log "\n";
    }
    my $final_sequence;
    if (defined($second_seq)) {
        $final_sequence = $first_seq . $second_seq;
    } else {
        $final_sequence = $first_seq;
    }

    ## Now write out the sequence, starting with 'first_seq', then iterate through everyone else
    my $input_id = $first_object->id();
    ## First write the reordered contig
    my $output_seq = Bio::Seq->new(-seq => $final_sequence,
                                   -id => $input_id);
    $out_assembly->write_seq($output_seq);

    ## The write all of the others.
    for my $obj (@other_objects) {
        $out_assembly->write_seq($obj);
    }
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    $log->close();
    return($result_data);
}

1;

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

=cut
