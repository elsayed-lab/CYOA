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

=head2 C<Get_DTR>

Extract the direct-terminal-repeats from a phageterm run.

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

=head2 C<Filter_Host_Kraken

=cut
sub Filter_Host_Kraken {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'input_fastq'],
        jdepends => '',
        jmem => 8,
        jprefix => '06',);
    my $comment = qq"## Use kraken results to choose a host species to filter against.\n";
    my $output_dir = qq"outputs/$options->{jprefix}filter_kraken_host";
    make_path($output_dir);
    my $out_r1_name = 'r1_host_filtered.fastq.xz';
    my $out_r2_name = 'r2_host_filtered.fastq.xz';
    my $output_files = qq"${output_dir}/${out_r1_name}:${output_dir}/${out_r2_name}";
    my $jstring = qq!
use Bio::Adventure;
use Bio::Adventure::Phage;
my \$result = Bio::Adventure::Phage::Get_Kraken_Host(\$h,
  comment => '${comment}',
  output => '${output_files}',
  output_dir => '${output_dir}',
  input => '$options->{input}',
  input_fastq => '$options->{input_fastq}',
  jdepends => '$options->{jdepends}',
  jname => 'kraken_host',
  jprefix => '$options->{jprefix}',);
!;
    my $host = $class->Submit(
        jdepends => $options->{jdepends},
        input => $options->{input},
        input_fastq => $options->{input_fastq},
        output => $output_files,
        output_dir => $output_dir,
        jmem => $options->{jmem},
        jname => 'hostfilter',
        jprefix => $options->{jprefix},
        jstring => $jstring,
        language => 'perl',
        shell => '/usr/bin/env perl',);
    return($host);
}

=head2 C<Get_Kraken_Host>

Read the report from kraken and hunt down the species name with the most reads.

=cut
sub Get_Kraken_Host {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'output'],
        jname => 'krakenfilter',
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
    my $out = FileHandle->new(">${output_dir}/kraken_filter.log");
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
    my $escaped_species = $most_species;
    $escaped_species =~ s/\s/\+/g;

    my $host_species;
    if (-f 'host_species.txt') {
        my $host_file = FileHandle->new("<host_species.txt");
        while (my $line = <$host_file>) {
            chomp $line;
            $host_species = $line;
        }
        $host_file->close();
        print $out "The host_species.txt file already exists and suggests using ${host_species}.\n";
    }

    if (!defined($host_species)) {
        my $search_url = qq"https://www.ncbi.nlm.nih.gov/assembly/?term=${escaped_species}";
        my $mech = WWW::Mechanize->new(autocheck => 1);
        print $out "Searching ${search_url} for appropriate download links.\n";
        my $search_data = $mech->get($search_url);
        my @search_links = $mech->find_all_links(
            tag => 'a', text_regex => qr/ASM/i);
        my $first_hit = $search_links[0];
        my ($assembly_link, $assembly_title) = @{$first_hit};
        my $accession = basename($assembly_link);
        my $downloaded_file = qq"$options->{libdir}/$options->{libtype}/${accession}.gbff.gz";
        if (-r $downloaded_file) {
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
        $host_species = $accession;

        ## Convert the assembly to fasta/gff/etc.
        print $out "Converting ${downloaded_file} assembly to fasta/gff.\n";
        my $converted = $class->Bio::Adventure::Convert::Gb2Gff(input => $downloaded_file);
    } ## End checking if the host_species was defined.

    my $cyoa_shell = Bio::Adventure->new(cluster => 0);
    my $index_input = qq"$options->{libdir}/$options->{libtype}/${host_species}.fasta";
    my $index_location = qq"$options->{libdir}/$options->{libtype}/indexes/${host_species}.1.ht2";
    if (! -f $index_location) {
        my $indexed = $cyoa_shell->Bio::Adventure::Index::Hisat2_Index(input => $index_input);
        ## Check to see if there is a text file containing the putative host species
        ## If it does not exist, write it with the species name.
    }
    ## Now perform the filter using our temp_cyoa.
    my $jprefix = $options->{jprefix} + 1;
    print $out "Filtering out reads which map to ${host_species}.\n";
    my $filter = $cyoa_shell->Bio::Adventure::Map::Hisat2(
        input => $options->{input_fastq},
        jprefix => $jprefix,
        do_htseq => 0,
        species => $host_species,);
    my $filtered_reads = $filter->{unaligned_comp};
    my ($in_r1, $in_r2) = split(/\:|\;|\,|\s+/, $filtered_reads);
    $in_r1 = File::Spec->rel2abs($in_r1);
    $in_r2 = File::Spec->rel2abs($in_r2);
    my ($out_r1, $out_r2) = split(/\:|\;|\,|\s+/, $options->{output});
    print $out "Symlinking final output files to $options->{output_dir}\n";
    if (! -f $out_r1) {
        my $s1 = symlink($in_r1, $out_r1);
    }
    if (! -f $out_r2) {
        my $s2 = symlink($in_r2, $out_r2);
    }
    return(%species_observed);
}

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

=head2 C<Classify_Phage>

Use tblastx against a database extracted from reference accessions
provided by the ICTV.  This should hopefully provide a reasonably
complete set of taxonomic classifications for viral assemblies.  This
function pretty much just calls Blast_Classify().

=cut
sub Classify_Phage {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        evalue => 0.01,
        blast_tool => 'tblastx',
        jmem => 12,
        jprefix => '18',
        library => 'ictv',
        modules => ['blastdb', 'blast'],
        topn => 5,);

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
my \$result = Bio::Adventure::Phage::Blast_Classify(\$h,
  comment => '${comment}',
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
        jdepends => $options->{jdepends},
        blast_tool => $options->{blast_tool},
        evalue => $options->{evalue},
        input => $options->{input},
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
        topn => $options->{topn},
        shell => '/usr/bin/env perl',);
    return($cjob);
}

=head2 C<Blast_Classify>

Does the actual work of attempting to classify a viral sequence
against the ICTV database.

=cut
sub Blast_Classify {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input',],
        evalue => 0.01,
        blast_tool => 'tblastx',
        jprefix => '18',
        jcpus => 6,
        library => 'ictv',
        modules => ['blast', 'blastdb'],
        output_log => 'classify.log',
        output_blast => 'ictv_hits.txt',
        output_dir => '.',
        output => 'ictv_filtered.tsv',
        score => 1000,
        topn => 5,);
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
    my $comment = qq!## This invokes blast and pulls the top entries in an attempt to classify a viral genome.
!;
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

    ##my ($stdout, $stderr, @returns)  = capture {
    ##    try {
    ##        @blast_output = $search->run($query);
    ##    } catch ($e) {
    ##        warn "An error occurred: $e";
    ##    }
    ##};

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
    my $index = qq"$options->{libdir}/codon_tables/$options->{species}.txt";
    if (!-r $index) {
        my $wrote_index = $class->Bio::Adventure::Phage::Make_Codon_Table(
            species => $options->{species});
    }
    my $test = ref($options->{suffixes});
    my @suffixes = split(/,/, $options->{suffixes});
    my $out_base = basename($options->{input}, @suffixes);
    my $jname = qq"${out_base}_vs_$options->{species}";
    my $output_dir = qq"outputs/$options->{jprefix}${jname}";
    make_path($output_dir);
    my $output_cai = qq"${output_dir}/${out_base}_cai.txt";
    my $random_sequences = qq"${output_dir}/${out_base}_random_sequences.txt";
    my $expected_cai = qq"${output_dir}/${out_base}_expected.txt";
    my $stderr = qq"${output_dir}/${out_base}.stderr";
    my $stdout = qq"${output_dir}/${out_base}.stdout";
    my $jstring = qq?mkdir -p ${output_dir}
caical -g 11 \\
  -f $options->{input} \\
  -h ${index} \\
  -o1 ${output_cai} \\
  -o2 ${random_sequences} \\
  -o3 ${expected_cai} \\
  2>${stderr} 1>${stdout}
?;

    my $job = $class->Submit(
        output => $output_cai,
        output_random => $random_sequences,
        output_expected => $expected_cai,
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
        jname => $jname,
        jprefix => $options->{jprefix},
        jstring => $jstring,);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload',);
    return($job);
}

=head2 C<Make_Codon_Table>

Given a reference assembly, print out a codon table for use with a codon adapatation calculator.
I wrote a little piece of code which works from a fasta/gff pair, but currently it is dumb.

=cut
sub Make_Codon_Table {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['species'],
        jprefix => '80',);
    my $out_table = qq"$options->{libdir}/codon_tables/$options->{species}.txt";

    ## I have a few suffixes for writing genbank files.
    my @potential_suffixes = ('gbff', 'gbk', 'gbf', 'gb', 'genbank');
    my @compressions = ('gz', 'xz', 'bz2');
    my $in_gbff = '';
  POTENTIAL: for my $potential (@potential_suffixes) {
      my $start = qq"$options->{libdir}/$options->{libtype}/$options->{species}.${potential}";
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

=head2 C<Phageterm>

Use phageterm on a viral assembly to look for terminal repeat
regions.  If they are found, rearrange the assembly to move them to
the ends.

=cut
sub Phageterm {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'library'],
        cpus => '8',
        modules => ['phageterm'],
        jmem => 12,
        jprefix => '14',);
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
  jname => 'phageterm_${job_name}',
  output_dir => '${output_dir}',
  output => '${output_file}',);
?;
    my $phageterm_run = $class->Submit(
        comment => $comment,
        input => $options->{input},
        output => $output_file,
        output_type => $dtr_type_file,
        output_dtr => $dtr_sequence_file,
        test_file => $dtr_test_file,
        shell => 'usr/bin/env perl',
        language => 'perl',
        jdepends => $options->{jdepends},
        jmem => $options->{jmem},
        jname => 'phageterm',
        jstring => $jstring,);
    $class->{language} = 'bash';
    $class->{shell} = '/usr/bin/env bash';
    return($phageterm_run);
}

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

sub Phageterm_old {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input', 'library'],
        cpus => '8',
        modules => ['phageterm'],
        jmem => 18,
        jprefix => '14',);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('PhageTerm.py');
    die('Could not find phageterm in your PATH.') unless($check);

    my $job_name = $class->Get_Job_Name();
    my $inputs = $class->Get_Paths($options->{input});
    my $cwd_name = basename(cwd());
    my $assembly_relative = $options->{library};
    my $assembly_full = abs_path($options->{library});
    my $assembly_name = basename(dirname($options->{library}));
    my $output_dir = qq"outputs/$options->{jprefix}phageterm_${assembly_name}";
    if (-d $output_dir) {
        my $removed = rmtree($output_dir);
    }

    my $uncompress_string = qq"";
    my $input_string = qq"";
    my $delete_string = qq"";
    if ($options->{input} =~ /\:|\;|\,|\s+/) {
        my @in = split(/\:|\;|\,|\s+/, $options->{input});
        ## I think abs_path only works on things which already exist.
        ## Which is a problem if we are running this in the middle of a pipeline
        my $in_dir = dirname($in[0]);
        make_path($in_dir);
        my $prefix = abs_path($in_dir);
        my $r1_dirname = dirname($in[0]);
        my $r2_dirname = dirname($in[1]);
        my $r1_filename = basename($in[0]);
        my $r2_filename = basename($in[1]);
        $uncompress_string = qq!
less \${start}/${r1_dirname}/${r1_filename} > r1.fastq && \\
  less \${start}/${r2_dirname}/${r2_filename} > r2.fastq
!;
        $input_string = qq! -f r1.fastq -p r2.fastq !;
        $delete_string = qq!rm r1.fastq && rm r2.fastq!;
    } else {
        my $r1_filename = basename($options->{input});
        my $r1_dirname = dirname($options->{input});
        my $r1 = abs_path($options->{input});
        $uncompress_string = qq!
less \${start}/${r1_dirname}/${r1_filename} > r1.fastq
!;
        $input_string = qq! -f r1.fastq !;
        $delete_string = qq!rm r1.fastq!;
    }

    my $comment = qq!## This is a script to run phageterm.
## Phageterm has some peculiarities which require one to be
## extra careful when running it.
!;
    my $test_file = qq"${output_dir}/${cwd_name}_direct-term-repeats.fasta";
    my $cwd_test_file = basename($test_file);
    my $output_file = "${output_dir}/${cwd_name}_sequence.fasta";
    my $jstring = qq?start=\$(pwd)
mkdir -p ${output_dir}
cd ${output_dir}
if [[ -f ${cwd_name}_sequence.fasta ]]; then
  rm -f ${cwd_name}_sequence.fasta
fi
if [[ -f ${cwd_name}_original_sequence.fasta ]]; then
  rm -f ${cwd_name}_original_sequence.fasta
fi
${uncompress_string}

PhageTerm.py ${input_string} \\
  -r \${start}/${assembly_relative} \\
  -c $options->{cpus} \\
  --report_title ${cwd_name} --nrt \\
  2>phageterm.err 1>phageterm.out

${delete_string}
ln -s \${start}/${assembly_relative} ${cwd_name}_original_sequence.fasta
## This is a little hack to make terminase reordering easier.
if [[ -f ${cwd_test_file} ]]; then
  ln -s ${cwd_test_file} direct-terminal-repeats.fasta
else
  ## If phageterm fails, it still writes the output file, but it is useless
  ## so just copy the assembly.
  rm ${cwd_name}_sequence.fasta
  ln -s \${start}/${assembly_relative} ${cwd_name}_sequence.fasta
fi

## I found another bug in phageterm, sometimes it finds a DTR, but leaves the genome blank.
## In this case, I can either use my reordering code to fix it, or just say eff it and
## remove the 'reordered' sequence.
reordered_lines=\$(wc -l ${cwd_name}_sequence.fasta | awk '{print \$1}')
if [[ \${reordered_lines} -eq '2' ]]; then
  rm ${cwd_name}_sequence.fasta
  ln -s \${start}/${assembly_relative} ${cwd_name}_sequence.fasta
fi

cd \${start}
?;

    my $phageterm = $class->Submit(
        cpus => $options->{cpus},
        comment => $comment,
        dtr_file => $test_file,
        jdepends => $options->{jdepends},
        jname => qq"phageterm_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => $options->{jmem},
        modules => $options->{modules},
        output => $output_file,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        test_file => $test_file,);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($phageterm);
}

=head2 C<Phastaf>

Invoke Torsten Seeman's phastaf tool in order to hunt for
phage-derived sequences in an assembly.

=cut
sub Phastaf {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        cpus => '8',
        modules => ['phastaf'],
        jmem => 12,
        jprefix => '14',);
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

    my $comment = qq!## This is a script to run phastaf.
!;
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
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => qq"phastaf_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => $options->{jmem},
        modules => $options->{modules},
        output => $output_file,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($phastaf);
}

=head2 C<Terminase_Reorder>

Reorder an assembly so that a terminase gene is at the beginning.
This will take an arbitrary assembly file, invoke prodigal on it to get ORFs,
pass them to Reorder_Search (above), and copy the assembly to a new file
with the most likely terminase sequence at the beginning.

We should potentially consider switching this from prodigal to phanotate.

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
        jmem => 8,
        jprefix => '15',
        species => 'phages',
        modules => ['fasta', 'blast', 'blastdb'],
        test_file => 'direct-term-repeats.fasta',);

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
    my $jstring = qq!
use Bio::Adventure;
use Bio::Adventure::Phage;
my \$result = Bio::Adventure::Phage::Terminase_ORF_Reorder_Worker(\$h,
  comment => '$comment',
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
        jdepends => $options->{jdepends},
        evalue => $options->{evalue},
        fasta_tool => 'fastx36',
        input => $options->{input},
        jmem => $options->{jmem},
        jname => 'terminase_reorder',
        jprefix => $options->{jprefix},
        jstring => $jstring,
        language => 'perl',
        library => $options->{library},
        modules => $options->{modules},
        query => $prodigal_cds,
        output => $final_output,
        output_dir => $output_dir,
        shell => '/usr/bin/env perl',
        test_file => $options->{test_file},);
    $tjob->{prodigal_job} = $term_prodigal;
    return($tjob);
}

=head2 C<Terminase_ORF_Reorder_Worker>

Reorder an assembly from the result of a blast(or presumbly other) search.

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
