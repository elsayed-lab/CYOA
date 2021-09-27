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

=head2 C<Get_DTR>

Extract the direct-terminal-repeats from a phageterm run.

=cut
sub Get_DTR {
    my ($class, %args) = @_;
    unless (-r $args{input}) {
        return(undef);
    }

    my $input_genome = Bio::SeqIO->new(-file => $args{input}, -format => 'Fasta');
    my $sequence = '';
    my $length = 0;
    my $id = '';
    ## I think there are usually 2 entries in these files, one for each of the methods
    ## phageterm uses to look for a repeat, but afaik they end up the same.
    my %dtr_features = ();
  INPUT: while (my $genome_seq = $input_genome->next_seq()) {
      next INPUT unless(defined($genome_seq->id));
      $id = $genome_seq->id;
      $sequence = $genome_seq->seq;
      $length = $genome_seq->length;
      my $dtr_feature = Bio::SeqFeature::Generic->new(
          -primary => 'misc_feature',
          -seq_id => $id,
          -source => 'PhageTerm',
          -start => 1,
          -end => $length,
          -strand => +1,
          -score => undef,
          -frame => 0,
          -tag => {
              'product' => 'Direct Terminal Repeat',
              'inference' => 'COORDINATES:profile:PhageTerm',
          },);
      $dtr_features{$id} = $dtr_feature;
  }
    return(%dtr_features);
}

sub Classify_Phage {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        evalue => 0.01,
        blast_tool => 'tblastx',
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
Bio::Adventure::Phage::Blast_Classify(\$h,
  comment => '${comment}',
  blast_tool => '$options->{blast_tool}',
  evalue => '$options->{evalue}',
  input => '$options->{input}',
  library => '$options->{library}',
  output => '${output_tsv}',
  output_blast => '${output_blast}',
  output_dir => '${output_dir}',
  output_log => '${output_log}',
  topn => '$options->{topn}',
);
?;

    my $cjob = $class->Submit(
        jdepends => $options->{jdepends},
        blast_tool => $options->{blast_tool},
        evalue => $options->{evalue},
        input => $options->{input},
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
    die("Could not find $options-{blast_tool} in your PATH.") unless($check);

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
    print $final_fh "contig\tquery_description\ttaxon\tname\tquery_length\thit_length\thit_accession\thit_description\thit_bit\thit_sig\thit_score\n";
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
    print $log_message;
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
          my $taxon;
          my $xref_found = 0;
          XREFLOOP: for my $hashref (@{$xref_aoh}) {
              my $first_name = $hashref->{virusgenbankaccession};
              my $second_name = $hashref->{virusrefseqaccession};
              my $hash_taxon = $hashref->{taxon};
              my $hash_longname = $hashref->{virusnames};
              if ($hit_name eq $first_name or $hit_name eq $second_name) {
                  $xref_found++;
                  print $log "A cross reference was found to the ICTV taxon: ${hash_taxon}.\n";
                  $longname = $hash_longname;
                  $taxon = $hash_taxon;
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
              taxon => $taxon,);
          push (@hit_lst, \%hit_datum);
          $result_data->{$query_name}->{hit_data}->{$hit_name} = \%hit_datum;
          $hit_count++;
      } ## End the hitloop.
  } ## End of the blast search

    print $log "The number of results recorded: ${result_count}.\n";
    ## Sort largest to smallest score.
    my @sorted = sort { $b->{score} <=> $a->{score} } @hit_lst;
    my $header = qq"query_name\tquery_description\ttaxon\tname\tlength\tquery_length\thit_acc\thit_description\thit_bit\thit_sig\thit_score\n";
    for my $d (@sorted) {
        my $hit_string = qq"$d->{query_name}\t$d->{query_descr}\t$d->{taxon}\t$d->{longname}\t$d->{length}\t$d->{query_length}\t$d->{acc}\t$d->{description}\t$d->{bit}\t$d->{sig}\t$d->{score}\n";
        print $final_fh $hit_string;
    }
    $final_fh->close();
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    $log->close();
    return($result_data);
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
        jprefix => '14',);
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('PhageTerm.py');
    die("Could not find phageterm in your PATH.") unless($check);
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
        my $r1_dirname = $inputs->{directory};
        my $r1_filename = basename($in[0]);
        my $r2_filename = basename($in[1]);
        $uncompress_string = qq!
less \${start}/${r1_dirname}/${r1_filename} > r1.fastq && \\
  less \${start}/${r1_dirname}/${r2_filename} > r2.fastq
!;
        $input_string = qq! -f r1.fastq -p r2.fastq !;
        $delete_string = qq!rm r1.fastq && rm r2.fastq!;
    } else {
        my $r1_filename = basename($options->{input});
        my $r1_dirname = dirname($options->{input});
        my $r1 = abs_path($options->{input});
        $uncompress_string = qq!
less \${start}/outputs/${r1_dirname}/${r1_filename} > r1.fastq
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
  --report_title ${cwd_name} \\
  2>phageterm.err 1>phageterm.out
sleep 5

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
        jname => "phageterm_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => 24,
        modules => $options->{modules},
        output => $output_file,
        prescript => $options->{prescript},
        postscript => $options->{postscript},
        test_file => $test_file,);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($phageterm);
}

sub Phastaf {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        cpus => '8',
        modules => ['phastaf'],
        jprefix => '14',
        );
    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('phastaf');
    die("Could not find phastaf in your PATH.") unless($check);

    my $job_name = $class->Get_Job_Name();
    my $input_paths = $class->Get_Paths($options->{input});
    my $input_full = $input_paths->{fullpath};
    my $output_dir = qq"outputs/$options->{jprefix}phastaf_$input_paths->{dirname}";
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
  2>${output_dir}/phastaf.err \\
  1>${output_dir}/phastaf.out
?;
    my $output_file = qq"${output_dir}/something.txt";
    my $phastaf = $class->Submit(
        cpus => $options->{cpus},
        comment => $comment,
        jdepends => $options->{jdepends},
        jname => "phastaf_${job_name}",
        jprefix => $options->{jprefix},
        jstring => $jstring,
        jmem => 12,
        modules => $options->{modules},
        output => $output_file,
        prescript => $options->{prescript},
        postscript => $options->{postscript},);
    $loaded = $class->Module_Loader(modules => $options->{modules},
                                    action => 'unload');
    return($phastaf);
}

=head2 C<Search_Reorder>

Reorder an assembly from the result of a blast(or presumbly other) search.

=cut
sub Search_Reorder {
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
        test_file => '',
        );

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
    print $log "Testing for $options->{test_file}\n";
    if (-r $options->{test_file}) {
        my $full_input = abs_path($options->{input});
        my $pwd = getcwd();
        my $full_new = qq"${pwd}/${new_filename}";
        print $log "The test file exists, do not reorder the genome.
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
    } else {
        print $log "The test file does not exist, performing terminase search.\n";
    }

    my $loaded = $class->Module_Loader(modules => $options->{modules});
    my $check = which('fastx36');
    die("Could not find fasta36 in your PATH.") unless($check);

    my $comment = qq!## This will reorder a sequence database from the results of a database search.
!;

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
        program => $options->{fasta_tool},
        );
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

    my $result_data = {};
    my $search_output = Bio::SearchIO->new(-file => $fasta_outfile, -format => 'fasta');
    my $element_count = 0;
    my $result_count = 0;
    while (my $result = $search_output->next_result) {
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
          push(@hit_data, $hit_datum);
          $result_data->{$query_name}->{hit_data} = \@hit_data;
          $hit_count++;
      } ## End the hitloop.
    } ## End of the individual fasta search  (maybe not needed?)

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
    ## From the 3 of terminase to the S of start (revcomp'd) as the first subsequence
    ## Then from the D of end to 3+1 as the second.
    my @other_objects = ();
    my $first_object;
    while (my $seqobj = $in_assembly->next_seq()) {
        my $current_id = $seqobj->id();
        if ($current_id eq $best_contig) {
            $first_object = $seqobj;
        } else {
            push(@other_objects, $seqobj);
        }
    }

    ## The last thing to remember is that prodigal encodes its position information
    ## as a set of # thing # thing # things, so lets grab out the information of interest.
    my @start_end_array = split(/# /, $best_description);
    my $best_seq_start = $start_end_array[1];
    $best_seq_start =~ s/ //g;
    my $best_seq_end = $start_end_array[2];
    $best_seq_end =~ s/ //g;
    my $best_seq_strand = $start_end_array[3];
    $best_seq_strand =~ s/ //g;

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
    } else {
        $first_start = 1;
        $first_end = $best_seq_end;
        $second_start = $best_seq_end + 1;
        $second_end = $sequence_end;
        ## yeah yeah I know one is supposed to use revcomp, but damn
        ## it annoys me when it yells at me for not using a sequence object.
        $first_seq = $first_object->subseq($first_start, $first_end);
        $first_seq = reverse($first_seq);
        $first_seq =~ tr/ATGCU/TACGA/;
        $second_seq = $first_object->subseq($second_start, $second_end);
        $second_seq = reverse($second_seq);
        $second_seq =~ tr/ATGCU/TACGA/;
    }
    my $final_sequence = $first_seq . $second_seq;

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


=head2 Terminase_Reorder

Reorder an assembly so that a terminase gene is at the beginning.
This will take an arbitrary assembly file, invoke prodigal on it to get ORFs,
pass them to Reorder_Search (above), and copy the assembly to a new file
with the most likely terminase sequence at the beginning.

=cut
sub Terminase_ORF_Reorder {
    my ($class, %args) = @_;
    my $options = $class->Get_Vars(
        args => \%args,
        required => ['input'],
        evalue => 0.01,
        fasta_tool => 'fastx36',
        gcode => '11',
        jprefix => '15',
        library => 'terminase',
        modules => ['fasta', 'blast', 'blastdb'],
        species => 'phages',
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

    my $term_prodigal = $class->Bio::Adventure::Annotation::Prodigal(
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
    my $jstring = qq?
use Bio::Adventure;
use Bio::Adventure::Phage;
Bio::Adventure::Phage::Search_Reorder(\$h,
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
?;
    my $tjob = $class->Submit(
        jdepends => $options->{jdepends},
        evalue => $options->{evalue},
        fasta_tool => 'fastx36',
        input => $options->{input},
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

1;
