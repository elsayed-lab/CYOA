# -*-Perl-*-
use Test::More qw"no_plan";
use Bio::Adventure;
use Cwd;
use File::Copy qw"cp";
use File::Path qw"remove_tree make_path rmtree";
use File::ShareDir qw"dist_file";
use String::Diff qw"diff";

my $start = getcwd();
my $new = 'test_output';
mkdir($new);
chdir($new);

## Copy the reads for running the tests.
my $input_r1 = dist_file('Bio-Adventure', 'r1.fastq.xz');
ok(cp($input_r1, 'r1.fastq.xz'), 'Copying r1.') if (!-r 'r1.fastq.xz');
my $input_r2 = dist_file('Bio-Adventure', 'r2.fastq.xz');
ok(cp($input_r2, 'r2.fastq.xz'), 'Copying r2.') if (!-r 'r2.fastq.xz');

## Invoke the pipeline, keep it within our test directory with basedir.
my $cyoa = Bio::Adventure->new(cluster => 0, basedir => cwd());
my $assemble = $cyoa->Bio::Adventure::Pipeline::Phage_Assemble(
    input => 'r1.fastq.xz:r2.fastq.xz',);

## Check the trimomatic output.
my $trimomatic_output = 'outputs/01trimomatic/r1-trimomatic.stdout';
ok(-f $trimomatic_output, 'The trimomatic output was created.');
my $expected = qq"Using PrefixPair: 'TACACTCTTTCCCTACACGACGCTCTTCCGATCT' and 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
Using Long Clipping Sequence: 'GGAATTCTCGGGTGCCAAGGAACTCCAA'
Using Long Clipping Sequence: 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA'
Using Long Clipping Sequence: 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
Using Long Clipping Sequence: 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
Using Long Clipping Sequence: 'TACACTCTTTCCCTACACGACGCTCTTCCGATCT'
Using Long Clipping Sequence: 'AGATCGGAAGAGCACACGTCTGAAC'
ILLUMINACLIP: Using 1 prefix pairs, 6 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences
Input Read Pairs: 108212 Both Surviving: 95613 (88.36%) Forward Only Surviving: 10173 (9.40%) Reverse Only Surviving: 865 (0.80%) Dropped: 1561 (1.44%)
TrimmomaticPE: Completed successfully
";
my $actual = qx"tail ${trimomatic_output}";
unless(ok($expected eq $actual, 'Did trimomatic return expected results?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

## Look at the fastqc outputs
my $fastqc_summary = 'outputs/02fastqc/r1_fastqc/summary.txt';
ok(-f $fastqc_summary, 'The trimomatic output was created.');
my $actual = qx"less ${fastqc_summary}";
my $expected = qq"PASS\tBasic Statistics\t63
WARN\tPer base sequence quality\t63
WARN\tPer tile sequence quality\t63
PASS\tPer sequence quality scores\t63
FAIL\tPer base sequence content\t63
FAIL\tPer sequence GC content\t63
PASS\tPer base N content\t63
WARN\tSequence Length Distribution\t63
WARN\tSequence Duplication Levels\t63
WARN\tOverrepresented sequences\t63
PASS\tAdapter Content\t63
";
unless(ok($expected eq $actual, 'Did fastqc return expected outputs?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

## Check out the results from RACER, this assumes the output compile flag is set.
my $racer_stdout = 'outputs/03racer/racer.stdout';
ok(-f $racer_stdout, 'The racer output was created.');
my $actual = qx"head ${racer_stdout}";
my $expected = qq"New peak Memory = 19 MB. goodRead array added 19 MB.
New peak Memory = 57 MB. readLocation array added 38 MB.
New peak Memory = 61 MB. charReads array added 4 MB.
New peak Memory = 67 MB. binReads array added 5 MB.
reading input file...done reading input file
Reads file processed. Successfully converted into binary file in 0 seconds.
Total reads with > 50% N's = 0
New peak Memory = 86 MB. counters array added 18 MB.
New peak Memory = 109 MB. witness array added 23 MB.
Genome length = 1000000
";
unless(ok($expected eq $actual, 'Did racer return expected outputs?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

## Look at the standard kraken report.
my $kraken_report = 'outputs/04kraken_standard/kraken_report.txt';
ok(-f $kraken_report, 'The kraken report was created.');
my $actual = qx"head ${kraken_report}";
my $expected = qq"d__Bacteria\t5003
d__Bacteria|p__Proteobacteria\t4884
d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria\t4748
d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales\t3523
d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae\t3096
d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Escherichia\t3023
d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Escherichia|s__Escherichia marmotae\t3023
d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Citrobacter\t72
d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Citrobacter|s__Citrobacter portucalensis\t72
d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Morganellaceae\t425
";
unless(ok($expected eq $actual, 'Did kraken provide the expected report?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

## See if the kraken-based filter worked.
my $kraken_filter = 'outputs/05filter_kraken_host/kraken_filter.log';
ok(-f $kraken_filter, 'The kraken filter log was created.');
my $actual = qx"head ${kraken_filter}";
my $expected = qq"Starting search for best kraken host strain.
Reading kraken report: outputs/04kraken_standard/kraken_report.txt.
species Microbulbifer sp. ALW1 was observed 71 times.
species Candidatus Solibacter usitatus was observed 118 times.
species Xylophilus rhododendri was observed 1 times.
species Proteus vulgaris was observed 1 times.
species Providencia rettgeri was observed 379 times.
species Pseudomonas sp. R76 was observed 1152 times.
species Citrobacter portucalensis was observed 72 times.
species Providencia huaxiensis was observed 33 times.
";
unless(ok($expected eq $actual, 'Did kraken provide the expected filter log?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

## Look at the viral kraken results.
my $kraken_viral = 'outputs/06kraken_viral/kraken_report.txt';
ok(-f $kraken_viral, 'The kraken viral report was created.');
my $actual = qx"head ${kraken_viral}";
my $expected = qq"d__Viruses\t90498
d__Viruses|k__Heunggongvirae\t90498
d__Viruses|k__Heunggongvirae|p__Uroviricota\t90498
d__Viruses|k__Heunggongvirae|p__Uroviricota|c__Caudoviricetes\t90498
d__Viruses|k__Heunggongvirae|p__Uroviricota|c__Caudoviricetes|o__Caudovirales\t90498
d__Viruses|k__Heunggongvirae|p__Uroviricota|c__Caudoviricetes|o__Caudovirales|f__Autographiviridae\t90477
d__Viruses|k__Heunggongvirae|p__Uroviricota|c__Caudoviricetes|o__Caudovirales|f__Autographiviridae|g__Kayfunavirus\t90464
d__Viruses|k__Heunggongvirae|p__Uroviricota|c__Caudoviricetes|o__Caudovirales|f__Autographiviridae|g__Kayfunavirus|s__Escherichia virus EcoDS1\t15752
d__Viruses|k__Heunggongvirae|p__Uroviricota|c__Caudoviricetes|o__Caudovirales|f__Autographiviridae|g__Kayfunavirus|s__Citrobacter virus CR44b\t10348
d__Viruses|k__Heunggongvirae|p__Uroviricota|c__Caudoviricetes|o__Caudovirales|f__Autographiviridae|g__Kayfunavirus|s__Escherichia virus ST31\t9695
";
unless(ok($expected eq $actual, 'Did kraken provide the expected filter log?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

## Check the unicycler outputs.
my $unicycler_assembly = 'outputs/07unicycler/test_output_final_assembly.fasta';
ok(-f $unicycler_assembly, 'The unicycler fasta output was created.');
my $actual = qx"head -n 1 ${unicycler_assembly}";
my $expected = qq">1 length=40082 depth=1.00x circular=true
";
unless(ok($expected eq $actual, 'Was the assembly created with expected size/depth?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

## Look at the filter_depth results.
my $depth_log = 'outputs/08filter_depth/filter_depth.log';
ok(-f $depth_log, 'The filter depth log was created.');
my $actual = qx"more ${depth_log}";
my $expected = qq"Starting depth coverage filter of outputs/07unicycler/test_output_final_assembly.fasta.
Any contigs with a coverage ratio vs. the highest coverage of < 0.2 will be dropped.
Writing filtered contigs to outputs/08filter_depth/final_assembly.fasta
The range of observed coverages is 1.00 <= x <= 1.00
Writing 1 with normalized coverage: 1
";
unless(ok($expected eq $actual, 'Is the depth log as expected?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

## Check the phastaf results.
my $phastaf_log = 'outputs/09phastaf_08filter_depth/diamond.coords';
ok(-f $phastaf_log, 'The phastaf coordinate file exists.');
my $actual = qx"head ${phastaf_log}";
my $expected = qq"1\t29601\t33485\tPHAGE_Cronob_Dev2_NC_023558-gi|651729639|ref|YP_009005149.1|\t91.8
1\t2474\t5152\tPHAGE_Cronob_Dev2_NC_023558-gi|651729605|ref|YP_009005115.1|\t94.0
1\t13177\t15342\tPHAGE_Cronob_Dev2_NC_023558-gi|651729621|ref|YP_009005131.1|\t90.6
1\t26882\t29161\tPHAGE_Cronob_Dev2_NC_023558-gi|651729638|ref|YP_009005148.1|\t88.7
1\t23392\t25746\tPHAGE_Citrob_CR44b_NC_023576-gi|589286983|ref|YP_009007177.1|\t75.8
1\t37788\t39527\tPHAGE_Escher_vB_EcoP_GA2A_NC_031943-gi|100054|ref|YP_009324988.1|\t97.4
1\t10828\t12501\tPHAGE_Entero_K1F_NC_007456-gi|77118185|ref|YP_338107.1|\t97.3
1\t18627\t20192\tPHAGE_Cronob_Dev2_NC_023558-gi|651729630|ref|YP_009005140.1|\t95.8
1\t21351\t22379\tPHAGE_Citrob_CR44b_NC_023576-gi|589286980|ref|YP_009007173.1|\t98.8
1\t16380\t17240\tPHAGE_Cronob_Dev2_NC_023558-gi|651729625|ref|YP_009005135.1|\t94.4
";
unless(ok($expected eq $actual, 'Is the depth log as expected?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}


chdir($start);
