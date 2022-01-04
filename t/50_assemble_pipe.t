# -*-Perl-*-
use Test::More qw"no_plan";
use Bio::Adventure;
use Cwd;
use File::Copy qw"cp";
use File::Path qw"remove_tree make_path rmtree";
use File::ShareDir qw"dist_file";
use String::Diff qw"diff";

my $start = getcwd();
my $a = 'test_output';
my $actual = '';
my $expected = '';
my $test_file = '';
mkdir($a);
chdir($a);

if (!defined($ENV{LESSOPEN})) {
    $ENV{LESSOPEN} = '| lesspipe %s';
}
if (!defined($ENV{LESS})) {
    $ENV{LESS} = '--buffers 0';
}

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
$test_file = 'outputs/01trimomatic/r1-trimomatic.stderr';
ok(-f $test_file, 'The trimomatic output was created.');
$expected = qq"Using PrefixPair: 'TACACTCTTTCCCTACACGACGCTCTTCCGATCT' and 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
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
$actual = qx"tail ${test_file}";
unless(ok($expected eq $actual, 'Did trimomatic return expected results?')) {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Look at the fastqc outputs
$test_file = 'outputs/02fastqc/r1_fastqc/summary.txt';
ok(-f $test_file, 'The trimomatic output was created.');
$actual = qx"less ${test_file}";
$expected = qq"PASS\tBasic Statistics\t63
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
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Check out the results from RACER, this assumes the output compile flag is set.
$test_file = 'outputs/03racer/racer.stdout';
ok(-f $test_file, 'The racer output was created.');
$actual = qx"grep changed ${test_file} | head";
$expected = qq"Number of changed positions\t\t62959
Number of changed positions after\t7704
Number of changed positions before\t55255
Number of changed positions\t\t937
Number of changed positions after\t934
Number of changed positions before\t3
Number of changed positions\t\t230
Number of changed positions after\t227
Number of changed positions before\t3
Number of changed positions\t\t104
";
unless(ok($expected eq $actual, 'Did racer return expected outputs?')) {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Look at the standard kraken report.
$test_file = 'outputs/04kraken_standard/kraken_report.txt';
ok(-f $test_file, 'The kraken report was created.');
$actual = qx"head ${test_file}";
$expected = qq"d__Bacteria\t5003
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
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## See if the kraken-based filter worked.
$test_file = 'outputs/05filter_kraken_host/kraken_filter.log';
ok(-f $test_file, 'The kraken filter log was created.');
$actual = qx"tail -n 2 ${test_file}";
$expected = qq"Filtering out reads which map to GCF_002900365.1.
Symlinking final output files to outputs/05filter_kraken_host
";
unless(ok($expected eq $actual, 'Did kraken provide the expected filter log?')) {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Look at the viral kraken results.
$test_file = 'outputs/06kraken_viral/kraken_report.txt';
ok(-f $test_file, 'The kraken viral report was created.');
$actual = qx"head ${test_file}";
$expected = qq"d__Viruses\t90498
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
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Check the unicycler outputs.
$test_file = 'outputs/07unicycler/test_output_final_assembly.fasta';
ok(-f $test_file, 'The unicycler fasta output was created.');
$actual = qx"head -n 1 ${test_file}";
$expected = qq">1 length=40082 depth=1.00x circular=true
";
unless(ok($expected eq $actual, 'Was the assembly created with expected size/depth?')) {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Look at the filter_depth results.
$test_file = 'outputs/08filter_depth/filter_depth.log';
ok(-f $test_file, 'The filter depth log was created.');
$actual = qx"more ${test_file}";
$expected = qq"Starting depth coverage filter of outputs/07unicycler/test_output_final_assembly.fasta.
Any contigs with a coverage ratio vs. the highest coverage of < 0.2 will be dropped.
Writing filtered contigs to outputs/08filter_depth/final_assembly.fasta
The range of observed coverages is 1.00 <= x <= 1.00
Writing 1 with normalized coverage: 1
";
unless(ok($expected eq $actual, 'Is the depth log as expected?')) {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Check the phastaf results.
$test_file = 'outputs/09phastaf_08filter_depth/diamond.coords';
ok(-f $test_file, 'The phastaf coordinate file exists.');
$actual = qx"head ${test_file}";
$expected = qq"1\t29601\t33485\tPHAGE_Cronob_Dev2_NC_023558-gi|651729639|ref|YP_009005149.1|\t91.8
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
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Check the ICTV results.
$test_file = 'outputs/10classify_08filter_depth/ictv_filtered.tsv';
ok(-f $test_file, 'The ICTV results exist.');
$actual = qx"head -n 3 ${test_file} | awk '{print \$2}'";
$expected = qq"query_description
length=40082
length=40082
";
unless(ok($expected eq $actual, 'Is the depth log as expected?')) {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Check the Rosalindplus results.
$test_file = 'outputs/11rosalindplus/final_assembly.log';
ok(-f $test_file, 'The Rosalindplus results exist.');
$actual = qx"more ${test_file}";
$expected = qq"Counting ORFs on the current Rosalind and Franklin strand.
If the Franklin strand is larger, flipping them.
1 was unchanged and has 49 plus and 0 minus ORFs.
";
unless(ok($expected eq $actual, 'Is the rosalindplus log as expected?')) {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Check the phageterm results.
$test_file = 'outputs/12phageterm_11rosalindplus/contig1_nrt.txt';
ok(-f $test_file, 'The phageterm results exist.');
$actual = qx"more ${test_file}";
$expected = qq"DTR (short)";
unless(ok($expected eq $actual, 'Are the phageterm results as expected?')) {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Check the coverage results.
$test_file = 'outputs/13assembly_coverage_test_output/coverage.tsv';
ok(-f $test_file, 'The coverage results exist.');
$actual = qx"more ${test_file}";
$expected = qq"#ID\tAvg_fold\tLength\tRef_GC\tCovered_percent\tCovered_bases\tPlus_reads\tMinus_reads\tRead_GC\tMedian_fold\tStd_Dev\tUnder_5/100
1\t1142.8501\t40082\t0.0000\t100.0000\t40082\t94321\t94610\t0.5010\t1117\t232.09\t0
";
unless(ok($expected eq $actual, 'Is the coverage tsv as expected?')) {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Check the terminase search results.
$test_file = 'outputs/14termreorder_12phageterm_11rosalindplus/terminase_summary.tsv';
ok(-f $test_file, 'The terminase search results exist.');
$actual = qx"head -n 3 ${test_file}";
$expected = qq"Name\tLength\tAccession\tDescription\tScore\tSignificance\tBit\tHitStrand\tQueryStrand
AZS06569.1\t551\tAZS06569\tterminase [Mycobacterium phage JacoRen57]\t26.4\t2.6\t26.4\t0\t1
AUV61411.1\t503\tAUV61411\tlarge terminase [Pontimonas phage phiPsal1]\t32.5\t0.24\t32.5\t0\t1
";
unless(ok($expected eq $actual, 'Did the terminase search return expected hits?')) {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Check the prokka run.
$test_file = 'outputs/15prokka/test_output.ffn';
ok(-f $test_file, 'The prokka results exist.');
$actual = qx"head ${test_file}";
$expected = qq">test_output_00001 hypothetical protein
ATGGAACGTAACGCTGACGCATACTATGAGCTGCTGAATGCAACCGTTAAAGCATTTAAC
GAGCGTGTTCAGTACGACGAAATAGCTAAAGGTGATGACTACCATGATGCGCTGCATGAA
GTCGTAGACGGTCAGGTTCCGCACTATTACCACGAGATCTTCACGGTGATGGCTGCTGAT
GGTATTGACATTGAGTTTGAAGACTCTGGGCTGATGCCTGAGACCAAGGACGTAACGCGC
ATACTGCAAGCTCGCATCTATGAGGCACTGTATAACGGCGTGTCTAATAGCTCGGATGTG
GTCTGGTTTGAGGCTGAAGAGAGCGACGAAGAGGGTAAGTATTGGGTAGTTGACGCTAAA
ACGGGACTATTCGCTGAGCAAGCTATACCTCTTGAGGTCGCTATTGCATCTGCCAAAGAC
CTCTATGCGGTAGGTCATCACATGAAAGTCGAAGACATTAACGATAACGTAGTGTTCGAC
CCTGCGGCTGAAGAGGACTGCGAGTGA
";
unless(ok($expected eq $actual, 'Did prokka return an expected contig?')) {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Check the prodigal run.
$test_file = 'outputs/16prodigal_test_output.fna/predicted_cds.fasta';
ok(-f $test_file, 'The prodigal results exist.');
$actual = qx"head ${test_file}";
$expected = qq">gnl|Prokka|test_output_1_1 # 1 # 117 # 1 # ID=1_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.504
TCTCACAGTTCAAGAACCTCAAGTCTCCCCATAGGCCCTCTTTCAGTCCAGACCAAAGGCCCTACCCCAG
TCTATCATAAGGTTGGACCGATGGTCAAGACTTCAGGTCAACGATAG
>gnl|Prokka|test_output_1_2 # 1267 # 1773 # 1 # ID=1_2;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.489
ATGGAACGTAACGCTGACGCATACTATGAGCTGCTGAATGCAACCGTTAAAGCATTTAACGAGCGTGTTC
AGTACGACGAAATAGCTAAAGGTGATGACTACCATGATGCGCTGCATGAAGTCGTAGACGGTCAGGTTCC
GCACTATTACCACGAGATCTTCACGGTGATGGCTGCTGATGGTATTGACATTGAGTTTGAAGACTCTGGG
CTGATGCCTGAGACCAAGGACGTAACGCGCATACTGCAAGCTCGCATCTATGAGGCACTGTATAACGGCG
TGTCTAATAGCTCGGATGTGGTCTGGTTTGAGGCTGAAGAGAGCGACGAAGAGGGTAAGTATTGGGTAGT
TGACGCTAAAACGGGACTATTCGCTGAGCAAGCTATACCTCTTGAGGTCGCTATTGCATCTGCCAAAGAC
";
unless(ok($expected eq $actual, 'Did prodigal return expected CDS predictions?')) {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Check the glimmer run.
$test_file = 'outputs/17glimmer/glimmer3.predict';
ok(-f $test_file, 'The glimmer results exist.');
$actual = qx"head ${test_file}";
$expected = qq">gnl|Prokka|test_output_1
orf00003      131      265  +2     0.84
orf00005      352      305  -2     1.29
orf00010     1052     1129  +2     7.87
orf00011     1178     1225  +2     9.71
orf00013     1267     1773  +1    14.08
orf00015     1821     1928  +3     2.10
orf00016     1931     2056  +2     3.05
orf00019     2262     2411  +3     4.05
orf00020     2447     2737  +2     2.27
";
unless(ok($expected eq $actual, 'Did glimmer return expected CDS predictions?')) {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Check the phanotate run.
$test_file = 'outputs/18phanotate/test_output_phanotate.tsv.xz';
ok(-f $test_file, 'The phanotate results exist.');
$actual = qx"less ${test_file} | head -n 5";
$expected = qq"#id:\tgnl|Prokka|test_output_1
#START\tSTOP\tFRAME\tCONTIG\tSCORE
1\t117\t+\tgnl|Prokka|test_output_1\t-1.248555528686707940866691777\t
183\t302\t+\tgnl|Prokka|test_output_1\t-0.2175130562377134455954126775\t
477\t617\t+\tgnl|Prokka|test_output_1\t-0.07018008835792925643848556090\t
";
unless(ok($expected eq $actual, 'Did phanotate return expected CDS predictions?')) {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Check the merge_cds.
$test_file = 'outputs/19merge_cds_predictions/test_output.tsv';
ok(-f $test_file, 'The merge_cds results exist.');
$actual = qx"head ${test_file}";
$expected = qq"locus_tag\tcontig\ttype\tsource\tstart\tend\tstrand\tcds_prediciton\taa_sequence
test_output_0001\ttest_output_1\tCDS\t\t1\t117\t1\tphanotate, score: -1.25\tSHSSRTSSLPIGPLSVQTKGPTPVYHKVGPMVKTSGQR
test_output_0002\ttest_output_1\tCDS\t\t131\t265\t1\tglimmer\tVVVETIGWDYWLSLSLLLAAGVTAGSQWVGWVETLVCSLVSQCN
test_output_0003\ttest_output_1\tCDS\t\t183\t302\t1\tphanotate, score: -0.218\tLLLALLLEVSGSGGSRLSYALWSLSVINAIMVTIHERKT
test_output_0004\ttest_output_1\tCDS\t\t305\t352\t-1\tglimmer\tLTTVAKVSRVASAMN
test_output_0005\ttest_output_1\tCDS\t\t477\t617\t1\tphanotate, score: -0.0702\tLDQKFETTSHSSRTSSLPIGPLSVQTKGPTPVYHKVGPMVKTSGQR
test_output_0006\ttest_output_1\tCDS\t\t1052\t1129\t1\tglimmer\tMRCKLYHIEATVATNTELVYLKSSD
test_output_0007\ttest_output_1\tCDS\t\t888\t1148\t-1\tphanotate, score: -0.107\tLLKSIPFSQRTSGRPVQCWSPPLLRCGTAYISSLLLVNYFLSSACCSYDLSGCLLNRDDPASSLSGCCRVVLTEAIKPQSRPIVNM
test_output_0008\ttest_output_1\tCDS\t\t1178\t1225\t1\tglimmer\tVINYRVFESTPEGPD
test_output_0009\ttest_output_1\tCDS\t\t1267\t1773\t1\tphanotate, score: -5930\tMERNADAYYELLNATVKAFNERVQYDEIAKGDDYHDALHEVVDGQVPHYYHEIFTVMAADGIDIEFEDSGLMPETKDVTRILQARIYEALYNGVSNSSDVVWFEAEESDEEGKYWVVDAKTGLFAEQAIPLEVAIASAKDLYAVGHHMKVEDINDNVVFDPAAEEDCE
";
unless(ok($expected eq $actual, 'Did we get expected merged CDS?')) {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## I messed up the jellyfish runs when standardizing the style.
## This test is therefore broken until I rerun.
## Check jellyfish
$test_file = 'outputs/20jellyfish_test_output/test_output_9.hist.xz';
ok(-f $test_file, 'The jellyfish results exist.');
$actual = qx"less ${test_file}";
$expected = qq"1 28869
2 4389
3 661
4 253
5 67
6 17
7 8
8 2
9 2
10 5
11 2
12 1
";
unless(ok($expected eq $actual, 'Did we get expected jellyfish output?')) {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Check aragorn
$test_file = 'outputs/21aragorn/aragorn.txt';
ok(-f $test_file, 'The aragorn results exist.');
$actual = qx"more ${test_file}";
$expected = qq">test_output_1
0 genes found
";
unless(ok($expected eq $actual, 'Did we get expected aragorn output?')) {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## I also broke the trnascan run...
$test_file = 'outputs/22trnascan/trnascan_relaxed.txt';
ok(-f $test_file, 'The trnascan results exist.');
$actual = qx"tail -n 4 ${test_file}";
$expected = qq"number of sequences= 1
number of bases tested (one strand)=41261
number of bases tested (both strands)= 82522
number of predicted tRNA=25
";
unless(ok($expected eq $actual, 'Did we get expected trnascan output?')) {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Trinotate
$test_file = 'outputs/23trinotate19merge_cds_predictions/test_output.ffn.tsv';
ok(-f $test_file, 'The trinotate results exist.');
$actual = qx"head -n 3 ${test_file}";
$expected = qq"#gene_id\ttranscript_id\tsprot_Top_BLASTX_hit\tRNAMMER\tprot_id\tprot_coords\tsprot_Top_BLASTP_hit\tphage_pep_BLASTX\tterminase_BLASTX\tphage_pep_BLASTP\tPfam\tSignalP\tTmHMM\teggnog\tKegg\tgene_ontology_BLASTX\tgene_ontology_BLASTP\tgene_ontology_Pfam\ttranscript\tpeptide
test_output_0001\ttest_output_0001\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.
test_output_0002\ttest_output_0002\t.\t.\t.\t.\t.\tCAJ29397.1^CAJ29397.1^Q:28-135,H:18-54^72.973%ID^E:1.74e-11^.^.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.
";
unless(ok($expected eq $actual, 'Did we get expected trinotate output?')) {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Abricate
$test_file = 'outputs/24abricate_19merge_cds_predictions/abricate_summary.txt';
ok(-f $test_file, 'The trinotate results exist.');
$actual = qx"more ${test_file}";
$expected = qq"#FILE\tNUM_FOUND
outputs/24abricate_19merge_cds_predictions/abricate_argannot.tsv\t0
outputs/24abricate_19merge_cds_predictions/abricate_card.tsv\t0
outputs/24abricate_19merge_cds_predictions/abricate_combined.tsv\t0
outputs/24abricate_19merge_cds_predictions/abricate_dbeth.tsv\t0
outputs/24abricate_19merge_cds_predictions/abricate_ecoh.tsv\t0
outputs/24abricate_19merge_cds_predictions/abricate_ecoli_vf.tsv\t0
outputs/24abricate_19merge_cds_predictions/abricate_megares.tsv\t0
outputs/24abricate_19merge_cds_predictions/abricate_mvir.tsv\t0
outputs/24abricate_19merge_cds_predictions/abricate_ncbi.tsv\t0
outputs/24abricate_19merge_cds_predictions/abricate_plasmidfinder.tsv\t0
outputs/24abricate_19merge_cds_predictions/abricate_resfinder.tsv\t0
outputs/24abricate_19merge_cds_predictions/abricate_vfdb.tsv\t0
";
unless(ok($expected eq $actual, 'Did we get expected abricate output?')) {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## It turns out that every invocation of interproscan is different in pretty much every file...
## interproscan
$test_file = 'outputs/25interproscan_19merge_cds_predictions/test_output.faa.gff3';
ok(-f $test_file, 'The interproscan results exist.');
$actual = qx"sort ${test_file} | head";
$expected = qq"AAESSAGNAKDSEDEARRIAESIKASGLIGYITRRSFENGFNVTTWNEVLLWEADGGYYR
AALQAKDFKAITDAIWNAKPFVPAGVVSAASLKDRTREAMLKAETEGLMFSSCTTLNAMT
AALSTVMSMGLAGIYYMA
AAPWLLAAV
ACP
ACPSFDLQRWLSTNELV
ADDAAFATVMRTMTDLSFFAKNAYMGVQNLTEIGGMLARGNVRALLHGVPMFRDLAFRNK
ADDAAFATVMRTMTDLSFFAKNAYMGVQNLTEIGGMLARGNVRALLHGVPMFRDLAFRNK
ADKFIPVEWLREATVRLPSGILIPKKGVK
ADKFIPVEWLREATVRLPSGILIPKKGVKK
";
unless(ok($expected eq $actual, 'Did we get expected interproscan output?')) {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## merge annotations 1
$test_file = 'outputs/26mergeannot/test_output_runlog.txt';
ok(-f $test_file, 'The merge annotations results exist.');
$actual = qx"head ${test_file}";
$expected = qq"Merging annotations and writing new output files:
gbf: outputs/26mergeannot/test_output.gbf, tbl: outputs/26mergeannot/test_output.tbl, xlsx: outputs/26mergeannot/test_output.xlsx.
Reading tsv data from outputs/19merge_cds_predictions/test_output.tsv to start.
Checking for ICTV classification data from outputs/10classify_08filter_depth/ictv_filtered.tsv.
Wrote outputs/26mergeannot/test_output.sbt with variables filled in.
Adding trinotate annotations from outputs/23trinotate19merge_cds_predictions/test_output.ffn.tsv.
Adding interproscan annotations from outputs/25interproscan_19merge_cds_predictions/interproscan.tsv.
Adding abricate annotations from outputs/24abricate_19merge_cds_predictions/abricate_combined.tsv.
Got DTR type: DTR (short).
Adding phageterm DTRs.
";
unless(ok($expected eq $actual, 'Did we get expected merge output?')) {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## merge annotations 2
$test_file = 'outputs/27mergeannot/test_output_stripped_runlog.txt';
ok(-f $test_file, 'The merge annotations results exist.');
$actual = qx"head ${test_file}";
$expected = qq"Merging annotations and writing new output files:
gbf: outputs/27mergeannot/test_output_stripped.gbf, tbl: outputs/27mergeannot/test_output_stripped.tbl, xlsx: outputs/27mergeannot/test_output_stripped.xlsx.
Reading tsv data from outputs/19merge_cds_predictions/test_output.tsv to start.
Checking for ICTV classification data from outputs/10classify_08filter_depth/ictv_filtered.tsv.
Wrote outputs/27mergeannot/test_output_stripped.sbt with variables filled in.
Adding trinotate annotations from outputs/23trinotate19merge_cds_predictions/test_output.ffn.tsv.
Adding interproscan annotations from outputs/25interproscan_19merge_cds_predictions/interproscan.tsv.
Adding abricate annotations from outputs/24abricate_19merge_cds_predictions/abricate_combined.tsv.
Got DTR type: DTR (short).
Adding phageterm DTRs.
";
unless(ok($expected eq $actual, 'Did we get expected merge output?')) {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Something something cgview...
$test_file = 'outputs/28cgview/test_output.xml';
ok(-f $test_file, 'The cgview input exists.');
$actual = qx"head -n 2 ${test_file}";
$expected = qq!<?xml version="1.0" encoding="ISO-8859-1"?>
<cgview backboneRadius="4000" backboneColor="rgb(102,102,102)" backboneThickness="40" featureSlotSpacing="10" labelLineLength="450" labelPlacementQuality="better" labelLineThickness="12" rulerPadding="130" tickThickness="18" shortTickThickness="18" arrowheadLength="60" rulerFont="SansSerif, plain, 130" rulerFontColor="rgb(0,0,0)" labelFont="SansSerif, plain, 130" isLinear="true" minimumFeatureLength="1.0" sequenceLength="41261" height="10000" width="10000" globalLabel="true" moveInnerLabelsToOuter="false" featureThickness="86.54" tickLength="45" useInnerLabels="true" shortTickColor="rgb(0,51,0)" longTickColor="rgb(0,51,0)" zeroTickColor="rgb(0,51,0)" showBorder="true" borderColor="black" backgroundColor="white" tickDensity="0.5">
!;
unless(ok($expected eq $actual, 'Did we get expected cgview output?')) {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## This is not working currently.
$test_file = 'outputs/29rnafold/test_output.fsa.tsv.xz';
ok(-f $test_file, 'The rnafold output exists.');
$actual = qx"less ${test_file} | head";
$expected = qq!contig	start	end	A	U	G	C	GC	paired	bp_percent	mfe	mfe_bp	mfe_gc	structure
test_output_1	-197	4	50	50	51	52	0.5124	73	0.3632	-63.90	-0.8753	-0.6204	..(((((....((((((((.((.....(((((...(((((...((((((.....((((.....((......)).....))))..((((((.((((((.(((....(((((.((.((...((((...)))).)).)).)))))))).))))))))).)))))))))...))))).)))))...))))))))))..)))))....
test_output_1	-194	7	50	50	51	52	0.5124	83	0.4129	-63.20	-0.7614	-0.6136	.........((((.((((((............(((((...((((((.....((((.....((......)).....))))..((((((.((((((.(((....(((((.((.((...((((...)))).)).)).)))))))).))))))))).)))))))))...)))))(((((........))))).))))))))))....
test_output_1	-191	10	48	52	50	53	0.5124	83	0.4129	-63.20	-0.7614	-0.6136	......((((.((((((............(((((...((((((.....((((.....((......)).....))))..((((((.((((((.(((....(((((.((.((...((((...)))).)).)).)))))))).))))))))).)))))))))...)))))(((((........))))).)))))))))).......
test_output_1	-188	13	50	51	51	51	0.5075	83	0.4129	-63.20	-0.7614	-0.6196	...((((.((((((............(((((...((((((.....((((.....((......)).....))))..((((((.((((((.(((....(((((.((.((...((((...)))).)).)).)))))))).))))))))).)))))))))...)))))(((((........))))).))))))))))..........
test_output_1	-185	16	51	50	50	52	0.5075	83	0.4129	-63.50	-0.7651	-0.6225	((((.((((((............(((((...((((((.....((((.....((......)).....))))..((((((.((((((.(((....(((((.((.((...((((...)))).)).)).)))))))).))))))))).)))))))))...)))))(((((........))))).)))))))))).............
test_output_1	-182	19	50	51	48	54	0.5075	91	0.4527	-58.40	-0.6418	-0.5725	..((((((............(((((...((((((.....((((.....((......)).....))))..((((((.((((((.(((....(((((.((.((...((((...)))).)).)).)))))))).))))))))).)))))))))...)))))(((((........))))).))))))....................
test_output_1	-179	22	51	51	48	53	0.5025	87	0.4328	-57.60	-0.6621	-0.5703	..........((((...(((((...((((((.....((((.....((......)).....))))..((((((.((((((.(((....(((((.((.((...((((...)))).)).)).)))))))).))))))))).)))))))))...)))))(((((........)))))((((....))))........))))......
test_output_1	-176	25	50	53	47	53	0.4975	81	0.4030	-58.90	-0.7272	-0.5890	..............(((((...((((((.....((((.....((......)).....))))..((((((.((((((.(((....(((((.((.((...((((...)))).)).)).)))))))).))))))))).)))))))))...))))).(((((...((((.(((.((((....)))).))....).))))...)))))
test_output_1	-173	28	50	51	47	55	0.5075	79	0.3930	-61.30	-0.7759	-0.6010	...........(((((...((((((.....((((.....((......)).....))))..((((((.((((((.(((....(((((.((.((...((((...)))).)).)).)))))))).))))))))).)))))))))...)))))((((((...((((.(((.((((....)))).))....).))))...))))))..
!;
unless(ok($expected eq $actual, 'Did we get expected rnafold output?')) {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Restriction catalog...
$test_file = 'outputs/30re_catalog/re_catalog.tsv';
ok(-f $test_file, 'The restriction endonuclease catalog exists.');
$actual = qx"head ${test_file}";
$expected = qq!RE	Site	Overhang	Cuts
AasI	GACNNNN^NNGTC	NN	30
AatI	AGG^CCT		1
AatII	GACGT^C	ACGT	13
AauI	T^GTACA	GTAC	9
Acc113I	AGT^ACT		1
Acc16I	TGC^GCA		10
Acc65I	G^GTACC	GTAC	1
AccB1I	G^GYRCC	GYRC	13
AccB7I	CCANNNN^NTGG	NNN	3
!;
unless(ok($expected eq $actual, 'Did we get expected restriction enzyme counts?')) {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## This is not working currently.
$test_file = 'outputs/31caical_test_output_vs_GCF_002900365.1/test_output_cai.txt';

## Phagepromoter
$test_file = 'outputs/32phagepromoter/output.fasta';
ok(-f $test_file, 'The phagepromoter output file exists.');
$actual = qx"head ${test_file}";
$expected = qq!>test_output_1:4118 host complement(69..98) score=0.716
TTGACCATCGGTCCAACCTTATGATAGACT
>test_output_1:4106 host complement(324..352) score=0.898
TTGACTACAGTAGCTAAGGTCAGTAGAGT
>test_output_1:4093 host complement(569..598) score=0.716
TTGACCATCGGTCCAACCTTATGATAGACT
>test_output_1:52 host (1036..1065) score=0.828
TTGACAAGCAGTAACGATGAGATGTAAGCT
>test_output_1:56 host (1134..1160) score=0.51
AATGCTCTTTAACAATCTGGATAAACT
!;
unless(ok($expected eq $actual, 'Did we get expected phage promoters?')) {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Rhotermination prediction
$test_file = 'outputs/33rhotermpredict_26mergeannot/predictions_coordinates_test_output_1.csv';
ok(-f $test_file, 'The rhotermpredict output file exists.');
$actual = qx"head ${test_file}";
$expected = qq!Region\tStart RUT\tEnd RUT\tStrand
T1\t0\t78\tplus
T2\t321\t399\tplus
T3\t700\t778\tplus
T4\t928\t1006\tplus
T5\t1631\t1709\tplus
T6\t2293\t2371\tplus
T7\t2854\t2932\tplus
T8\t3125\t3203\tplus
T9\t3622\t3700\tplus
!;
unless(ok($expected eq $actual, 'Did we get expected rhotermpredict output?')) {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

chdir($start);