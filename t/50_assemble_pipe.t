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
my $actual = '';
my $expected = '';
my $test_file = '';
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
$test_file = 'outputs/01trimomatic/r1-trimomatic.stdout';
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
    my ($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
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
    my ($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

## Check out the results from RACER, this assumes the output compile flag is set.
$test_file = 'outputs/03racer/racer.stdout';
ok(-f $test_file, 'The racer output was created.');
$actual = qx"head ${test_file}";
$expected = qq"New peak Memory = 19 MB. goodRead array added 19 MB.
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
    my ($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
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
    my ($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

## See if the kraken-based filter worked.
$test_file = 'outputs/05filter_kraken_host/kraken_filter.log';
ok(-f $test_file, 'The kraken filter log was created.');
$actual = qx"head ${test_file}";
$expected = qq"Starting search for best kraken host strain.
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
    my ($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
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
    my ($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

## Check the unicycler outputs.
$test_file = 'outputs/07unicycler/test_output_final_assembly.fasta';
ok(-f $test_file, 'The unicycler fasta output was created.');
$actual = qx"head -n 1 ${test_file}";
$expected = qq">1 length=40082 depth=1.00x circular=true
";
unless(ok($expected eq $actual, 'Was the assembly created with expected size/depth?')) {
    my ($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
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
    my ($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
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
    my ($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

## Check the ICTV results.
$test_file = 'outputs/10classify_08filter_depth/ictv_filtered.tsv';
ok(-f $test_file, 'The ICTV results exist.');
$actual = qx"head -n 3 ${test_file}";
$expected = qq"contig\tquery_description\ttaxon\tname\tquery_length\thit_length\thit_accession\thit_description\thit_bit\thit_sig\thit_score\thit_family\thit_genus
1\tlength=40082 coverage=1x circular=true\tDuplodnaviria Heunggongvirae Uroviricota Caudoviricetes Caudovirales Autographiviridae Kayfunavirus Citrobacter virus CR44b\tCitrobacter phage CR44b\t39207\t40082\tHG818823\tCitrobacter phage CR44b complete genome.\t3223\t0.0\t7029\tAutographiviridae\tKayfunavirus
1\tlength=40082 coverage=1x circular=true\tDuplodnaviria Heunggongvirae Uroviricota Caudoviricetes Caudovirales Autographiviridae Kayfunavirus Escherichia virus ST31\tEscherichia phage ST31\t39693\t40082\tKY962008\tEscherichia phage ST31, complete genome.\t3013\t0.0\t6570\tAutographiviridae\tKayfunavirus
";
unless(ok($expected eq $actual, 'Is the depth log as expected?')) {
    my ($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
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
    my ($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

## Check the phageterm results.
$test_file = 'outputs/12phageterm_11rosalindplus/contig1_nrt.txt';
ok(-f $test_file, 'The phageterm results exist.');
$actual = qx"more ${test_file}";
$expected = qq"DTR (short)";
unless(ok($expected eq $actual, 'Are the phageterm results as expected?')) {
    my ($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

## Check the coverage results.
$test_file = 'outputs/13assembly_coverage_test_output/coverage.tsv';
ok(-f $test_file, 'The coverage results exist.');
$actual = qx"more ${test_file}";
$expected = qq"#ID\tAvg_fold\tLength\tRef_GC\tCovered_percent\tCovered_bases\tPlus_reads\tMinus_reads\tRead_GC\tMedian_fold\tStd_Dev\tUnder_5/100
1\t1142.8501\t40082\t0.0000\t100.0000\t40082\t94321\t94610\t0.5010\t1117\t232.09\t0
";
unless(ok($expected eq $actual, 'Is the coverage tsv as expected?')) {
    my ($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

## Check the terminase search results.
$test_file = 'outputs/14termreorder_12phageterm_11rosalindplus/terminase_summary.tsv';
ok(-f $test_file, 'The terminase search results exist.');
$actual = qx"head ${test_file}";
$expected = qq"Name\tLength\tAccession\tDescription\tScore\tSignificance\tBit\tHitStrand\tQueryStrand
AZS06569.1\t551\tAZS06569\tterminase [Mycobacterium phage JacoRen57]\t26.4\t2.6\t26.4\t0\t1
AUV61411.1\t503\tAUV61411\tlarge terminase [Pontimonas phage phiPsal1]\t32.5\t0.24\t32.5\t0\t1
QIN97975.1\t423\tQIN97975\tterminase [Salmonella phage pink]\t26.7\t4.8\t26.7\t0\t-1
ATN92724.1\t499\tATN92724\tterminase [Escherichia phage APC_JM3.2]\t28.3\t1.9\t28.3\t0\t-1
YP_009614610.1\t839\tYP_009614610\tterminase [Aeromonas phage pIS4-A]\t30.4\t0.72\t30.4\t0\t1
QCW19659.1\t696\tQCW19659\tterminase [Vibrio phage Va_90-11-286_p16]\t27.6\t1.9\t27.6\t0\t1
BAV81144.1\t505\tBAV81144\tterminase [Vibrio phage CKB-S2]\t26.7\t4.6\t26.7\t0\t1
ASZ73214.1\t492\tASZ73214\tterminase [Arthrobacter phage JayCookie]\t30.3\t1.6\t30.3\t0\t-1
ATS92560.1\t186\tATS92560\tterminase [Klebsiella phage 1611E-K2-1]\t26.7\t1.3\t26.7\t0\t-1
";
unless(ok($expected eq $actual, 'Did the terminase search return expected hits?')) {
    my ($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
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
    my ($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
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
    my ($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
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
    my ($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
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
    my ($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
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
    my ($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

## I messed up the jellyfish runs when standardizing the style.
## This test is therefore broken until I rerun.
## Check jellyfish
$test_file = 'outputs/20jellyfish_test_output/test_output_9.stdout';
ok(-f $test_file, 'The jellyfish results exist.');
$actual = qx"tail ${test_file}";
$expected = qq">1
ATACCTGCG
>1
AAGCCAGAA
>1
CCGTCCCAT
>1
TAGGAATCA
>1
AGCGAGAGT
";
unless(ok($expected eq $actual, 'Did we get expected jellyfish output?')) {
    my ($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

## Check aragorn
$test_file = 'outputs/21aragorn/aragorn.txt';
ok(-f $test_file, 'The aragorn results exist.');
$actual = qx"more ${test_file}";
$expected = qq">test_output_1
0 genes found
";
unless(ok($expected eq $actual, 'Did we get expected aragorn output?')) {
    my ($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

## I also broke the trnascan run...
$test_file = 'outputs/22trnascan/aragorn.txt';
ok(-f $test_file, 'The trnascan results exist.');
$actual = qx"more ${test_file}";
$expected = qq">test_output_1
0 genes found
";
unless(ok($expected eq $actual, 'Did we get expected trnascan output?')) {
    my ($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
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
    my ($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
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
    my ($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

## interproscan
$test_file = 'outputs/25interproscan_19merge_cds_predictions/interproscan.tsv';
ok(-f $test_file, 'The interproscan results exist.');
$actual = qx"head -n 3 ${test_file}";
$expected = qq"test_output_0061\t4eb154f61b29dbd33cf415f020f10f23\t87\tPfam\tPF11123\tDNA packaging protein\t6\t72\t3.1E-16\tT\t21-12-2021\tIPR024345\tDNA maturase, bacteriophage T7-like
test_output_0051\tf8d393226c15f8244492a90ec5880836\t188\tPfam\tPF17212\tTail tubular protein\t13\t179\t7.1E-45\tT\t21-12-2021\tIPR033767\tTail tubular protein Gp11
test_output_0028\t0994f22652375b1bbe2ca5661ffaeffb\t150\tPfam\tPF05367\tPhage endonuclease I\t1\t149\t1.6E-81\tT\t21-12-2021\tIPR008029\tBacteriophage T7, Gp3, endodeoxynuclease I
";
unless(ok($expected eq $actual, 'Did we get expected interproscan output?')) {
    my ($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
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
unless(ok($expected eq $actual, 'Did we get expected interproscan output?')) {
    my ($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
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
unless(ok($expected eq $actual, 'Did we get expected interproscan output?')) {
    my ($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

## Something something cgview...
$test_file = 'outputs/28cgview/test_output.xml';
ok(-f $test_file, 'The cgview input exists.');
$actual = qx"head -n 2 ${test_file}";
$expected = qq!<?xml version="1.0" encoding="ISO-8859-1"?>
<cgview backboneRadius="4000" backboneColor="rgb(102,102,102)" backboneThickness="40" featureSlotSpacing="10" labelLineLength="450" labelPlacementQuality="better" labelLineThickness="12" rulerPadding="130" tickThickness="18" shortTickThickness="18" arrowheadLength="60" rulerFont="SansSerif, plain, 130" rulerFontColor="rgb(0,0,0)" labelFont="SansSerif, plain, 130" isLinear="true" minimumFeatureLength="1.0" sequenceLength="41261" height="10000" width="10000" globalLabel="true" moveInnerLabelsToOuter="false" featureThickness="86.54" tickLength="45" useInnerLabels="true" shortTickColor="rgb(0,51,0)" longTickColor="rgb(0,51,0)" zeroTickColor="rgb(0,51,0)" showBorder="true" borderColor="black" backgroundColor="white" tickDensity="0.5">
!;
unless(ok($expected eq $actual, 'Did we get expected interproscan output?')) {
    my ($old, $new) = diff($expected, $actual);
    diag("--\n${old}\n--\n${new}\n");
}

chdir($start);
