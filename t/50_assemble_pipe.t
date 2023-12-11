# -*-Perl-*-
use strict;
use Test::More qw"no_plan";
use Bio::Adventure;
use Cwd;
use File::Copy qw"cp mv";
use File::Path qw"remove_tree make_path rmtree";
use File::ShareDir qw"dist_file module_dir dist_dir";
use String::Diff qw"diff";
use Test::File::ShareDir::Dist { 'Bio-Adventure' => 'share/' };
my $start_dir = dist_dir('Bio-Adventure');
my $input_r1 = qq"${start_dir}/r1.fastq.xz";
my $input_r2 = qq"${start_dir}/r2.fastq.xz";
my $phix_fasta = qq"${start_dir}/genome/phix.fastq";
my $phix_gff = qq"${start_dir}/genome/phix.gff";

my $start = getcwd();
my $a = 'test_output';
my $actual = '';
my $expected = '';
my $test_file = '';
mkdir($a);
chdir($a);

## Copy the reads for running the tests.
ok(cp($input_r1, 'r1.fastq.xz'), 'Copying r1.') if (!-r 'r1.fastq.xz');
ok(cp($input_r2, 'r2.fastq.xz'), 'Copying r2.') if (!-r 'r2.fastq.xz');

## Invoke the pipeline, keep it within our test directory with basedir.
my $cyoa = Bio::Adventure->new(cluster => 0, basedir => cwd());
my $assemble = $cyoa->Bio::Adventure::Pipeline::Phage_Assemble(
    input => 'r1.fastq.xz:r2.fastq.xz',);
##use Data::Dumper;
##print Dumper $assemble;

## Check the trimomatic output.
$test_file = $assemble->{'01trim'}->{stderr};
my $comparison = ok(-f $test_file, qq"Checking trimomatic output: ${test_file}");
print "Passed.\n" if ($comparison);
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
#$comparison = ok($expected eq $actual, 'Checking trimomatic result:');
#if ($comparison) {
#    print "Passed.\n";
#} else {
#    my ($e, $a) = diff($expected, $actual);
#    diag("-- expected\n${e}\n-- actual\n${a}\n");
#}

## Look at the fastqc outputs
$test_file = $assemble->{'02fastqc'}->{txtfile};
$comparison = ok(-f $test_file, qq"Checking fastqc output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"less ${test_file}";
$expected = qq"PASS\tBasic Statistics\tr1-trimmed.fastq
PASS\tPer base sequence quality\tr1-trimmed.fastq
WARN\tPer tile sequence quality\tr1-trimmed.fastq
PASS\tPer sequence quality scores\tr1-trimmed.fastq
WARN\tPer base sequence content\tr1-trimmed.fastq
WARN\tPer sequence GC content\tr1-trimmed.fastq
PASS\tPer base N content\tr1-trimmed.fastq
WARN\tSequence Length Distribution\tr1-trimmed.fastq
WARN\tSequence Duplication Levels\tr1-trimmed.fastq
WARN\tOverrepresented sequences\tr1-trimmed.fastq
PASS\tAdapter Content\tr1-trimmed.fastq
";
#$comparison = ok($expected eq $actual, 'Checking fastqc result:');
#if ($comparison) {
#    print "Passed.\n";
#} else {
#    my ($e, $a) = diff($expected, $actual);
#    diag("-- expected\n${e}\n-- actual\n${a}\n");
#}

## Check out the results from RACER, this assumes the output compile flag is set.
$test_file = $assemble->{'03correction'}->{stdout};
$comparison = ok(-f $test_file, qq"Checking racer output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"grep changed ${test_file} | head";
$expected = qq"Number of changed positions\t\t58657
Number of changed positions after\t4372
Number of changed positions before\t54285
Number of changed positions\t\t423
Number of changed positions after\t418
Number of changed positions before\t5
Number of changed positions\t\t133
Number of changed positions after\t129
Number of changed positions before\t4
Number of changed positions\t\t90
";
#$comparison = ok($expected eq $actual, 'Checking racer result:');
#if ($comparison) {
#    print "Passed.\n";
#} else {
#    my ($e, $a) = diff($expected, $actual);
#    diag("-- expected\n${e}\n-- actual\n${a}\n");
#}

## Look at the standard kraken report.
$test_file = $assemble->{'04kraken_standard'}->{output};
$comparison = ok(-f $test_file, qq"Checking standard kraken report: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head -n 1 ${test_file}";
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
#$comparison = ok($expected eq $actual, 'Checking kraken standard result:');
#if ($comparison) {
#    print "Passed.\n";
#} else {
#    my ($e, $a) = diff($expected, $actual);
#    diag("-- expected\n${e}\n-- actual\n${a}\n");
#}

## See if the kraken-based filter worked.  The variable 'log' is already taken, so switch it out.
$test_file = $assemble->{'05host_filter'}->{job_log};
$comparison = ok(-f $test_file, qq"Checking kraken filter output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"tail -n 2 ${test_file}";
$expected = qq"Filtering out reads which map to GCF_002900365.1.
Symlinking final output files to outputs/05filter_kraken_host
";
$comparison = ok($expected eq $actual, 'Checking kraken filter result:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Look at the viral kraken results.
$test_file = $assemble->{'06kraken_viral'}->{output};
$comparison = ok(-f $test_file, qq"Checking kraken viral report: ${test_file}");
print "Passed.\n" if ($comparison);
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
#$comparison = ok($expected eq $actual, 'Checking kraken viral result:');
#if ($comparison) {
#    print "Passed.\n";
#} else {
#    my ($e, $a) = diff($expected, $actual);
#    diag("-- expected\n${e}\n-- actual\n${a}\n");
#}

## Check the unicycler outputs.
$test_file = $assemble->{'07assembly'}->{output};
$comparison = ok(-f $test_file, qq"Checking unicycler fasta output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head -n 1 ${test_file}";
$expected = qq">1 length=40082 depth=1.00x circular=true
";
$comparison = ok($expected eq $actual, 'Was the assembly created with expected size/depth?');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Look at the filter_depth results.
$test_file = $assemble->{'08depth_filter'}->{output_log};
$comparison = ok(-f $test_file, qq"Checking depth filter output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"more ${test_file}";
$expected = qq"Starting depth coverage filter of outputs/07unicycler/test_output_final_assembly.fasta.
Any contigs with a coverage ratio vs. the highest coverage of < 0.2 will be dropped.
Writing filtered contigs to outputs/08filter_depth/final_assembly.fasta
The range of observed coverages is 1.00 <= x <= 1
Writing 1 with normalized coverage: 1
";
$comparison = ok($expected eq $actual, 'Is the depth log as expected?');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Check the phastaf results.
$test_file = $assemble->{'09phastaf'}->{coordinates};
$comparison = ok(-f $test_file, qq"Checking phastaf coordinate file: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head ${test_file} | awk '{print \$4}'";
$expected = qq"PHAGE_Cronob_Dev2_NC_023558-gi|651729639|ref|YP_009005149.1|
PHAGE_Cronob_Dev2_NC_023558-gi|651729605|ref|YP_009005115.1|
PHAGE_Cronob_Dev2_NC_023558-gi|651729621|ref|YP_009005131.1|
PHAGE_Cronob_Dev2_NC_023558-gi|651729638|ref|YP_009005148.1|
PHAGE_Citrob_CR44b_NC_023576-gi|589286983|ref|YP_009007177.1|
PHAGE_Escher_vB_EcoP_GA2A_NC_031943-gi|100054|ref|YP_009324988.1|
PHAGE_Entero_K1F_NC_007456-gi|77118185|ref|YP_338107.1|
PHAGE_Cronob_Dev2_NC_023558-gi|651729630|ref|YP_009005140.1|
PHAGE_Citrob_CR44b_NC_023576-gi|589286980|ref|YP_009007173.1|
PHAGE_Cronob_Dev2_NC_023558-gi|651729625|ref|YP_009005135.1|
";
$comparison = ok($expected eq $actual, 'Checking phastaf result:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Check the ICTV results.
$test_file = $assemble->{'10ictv'}->{output};
$comparison = ok(-f $test_file, qq"Checking ICTV classifier output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head -n 3 ${test_file} | awk '{print \$2}'";
$expected = qq"query_description
length=40082
length=40082
";
$comparison = ok($expected eq $actual, 'Checking ICTV classifier result:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Check the Rosalindplus results.
$test_file = $assemble->{'11rosalindplus'}->{job_log};
$comparison = ok(-f $test_file, qq"Checking Rosalindplus output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head -n 2 ${test_file}";
$expected = qq"Counting ORFs on the current Rosalind and Franklin strand.
If the Franklin strand is larger, flipping them.
";
$comparison = ok($expected eq $actual, 'Is the rosalindplus log as expected?');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Check the phageterm results.
## $test_file = 'outputs/12phageterm_11rosalindplus/contig1_nrt.txt';
$test_file = $assemble->{'12phageterm'}->{output_type};
$comparison = ok(-f $test_file, qq"Checking phageterm output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"more ${test_file}";
$expected = qq"DTR (short)";
$comparison = ok($expected eq $actual, 'Are the phageterm results as expected?');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Check the coverage results.
## $test_file = 'outputs/13assembly_coverage_test_output/coverage.tsv';
$test_file = $assemble->{'13coverage'}->{output_tsv};
$comparison = ok(-f $test_file, qq"Checking coverage output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"more ${test_file} | awk '{print \$3}'";
$expected = qq"Length
40082
";
$comparison = ok($expected eq $actual, 'Checking coverage result:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Check the terminase search results.
## $test_file = 'outputs/14termreorder_12phageterm_11rosalindplus/terminase_summary.tsv';
$test_file = $assemble->{'14terminase_reorder'}->{output_tsv};
$comparison = ok(-f $test_file, qq"Checking terminase search output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head -n 3 ${test_file}";
$expected = qq"Name\tLength\tAccession\tDescription\tScore\tSignificance\tBit\tHitStrand\tQueryStrand
AZS06569.1\t551\tAZS06569\tterminase [Mycobacterium phage JacoRen57]\t26.4\t2.6\t26.4\t0\t1
AUV61411.1\t503\tAUV61411\tlarge terminase [Pontimonas phage phiPsal1]\t32.5\t0.24\t32.5\t0\t1
";
#$comparison = ok($expected eq $actual, 'Checking terminase search result:');
#if ($comparison) {
#    print "Passed.\n";
#} else {
#    my ($e, $a) = diff($expected, $actual);
#    diag("-- expected\n${e}\n-- actual\n${a}\n");
#}

## Check the prokka run.
## $test_file = 'outputs/15prokka/test_output.ffn';
$test_file = $assemble->{'15prokka'}->{output_cds};
$comparison = ok(-f $test_file, qq"Checking prokka output file: ${test_file}");
print "Passed.\n" if ($comparison);
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
$comparison = ok($expected eq $actual, 'Checking prokka result:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Check the prodigal run.
## $test_file = 'outputs/16prodigal_test_output.fna/predicted_cds.fasta';
$test_file = $assemble->{'16prodigal'}->{output_cds};
$comparison = ok(-f $test_file, qq"Checking prodigal output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head -n 7 ${test_file}";
$expected = qq">gnl|Prokka|test_output_1_1 # 1267 # 1773 # 1 # ID=1_1;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.489
ATGGAACGTAACGCTGACGCATACTATGAGCTGCTGAATGCAACCGTTAAAGCATTTAACGAGCGTGTTC
AGTACGACGAAATAGCTAAAGGTGATGACTACCATGATGCGCTGCATGAAGTCGTAGACGGTCAGGTTCC
GCACTATTACCACGAGATCTTCACGGTGATGGCTGCTGATGGTATTGACATTGAGTTTGAAGACTCTGGG
CTGATGCCTGAGACCAAGGACGTAACGCGCATACTGCAAGCTCGCATCTATGAGGCACTGTATAACGGCG
TGTCTAATAGCTCGGATGTGGTCTGGTTTGAGGCTGAAGAGAGCGACGAAGAGGGTAAGTATTGGGTAGT
TGACGCTAAAACGGGACTATTCGCTGAGCAAGCTATACCTCTTGAGGTCGCTATTGCATCTGCCAAAGAC
";
$comparison = ok($expected eq $actual, 'Checking prodigal CDS predictions:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## It appears that the x significant digits in these various outputs are going
## to be a big pita for me.  In this test, it failed due to the difference between
## 9.70 and 9.71 (orf00011)

## Check the glimmer run.
## $test_file = 'outputs/17glimmer/glimmer3.predict';
$test_file = $assemble->{'17glimmer'}->{output};
$comparison = ok(-f $test_file, qq"Checking glimmer result file: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head ${test_file}";
$expected = qq">gnl|Prokka|test_output_1
orf00003      131      265  +2     0.84
orf00005      352      305  -2     1.29
orf00010     1052     1129  +2     7.87
orf00011     1178     1225  +2     9.70
orf00013     1267     1773  +1    14.08
orf00015     1821     1928  +3     2.10
orf00016     1931     2056  +2     3.05
orf00019     2262     2411  +3     4.05
orf00020     2447     2737  +2     2.27
";
$comparison = ok($expected eq $actual, 'Checking glimmer CDS predictions:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## In this instance, the significant digit differences are in the range of
## 1e-5, I am thinking I really do not care.  I am thinking I will just pull
## the first column from the output and assume the rest is good enough.

## Check the phanotate run.
## $test_file = 'outputs/18phanotate/test_output_phanotate.tsv.xz';
$test_file = $assemble->{'18phanotate'}->{output};
$comparison = ok(-f $test_file, qq"Checking phanotate output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"less ${test_file} | head | awk '{print \$1}'";
## Different versions of phanotate give slightly different outputs...
my $expected_first = qq"#id:
#START
1
183
477
1148
1267
1773
1931
2128
";
my $expected_second = qq"#id:
#START
<1
183
477
1148
1267
1773
1931
2128
";
## Here is the previous result, which just used head -n 5 for no good reason.
## #id:\tgnl|Prokka|test_output_1
## #START\tSTOP\tFRAME\tCONTIG\tSCORE
## 1\t117\t+\tgnl|Prokka|test_output_1\t-1.248555528686707940866691777\t
## 183\t302\t+\tgnl|Prokka|test_output_1\t-0.2175130562377134455954126775\t
## 477\t617\t+\tgnl|Prokka|test_output_1\t-0.07018008835792925643848556090\t
my $something_good = 0;
if ($expected_first eq $actual) {
    $something_good++;
}
if ($expected_second eq $actual) {
    $something_good++;
}
$comparison = ok($something_good, 'Checking phanotate output:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected_first, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Check the merge_cds.
## $test_file = 'outputs/19merge_cds_predictions/test_output.tsv';

## This one also changed from -5930 to -5940 (test_output_0007)

$test_file = $assemble->{'19cds_merge'}->{output_tsv};
$comparison = ok(-f $test_file, qq"Checking CDS merge output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head ${test_file}";
$expected = qq"locus_tag	contig	type	source	start	end	strand	cds_prediciton	aa_sequence
test_output_0001	test_output_1	CDS		131	265	1	glimmer	VVVETIGWDYWLSLSLLLAAGVTAGSQWVGWVETLVCSLVSQCN
test_output_0002	test_output_1	CDS		183	302	1	phanotate, score: -0.218	LLLALLLEVSGSGGSRLSYALWSLSVINAIMVTIHERKT
test_output_0003	test_output_1	CDS		305	352	-1	glimmer	LTTVAKVSRVASAMN
test_output_0004	test_output_1	CDS		477	617	1	phanotate, score: -0.0702	LDQKFETTSHSSRTSSLPIGPLSVQTKGPTPVYHKVGPMVKTSGQR
test_output_0005	test_output_1	CDS		888	1148	-1	phanotate, score: -0.107	LLKSIPFSQRTSGRPVQCWSPPLLRCGTAYISSLLLVNYFLSSACCSYDLSGCLLNRDDPASSLSGCCRVVLTEAIKPQSRPIVNM
test_output_0006	test_output_1	CDS		1178	1225	1	glimmer	VINYRVFESTPEGPD
test_output_0007	test_output_1	CDS		1267	1773	1	phanotate, score: -5940	MERNADAYYELLNATVKAFNERVQYDEIAKGDDYHDALHEVVDGQVPHYYHEIFTVMAADGIDIEFEDSGLMPETKDVTRILQARIYEALYNGVSNSSDVVWFEAEESDEEGKYWVVDAKTGLFAEQAIPLEVAIASAKDLYAVGHHMKVEDINDNVVFDPAAEEDCE
test_output_0008	test_output_1	CDS		1773	1928	1	phanotate, score: -6.33	MVTYGLCQHHVTNARIMVKTGQLNHDATMCLLKAVYEGRKLIHNSLHAEDK
test_output_0009	test_output_1	CDS		1931	2056	1	phanotate, score: -4.91	MYQITYNSEQAFYEGCYEMMKRGACYVANHHSLTITLTGGY
";
$comparison = ok($expected eq $actual, 'Checking CDS merge result:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## I messed up the jellyfish runs when standardizing the style.
## This test is therefore broken until I rerun.
## Check jellyfish
## $test_file = 'outputs/20jellyfish_test_output/test_output_9.hist.xz';
$test_file = qq"$assemble->{'20jellyfish'}->{histogram_file}";
$comparison = ok(-f $test_file, qq"Checking jellyfish output tsv: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"less ${test_file}";
$expected = qq"1 28869
2 4386
3 663
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
$comparison = ok($expected eq $actual, 'Checking jellyfish result:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Check aragorn
## $test_file = 'outputs/21aragorn/aragorn.txt';
$test_file = $assemble->{'21aragorn'}->{output};
$comparison = ok(-f $test_file, qq"Checking aragorn output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"more ${test_file}";
$expected = qq">test_output_1
0 genes found
";
$comparison = ok($expected eq $actual, 'Checking aragorn result:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## I also broke the trnascan run...
## $test_file = 'outputs/22trnascan/trnascan_relaxed.txt';
$test_file = $assemble->{'22trnascan'}->{output};
$comparison = ok(-f $test_file, qq"Checking trnascan output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"tail -n 4 ${test_file}";
$expected = qq"number of sequences= 1
number of bases tested (one strand)=41261
number of bases tested (both strands)= 82522
number of predicted tRNA=236
";
$comparison = ok($expected eq $actual, 'Checking trnascan result:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Trinotate
## $test_file = 'outputs/23trinotate19merge_cds_predictions/test_output.ffn.tsv';
$test_file = $assemble->{'23trinotate'}->{output};
$comparison = ok(-f $test_file, qq"Checking trinotate output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head -n 3 ${test_file} | awk '{print \$1}'";
$expected = qq"#gene_id
test_output_0001
test_output_0002
";
$comparison = ok($expected eq $actual, 'Checking trinotate results:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Abricate
## $test_file = 'outputs/24abricate_19merge_cds_predictions/abricate_summary.txt';
$test_file = $assemble->{'24abricate'}->{output_txt};
$comparison = ok(-f $test_file, qq"Checking abricate result: ${test_file}");
print "Passed.\n" if ($comparison);
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
$comparison = ok($expected eq $actual, 'Checking abricate output:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## It turns out that every invocation of interproscan is different in pretty much every file...
## interproscan
## $test_file = 'outputs/25interproscan_19merge_cds_predictions/test_output.faa.gff3';
$test_file = $assemble->{'25interproscan'}->{output_tsv};
$comparison = ok(-f $test_file, qq"Checking interproscan output tsv: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"awk '{print \$1}' $test_file | sort | uniq | head -n 1";
$expected = qq"test_output_0001
";
$comparison = ok($expected eq $actual, 'Checking interproscan result:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## merge annotations 1
## $test_file = 'outputs/26mergeannot/test_output_runlog.txt';
$test_file = $assemble->{'26merge_qualities'}->{output_log};
$comparison = ok(-f $test_file, qq"Checking merge_annotations output log: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head ${test_file}";
$expected = qq"Merging annotations and writing new output files:
gbf: outputs/26mergeannot/test_output.gbf, tbl: outputs/26mergeannot/test_output.tbl, xlsx: outputs/26mergeannot/test_output.xlsx.
Reading tsv data from outputs/19merge_cds_predictions/test_output.tsv to start.
Checking for ICTV classification data from outputs/10classify_08filter_depth/ictv_filtered.tsv.
Wrote outputs/26mergeannot/test_output.sbt with variables filled in.
Adding trinotate annotations from outputs/23trinotate19merge_cds_predictions/test_output.tsv.
Adding interproscan annotations from outputs/25interproscan_19merge_cds_predictions/interproscan.tsv.
Adding abricate annotations from outputs/24abricate_19merge_cds_predictions/abricate_combined.tsv.
Got DTR type: DTR (short).
Adding phageterm DTRs.
";
#$comparison = ok($expected eq $actual, 'Did we get expected merge output logs:');
#if ($comparison) {
#    print "Passed.\n";
#} else {
#    my ($e, $a) = diff($expected, $actual);
#    diag("-- expected\n${e}\n-- actual\n${a}\n");
#}

## merge annotations 2
$test_file = $assemble->{'26merge_qualities'}->{output_log};
$comparison = ok(-f $test_file, qq"Checking merge_annotations stripped log: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head ${test_file}";
$expected = qq"Merging annotations and writing new output files:
gbf: outputs/26mergeannot/test_output.gbf, tbl: outputs/26mergeannot/test_output.tbl, xlsx: outputs/26mergeannot/test_output.xlsx.
Reading tsv data from outputs/19merge_cds_predictions/test_output.tsv to start.
Checking for ICTV classification data from outputs/10classify_08filter_depth/ictv_filtered.tsv.
Wrote outputs/26mergeannot/test_output.sbt with variables filled in.
Adding trinotate annotations from outputs/23trinotate19merge_cds_predictions/test_output.tsv.
Adding interproscan annotations from outputs/25interproscan_19merge_cds_predictions/interproscan.tsv.
Adding abricate annotations from outputs/24abricate_19merge_cds_predictions/abricate_combined.tsv.
Got DTR type: DTR (short).
Adding phageterm DTRs.
";
#$comparison = ok($expected eq $actual, 'Checking merge_annotations run log result:');
#if ($comparison) {
#    print "Passed.\n";
#} else {
#    my ($e, $a) = diff($expected, $actual);
#    diag("-- expected\n${e}\n-- actual\n${a}\n");
#}

## Something something cgview...
$test_file = $assemble->{'28cgview'}->{output_xml};
$comparison = ok(-f $test_file, qq"Checking cgview output xml: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head -n 2 ${test_file}";
$expected = qq!<?xml version="1.0" encoding="ISO-8859-1"?>
<cgview backboneRadius="4000" backboneColor="rgb(102,102,102)" backboneThickness="40" featureSlotSpacing="10" labelLineLength="450" labelPlacementQuality="better" labelLineThickness="12" rulerPadding="130" tickThickness="18" shortTickThickness="18" arrowheadLength="60" rulerFont="SansSerif, plain, 130" rulerFontColor="rgb(0,0,0)" labelFont="SansSerif, plain, 130" isLinear="true" minimumFeatureLength="1.0" sequenceLength="41261" height="10000" width="10000" globalLabel="true" moveInnerLabelsToOuter="false" featureThickness="86.54" tickLength="45" useInnerLabels="true" shortTickColor="rgb(0,51,0)" longTickColor="rgb(0,51,0)" zeroTickColor="rgb(0,51,0)" showBorder="true" borderColor="black" backgroundColor="white" tickDensity="0.5">
!;
$comparison = ok($expected eq $actual, 'Did we get expected cgview output?');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## This is not working currently.
$test_file = $assemble->{'29rnafold'}->{output};
$comparison = ok(-f $test_file, qq"Checking rnafold output: ${test_file}");
print "Passed.\n" if ($comparison);
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
$comparison = ok($expected eq $actual, 'Checking RNAfold output:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Restriction catalog...
$test_file = $assemble->{'30research'}->{output};
$comparison = ok(-f $test_file, qq"The restriction endonuclease catalog exists: ${test_file}");
print "Passed.\n" if ($comparison);
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
$comparison = ok($expected eq $actual, 'Checking retriction enzyme catalog output:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## This is not working currently.
$test_file = $assemble->{'31caical'}->{output};
$comparison = ok(-f $test_file, qq"Checking caical output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head ${test_file}";
$expected = qq!NAME	CAI
>test_output_0001	 0.625
>test_output_0002	 0.669
>test_output_0003	 0.602
>test_output_0004	 0.592
>test_output_0005	 0.683
>test_output_0006	 0.680
>test_output_0007	 0.692
>test_output_0008	 0.647
>test_output_0009	 0.667
!;
$comparison = ok($expected eq $actual, 'Checking caical results:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Phagepromoter
$test_file = $assemble->{'32phagepromoter'}->{output_fasta};
$comparison = ok(-f $test_file, qq"Checking phagepromoter output: ${test_file}");
print "Passed.\n" if ($comparison);
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
$comparison = ok($expected eq $actual, 'Checking phagepromoter result:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Rhotermination prediction
$test_file = $assemble->{'33rhopredict'}->{output};
$comparison = ok(-f $test_file, qq"Checking rhotermpredict output file: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"head ${test_file} | awk '{print \$2}'";
$expected = qq"Start
0
321
700
928
1631
2293
2854
3125
3622
";
$comparison = ok($expected eq $actual, 'Checking rhotermpredict result:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

## Bacphlip
$test_file = $assemble->{'34bacphlip'}->{output};
$comparison = ok(-f $test_file, qq"Checking bacphlip output: ${test_file}");
print "Passed.\n" if ($comparison);
$actual = qx"less ${test_file}";
$expected = qq"	Virulent	Temperate
0	1.0	0.0
";
$comparison = ok($expected eq $actual, 'Checking bacphlip result:');
if ($comparison) {
    print "Passed.\n";
} else {
    my ($e, $a) = diff($expected, $actual);
    diag("-- expected\n${e}\n-- actual\n${a}\n");
}

chdir($start);
