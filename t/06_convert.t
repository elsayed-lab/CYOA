# -*-Perl-*-
BEGIN {
    use Test::More qw"no_plan";
    use Bio::Adventure;
    use File::Copy qw"cp";
    use Archive::Extract;
    use String::Diff qw( diff_fully diff diff_merge diff_regexp );
}

my $cyoa = Bio::Adventure->new();

diag("Copying data file to cwd().");
ok(cp("t/data/genome/mgas_5005.gb.xz", "mgas_5005.gb.xz"), '01 copy.');

my $ae = Archive::Extract->new(archive => "mgas_5005.gb.xz");
diag("Can I extract the compressed genbank archive?");
ok($ae->extract(), '02 extract.');

diag("Can I run Gb2Gff()?");
my $gff = Bio::Adventure::Convert::Gb2Gff($cyoa, input => qq"mgas_5005.gb");
diag("Do I get the expected number of sequences, nucleotides, and features?");
ok($gff->{num_sequences} eq '1', '03 number sequences.');
ok($gff->{total_nt} eq '1838562', '04 number nucleotides.');
ok($gff->{num_features} eq '3858', '05 number features.');
diag("Do the output files have the expected number of characters/lines?");

my ($expected, $actual);
diag("Check the gff file of all features.");
$actual = qx"wc mgas_5005_all.gff";
$expected = qq"   3862   38462 1211002 mgas_5005_all.gff\n";
ok($actual eq $expected, 'wc the all gff.');

diag("Check the fasta file of CDS features.");
$actual = qx"wc mgas_5005_cds.fasta";
$expected = qq"  29261   39835 1727732 mgas_5005_cds.fasta\n";
ok($actual eq $expected, 'wc the cds fasta.');

diag("Check the gff file of CDS features.");
$actual = qx"wc mgas_5005_cds.gff";
$expected = qq"  1842  20186 977161 mgas_5005_cds.gff\n";
ok($actual eq $expected, 'wc the cds gff.');

diag("Check the fasta file of the entire genome.");
$actual = qx"wc mgas_5005.fasta";
$expected = qq"  30644   30649 1869265 mgas_5005.fasta\n";
ok($actual eq $expected, 'wc the fasta.');

diag("Check the gff file of the entire genome.");
$actual = qx"wc mgas_5005.gff";
$expected = qq"  1927  17336 221048 mgas_5005.gff\n";
ok($actual eq $expected, 'wc the gff.');

diag("Check the gff file of the interCDS regions.");
$actual = qx"wc mgas_5005_interCDS.gff";
$expected = qq"  1842  20186 977155 mgas_5005_interCDS.gff\n";
ok($actual eq $expected, 'wc the intercds gff.');

diag("Check the peptide fasta file.");
$actual = qx"wc mgas_5005_pep.fasta";
$expected = qq" 11580  22154 648129 mgas_5005_pep.fasta\n";
ok($actual eq $expected, 'wc the peptide fasta.') ||
    (diag explain($actual, $expected));

my @files = ('mgas_5005_cds.gff', 'mgas_5005.gb.xz', 'mgas_5005_pep.fasta', 'mgas_5005_all.gff',
             'mgas_5005.fasta', 'mgas_5005.gff', 'mgas_5005_cds.fasta', 'mgas_5005.gb',
             'mgas_5005_interCDS.gff');

for my $t (@files) {
    diag("delete ${t}");
    ok(unlink($t));
}

##my($old, $new) = String::Diff::diff($expected_mgas_5005_all, $mgas_5005_all);
##print STDERR "$old\n";
##print STDERR "---\n";
##print STDERR "$new\n";
##print STDERR "---\n";
