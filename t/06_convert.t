# -*-Perl-*-
BEGIN {
    use Test::More qw"no_plan";
    use HPGL;
    use File::Copy qw"cp";
    use Archive::Extract;
    use String::Diff qw( diff_fully diff diff_merge diff_regexp );
}
diag("Copying data file to cwd().");
ok(cp("t/data/genome/mgas_5005.gb.xz", "mgas_5005.gb.xz"));
my $ae = new Archive::Extract(archive => "mgas_5005.gb.xz");
diag("Can I extract the compressed genbank archive?");
ok($ae->extract());
diag("Can I create a new HPGL using the genbank file?");
my $hpgl = new HPGL(input => qq"mgas_5005.gb", pbs => 0);
diag("Can I run Gb2Gff()?");
my $gff = $hpgl->Gb2Gff();
diag("Do I get the expected number of sequences, nucleotides, and features?");
ok($gff->{num_sequences} eq '1');
ok($gff->{total_nt} eq '1838562');
ok($gff->{num_features} eq '3858');
diag("Do the output files have the expected number of characters/lines?");

diag("Check the gff file of all features.");
my $mgas_5005_all = qx"wc mgas_5005_all.gff";
my $expected_mgas_5005_all = qq"   3862   38462 1211002 mgas_5005_all.gff\n";
ok($mgas_5005_all eq $expected_mgas_5005_all);

diag("Check the fasta file of CDS features.");
my $mgas_5005_cds = qx"wc mgas_5005_cds.fasta";
my $expected_mgas_5005_cds = qq"  29261   29261 1647918 mgas_5005_cds.fasta\n";
ok($mgas_5005_cds eq $expected_mgas_5005_cds);

diag("Check the gff file of CDS features.");
my $mgas_5005_cds_gff = qx"wc mgas_5005_cds.gff";
my $expected_mgas_5005_cds_gff = qq"  1842  20186 977161 mgas_5005_cds.gff\n";
ok($mgas_5005_cds_gff eq $expected_mgas_5005_cds_gff);

diag("Check the fasta file of the entire genome.");
my $mgas_5005_fasta = qx"wc mgas_5005.fasta";
my $expected_mgas_5005_fasta = qq"  30644   30649 1869265 mgas_5005.fasta\n";
ok($mgas_5005_fasta eq $expected_mgas_5005_fasta);

diag("Check the gff file of the entire genome.");
my $mgas_5005_gff = qx"wc mgas_5005.gff";
my $expected_mgas_5005_gff = qq"  1927  17336 221048 mgas_5005.gff\n";
ok($mgas_5005_gff eq $expected_mgas_5005_gff);

diag("Check the gff file of the interCDS regions.");
my $mgas_5005_intercds_gff = qx"wc mgas_5005_interCDS.gff";
my $expected_mgas_5005_intercds_gff = qq"  1842  20186 977155 mgas_5005_interCDS.gff\n";
ok($mgas_5005_intercds_gff eq $expected_mgas_5005_intercds_gff);

diag("Check the peptide fasta file.");
my $mgas_5005_peptide_fasta = qx"wc mgas_5005_pep.fasta";
my $expected_mgas_5005_peptide_fasta = qq" 11580  11580 568315 mgas_5005_pep.fasta\n";
ok($mgas_5005_peptide_fasta eq $expected_mgas_5005_peptide_fasta);

diag("Can I clean up the mess?");
ok(unlink('mgas_5005_all.gff') and unlink('mgas_5005_cds.fasta') and unlink('mgas_5005_cds.gff') and unlink('mgas_5005.fasta') and unlink('mgas_5005.gb') and unlink('mgas_5005.gb.xz') and unlink('mgas_5005.gff') and unlink('mgas_5005_interCDS.gff') and unlink('mgas_5005_pep.fasta'));

##my($old, $new) = String::Diff::diff($expected_mgas_5005_all, $mgas_5005_all);
##print STDERR "$old\n";
##print STDERR "---\n";
##print STDERR "$new\n";
##print STDERR "---\n";

