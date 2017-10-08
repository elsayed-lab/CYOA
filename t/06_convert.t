# -*-Perl-*-
use Test::More qw"no_plan";
use Bio::Adventure;
use File::Copy qw"cp";
use Archive::Extract;
use String::Diff qw"diff";

my $cyoa = Bio::Adventure->new();

ok(cp('share/genome/mgas_5005.gb.xz', 'mgas_5005.gb.xz'),
   'Copy data files.');

my $ae = Archive::Extract->new(archive => 'mgas_5005.gb.xz');
ok($ae->extract(),
   '02 extract.');

my $gff = Bio::Adventure::Convert::Gb2Gff($cyoa, input => qq'mgas_5005.gb');
ok($gff->{num_sequences} eq '1',
   '03 number sequences.');
ok($gff->{total_nt} eq '1838562',
   '04 number nucleotides.');
ok($gff->{num_features} eq '3858',
   '05 number features.');

my ($expected, $actual);
$actual = qx"wc mgas_5005_all.gff";
$expected = qq"   3862   38462 1211002 mgas_5005_all.gff\n";
ok($actual eq $expected,
   'wc the all gff.');

$actual = qx"wc mgas_5005_cds.fasta";
$expected = qq"  29261   39835 1727732 mgas_5005_cds.fasta\n";
ok($actual eq $expected,
   'wc the cds fasta.');

$actual = qx"wc mgas_5005_cds.gff";
$expected = qq"  1842  20186 977161 mgas_5005_cds.gff\n";
ok($actual eq $expected,
   'wc the cds gff.');

$actual = qx"wc mgas_5005.fasta";
$expected = qq"  30644   30649 1869265 mgas_5005.fasta\n";
ok($actual eq $expected,
   'wc the fasta.');

$actual = qx"wc mgas_5005.gff";
$expected = qq"  1927  17336 221048 mgas_5005.gff\n";
ok($actual eq $expected,
   'wc the gff.');

$actual = qx"wc mgas_5005_interCDS.gff";
$expected = qq"  1842  20186 977155 mgas_5005_interCDS.gff\n";
ok($actual eq $expected,
   'wc the intercds gff.');

$actual = qx"wc mgas_5005_pep.fasta";
$expected = qq" 11580  22154 648129 mgas_5005_pep.fasta\n";
ok($actual eq $expected,
   'wc the peptide fasta.');

my @files = ('mgas_5005_cds.gff', 'mgas_5005.gb.xz', 'mgas_5005_pep.fasta', 'mgas_5005_all.gff',
             'mgas_5005.fasta', 'mgas_5005.gff', 'mgas_5005_cds.fasta', 'mgas_5005.gb',
             'mgas_5005_interCDS.gff');

for my $t (@files) {
    ok(unlink($t),
       "Delete $t.");
}
