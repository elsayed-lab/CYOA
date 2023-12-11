# -*-Perl-*-
use strict;
use Test::More qw"no_plan";
use Bio::Adventure;
use Cwd;
use File::Basename qw"basename dirname";
use File::Copy qw"cp mv";
use File::Path qw"remove_tree make_path rmtree";
use File::ShareDir qw"dist_file module_dir dist_dir";
use String::Diff qw"diff";
use Test::File::ShareDir::Dist { 'Bio-Adventure' => 'share/' };

my $start = getcwd();
my $new = 'test_output';
mkdir($new);
chdir($new);
my $start_dir = dist_dir('Bio-Adventure');

my $groupa_strep = qq"${start_dir}/genome/mgas_5005.gb.xz";
my $groupa_local = basename($groupa_strep);
ok(cp($groupa_strep, $groupa_local), 'Copy the group A file.');
ok(-r $groupa_local, 'Group A streptococcus file exists.');
my $groupb_strep = qq"${start_dir}/genome/sagalactiae_cjb111.gb.xz";
my $groupb_local = basename($groupb_strep);
ok(cp($groupb_strep, $groupb_local), 'Copy the group B file.');
ok(-r $groupb_local, 'Group B streptococcus file exists.');
my $cyoa = Bio::Adventure->new(
    basedir => cwd(),
    libdir => cwd());

## First extract the amino acid sequences of our two streptococci genbank files.
my $groupa_convert = $cyoa->Bio::Adventure::Convert::Gb2Gff(
    input => $groupa_local, output_dir => '.',);
my $status = $cyoa->Wait(job => $groupa_convert);
ok($groupa_convert, 'Converted the group A strep genbank to fasta/gff/etc.');
ok(-r $groupa_convert->{output_pep_fasta}, 'Created output peptide Group A fasta file.');
my $groupb_convert = $cyoa->Bio::Adventure::Convert::Gb2Gff(
    input => $groupb_local, output_dir => '.',);
$status = $cyoa->Wait(job => $groupa_convert);
ok($groupb_convert, 'Converted the group B strep genbank to fasta/gff/etc.');
ok(-r $groupb_convert->{output_pep_fasta}, 'Created output peptide Group B fasta file.');

## Now invoke orthofinder
my $orthos = $cyoa->Bio::Adventure::Align::OrthoFinder(
    input => qq"$groupa_convert->{output_pep_fasta}:$groupb_convert->{output_pep_fasta}");
$status = $cyoa->Wait(job => $orthos);

ok(-r $orthos->{named_out}, 'Created tsv of named orthologs.');
my $expected = qq!Orthogroup
"OG0000000"	"AAZ50854.1 GI_71852831 ; amino acid transport ATP-binding protein, AAZ51149.1 ftsE ; GI_71853126 ; cell division ATP-binding protein, AAZ51186.1 sagG ; GI_71853163 ; streptolysin S export ATP-binding protein, AAZ51426.1 srtF ; GI_71853403 ; lantibiotic transport ATP-binding protein, AAZ51473.1 proV ; GI_71853450 ; glycine betaine transport ATP-binding protein, AAZ51585.1 GI_71853562 ; ABC transporter ATP-binding protein, AAZ51601.1 GI_71853578 ; histidine transport ATP-binding protein, AAZ51695.1 glnQ.2 ; GI_71853672 ; glutamine transport ATP-binding protein, AAZ51855.1 artP ; GI_71853832 ; arginine transport ATP-binding protein, AAZ51980.1 GI_71853957 ; transporter, AAZ52136.1 GI_71854113 ; transporter, AAZ52139.1 GI_71854116 ; cobalt transport ATP-binding protein cbiO, AAZ52143.1 GI_71854120 ; ABC transporter ATP-binding protein, AAZ52247.1 salX ; GI_71854224 ; lantibiotic transport ATP-binding protein"	"AAZ50854.1 GI_71852831 ; amino acid transport ATP-binding protein"	""	"QOW75830.1 ATP-binding cassette domain-containing protein, QOW75889.1 ATP-binding cassette domain-containing protein, QOW75941.1 amino acid ABC transporter ATP-binding protein, QOW75969.1 ATP-binding cassette domain-containing protein, QOW75973.1 ABC transporter ATP-binding protein, QOW75989.1 ABC transporter ATP-binding protein, QOW76092.1 ABC transporter ATP-binding protein, QOW76121.1 ftsE ; cell division ATP-binding protein FtsE, QOW76148.1 amino acid ABC transporter ATP-binding protein, QOW76318.1 ABC transporter ATP-binding protein, QOW76332.1 amino acid ABC transporter ATP-binding protein, QOW76362.1 ABC transporter ATP-binding protein, QOW76375.1 betaine/proline/choline family ABC transporter ATP-binding protein, QOW76483.1 amino acid ABC transporter ATP-binding protein, QOW76495.1 ABC transporter ATP-binding protein, QOW76500.1 sugar ABC transporter ATP-binding protein, QOW76676.1 ABC transporter ATP-binding protein, QOW76764.1 ABC transporter ATP-binding protein, QOW77065.1 ABC transporter ATP-binding protein, QOW77115.1 ABC transporter ATP-binding protein, QOW77116.1 ABC transporter ATP-binding protein, QOW77166.1 ABC transporter ATP-binding protein, QOW77167.1 ABC transporter ATP-binding protein, QOW77187.1 ABC transporter ATP-binding protein, QOW77207.1 amino acid ABC transporter ATP-binding protein, QOW77578.1 ATP-binding cassette domain-containing protein, QOW77582.1 ABC transporter ATP-binding protein, QOW77612.1 ABC transporter ATP-binding protein, QOW77645.1 amino acid ABC transporter ATP-binding protein, QOW77719.1 ABC transporter ATP-binding protein"	"QOW75830.1 ATP-binding cassette domain-containing protein"	""
!;
my $actual = qx"head -n 2 $orthos->{named_out}";
unless(ok($expected eq $actual, 'Is the resulting set of named orthologs as expected?')) {
    my($old, $new) = diff($expected, $actual);
    diag("--Expected--\n${old}\n--Actual--\n${new}\n");
}
chdir($start);
