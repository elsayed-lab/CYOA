# -*-Perl-*-
use Test::More qw"no_plan";
use Bio::Adventure;
use File::Path qw"remove_tree";
use File::Copy qw"cp";
use String::Diff qw( diff_fully diff diff_merge diff_regexp );

my $cyoa = new Bio::Adventure(pbs => 0, species => 'phix', libdir => 't/data');

ok(cp("t/data/genome/phix.fasta", "phix.fasta"),
   'Copied phix.fasta');
ok(cp("t/data/genome/phix.gff", "phix.gff"),
   'Copied phix.gff');

ok(Bio::Adventure::Convert::Gff2Fasta($cyoa, genome => 'phix.fasta', gff => 'phix.gff'),
   'Run gff2fasta.');

ok(Bio::Adventure::Align_Fasta::Split_Align_Fasta($cyoa, query => 'phix_cds_nt.fasta',
                                                  library => 'phix_cds_nt.fasta', number => 1,
                                                  parse => 0),
   'Run Split_Align_Fasta.');
