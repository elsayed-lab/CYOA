# -*-Perl-*-
BEGIN {
    use Test::More qw"no_plan";
    use HPGL;
    use File::Copy qw"cp";
    use Archive::Extract;
    use String::Diff qw( diff_fully diff diff_merge diff_regexp );
    use LWP;
}

my $hpgl = new HPGL;
diag("Downloading the newest tritrypdb text file for lmajor.");
my $txt_file = $hpgl->TriTryp_DL_Text();
print STDERR "The actual file downloaded is: $txt_file\n";
ok($txt_file =~ m/^TriTrypDB.*Gene\.txt$/);
$hpgl->{input} = $txt_file;
$hpgl->{species} = 'lmajor';
diag("Can I translate the TriTrypDB file to a simplified table?");
my $counts = $hpgl->TriTryp2Text(ortho => 1, go => 1);
ok($counts);
diag("Did I get more than 9370 entries in the master file?");
print STDERR "HOW MANY? $counts->{master_lines}\n";
ok($counts->{master_lines} >= 9370);
diag("Did I get >= 325000 lines in the orthology table?");
ok($counts->{ortho_lines} >= 325000);
diag("Did I get >= 20000 lines in the go table?");
ok($counts->{go_lines} >= 20000);
diag("Can I clean up the mess?");
ok(unlink("lmajor_master.tab") and unlink("lmajor_ortho.txt") and unlink("lmajor_go.txt") and unlink("lmajor_go.gaf"));
my $tritryp_remove = qx"rm TriTrypDB*.txt";

