BEGIN {
    use Test::More qw"no_plan";
    use HPGL;
    use File::Path qw"remove_tree";
}

diag("Can I clean up the mess?");
my $hpgl = new HPGL(input => qq"test_forward.fastq.gz", pbs => 0);
ok($hpgl->Cleanup());
ok(remove_tree(scripts));
ok(remove_tree("t/data/genome/indexes"));
