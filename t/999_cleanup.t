BEGIN {
    use Test::More qw"no_plan";
    use CYOA;
    use File::Path qw"remove_tree";
}

diag("Can I clean up the mess?");
my $cyoa = new CYOA(input => qq"test_forward.fastq.gz", pbs => 0);
ok($cyoa->Cleanup());
#ok(remove_tree('scripts'));
#ok(remove_tree("t/data/genome/indexes"));
