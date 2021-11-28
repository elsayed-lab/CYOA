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
mkdir($new);
chdir($new);

make_path('genome/indexes'); ## Make a directory for the phix indexes.
my $input_fasta = dist_file('Bio-Adventure', 'genome/phix.fasta');
if (!-r 'genome/phix.fasta') {
    ok(cp($input_fasta, 'genome/phix.fasta'), 'Copying phix fasta file.');
}
my $cyoa = Bio::Adventure->new(cluster => 0, basedir => cwd());
my $jelly = $cyoa->Bio::Adventure::Count::Jellyfish(
    input => $input_fasta,
    jprefix => 40,);

my $first_job_outputs = $jelly->{compression}->{output};
my @output_files = split(/:/, $first_job_outputs);
for my $o (@output_files) {
    ok (-f $o, "Output file create: ${o}");
}
