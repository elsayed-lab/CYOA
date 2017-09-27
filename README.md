Adventure: Choose your own Adventure in Sequence processing.
============================================================

Over the few years I have been playing in Dr. El-Sayed's lab, I found myself copy-pasting snippets
of Perl for common tasks like sequence alignment, setting up PBS jobs, counting reads, whatever.
Eventually, having a large number of these repetitive scripts came to seem stupid.  Therefore I
combined them into a single Perl package and attempted to give them a single interface: cyoa.

# Rebuilding the distribution

> dzil build

# Installation

> git clone git@github.com:abelew/CYOA.git && cd CYOA && perl Makefile.PL && make && make install

If you have some free time:
> make test

That will test (currently) the ability to:
* load HPGL.pm
* look for perl requirements
* run fastqc on a small dataset
* run trimomatic on the same data
* run bowtie, convert the sam->bam, and count reads of phix
* test parsing of a tritryp text file
* run tophat for phix
* run bwa for phix
* perform a fasta36 search for some leishmania data
* clean up the mess

Final installation:
> make install

# Usage

I have recently been using it via a small script 'cyoa' in bin/, for example:

> cyoa

brings up a readline menu interface asking what you want to do.

> cyoa --query lmajor_cds.fasta --library lmajor_cds.fasta --task align --method fastasplit

Splits the set of annotated coding sequences from lmajor into 200 pieces and aligns them against the
library of lmajor_cds sequences.  Upon completion, it counts the number of hits and therefore provides a
rough metric of the set of multi-copy genes.

> cyoa --input test-forward.fastq:test-reverse.fastq --task rnaseq --method fastqc

Runs fastqc

> cyoa --input test-forward.fastq:test-reverse.fastq --task rnaseq --method tophat --species lmajor

Runs tophat using some generic options with databases specific for L.major.

> cyoa --query test.fasta --task alignment --method blastsplit --blast_tool blastp --library nr

Splits test.fasta into a bunch of pieces, aligns them against nr using blastp

> cyoa --input test.fastq --task tnseq --method sort --indexfile indexes.txt

Sorts test.fastq using a set of index->samples in indexes.txt

> cyoa --input test_forward.fastq:test_reverse.fastq --task pipeline --method priboseq --species lmajor

Performs a ribosome profiling pipeline of steps to fastqc the data,
trim it, graph the qualities etc, rRNA search the data, align the
non-rRNA against the lmajor genome, convert to bam, and count. It
might also try and count the A/P/E sites, I forget if I added that.
crap I have to go collect rna.
