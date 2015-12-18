# HPGL

Some simple Perl wrappers for (pre)processing our RNAseq data.

## Installation

> perl Makefile.PL

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

## Usage

I have recently been using it via a small script 'cyoa' in bin/, for example:
> cyoa
brings up a readline menu interface asking what you want to do.
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
