#!/usr/bin/env Rscript
library(hpgltools)
startdir <- "/home/trey/CYOA/test_output/"
input_tsv <- "outputs/26mergeannot/test_output.tsv"
input_fsa <- "outputs/27mergeannot/test_output_stripped.fsa"
setwd(startdir)

annot <- as.data.frame(readr::read_tsv(file=input_tsv))
annot[["cog"]] <- "Z"

## start, end, width, strand
circos_cfg <- circos_prefix(annot, name="phage", chr_column="contig", id_column="locus_tag",
                            cog_column = "cog")
circos_kary <- circos_karyotype(circos_cfg, fasta=input_fsa)
circos_plus <- circos_plus_minus(circos_cfg)
circos_suffix <- circos_suffix(circos_cfg)
circos_made <- sm(circos_make(circos_cfg, target = "phage"))
