#!/usr/bin/env Rscript
library(hpgltools)
startdir <- "/home/trey/CYOA/test_output/"
input_tsv <- "outputs/26mergeannot/test_output.tsv"
input_fsa <- "outputs/27mergeannot/test_output_stripped.fsa"
setwd(startdir)

annot <- as.data.frame(readr::read_tsv(file=input_tsv))
annot[["cog"]] <- "Z"

init <- annot[, c("locus_tag", "start", "end")]
colnames(init) <- c("locus_tag", "x", "y")

library(circlize)
circos.par("track.height" = 0.1)
circos.initializeWithIdeogram(chromosome.index = "chr1")

circos.genomicInitialize(init)
circos.track(ylim = c(0, 1))
