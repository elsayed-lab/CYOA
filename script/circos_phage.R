#!/usr/bin/env Rscript
library(hpgltools)
startdir <- "/bio/getnet/preprocessing/test_assemblies/EPr2v2/outputs/23mergeannot"
input_gbk <- "EPr2v2.gbk"
input_fsa <- "EPr2v2.fsa"
setwd(startdir)
annot <- as.data.frame(sm(load_genbank_annotations(file=input)[["cds"]]))
annot[["seqnames"]] <- "EPr2v2"

## start, end, width, strand
circos_cfg <- circos_prefix(annot, name="phage", id_column="gene_id")
circos_kary <- circos_karyotype(circos_cfg, fasta=input_fsa)
circos_plus <- circos_plus_minus(circos_cfg)
circos_suffix <- circos_suffix(circos_cfg)
circos_made <- sm(circos_make(circos_cfg, target = "phage"))
