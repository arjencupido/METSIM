library(tidyverse)
load(file = 'C:/Users/arjen/Documents/research/UCLA/METSIM/data/R/FINAL DATA/cleaned datasets/metsim_microbiome_pso_updated.rda')
setwd('C:/Users/arjen/Documents/research/UCLA/METSIM/data/R/FINAL DATA/cleaned datasets/')
ps3
otus <- data.frame((otu_table(ps3)))
otus$Sample_ID <- rownames(otus)
otus <- subset(otus, select=c(23920,1:23919))
write_delim(otus, file= "OTUS.txt")


samples <- data.frame(sample_data(ps3))
samples <- subset(samples, select=c(456,1:455))
write_delim(samples, file= "SAMPLES.txt")

taxatable <- data.frame(tax_table(ps3))
taxatable$ASV <- rownames(taxatable)
taxatable <- subset(taxatable, select=c(8, 1:7))
write_delim(taxatable, file = "TAXTABLE.txt")


refseqences <- data.frame(refseq(ps3))
refseqences$ASV <- rownames(refseqences)
write_delim(refseqences, file = "REFSEQUENCES.txt")
