#!/usr/bin/env Rscript

dir = "."
samples = read.table(file.path(dir, 'samples.csv'))
files = file.path(dir, samples[,1])

ciriq_index <- grep("ciriquant", files)
ciriquant <- read.table(files[ciriq_index], header = F, sep ="\t")

circrna_finder_index <- grep("circrnafinder", files)
circrna_finder <- read.table(files[circrna_finder_index], header = F, sep="\t")

circexplorer2_index = grep("circexplorer", files)
circexplorer2 <- read.table(files[circexplorer2_index], header = F, sep ="\t")

dcc_index = grep("dcc", files)
dcc <- read.table(files[dcc_index], header = F, sep="\t")

mapsplice_index = grep("mapsplice", files)
mapsplice <- read.table(files[mapsplice_index], header = F, sep="\t")

find_circ_index <- grep("find_circ", files)
find_circ <- read.table(files[find_circ_index], header = F, sep="\t")

mat <- rbind(ciriquant, circexplorer2, circrna_finder, find_circ, dcc, mapsplice)
mat$id <- with(mat, paste0(V1, sep="_", V2, sep="_", V3, sep="_", V4))
mat1 <- mat[order(mat[,6], -abs(mat[,5])),]
mat1 <- mat1[!duplicated(mat1$id),]

mat1 <- mat1[,1:5]

write.table(mat1, "combined_counts.bed", sep="\t", row.names = F, col.names = F, quote = F)
