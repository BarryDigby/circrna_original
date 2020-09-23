#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

get_upregulated <- function(df){
    key <- intersect(rownames(df)[which(df$log2FoldChange>=1)],
              rownames(df)[which(df$pvalue<=0.05)])
    
    results <- as.data.frame((df)[which(rownames(df) %in% key),])
    return(results)
  }

get_downregulated <- function(df){
  key <- intersect(rownames(df)[which(df$log2FoldChange<=-1)],
            rownames(df)[which(df$pvalue<=0.05)])
  
  results <- as.data.frame((df)[which(rownames(df) %in% key),])
  return(results)
}

library(DESeq2)
library(biomaRt)

## read in RNA-Seq and calc size factors
countData <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id", check.names=F))
colData <- read.csv(args[1], sep="\t", row.names=1)

dds <- DESeqDataSetFromMatrix(countData, colData, design = ~ condition)
dds <- DESeq(dds)

size_factors <- sizeFactors(dds)


## read in circRNA mat
circ <- read.table("circRNA_matrix.txt", sep ="\t", header = T)

## merge coordinates for row names
circ$circ <- with(circ, paste0(Chr, sep="_", Start, sep="_", Stop, sep="_", Strand))
rownames(circ) <- circ$circ
circ <- subset(circ, select=-c(Chr, Start, Stop, Strand, circ))

## make sure colnames == rownames(pheno)
reorder <- rownames(colData)
circ <- circ[, reorder]

## now have circ mat

dds_circ <- DESeqDataSetFromMatrix(circ, colData, design = ~ condition)
sizeFactors(dds_circ) = c(size_factors)
dds_circ <- DESeq(dds_circ)
res <- results(dds_circ)

up_regulated <- get_upregutaled(res)
down_regulated <- get_downregulated(res)

write.table(up_regulated, "up_reg.txt", sep="\t", quote=F)
write.table(down_regulated, "down_reg.txt", sep="\t", quote=F)

