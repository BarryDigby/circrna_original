#!/usr/bin/Rscript

library(biomaRt)

x <- read.table("merged_reports.txt", sep="\t", header=T)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
results <- getBM(attributes=c("entrezgene_description", "hgnc_symbol"), mart=mart)
colnames(results) <- c("Description", "Parent_Gene")

x <- merge(x, results, by="Parent_Gene")
x <- subset(x, select=c(circRNA_ID, Parent_Gene, Description, Strand, Log2FC, pvalue, Adjusted_pvalue))

write.table(x, "DE_circRNA_Report.txt", quote=F, sep="\t", row.names=F)
