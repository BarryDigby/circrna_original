# circseq
Nextflow pipeline to scan high throughput sequencing data for circRNAs, perform differential expression and generate circRNA-miRNA targets.

# Notes
The pipeline uses RNA-Seq data to scan for the presence of circRNAs using a combination of `CIRCexplorer2`, `CIRIquant`, `DCC`, `find_circ`, `circrna_finder`, `mapsplice`. There is an option to use all 6 tools `--tool combine` or to use a tool on its own `--tool CIRIquant`. Currently developing the tool using a small test dataset, which works for the combine flag however using just one tool on the toy dataset can result in very few circRNAs being called which is an issue for DESeq2.

***

Differential expression -- the contrast of interest must be labelled as `condition` as this is hardcoded into the automated DE analysis script. The reference level, i.e the wild type condition to compare against (tumor vs **normal**, treated vs **untreated**) **must be called normal**. 

For example:

| sample   	| condition 	|
|----------	|-----------	|
| rep1_trt 	| treated   	|
| rep2_trt 	| treated   	|
| rep3_trt 	| treated   	|
| rep1     	| normal    	|
| rep2     	| normal    	|
| rep3     	| normal    	|

***

# TO-DO
Facilitate single end data. 
Incorporate 'stranded' datasets.

# TOOLS
The container hosts the following circRNA tools:

- CIRIquant (CIRI2)
- CIRCexplorer2
- circrna_finder
- MapSplice
- find_circ
- DCC
- UROBORUS **still under testing

## Gotchas
Thank you Simone Coughlan for the catch :) 
```
*--params.star_index:* "/data/bdigby/grch38/index/star_index"

*--params.bwa_index:* "/data/bdigby/grch38/index/bwa"

*--params.hisat2_index:* "/data/bdigby/grch38/index/hisat2"

*--params.bowtie2_index:* "/data/bdigby/grch38/index/bowtie2/*" (need to point to files to collect for find_circ)

*--params.bowtie_index:* "/data/bdigby/grch38/index/bowtie/*" 
```

