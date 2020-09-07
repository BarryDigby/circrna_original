# circRNA
Nextflow pipeline to scan high throughput sequencing data for circRNAs, perform differential expression and generate circRNA-miRNA targets

# TOOLS
The container hosts the following circRNA tools:

- CIRIquant (CIRI2)
- CIRCexplorer2
- MapSplice
- find_circ
- UROBORUS

## Gotchas
Thank you Simone Coughlan for the catch :) 
```
*--params.star_index:* Must put provide full path 

*--params.bwa_index:* "/data/bdigby/grch38/index/bwa"

*--params.hisat2_index:* "/data/bdigby/grch38/index/hisat2"

*--params.bowtie2_index:* "/data/bdigby/grch38/index/bowtie2/*" (need to point to files to collect for find_circ)

*--params.bowtie_index:* "/data/bdigby/grch38/index/bowtie" just path for mapsplice 
```
