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

*--params.bwa_index:*  Must put in quotes! probaby because of the glob pattern, safe to put both these paths in "" 
```
