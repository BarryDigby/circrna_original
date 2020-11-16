# nf-circrna
Nextflow pipeline to scan high throughput sequencing data for circRNAs, perform differential expression and generate circRNA-miRNA targets.

pipeline currently undergoing linting at https://github.com/BarryDigby/test. Test dataset currently being tested. 
Once `BarryDigby/test` passes `nf-core` specifications, it will be pushed to this repo (or renamed `circrna`). 

Completed tests:

- [ ] nextflow_schema.json. Added pipeline parameters but need to add all relevant information for each flag. 
- [x] environment.yml file pinned version numbers. 
- [ ] test dataset

**Test Data:** circRNA data for Chromosome 1 in fastq files (1 set of fastq pairs). map the fastq files, split chromosome 'in half' and use for DE (must have 3 replicates). Fastq files are suitably small ~2Mb. Figure out how to store the data on github + bwa + hisat2 index files for minimal test run. use ciriquant for the analysis. do not use all 6 tools for test dataset.

Simulated reads generated using: 

```
#!/usr/bin/bash

list=("normal_rep1" "normal_rep2" "normal_rep3")

for i in ${list[@]}; do 

	perl CIRI_simulator.pl -O ${i} -G arm2.gtf -C 10 -LC 10 -R 2 -LR 2 -L 100 -E 1 -D ./ -CHR1 1 -M 320 -M2 550 -PM 15 -S 70 -S2 70 -SE 0 -PSI 0
done

list1=("tumor_rep1" "tumor_rep2" "tumor_rep3")

for i in ${list1[@]}; do

	perl CIRI_simulator.pl -O ${i} -G arm1.gtf -C 10 -LC 10 -R 2 -LR 2 -L 100 -E 1 -D ./ -CHR1 1 -M 320 -M2 550 -PM 15 -S 70 -S2 70 -SE 0 -PSI 0

done

mv *.fq fastq/

gzip fastq/*.*
```

After running the script, rep2 for each sample was deleted and duplicated using rep1. This was to satisfy DESeq2 which complains about too many zeros if each replicate is completely randomn and share no circRNAs / RNA. 


***

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

