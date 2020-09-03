# circRNA
Nextflow pipeline to scan high throughput sequencing data for circRNAs, perform differential expression and generate circRNA-miRNA targets

### To Do
Currently I have tested bwa and bams. need to test bwa and fastq, and then both with STAR

Make sure every possible combination of flags works!! 

circRNA step should be straight forward. 

have scripts to pull sequences and make circRNA files.

### After the fact
Once you have the script at the stage where it conducts every step up to circRNA prediction, 
go back and add customised memory options to steps that need it (STAR). 
do this at intervals so it doesnt pile up on you 

## Gotchas

```
--params.star_index 
