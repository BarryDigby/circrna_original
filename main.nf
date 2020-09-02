#!/usr/bin/env nextflow

/*
================================================================================
                                circRNA analysis
================================================================================
Started August 2020
--------------------------------------------------------------------------------
Description:
  (To my knowledge) the first circRNA tool to scan RNA-Seq data for circRNAs,
  perform differential gene expression AND conduct circRNA-miRNA predictions
--------------------------------------------------------------------------------
 @Homepage
 https://github.com/BarryDigby/circRNA
 --------------------------------------------------------------------------------
 @Documentation
 Work in progress
--------------------------------------------------------------------------------
*/

/* 
================================================================================
                                Workflow Plan
================================================================================
--------------------------------------------------------------------------------
1. Downloading Reference Files
   Provide option for GRCh37 or GRCh38 from gencode
   Process GTF for CIRCexplorer2 compatability
   
   Parameters:
     params.version = GRCh37 | GRCh38
     params.fasta = 'path/to/genome.fa'
     params.gencode_annotation_gtf = 'path/to/genome.gtf'
     params.circexplorer2_ref = 'path/to/genome.txt'
--------------------------------------------------------------------------------
2. Create Genome Index
   Provide option for STAR, BWA, hisat2(single end, TO DO)
   
   Parameters:
    params.aligner = star, bwa 
    params.star_idx = 'star_index/path/'
    params.bwa_idx = 'bwa_index/path'
*/


//pre-defined functions for render command
//==============================================================================
ANSI_RESET = "\u001B[0m";
ANSI_BLACK = "\u001B[30m";
ANSI_RED = "\u001B[31m";
ANSI_GREEN = "\u001B[32m";
ANSI_YELLOW = "\u001B[33m";
ANSI_BLUE = "\u001B[34m";
ANSI_PURPLE = "\u001B[35m";
ANSI_CYAN = "\u001B[36m";
ANSI_WHITE = "\u001B[37m";

def print_red = {  str -> ANSI_RED + str + ANSI_RESET }
def print_black = {  str -> ANSI_BLACK + str + ANSI_RESET }
def print_green = {  str -> ANSI_GREEN + str + ANSI_RESET }
def print_yellow = {  str -> ANSI_YELLOW + str + ANSI_RESET }
def print_blue = {  str -> ANSI_BLUE + str + ANSI_RESET }
def print_cyan = {  str -> ANSI_CYAN + str + ANSI_RESET }
def print_purple = {  str -> ANSI_PURPLE + str + ANSI_RESET }
def print_white = {  str -> ANSI_WHITE + str + ANSI_RESET }

//help information
params.help = null
if (params.help) {
    log.info ''
    log.info print_purple('-------------------------------------------------------------------------')
    log.info "circRNA_tool: a Nextflow-based circRNA analysis pipeline"
    log.info "Utilises NGS tools to scan RNA-Seq data for the presence of circRNAs before conducting"
    log.info "differential expression analysis and circRNA - miRNA interaction prediciton." 
    log.info "To run the pipeline, users need only install nextflow and singularity as the tools are"
    log.info "self contained within the Docker container."
    log.info print_purple('-------------------------------------------------------------------------')
    log.info ''
    log.info print_yellow('Usage: ')
    log.info print_yellow('   The typical commands for running the pipeline are as follows:\n') +
            print_purple('       Nextflow run LncRNAanalysisPipe.nf \n') +

            print_yellow('    General arguments:             Input and output setting\n') +
            print_cyan('      --inputdir <path>         ') + print_green('Path to input data(optional), current path default\n') +
            print_cyan('      --input_file              ') + print_green('Input file type. fasta (default) or bam ') +
            '\n' +
            print_yellow('    References:                      If not specified in the configuration file or you wish to overwrite any of the references.\n') +
            print_cyan('      --fasta                       ') + print_green('Path to Fasta reference(required)\n') +
            print_cyan('      --gencode_annotation_gtf      ') + print_green('An annotation file from GENCODE database in GTF format (required)\n') +
            print_cyan('      --lncipedia_gtf               ') + print_green('An annotation file created for CIRCexplorer2 analysis (required)\n') 

            log.info '------------------------------------------------------------------------'
    log.info print_yellow('Contact information: b.digby1@nuigalway.ie')
    log.info '------------------------------------------------------------------------'
    exit 0
}

//Step 1
params.version = 'GRCh38'
//Step 2
params.star_overhang = '49' 



/*
 * Step 1: Downlaod Reference Files
 */

process download_genome {

        publishDir "$params.outdir/reference", mode: 'copy'

        output:
        file('*.fa') into fasta_downloaded
        file('*.txt') into gene_annotation_created
        file('*.gtf') into gencode_gtf_downloaded
        
        shell:
        if(params.version == 'GRCh37'){
          $/
          wget --no-check-certificate ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/gencode.v34lift37.annotation.gtf.gz
          wget --no-check-certificate ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz
          gunzip gencode.v34lift37.annotation.gtf.gz
          gunzip GRCh37.primary_assembly.genome.fa.gz
          mv gencode.v34.primary_assembly.annotation.gtf GRCh37.gtf
          mv GRCh37.primary_assembly.genome.fa GRCh37.fa.tmp
          sed 's/\s.*$//' GRCh37.fa.tmp > GRCh37.fa
          gtfToGenePred -genePredExt -geneNameAsName2 GRCh37.gtf GRCh37.genepred
          perl -alne '$"="\t";print "@F[11,0..9]"' GRCh37.genepred > GRCh37.txt
          rm GRCh37.fa.tmp
          /$
          
        }else if(params.version == 'development'){
          $/
          wget --no-check-certificate ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.primary_assembly.annotation.gtf.gz
          wget --no-check-certificate ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz         
          gunzip gencode.v34.primary_assembly.annotation.gtf.gz
          gunzip Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz
          mv gencode.v34.primary_assembly.annotation.gtf GRCh38.gtf
          mv Homo_sapiens.GRCh38.dna.chromosome.20.fa chr20.fa
          gtfToGenePred -genePredExt -geneNameAsName2 GRCh38.gtf GRCh38.genepred
          perl -alne '$"="\t";print "@F[11,0..9]"' GRCh38.genepred > GRCh38.txt          
          /$
          
        }else if(params.version == 'GRCh38'){
          $/
          wget --no-check-certificate ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.primary_assembly.annotation.gtf.gz
          wget --no-check-certificate ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz
          gunzip gencode.v34.primary_assembly.annotation.gtf.gz
          gunzip GRCh38.primary_assembly.genome.fa.gz
          mv gencode.v34.primary_assembly.annotation.gtf GRCh38.gtf
          mv GRCh38.primary_assembly.genome.fa GRCh38.fa.tmp
          sed 's/\s.*$//' GRCh38.fa.tmp > GRCh38.fa
          gtfToGenePred -genePredExt -geneNameAsName2 GRCh38.gtf GRCh38.genepred
          perl -alne '$"="\t";print "@F[11,0..9]"' GRCh38.genepred > GRCh38.txt
          rm GRCh38.fa.tmp
          /$
          }
}

ch_fasta = params.fasta ? Channel.value(file(params.fasta)) : fasta_downloaded
ch_gene_annotation = params.gene_annotation ? Channel.value(file(params.gene_annotation)) : gene_annotation_created
ch_gencode_gtf = params.gencode_gtf ? Channel.value(file(params.gencode_gtf)) : gencode_gtf_downloaded

/*
 * Step 2: Create Genome Index
 */ 
 
if(params.aligner == 'star' && params.star_index == false){
    process star_index {
    
          publishDir "params.outdir/index", mode:'copy'
          
          input:
              file(fatsa) from ch_fasta
              file(gtf) from ch_gencode_gtf
              
          output:
              file("star_index") into star_built
              
          script:
          """
          mkdir star_index
          
          STAR \
          --runMode genomeGenerate \
          --runThreadN 8 \
          --sjdbGTFfile $gtf \
          --sjdbOverhang $star_overhang \
          --genomeDir star_index/ \
          --genomeFastaFiles $fasta
          """
          }
} else if(params.aligner == 'bwa' && params.bwa_index == false){
    process bwa_index {
    
        publishDir "params.outdir/index/bwa", mode:'copy'
        
        input:
            file(fasta) from fasta_ch
            
        output:
            file("${fasta}.*") into bwa_built
            
        script:
        """
        bwa index ${fasta}
        """
        }
 }
