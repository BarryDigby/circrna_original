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
*/

//help information
params.help = null
if (params.help) {
    log.info ''
    log.info print_purple('------------------------------------------------------------------------')
    log.info "circRNA_tool: a Nextflow-based circRNA analysis pipeline"
    log.info "circRNA_tool utilises NGS tools to scan RNA-Seq data for the presence of circRNAs before conducting differential"
    log.info "expression analysis and circRNA - miRNA interaction prediciton. To run the pipeline, users need only install nextflow and "
    log.info "singularity as the tools are self contained within the container."
    log.info print_purple('------------------------------------------------------------------------')
    log.info ''
    log.info print_yellow('Usage: ')
    log.info print_yellow('    The typical command for running the pipeline is as follows (we do not recommend users passing configuration parameters through command line, please modify the config.file instead):\n') +
            print_purple('       Nextflow run LncRNAanalysisPipe.nf \n') +

            print_yellow('    General arguments:             Input and output setting\n') +
            print_cyan('      --inputdir <path>         ') + print_green('Path to input data(optional), current path default\n') +
            print_cyan('      --input_file              ') + print_green('Input file type. fasta (default) or bam ') +
            '\n' +
            print_yellow('    Options:                         General options for run this pipeline\n') +
            print_cyan('      --merged_gtf <gtffile>        ') + print_green('Start analysis with assemblies already produced and skip fastqc/alignment step, DEFAOUL NULL\n') +
            print_cyan('      --design <file>               ') + print_green('A flat file stored the experimental design information ( required when perform differential expression analysis)\n') +
            print_cyan('      --singleEnd                   ') + print_green('Reads type, True for single ended \n') +
            print_cyan('      --unstrand                    ') + print_green('RNA library construction strategy, specified for \'unstranded\' library \n') +
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

process download_genome {

        publishDir "$params.outdir/Reference", mode: 'copy'

        output:
        file('*.fa') into genome_ch
        file('*.txt') into txt_ref_ch
        file('*.gtf') into gtf_ch

        when: 'download' in step
        
        shell:
        if( params.version == 'GRCh37' )
          '''
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
          '''
        
        else
          '''
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
          '''
}

