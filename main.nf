#!/usr/bin/env nextflow

/*
================================================================================
                                  deNovo-circRNA
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

if(params.outDir){
  params.refDir = "${params.outDir}"
}
if(!params.outDir){
  params.refDir = "${workflow.launchDir}/${params.version}"
}


stepList = defineStepList()
step = params.step ? params.step.toLowerCase().replaceAll('-', '').replaceAll('_', '') : ''

if (step.contains(',')) exit 1, 'You can choose only one step, see --help for more information'
if (!checkParameterExistence(step, stepList)) exit 1, "Unknown step ${step}, see --help for more information"


params.threads = 8
params.readlength = 44


/*
================================================================================
                                  Preprocessing
================================================================================
*/

process Download_Genome {

        publishDir "$params.refDir/Reference", mode: 'copy'

        output:
        file('*.fa') into genome_ch
        file('*.txt') into txt_ann_ch
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

process STAR_Index {

        publishDir "$params.refDir/Index/STAR", mode: 'copy'
        
        input:
        file(a) from genome_ch
        file(gtf) from gtf_ch

        output:
        file ('*')

        script:
        """
        mkdir -p "$params.refDir/Index/STAR"

        STAR --runMode genomeGenerate \
        --genomeDir "$params.refDir/Index/STAR" \
        --genomeFastaFiles $a \
        --sjdbGTFfile $gtf \
        --sjdbOverhang ${params.readlength} \
        --runThreadN ${params.threads} \
        --genomeSAindexNbases 12
        """
}

