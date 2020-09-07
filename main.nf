#!/usr/bin/env nextflow

/*
================================================================================
                                circRNA analysis
================================================================================
Started August 2020
--------------------------------------------------------------------------------
Description:
  (To my knowledge) the first circRNA tool to scan RNA-Seq data for circRNAs and
  conduct circRNA-miRNA predictions
--------------------------------------------------------------------------------
 @Homepage
 https://github.com/BarryDigby/circRNA
 --------------------------------------------------------------------------------
 @Documentation
 Work in progress
--------------------------------------------------------------------------------
*/

/*
 * Parameters
 */
params.outdir ='./'
params.fasta = ''
params.gencode_gtf = ''
params.gene_annotation = ''
params.version = ''
params.tool = ''
params.fasta_fai = ''
params.bwa_index = ''
params.star_index = ''


toolList = defineToolList()
tool = params.tool ? params.tool.toLowerCase().replaceAll('-', '').replaceAll('_', '') : ''

/*
 * Step 1:
 * Download Gencode Reference Files
 */
 
 
process download_genome {
        
        publishDir "$params.outdir/reference", mode: 'copy'

        output:
        file('*.fa') into fasta_downloaded
        file('*.txt') into gene_annotation_created
        file('*.gtf') into gencode_gtf_downloaded
        
        when: !(params.fasta) && !(params.gencode_gtf) && !(params.gene_annotation)
        
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

// ternary operators: result = condition ? value_if_true : value_if_false
ch_fasta = params.fasta ? Channel.value(file(params.fasta)) : fasta_downloaded
ch_gene_annotation = params.gene_annotation ? Channel.value(file(params.gene_annotation)) : gene_annotation_created
ch_gencode_gtf = params.gencode_gtf ? Channel.value(file(params.gencode_gtf)) : gencode_gtf_downloaded


ch_fasta.view()
ch_gene_annotation.view()
ch_gencode_gtf.view()

/*
 * Step 2:
 * Create Genome Index files
 * Samtools needed for miRNA prediction - must be made
 */

process samtools_index{

        publishDir "$params.outdir/reference", mode:'copy'
        
        input:
            file(fasta) from ch_fasta
            
        output:
            file("${fasta}.fai") into fasta_fai_built
        
        when: !(params.fasta_fai)
        
        script:
        """
        samtools faidx $fasta
        """
}
        
ch_fai = params.fasta_fai ? Channel.value(file(params.fasta_fai)) : fasta_fai_built

process bwa_index{

        publishDir "$params.outdir/index/bwa", mode:'copy'
        
        input:
            file(fasta) from ch_fasta
            
        output:
            file("${fasta.baseName}.*") into bwa_built
            val("$params.outdir/index/bwa") into bwa_path
        
        when: !(params.bwa_index) && 'ciriquant' in tool 
        
        script:
        """
        bwa index ${fasta} -p ${fasta.baseName}
        """
}

ch_bwa_index = params.bwa_index ? Channel.value(params.bwa_index) : bwa_path 

ch_bwa_index.view()

process star_index{

        publishDir "$params.outdir/index", mode:'copy'
          
        input:
            file(fasta) from ch_fasta
            file(gtf) from ch_gencode_gtf
              
        output:
            file("star_index") into star_built
              
        when: !(params.star_index) && 'circexplorer2' in tool
        
        script:
        """
        mkdir star_index
          
        STAR \
        --runMode genomeGenerate \
        --runThreadN 8 \
        --sjdbGTFfile $gtf \
        --genomeDir star_index/ \
        --genomeFastaFiles $fasta
        """
}

ch_star_index = params.star_index ? Channel.value(file(params.star_index)) : star_built
ch_star_index.view()







// Define list of available tools
def defineToolList() {
    return [
        'ciriquant',
        'circexplorer2',
        'find_circ',
        'mapsplice',
        'uroborus'
        ]
}
