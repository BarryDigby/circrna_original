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
params.hisat2_index = ''
params.bowtie_index = ''
params.bowtie2_index = ''
params.mapsplice_ref = ''
params.ciriquant_yml = ''

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
        bwa index -a bwtsw $fasta -p ${fasta.baseName}
        """
}

ch_bwa_index = params.bwa_index ? Channel.value(params.bwa_index) : bwa_path 

ch_bwa_index.view()

process hisat2_index{

        publishDir "$params.outdir/index/hisat2", mode: 'copy'
        
        input:
            file(fasta) from ch_fasta
         
        output:
            file("${fasta.baseName}.*.ht2") into hisat2_built
            val("$params.outdir/index/hisat2") into hisat2_path
            
        when: !(params.hisat2_index) && 'ciriquant' in tool
        
        script:
        """
        hisat2-build $fasta ${fasta.baseName}
        """
}

ch_hisat2_index = params.hisat2_index ? Channel.value(params.hisat2_index) : hisat2_path
ch_hisat2_index.view()

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

process bowtie_index{

        publishDir "$params.outdir/index/bowtie", mode:'copy'
        
        input:
            file(fasta) from ch_fasta
            
        output:
            file ('"${fasta.baseName.*') into bowtie_built
        
        when: !(params.bowtie_index) && 'mapsplice' in tool

        script:
        """
        bowtie-build $fasta ${fasta.baseName}
        """
}

ch_bowtie_index = params.bowtie_index ? Channel.value(file(params.bowtie_index)) : bowtie_built
ch_bowtie_index.view() 
 
process bowtie2_index{

        publishDir "$params.outdir/index/bowtie2", mode:'copy'
        
        input:
            file(fasta) from ch_fasta
            
        output:
            file ('"${fasta.baseName.*') into bowtie2_built
        
        when: !(params.bowtie2_index) && ('find_circ' in tool || 'uroborus' in tool)

        script:
        """
        bowtie2-build $fasta ${fasta.baseName}
        """
}

ch_bowtie2_index = params.bowtie2_index ? Channel.value(file(params.bowtie2_index)) : bowtie2_built
ch_bowtie2_index.view()


/*
 * Step 3:
 * miscellaneous circRNA tool requirements
 */
 
 
process split_fasta{

        publishDir "$params.outdir/index/mapsplice", mode:'copy'
        
        input:
            file(fasta) from ch_fasta
            
        output:
             file("*.fa") into split_fasta
             
        when 'mapsplice' in tool
        
        shell:
        '''
        awk '$0 ~ "^>" { match($1, /^>([^:]+)/, id); filename=id[1]} {print >> filename".fa"}' !{fasta}
        '''
}

ch_mapsplice_ref = params.mapsplice_ref ? Channel.value(params.mapsplice_ref) : split_fasta
ch_mapsplice_ref.view()

process ciriquant_yml{

      publishDir "$params.outdir", mode:'copy'
      
      input:
          val(gencode_gtf_path) from ch_gencode_gtf
          val(fasta_path) from ch_fasta
          val(bwa_path) from ch_bwa_index
          val(hisat2_path) from ch_hisat2_index
          
      output:
          file("travis.yml") into travis_built
          
      when: !(params.ciriquant_yml) && 'ciriquant' in tool
      
      script:
      """
      export bwa=`whereis bwa | cut -f2 -d':'`
      export hisat2=`whereis hisat2 | cut -f2 -d':'`
      export stringtie=`whereis stringtie | cut -f2 -d':'`
      export samtools=`whereis samtools | cut -f2 -d':'`

      touch travis.yml
      printf "name: ciriquant\n\
      tools:\n\
       bwa: \$bwa\n\
       hisat2: \$hisat2\n\
       stringtie: \$stringtie\n\
       samtools: \$samtools\n\n\
      reference:\n\
       fasta: ${fasta_path}\n\
       gtf: ${gencode_gtf_path}\n\
       bwa_index: ${bwa_path}\n\
       hisat_index: ${hisat2_path}" >> travis.yml
      """
}

ch_ciriquant_yml = params.ciriquant_yml ? Channel.value(file(params.ciriquant_yml)) : travis_built
      
      


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
