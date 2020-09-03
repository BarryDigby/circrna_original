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
   Provide option for STAR, BWA, hisat2(single end, TO DO list (:)
   
   Parameters:
    params.aligner = star, bwa 
    params.star_idx = 'star_index/path/'
    params.bwa_idx = 'bwa_index/path'
--------------------------------------------------------------------------------
3. Pre process Reads
   Convert bam to fastq
   Trim the reads (conservatively)
   run fastqc - multiqc on the reads
   
   Parameters:
    params.input_type = fastq | bam
--------------------------------------------------------------------------------
4. Align Reads
   Align using BWA 
   Align using STAR
   Align using hisat2 (TO DO)
   
   Parameters
    aligner (alread init)
*/


/*
 * PARAMETERS
*/


//Step 1
params.outdir = './'
params.version = 'GRCh38'
params.fasta = '/data/bdigby/grch38/reference/GRCh38.fa'
params.gene_annotation = '/data/bdigby/grch38/reference/GRCh38.txt'
params.gencode_gtf = '/data/bdigby/grch38/reference/GRCh38.gtf'
//Step 2
params.star_overhang = '49' 
params.star_index = '/data/bdigby/grch38/index/star_index'
params.bwa_index = ''
params.aligner = 'bwa'
//Step 3
params.inputdir = '/data/bdigby/circTCGA/fastq/'
params.input_type = 'bam'
params.fastq_glob = '*_R{1,2}.fastq.gz'
params.bam_glob = '*.bam'
params.reads = '' // leave empty 
//Step 4
params.adapters = '/data/bdigby/grch38/adapters.fa'
//Step 5


/*
 * Step 1: Download Reference Files
 */

process download_genome {
        
        publishDir "$params.outdir/reference", mode: 'copy'

        output:
        file('*.fa') into fasta_downloaded
        file('*.txt') into gene_annotation_created
        file('*.gtf') into gencode_gtf_downloaded
        
        shell:
        if(params.version == 'GRCh37' && !(params.fasta) && !(params.gencode_gtf) && !(params.gene_annotation)){
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
          
        }else if(params.version == 'development' && !(params.fasta) && !(params.gencode_gtf) && !(params.gene_annotation)){
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
          
        }else if(params.version == 'GRCh38' && !(params.fasta) && !(params.gencode_gtf) && !(params.gene_annotation)){
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
 
if(params.aligner == 'star' && !(params.star_index)){
    process star_index {
    
          publishDir "$params.outdir/index", mode:'copy'
          
          input:
              file(fasta) from ch_fasta
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
          --sjdbOverhang $params.star_overhang \
          --genomeDir star_index/ \
          --genomeFastaFiles $fasta
          """
          }
          ch_star_index = params.star_index ? Channel.value(file(params.star_index)) : star_built
} else if(params.aligner == 'bwa' && !(params.bwa_index)){
    process bwa_index {
    
        publishDir "$params.outdir/index/bwa", mode:'copy'
        
        input:
            file(fasta) from ch_fasta
            
        output:
            file("${fasta}.*") into bwa_built
            
        script:
        """
        bwa index ${fasta}
        """
        }
        ch_bwa_index = params.bwa_index ? Channel.value(file(params.bwa_index)) : bwa_built
 }

// if paths provided to index (skips above step) then must assign to ch
ch_star_index = Channel.value(file(params.star_index))
ch_bwa_index = Channel.value(file(params.bwa_index))

/*
 * Step3: Stage Fastq files
 */
 
 // stage bam files
bam_files = params.inputdir + params.bam_glob

if(params.input_type == 'bam'){
   ch_bam = Channel.fromPath( bam_files )
                   .map{ file -> [file.baseName, file]}
      process bam_to_fq{

          input:
              tuple val(base), file(bam) from ch_bam

          output:
              tuple val(base), file('*.fastq.gz') into fastq_built

          script:
          """
          picard -Xmx8g \
          SamToFastq \
          I=$bam \
          F=${base}_R1.fastq.gz F2=${base}_R2.fastq.gz \
          VALIDATION_STRINGENCY=LENIENT
          """
        }
    }else if(params.input_type == 'fastq'){
          fastq_build = params.inputdir + params.fastq_glob
          Channel.fromFilePairs( fastq_build )
                 .set{ fastq_built }
}
    
ch_reads = params.reads ? Channel.value(file(params.reads)) : fastq_built
 

/*
 * Step 4: Trim
 */
 
 process bbduk {
    
        publishDir "$params.outdir/trimmed_reads", mode:'copy'
    
        input:
            tuple val(base), file(fastq) from ch_reads
            path adapters from params.adapters
        output:
            tuple val(base), file('*.fastq.gz') into trim_reads_built
            
        script:
        """
        bbduk.sh -Xmx4g \
        in1=${fastq[0]} \
        in2=${fastq[1]} \
        out1=${base}_1.fastq.gz \
        out2=${base}_2.fastq.gz \
        ref=$adapters \
        minlen=30 \
        ktrim=r \
        k=12 \
        qtrim=r \
        trimq=20
        """
}


// QC left out for now. collect as you go for QC and run at the end
// MultiQC does not work with python 2.7, will have to think of workaround. 


if(params.aligner == 'bwa'){
    process bwa_align{
    
        publishDir "$params.outdir/bwa_alignment", mode:'copy', overwrite: true
    
          input:
              tuple val(base), file(fastq) from trim_reads_built
              file(index) from ch_bwa_index.collect()
              
          output:
              tuple val(base), file(sam) into circexplorer2_input
              
          script:
          """
          bwa mem -T 19 -t 8 $index ${fastq[0]} ${fastq[1]} > ${base}.sam
          """
          }
} else if(params.aligner == 'star'){
    process star_align{
    
        publishDir "$params.outdir/star_alignment", mode:'copy', overwrite: true
    
        input:
            tuple val(base), file(fastq) from trim_reads_built
            file(gtf) from ch_gencode_gtf
            val(star_idx) from ch_star_index
            
        output:
            tuple val(base), file("${base}.Chimeric.out.junction") into circexplorer2_input
            
        script:
        """
        STAR    \
        --runThreadN 8 \
        --twopassMode Basic \
        --twopass1readsN -1 \
        --genomeLoad NoSharedMemory \
        --genomeDir $star_idx \
        --readFilesIn ${fastq[0]},${fastq[1]} \
        --readFilesCommand zcat \
        --outFileNamePrefix ${base}. \
        --outSJfilterOverhangMin 15 15 15 15 \
        --outFilterMultimapNmax 1 \
        --outFilterMultimapScoreRange 1 \
        --outFilterScoreMin 1 \
        --outFilterMatchNminOverLread 0.33 \
        --outFilterMismatchNmax 10 \
        --outFilterMismatchNoverLmax 0.05 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --alignSJoverhangMin 1 \
        --alignSJDBoverhangMin 1 \
        --alignSoftClipAtReferenceEnds No \
        --chimSegmentMin 10 \
        --chimScoreMin 15 \
        --chimScoreSeparation 10 \
        --chimJunctionOverhangMin 15 \
        --sjdbGTFfile $gtf  \
        --sjdbOverhang $params.star_overhang \
        --sjdbScore 2 \
        --chimOutType Junctions \
        --outSAMtype BAM SortedByCoordinate
        """
        }
}
            
