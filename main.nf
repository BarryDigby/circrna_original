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
    log.info print_purple('------------------------------------------------------------------------')
    log.info "tool_name: A Nextflow based circular RNA analsis pipeline"
    log.info "tool_name integrates several NGS processing tools to identify novel circRNAs from "
    log.info "un-processed RNA sequencing data. To run this pipeline users need to install nextflow"
    log.info "and singularity. "
    print_purple('------------------------------------------------------------------------')
    log.info ''
    log.info print_yellow('Usage: ') +
    
            print_purple('Nextflow run BarryDigby/circRNA --profile singularity, standard <options> \n') +

            print_yellow('    General arguments:             Input and output setting\n') +
            print_cyan('      --inputdir <path>            ') + print_green('Path to input data\n') +
            print_cyan('      --input_type <str>           ') + print_green('Input data type. Supported: \'fastq\', \'bam\'\n') +
            print_cyan('      --fastq_glob <\"str\">         ') + print_green('Glob pattern of fastq files e.g: \'_R{1,2}.fastq.gz\'\n') +
            print_cyan('      --bam_glob <\"str\">           ') + print_green('Glob pattern of bam files expected: \'*.bam\'\n') +
            print_cyan('      --aligner <str>              ') + print_green('Aligner to use for analysis. Supported: \'bwa\', \'star\'\n') +
            '\n' +
            print_yellow('    Input Files:\n') +
            print_cyan('      --fasta <path>               ') + print_green('Path to genome fasta file\n') +
            print_cyan('      --fasta_fai <path>           ') + print_green('Path to genome fasta fai file\n') +
            print_cyan('      --gencode_gtf <path>         ') + print_green('Path to genocde gtf file\n') + 
            print_cyan('      --gene_annotation <path>     ') + print_green('Path to gene annotation file \n') + 
            print_cyan('      --star_index <\"str\">         ') + print_green('Path to STAR index\n') +
            print_cyan('      --bwa_index <\"str/*glob\">    ') + print_green('Path to BWA index\n') +
            print_cyan('      --adapters <path>            ') + print_green('Fasta file containing adapters to trim\n') +
            print_cyan('      --mirna_database <path>      ') + print_green('Fasta file containing mature miRNA sequences\n') +
            '\n' + 
            print_yellow('    Paramters:\n') +
            print_cyan('      --star_overhang <int>        ') + print_green('Parameter for STAR (Read length - 1)\n') +

            log.info ('------------------------------------------------------------------------')
            log.info print_yellow('Contact information: b.digby237@gmail.com') 
            log.info print_yellow('O\'Broin Lab, National University of Ireland Galway')
    log.info ('------------------------------------------------------------------------')
    exit 0
}




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
--------------------------------------------------------------------------------
5. circRNA Discovery
   CIRCexplorer2 
   Need to configure container(s) for more tools. 
   
   Parameters
    even though only 1 for now, 
    code if else options. 
--------------------------------------------------------------------------------
6. Filter circRNA results
   remove low read counts
--------------------------------------------------------------------------------
7. predict circRNA - miRNA interactions
   make sure each circRNA in each sample is 
   done in parallel! might have to stage
   circRNAs in script prior. 
*/


/*
 * PARAMETERS
*/


//Step 1
params.outdir = './'
params.version = ''
params.fasta = ''
params.fasta_fai = ''
params.gene_annotation = ''
params.gencode_gtf = ''
//Step 2
params.star_overhang = '49' 
params.star_index = ''
params.bwa_index = ''
params.aligner = ''
//Step 3
params.inputdir = '/data/bdigby/circTCGA/fastq/'
params.input_type = 'bam'
params.fastq_glob = '*_R{1,2}.fastq.gz'
params.bam_glob = '*.bam'
//Step 4
params.adapters = '/data/bdigby/grch38/adapters.fa'
//Step 5
params.circrna_tool = ''
//Step 6
// samtools_index
//Step 7
params.mirna_database = ''
params.circrna_discovered = '' // leave empty

/*
 * Step 1: Download Reference Files
 * LOGIC:
 * If no reference files provided, download them
 * If provided, assign to channel 
 */

if(!(params.fasta) && !(params.gencode_gtf) && !(params.gene_annotation)){
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
}else{
          ch_fasta = Channel.value(file(params.fasta)) 
          ch_gene_annotation = Channel.value(file(params.gene_annotation))
          ch_gencode_gtf =Channel.value(file(params.gencode_gtf)) 
}

// ternary operators: result = condition ? value_if_true : value_if_false
ch_fasta = params.fasta ? Channel.value(file(params.fasta)) : fasta_downloaded
ch_gene_annotation = params.gene_annotation ? Channel.value(file(params.gene_annotation)) : gene_annotation_created
ch_gencode_gtf = params.gencode_gtf ? Channel.value(file(params.gencode_gtf)) : gencode_gtf_downloaded



/*
 * Step 2: Create Genome Index
 * LOGIC
 * if not made, make them
 * if made, assign to channel
 * need samtools index for generating
 * inputs to miRNA
 */ 
 
if(!(params.fasta_fai)){
  process samtools_index{

        publishDir "$params.outdir/reference", mode:'copy'
        
        input:
            file(fasta) from ch_fasta
            
        output:
            file("${fasta}.fai") into fasta_fai_built
            
        script:
        """
        samtools faidx $fasta
        """
        }
        
        ch_fai = params.fasta_fai ? Channel.value(file(params.fasta_fai)) : fasta_fai_built
} else{

        ch_fai = Channel.value(file(params.fasta_fai))
}
 
 
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

}else if(params.aligner == 'star' && params.star_index){
          
        ch_star_index = params.star_index
          
}else if(params.aligner == 'bwa' && !(params.bwa_index)){
    process bwa_index {
    
        publishDir "$params.outdir/index/bwa", mode:'copy'
        
        input:
            file(fasta) from ch_fasta
            
        output:
            file("${fasta.baseName}.*") into bwa_built
            
        script:
        """
        bwa index ${fasta} -p ${fasta.baseName}
        """
        }
        
        ch_bwa_index = params.bwa_index ? Channel.value(file(params.bwa_index)) : bwa_built
 
 } else if(params.aligner == 'bwa' && params.bwa_index){
 
        ch_bwa_index = Channel.fromPath(params.bwa_index)                   
}



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
        F=${base}_R1.fastq.gz \
        F2=${base}_R2.fastq.gz \
        VALIDATION_STRINGENCY=LENIENT
        """
        }
}else if(params.input_type == 'fastq'){
         fastq_build = params.inputdir + params.fastq_glob
         Channel.fromFilePairs( fastq_build )
                .set{ fastq_built }
}
    
ch_reads = fastq_built
 

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
              file(fasta) from ch_fasta
              
          output:
              tuple val(base), file("${base}.sam") into circexplorer2_input
              
          script:
          """
          bwa mem -T 19 -t 8 ${fasta.baseName} $fastq > ${base}.sam
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
            

/*
 * Step 5 circRNA discovery
*/

if(params.circrna_tool == 'circexplorer2' && params.aligner == 'bwa'){
    process circexplorer2_bwa{
    
        publishDir "$params.outdir/circexplorer2", mode:'copy'
        
        input:
            tuple val(base) file(sam) from circexplorer2_input
            file(fasta) from ch_fasta
            file(gene_annotation) from ch_gene_annotation
            
        output:
            tuple val(base), file("${base}.BWA.circRNA.txt") into circrna_discovered
            
        script:
        """
        CIRCexplorer2 parse -t BWA $sam -b ${base}.BWA.junction.bed

        CIRCexplorer2 annotate -r $gene_annotation -g $fasta -b ${base}.BWA.junction.bed -o ${base}.BWA.circRNA.txt
        """
        }  
}else if(params.circrna_tool == 'circexplorer2' && params.aligner == 'star'){
    process circexplorer2_star{
    
        publishDir "$params.outdir/circexplorer2", mode:'copy'
        
        input:
            tuple val(base), file(chimeric_reads) from circexplorer2_input
            file(fasta) from ch_fasta
            file(gene_annotation) from ch_gene_annotation
            
        output:
            tuple val(base), file("${base}.STAR.circRNA.txt") into circrna_discovered
            
        script:
        """
        CIRCexplorer2 parse -t STAR $chimeric_reads -b ${base}.STAR.junction.bed

        CIRCexplorer2 annotate -r $gene_annotation -g $fasta -b ${base}.STAR.junction.bed -o ${base}.STAR.circRNA.txt
        """
        }
}


/*
 * Step 6
 * Filter low count circRNAs
 * Generate fasta file of sequences
 */
 
 
process get_sequences{
    
        publishDir "$params.outdir/circexplorer2", mode:'copy'
        
        input:
            tuple val(base), file(circrna_discovered) from circrna_discovered
            file(fasta) from ch_fasta
            file(fai) from ch_fai
            
        output:
            tuple val(base), file("${base}.fa") into circrna_fasta
            
        script:
        """
        bash filter_circexplorer2.sh $fasta $circrna_discovered
        """
}


/*
 * Step 7 
 * miRanda analysis
 *
 */
 
process miRanda{
 
        publishDir "$params.outdir/mirna_bindsites", mode:'copy'
        
        input:
            tuple val(base), file(circrna_fasta) from circrna_fasta
            file(miRs) from params.mirna_database
            
        output:
            tuple val(base), file("${base}.mirna.bindsites.txt") into mirna_predictions
            
        script:
        """
        miranda $miRs $circrna_fasta -out ${base}.mirna.bindsites.out -quiet
        
        echo "miRNA Target  Score Energy-Kcal/Mol Query-Aln(start-end) Subject-Al(Start-End) Al-Len Subject-Identity Query-Identity" > ${base}.mirna.bindsites.txt
        
        grep -A 1 "Scores for this hit:" ${base}.mirna.bindsites.out | sort | grep ">" | cut -c 2- >> ${base}.mirna.bindsites.txt
        """
  }
