#!/usr/bin/env nextflow

/*
================================================================================
                                circRNA analysis
================================================================================
Started August 2020
--------------------------------------------------------------------------------
Description:
  (To my knowledge) the first circRNA pipeline to scan RNA-Seq data for circRNAs and
  conduct differential expression analysis + circRNA-miRNA predictions
--------------------------------------------------------------------------------
 @Homepage
 https://github.com/BarryDigby/circRNA
 -------------------------------------------------------------------------------
 @Documentation
 Work in progress
--------------------------------------------------------------------------------
*/

/*
================================================================================
                                  Help Flags
================================================================================
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

            print_yellow('      Mandatory arguments:\n') +
            print_cyan('      --inputdir <path>            ') + print_green('Path to input data\n') +
            print_cyan('      --input_type <str>           ') + print_green('Input data type. Supported: fastq, bam\n') +
            print_cyan('      --fastq_glob <str>           ') + print_green('Glob pattern of fastq files e.g: \'_R{1,2}.fastq.gz\'\n') +
            print_cyan('      --bam_glob <str>             ') + print_green('Glob pattern of bam files. Expected: \'*.bam\'\n') +
            print_cyan('      --tool <str>                 ') + print_green('circRNA tool to use for analysis. \n') +
            print_green('                                   Supported: CIRCexplorer2, CIRIquant, find_circ\n') +
            print_green('                                   UROBORUS, mapsplice, DCC, circRNA_finder\n') +
            print_cyan('      --version <str>              ') + print_green('Genome version. Supported: GRCh37, GRCh38\n') +
            '\n' +
            print_yellow('    Input Files:            if left empty will be generated\n') +
            print_cyan('      --fasta <path>               ') + print_green('Path to genome fasta file\n') +
            print_cyan('      --fasta_fai <path>           ') + print_green('Path to genome fasta fai file\n') +
            print_cyan('      --gencode_gtf <path>         ') + print_green('Path to genocde gtf file\n') + 
            print_cyan('      --gene_annotation <path>     ') + print_green('Path to gene annotation file \n') + 
            print_cyan('      --star_index <str>           ') + print_green('Path to STAR index\n') +
            print_cyan('      --bwa_index <str>            ') + print_green('Path to BWA index\n') +
            print_cyan('      --bowtie_index <str>         ') + print_green('Path to Bowtie index (must include glob for files)\n') +
            print_cyan('      --bowtie2_index <str>        ') + print_green('Path to Bowtie2 index (must include glob for files)\n') +
            print_cyan('      --hisat2_index <str>         ') + print_green('Path to Hisat2 index\n') +
            print_cyan('      --ciriquant_yml <str>        ') + print_green('Path to CIRIquant yml configuration file\n') +
            print_cyan('      --adapters <path>            ') + print_green('Fasta file containing adapters to trim\n') +
            print_cyan('      --mirna_database <path>      ') + print_green('Fasta file containing mature miRNA sequences\n') +


            log.info ('------------------------------------------------------------------------')
            log.info print_yellow('Contact information: b.digby237@gmail.com') 
            log.info print_yellow('O\'Broin Lab, National University of Ireland Galway')
    log.info ('------------------------------------------------------------------------')
    exit 0
}


/*
================================================================================
                                  Paramaters
================================================================================
*/


params.outdir = null
params.fasta = null
params.gencode_gtf = null
params.gene_annotation = null
params.version = null
params.tool = null
params.fasta_fai = null
params.bwa_index = null
params.star_index = null
params.hisat2_index = null
params.bowtie_index = null
params.bowtie2_index = null
params.fasta_chr = null
params.ciriquant_yml = null
params.inputdir = null
params.input_type = null
params.fastq_glob = null
params.bam_glob = null
params.adapters = null
params.phenotype = null

toolList = defineToolList()
tool = params.tool ? params.tool.split(',').collect{it.trim().toLowerCase()} : []
if (!checkParameterList(tool, toolList)) exit 1, 'Unknown tool, see --help for more information'


/*
================================================================================
                          Download Files
================================================================================
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

process download_mirbase{
	errorStrategy 'retry'
   	maxRetries 10
	
	publishDir "$params.outdir/assets", mode:'copy'
	
	output:
		file("hsa_mature.fa") into miranda_miRs
		
	script:
	"""
	wget --no-check-certificate ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz
	gunzip mature.fa.gz
	grep "sapiens" -A1 mature.fa | awk '!/--/' > hsa_mature.fa
	"""
}

// TO DO: add a retry attempt for process below (it sometimes fails to resolve the link)

process download_targetscan{
	errorStrategy 'retry'
   	maxRetries 10	
	
	publishDir "$params.outdir/assets", mode:'copy'

	output:
	file("hsa_miR.txt") into targetscan_miRs
	file("hsa_miR_for_context_scores.txt") into targetscan_miRs_context_scores

	script:
	"""
	wget --no-check-certificate http://www.targetscan.org/vert_72/vert_72_data_download/miR_Family_Info.txt.zip
	jar xvf miR_Family_Info.txt.zip
	grep 9606 miR_Family_Info.txt > hsa_miR_Family_Info.txt
	awk -v OFS="\t" '{print \$1, \$2, \$3}' hsa_miR_Family_Info.txt > hsa_miR.txt
	awk -v OFS="\t" '{print \$1, \$3, \$4, \$5}' hsa_miR.txt > hsa_miR_for_context_scores.txt
	"""
}


/*
================================================================================
                          Create Genome Index
================================================================================
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
            val("$launchDir/index/bwa") into bwa_path
        
        when: !(params.bwa_index) && ('ciriquant' in tool || 'combine' in tool)
        
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
            val("$launchDir/index/hisat2") into hisat2_path
            
        when: !(params.hisat2_index) && ('ciriquant' in tool || 'combine' in tool)
        
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
              
        when: !(params.star_index) && ('circexplorer2' in tool || 'circrna_finder' in tool || 'dcc' in tool || 'combine' in tool)
        
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
            file ("${fasta.baseName}.*") into bowtie_built
            val("$launchDir/index/bowtie") into bowtie_path
            
        when: !(params.bowtie_index) && ('mapsplice' in tool || 'uroborus' in tool || 'combine' in tool)

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
            file ("${fasta.baseName}.*") into bowtie2_built
            
        when: !(params.bowtie2_index) && ('find_circ' in tool || 'uroborus' in tool || 'combine' in tool)

        script:
        """
        bowtie2-build $fasta ${fasta.baseName}
        """
}

ch_bowtie2_index = params.bowtie2_index ? Channel.value(file(params.bowtie2_index)) : bowtie2_built
ch_bowtie2_index.view()


/*
================================================================================
                       Misc. circRNA requirements
================================================================================
*/
 
 
process split_fasta{

        publishDir "$params.outdir/index/chromosomes", mode:'copy'
        
        input:
            file(fasta) from ch_fasta
            
        output:
             file("*.fa") into split_fasta
             val("$launchDir/index/chromosomes") into split_fasta_path
             
        when: ('mapsplice' in tool || 'find_circ' in tool || 'combine' in tool)
        
        shell:
        '''
	## Add catch for test data (uses only 1 chr, no action needed)
	n_chr=$(grep '>' !{fasta} | wc -l)
	if [[ $n_chr -gt 1 ]];
	then
        	awk '/^>/ {F=substr($0, 2, length($0))".fa"; print >F;next;} {print >> F;}' < !{fasta}
		rm !{fasta}
	else
		:
	fi
        '''
}

ch_fasta_chr = params.fasta_chr ? Channel.value(params.fasta_chr) : split_fasta_path
ch_fasta_chr.view()

process ciriquant_yml{
        
        publishDir "$params.outdir", mode:'copy'
      
        input:
            file(fasta) from ch_fasta
            val(gencode_gtf_path) from ch_gencode_gtf
            val(fasta_path) from ch_fasta
            val(bwa_path) from ch_bwa_index
            val(hisat2_path) from ch_hisat2_index

        output:
            file("travis.yml") into yml_built

        when: !(params.ciriquant_yml) && ('ciriquant' in tool || 'combine' in tool)

        script:
        index_prefix = fasta.toString() - ~/.fa/
        """
        export bwa=`whereis bwa | cut -f2 -d':'`
        export hisat2=`whereis hisat2 | cut -f2 -d':'`
        export stringtie=`whereis stringtie | cut -f2 -d':'`
        export samtools=`whereis samtools | cut -f2 -d':' | awk '{print \$1}'`

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
         bwa_index: ${bwa_path}/${index_prefix}\n\
         hisat_index: ${hisat2_path}/${index_prefix}" >> travis.yml
        """
}

ch_ciriquant_yml = params.ciriquant_yml ? Channel.value(file(params.ciriquant_yml)) : yml_built
      
/*
================================================================================
                          Process Input Data
================================================================================
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

// Test dataset does not like being trimmed, omit for testing. 

//process bbduk {
//    
//        publishDir "$params.outdir/trimmed_reads", mode:'copy'
//    
//        input:
//            tuple val(base), file(fastq) from ch_reads
//            path adapters from params.adapters
//            
//        output:
//            tuple val(base), file('*.fastq.gz') into trim_reads_built
//            
//        script:
//        """
//        bbduk.sh -Xmx4g \
//        in1=${fastq[0]} \
//        in2=${fastq[1]} \
//        out1=${base}_1.fastq.gz \
//        out2=${base}_2.fastq.gz \
//        ref=$adapters \
//        minlen=30 \
//        ktrim=r \
//        k=12 \
//        qtrim=r \
//        trimq=20
//        """
//}

(circexplorer2_reads, find_circ_reads, ciriquant_reads, mapsplice_reads, uroborus_reads, circrna_finder_reads, dcc_reads, dcc_reads_mate1, dcc_reads_mate2, hisat2_reads) = ch_reads.into(10)

// fastqc , multiqc omitted until I figure out how to incorporate multiQC separately (multiqc is python3, container is python2). 
// perhaps use multiqc container for this. 

/*
================================================================================
                             circRNA Discovery
================================================================================
*/
 
// CIRCexplorer2

process star_align{

        publishDir "$params.outdir/star_alignment", mode:'copy', overwrite: true
    
        input:
            tuple val(base), file(fastq) from circexplorer2_reads
            file(gtf) from ch_gencode_gtf
            val(star_idx) from ch_star_index
            
        output:
            tuple val(base), file("${base}.Chimeric.out.junction") into circexplorer2_input
         
        when: ('circexplorer2' in tool || 'combine' in tool)
        
        script:
        """
        STAR    \
        --runThreadN 16 \
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
        --sjdbScore 2 \
        --chimOutType Junctions \
        --outSAMtype BAM SortedByCoordinate
        """
}


process circexplorer2_star{
    
        publishDir "$params.outdir/circrna_discovery/circexplorer2", mode:'copy'
        
        input:
            tuple val(base), file(chimeric_reads) from circexplorer2_input
            file(fasta) from ch_fasta
            file(gene_annotation) from ch_gene_annotation
            
        output:
            tuple val(base), file("${base}_circexplorer2.bed") into circexplorer2_results
        
        when: ('circexplorer2' in tool || 'combine' in tool)
        
        script:
        """
        CIRCexplorer2 parse -t STAR $chimeric_reads -b ${base}.STAR.junction.bed
        CIRCexplorer2 annotate -r $gene_annotation -g $fasta -b ${base}.STAR.junction.bed -o ${base}.txt
        
        awk '{if(\$13 > 1) print \$0}' ${base}.txt | awk -v OFS="\t" '{print \$1,\$2,\$3,\$6,\$13}' > ${base}_circexplorer2.bed
        """
}


// find_circ

process find_anchors{

        input:
            tuple val(base), file(fastq) from find_circ_reads
            file(fasta) from ch_fasta
            file(bowtie2_index) from ch_bowtie2_index.collect()
        
        output:
            tuple val(base), file("${base}_anchors.qfa.gz") into ch_anchors
            
        when: ('find_circ' in tool || 'combine' in tool)
        
        script:
        """
        bowtie2 -p 16 --very-sensitive --mm -D 20 --score-min=C,-15,0 \
        -x ${fasta.baseName} -q -1 ${fastq[0]} -2 ${fastq[1]} \
        | samtools view -hbuS - | samtools sort --threads 16 -m 2G - > ${base}.bam

        samtools view -hf 4 ${base}.bam | samtools view -Sb - > ${base}_unmapped.bam

        unmapped2anchors.py ${base}_unmapped.bam | gzip > ${base}_anchors.qfa.gz
        """
}


process find_circ{

        publishDir "$params.outdir/circrna_discovery/find_circ", mode:'copy'
        
        input:
            tuple val(base), file(anchors) from ch_anchors
            file(bowtie2_index) from ch_bowtie2_index.collect()
            file(fasta) from ch_fasta
            val(fasta_chr_path) from ch_fasta_chr
        
        output:
            tuple val(base), file("${base}_find_circ.bed") into find_circ_results
         
        when: ('find_circ' in tool || 'combine' in tool)
        
        script:
        """
        bowtie2 -p 16 --reorder --mm -D 20 --score-min=C,-15,0 -q -x ${fasta.baseName} \
        -U $anchors | python /opt/conda/envs/circrna/bin/find_circ.py -G $fasta_chr_path -p ${base} -s ${base}.sites.log > ${base}.sites.bed 2> ${base}.sites.reads

        echo "# chrom:start:end:name:n_reads:strand:n_uniq:best_qual_A:best_qual_B:spliced_at_begin:spliced_at_end:tissues:tiss_counts:edits:anchor_overlap:breakpoints" > tmp.txt

        cat tmp.txt | tr ':' '\t' > ${base}.bed

        grep circ ${base}.sites.bed | grep -v chrM | python /opt/conda/envs/circrna/bin/sum.py -2,3 | python /opt/conda/envs/circrna/bin/scorethresh.py -16 1 | python /opt/conda/envs/circrna/bin/scorethresh.py -15 2 | python /opt/conda/envs/circrna/bin/scorethresh.py -14 2 | python /opt/conda/envs/circrna/bin/scorethresh.py 7 2 | python /opt/conda/envs/circrna/bin/scorethresh.py 8,9 35 | python /opt/conda/envs/circrna/bin/scorethresh.py -17 100000 >> ${base}.txt
        
	tail -n +2 ${base}.txt | awk -v OFS="\t" '{print \$1,\$2,\$3,\$6,\$5}' > ${base}_find_circ.bed
	"""
}


// circRNA_finder

process circrna_finder_star{

        input:
            tuple val(base), file(fastq) from circrna_finder_reads
            val(star_index) from ch_star_index
            
        output:
            tuple val(base), file("${base}") into circrna_finder_star
            
        when: ('circrna_finder' in tool || 'combine' in tool)
        
        script:
        """      
        STAR \
        --genomeDir $star_index \
        --readFilesIn ${fastq[0]} ${fastq[1]} \
        --readFilesCommand zcat \
        --runThreadN 16 \
        --chimSegmentMin 20 \
        --chimScoreMin 1 \
        --chimOutType Junctions SeparateSAMold \
        --alignIntronMax 100000 \
        --outFilterMismatchNmax 4 \
        --alignTranscriptsPerReadNmax 100000 \
        --outFilterMultimapNmax 2 \
        --outFileNamePrefix ${base}/${base}.
        """
}


process circrna_finder{

        publishDir "$params.outdir/circrna_discovery/circrna_finder", mode:'copy'
        
        input:
            tuple val(base), file(star_dir) from circrna_finder_star
         
        output:
            tuple val(base), file("${base}_circrna_finder.bed") into circrna_finder_results
            
        when: ('circrna_finder' in tool || 'combine' in tool)
        
        script:
        """
        postProcessStarAlignment.pl --starDir ${star_dir}/ --outDir ./

	tail -n +2 ${base}.filteredJunctions.bed | awk '{if(\$5 > 1) print \$0}' | awk  -v OFS="\t" -F"\t" '{print \$1,\$2,\$3,\$6,\$5}' > ${base}_circrna_finder.bed
        """
}

// DCC

process dcc_pair{

        input:
            tuple val(base), file(fastq) from dcc_reads
            val(star_index) from ch_star_index

        output:
            tuple val(base), file("samples") into dcc_samples
            
        when: ('dcc' in tool || 'combine' in tool)
        
        script:
        """
        STAR \
        --runThreadN 16 \
        --genomeDir $star_index \
        --outSAMtype BAM SortedByCoordinate \
        --readFilesIn ${fastq[0]} ${fastq[1]} \
        --readFilesCommand zcat \
        --outFileNamePrefix samples/${base}. \
        --outReadsUnmapped Fastx \
        --outSJfilterOverhangMin 15 15 15 15 \
        --alignSJoverhangMin 15 \
        --alignSJDBoverhangMin 15 \
        --outFilterMultimapNmax 20 \
        --outFilterScoreMin 1 \
        --outFilterMatchNmin 1 \
        --outFilterMismatchNmax 2 \
        --chimSegmentMin 15 \
        --chimScoreMin 15 \
        --chimScoreSeparation 10 \
        --chimJunctionOverhangMin 15 
        """
}

process dcc_1{

        input:
            tuple val(base), file(fastq) from dcc_reads_mate1
            val(star_index) from ch_star_index

        output:
            tuple val(base), file("mate1") into dcc_mate1
            
       when: ('dcc' in tool || 'combine' in tool)
        
        script:
        """
        STAR \
        --runThreadN 16 \
        --genomeDir $star_index \
        --outSAMtype None \
        --readFilesIn ${fastq[0]} \
        --readFilesCommand zcat \
        --outFileNamePrefix mate1/${base}. \
        --outReadsUnmapped Fastx \
        --outSJfilterOverhangMin 15 15 15 15 \
        --alignSJoverhangMin 15 \
        --alignSJDBoverhangMin 15 \
        --seedSearchStartLmax 30 \
        --outFilterMultimapNmax 20 \
        --outFilterScoreMin 1 \
        --outFilterMatchNmin 1 \
        --outFilterMismatchNmax 2 \
        --chimSegmentMin 15 \
        --chimScoreMin 15 \
        --chimScoreSeparation 10 \
        --chimJunctionOverhangMin 15
        """
}

process dcc_2{

        input:
            tuple val(base), file(fastq) from dcc_reads_mate2
            val(star_index) from ch_star_index

        output:
            tuple val(base), file("mate2") into dcc_mate2
            
        when: ('dcc' in tool || 'combine' in tool)
	
        script:
        """
        STAR \
        --runThreadN 16 \
        --genomeDir $star_index \
        --outSAMtype None \
        --readFilesIn ${fastq[1]} \
        --readFilesCommand zcat \
        --outFileNamePrefix mate2/${base}. \
        --outReadsUnmapped Fastx \
        --outSJfilterOverhangMin 15 15 15 15 \
        --alignSJoverhangMin 15 \
        --alignSJDBoverhangMin 15 \
        --seedSearchStartLmax 30 \
        --outFilterMultimapNmax 20 \
        --outFilterScoreMin 1 \
        --outFilterMatchNmin 1 \
        --outFilterMismatchNmax 2 \
        --chimSegmentMin 15 \
        --chimScoreMin 15 \
        --chimScoreSeparation 10 \
        --chimJunctionOverhangMin 15
        """
}

// collect runs according to val(base) in tuple
ch_dcc_dirs = dcc_samples.join(dcc_mate1).join(dcc_mate2)

process dcc{

        publishDir "$params.outdir/circrna_discovery/dcc", mode:'copy'

        input:
            tuple val(base), file(samples), file(mate1), file(mate2) from ch_dcc_dirs
            file(gtf) from ch_gencode_gtf
            file(fasta) from ch_fasta

        output:
            tuple val(base), file("${base}_dcc.bed") into dcc_results

       when: ('dcc' in tool || 'combine' in tool)
        
        script:
        COJ="Chimeric.out.junction"
        """
        sed -i 's/^chr//g' $gtf

        printf "samples/${base}.${COJ}" > samplesheet
        printf "mate1/${base}.${COJ}" > mate1file
        printf "mate2/${base}.${COJ}" > mate2file

        DCC @samplesheet -mt1 @mate1file -mt2 @mate2file -D -an $gtf -Pi -ss -F -M -Nr 1 1 -fg -A $fasta -N -T 8
        
        awk '{print \$6}' CircCoordinates >> strand
        paste CircRNACount strand | tail -n +2 | awk -v OFS="\t" '{print \$1,\$2,\$3,\$5,\$4}' >> ${base}_dcc.txt
	bash filter_DCC.sh ${base}_dcc.txt
        """
}

// CIRIquant


process ciriquant{

        publishDir "$params.outdir/circrna_discovery/ciriquant", mode:'copy'
        
        input:
            tuple val(base), file(fastq) from ciriquant_reads
            file(ciriquant_yml) from ch_ciriquant_yml

        output:
            tuple val(base), file("${base}_ciriquant.bed") into ciriquant_results
            
        when: ('ciriquant' in tool || 'combine' in tool)
        
        script:
        """
        CIRIquant -t 16 \
        -1 ${fastq[0]} \
        -2 ${fastq[1]} \
        --config $ciriquant_yml \
        --no-gene \
        -o ${base} \
        -p ${base}
        
        mv ${base}/${base}.gtf ${base}_ciriquant.gtf
	
	bash filter_CIRIquant.sh ${base}_ciriquant.gtf
        """
}
      
      
// mapsplice

process mapsplice_align{

        publishDir "$params.outdir/mapsplice", mode:'copy'
        
        input:
            tuple val(base), file(fastq) from mapsplice_reads
            val(mapsplice_ref) from ch_fasta_chr
            file(bowtie_index) from ch_bowtie_index.collect()
            file(gtf) from ch_gencode_gtf

        output:
            tuple val(base), file("${base}/fusions_raw.txt") into mapsplice_fusion

        when: ('mapsplice' in tool || 'combine' in tool)

        script:
        prefix = gtf.toString() - ~/.gtf/
        """
        gzip -d --force ${fastq[0]}
        gzip -d --force ${fastq[1]}
        
        mapsplice.py \
        -c $mapsplice_ref \
        -x $prefix \
        -1 ${base}_r1.fastq \
        -2 ${base}_r2.fastq \
        -p 8 \
        --bam \
        --seglen 25 \
        --min-map-len 40 \
        --fusion-non-canonical \
        --min-fusion-distance 200 \
        --gene-gtf $gtf \
        -o $base
        """
}


process mapsplice_parse{

        publishDir "$params.outdir/circrna_discovery/mapsplice", mode:'copy'
        
        input:
            tuple val(base), file(raw_fusion) from mapsplice_fusion
            file(fasta) from ch_fasta
            file(gene_annotation) from ch_gene_annotation
            
        output:
            tuple val(base), file("${base}_mapsplice.bed") into mapsplice_results
        
        when: ('mapsplice' in tool || 'combine' in tool)
        
        script:
        """
        CIRCexplorer2 parse -t MapSplice $raw_fusion -b ${base}.mapsplice.junction.bed

        CIRCexplorer2 annotate -r $gene_annotation -g $fasta -b ${base}.mapsplice.junction.bed -o ${base}.txt
	
	awk '{if(\$13 > 1) print \$0}' ${base}.txt | awk -v OFS="\t" '{print \$1,\$2,\$3,\$6,\$13}' > ${base}_mapsplice.bed
        """
}


// UROBORUS

process tophat_align{

        input:
            tuple val(base), file(fastq) from uroborus_reads
            file(bowtie2_index) from ch_bowtie2_index.collect()
            file(fasta) from ch_fasta
            
        output:
            tuple val(base), file("unmapped.bam") into tophat_unmapped_bam
            tuple val(base), file("accepted_hits.bam") into tophat_accepted_hits
        
        when: 'uroborus' in tool
        
        script:
        """
        tophat -p 8 -o ${base} ${fasta.baseName} ${fastq[0]} ${fastq[1]}
        mv ${base}/unmapped.bam ./ 
        mv ${base}/accepted_hits.bam ./
        """
}
            
            
process uroborus{

        publishDir "$params.outdir/circrna_discovery/uroborus", mode:'copy'
        
        input:
            tuple val(base), file(unmapped_bam) from tophat_unmapped_bam
            tuple val(base), file(accepted_hits) from tophat_accepted_hits
            file(bowtie_index) from ch_bowtie_index.collect()
            file(gtf) from ch_gencode_gtf
            val(uroborus_ref) from ch_fasta_chr
            file(fasta) from ch_fasta
            
        output:
            file("${base}.txt") into uroborus_results
            
        when: 'uroborus' in tool
        
        script:
        """
        samtools view $unmapped_bam > unmapped.sam
        
        perl /opt/conda/envs/circrna/bin/UROBORUS.pl \
        -index ${fasta.baseName} \
        -gtf $gtf \
        -fasta $uroborus_ref \
        unmapped.sam $accepted_hits &> uroborus_logs.txt
        
        mv circRNA_list.txt ${base}.txt
        """
}



/*
================================================================================
                                  RNA-Seq
================================================================================
*/

hisat2_files = params.hisat2_index + "/*"
ch_hisat2_index_files = Channel.fromPath( hisat2_files )

process Hisat2_align{

        input:
                tuple val(base), file(fastq) from hisat2_reads
                file(hisat2_index) from ch_hisat2_index_files.collect()
                file(fasta) from ch_fasta

        output:
                tuple val(base), file("${base}.bam") into hisat2_bam

        script:
        """
        hisat2 -p 16 --dta -q -x ${fasta.baseName} -1 ${fastq[0]} -2 ${fastq[1]} -t | samtools view -bS - | samtools sort --threads 16 -m 2G - > ${base}.bam
        """
}


process StringTie{

        publishDir "$params.outdir/rna-seq", mode:'copy'

        input:
                tuple val(base), file(bam) from hisat2_bam
                file(gtf) from ch_gencode_gtf
        output:
                file("${base}") into stringtie_dir

        script:
        """
        mkdir ${base}/
        stringtie $bam -e -G $gtf -C ${base}/${base}_cov.gtf -p 16 -o ${base}/${base}.gtf -A ${base}/${base}_genes.list
        """
}
 
/*
================================================================================
                         circRNA Differential Expression
================================================================================
*/


// need to consolidate circRNA results into single channel.
// apply filtering steps before this step! 
// simple awk|grep one liners will do after tool generates results. 

// CONSOLIDATION

if('combine' in tool){

	combined_tool = ciriquant_results.join(circexplorer2_results).join(dcc_results).join(circrna_finder_results).join(find_circ_results).join(mapsplice_results)

        process consolidate_algorithms{
                        echo true
                        publishDir "$params.outdir/circrna_discovery/matrix", mode:'copy'

                        input:
                                tuple val(base), file(ciriquant), file(circexplorer2), file(dcc), file(circrna_finder), file(find_circ), file(mapsplice) from combined_tool

                        output:
                                file("${base}.bed") into sample_counts

                        script:
                        """
			## make tool output csv file
                        files=\$(ls *.bed)

                        for i in \$files; do
                                printf "\$i\n" >> samples.csv
                        done
			
			## Add catch for empty file in tool output
			bash "$projectDir"/bin/check_empty.sh
			
			## Bring forward circRNAs called by at least 2 tools
                        Rscript "$projectDir"/bin/consolidate_algorithms.R samples.csv
			
                        mv combined_counts.bed ${base}.bed
                        """
                        }
			
	process get_counts_combined{
			publishDir "$params.outdir/circrna_discovery/matrix", mode:'copy'
			
			input:
				file(bed) from sample_counts.collect()
				
			output:
				file("circRNA_matrix.txt") into circRNA_counts
				
			script:
			"""
			python "$projectDir"/bin/circRNA_counts_matrix.py > circRNA_matrix.txt
			"""
			}

} else{

        single_tool = ciriquant_results.mix(circexplorer2_results, dcc_results, circrna_finder_results, find_circ_results)

        process get_counts_single{

                        echo true
                        publishDir "$params.outdir/circrna_discovery/matrix", mode:'copy'


                        input:
                                file(bed) from single_tool.collect()
				val(tool) from params.tool
                        output:
                                file("circRNA_matrix.txt") into circRNA_counts

                        script:
                        """
                        for b in *.bed; do 
				foo=\${b%".bed"}; 
				bar=\${foo%"_${tool}"}; 
				mv \$b \${bar}.bed
			done
			
			python "$projectDir"/bin/circRNA_counts_matrix.py > circRNA_matrix.txt
                        """
                        }
}

ch_phenotype = file(params.phenotype)

process diff_exp{

	publishDir "$params.outdir/Differential_Expression", mode:'copy'
	
	input:
		file(gtf_dir) from stringtie_dir.collect()
		file(circ_matrix) from circRNA_counts
		file(phenotype) from ch_phenotype
		
	output:
		file("RNA-Seq") into rnaseq_dir
		file("circRNA") into circrna_dir
		
	script:
	"""
	for i in \$(ls -d */); do sample=\${i%"/"}; file=\${sample}.gtf; touch samples.txt; printf "\$sample\t\${i}\${file}\n" >> samples.txt; done
	
	prepDE.py -i samples.txt
	
	Rscript "$projectDir"/bin/DEA.R gene_count_matrix.csv $phenotype $circ_matrix
	"""
}

(circrna_dir_mature_seq, circrna_dir_parent_gene, circrna_dir_report) = circrna_dir.into(3)
(rnaseq_dir_parent_gene, rnaseq_dir_report) = rnaseq_dir.into(2)

/*
================================================================================
                         circRNA - miRNA prediction
================================================================================
*/

//Create filtered gtf in its own channel so sym link used,
// saves space in work dir. 


process remove_unwanted_biotypes{

	input:
		file(gtf) from ch_gencode_gtf
	
	output:
		file("filt.gtf") into ch_gtf_filtered
		
	script:
	"""
	cp "$projectDir"/bin/unwanted_biotypes.txt ./
	
	grep -vf unwanted_biotypes.txt $gtf > filt.gtf
	"""
}


process get_mature_seq{

	publishDir "$params.outdir", mode:'copy', pattern: 'bed12/*.bed'
	
	input:
		file(fasta) from ch_fasta
		file(fai) from ch_fai
		file(gtf) from ch_gtf_filtered
		file(circRNA) from circrna_dir_mature_seq
		
	output:
		file("miranda/*.fa") into miranda_sequences
		file("targetscan/*.txt") into targetscan_sequences
		file("bed12/*.bed") into bed_files
		
	script:
	up_reg = "${circRNA}/*up_regulated_differential_expression.txt"
	down_reg = "${circRNA}/*down_regulated_differential_expression.txt"
	"""
	# Extract circRNA ID's from DESeq2 DECs. 
	awk '{print \$1}' $up_reg | tail -n +2 > up_reg_circ.txt
	awk '{print \$1}' $down_reg | tail -n +2 > down_reg_circ.txt
	
	# Split ID to BED file
	bash "$projectDir"/bin/ID_to_BED.sh up_reg_circ.txt
	bash "$projectDir"/bin/ID_to_BED.sh down_reg_circ.txt

	# Consolidate DEC BED files
	cat *.bed > de_circ.bed
	
	# Create BED12 files
	bash "$projectDir"/bin/get_mature_seq.sh 
	
	# Create miRanda inputs
	bedtools getfasta -fi $fasta -bed de_circ_exon_annotated.bed -s -split -name > de_circ_sequences.fa_tmp
	grep -A 1 '>' de_circ_sequences.fa_tmp | cut -d: -f1,2,3 > de_circ_sequences.fa && rm de_circ_sequences.fa_tmp
	mkdir -p miranda
	awk -F '>' '/^>/ {F=sprintf("miranda/%s.fa",\$2); print > F;next;} {print >> F;}' < de_circ_sequences.fa
	
	# Create TargetScan inputs
	bedtools getfasta -fi $fasta -bed de_circ_exon_annotated.bed -s -split -tab | sed 's/(/:/g' | sed 's/)//g' > de_circ_seq_tab.txt_tmp
	awk -v OFS="\t" '{print \$1, 9606, \$2}' de_circ_seq_tab.txt_tmp > de_circ_seq_tab.txt && rm de_circ_seq_tab.txt_tmp
	mkdir -p targetscan
	while IFS='' read -r line; do name=\$(echo \$line | awk '{print \$1}'); echo \$line | sed 's/ /\t/g' >> targetscan/\${name}.txt; done < de_circ_seq_tab.txt
	"""
}

process miRanda{

	publishDir "$params.outdir/miRanda", mode:'copy'
	
	input:
		file(mirbase) from miranda_miRs
		file(miranda) from miranda_sequences.flatten()
	
	output:
		file("*.miRanda.txt") into miranda_out
		file("*.mature_len.txt") into mature_len
	script:
	prefix = miranda.toString() - ~/.fa/
	"""
	grep -v '>' $miranda | wc -c > ${prefix}.mature_len.txt
	miranda $mirbase $miranda -strict -out ${prefix}.bindsites.out -quiet
        echo "miRNA Target Score Energy_KcalMol Query_Start Query_End Subject_Start Subject_End Aln_len Subject_Identity Query_Identity" | tr ' ' '\t' > ${prefix}.miRanda.txt
        grep -A 1 "Scores for this hit:" ${prefix}.bindsites.out | sort | grep ">" | cut -c 2- | tr ' ' '\t' >> ${prefix}.miRanda.txt
	"""
}

process targetscan{

	publishDir "$params.outdir/TargetScan", mode:'copy'
	
	input:
		file(miR) from targetscan_miRs
		file(circ) from targetscan_sequences.flatten()
		
	output:
		file("*.targetscan.txt") into targetscan_out 
		
	script:
	prefix = circ.toString() - ~/.txt/
	"""
	targetscan_70.pl $miR $circ ${prefix}.targetscan.txt
	"""
}


process get_parent_gene{

	input:
		file(gtf) from ch_gtf_filtered
		file(circRNA) from circrna_dir_parent_gene
		
	output:
		file("parent_genes/*.txt") into parent_genes
		
	script:
	up_reg = "${circRNA}/*up_regulated_differential_expression.txt"
	down_reg = "${circRNA}/*down_regulated_differential_expression.txt"
	"""
	# Extract circRNA ID's from DESeq2 DECs. 
	awk '{print \$1}' $up_reg | tail -n +2 > up_reg_circ.txt
	awk '{print \$1}' $down_reg | tail -n +2 > down_reg_circ.txt
	
	# Split ID to BED file
	bash "$projectDir"/bin/ID_to_BED.sh up_reg_circ.txt
	bash "$projectDir"/bin/ID_to_BED.sh down_reg_circ.txt

	# Consolidate DEC BED files
	cat *.bed > de_circ.bed
	
	bash "$projectDir"/bin/get_parent_genes.sh
	"""
}


/*
================================================================================
                         circRNA Plots
================================================================================
*/


// Create tuples, merge channels by simpleName for report. 
ch_mature_len = mature_len.map{ file -> [file.simpleName, file]}
ch_parent_genes_tmp = parent_genes.flatten()
ch_parent_genes = ch_parent_genes_tmp.map{ file -> [file.simpleName, file]}
ch_targetscan = targetscan_out.map{ file -> [file.simpleName, file]}
ch_miranda = miranda_out.map{ file -> [file.simpleName, file]}
ch_bed_tmp = bed_files.flatten()
ch_bed = ch_bed_tmp.map{ file -> [file.simpleName, file]}


ch_report = ch_targetscan.join(ch_miranda).join(ch_bed).join(ch_parent_genes).join(ch_mature_len)

// must combine folders here or else process uses once then exits. 
ch_DESeq2_dirs = circrna_dir_report.combine(rnaseq_dir_report)

process make_circRNA_plots{
	publishDir "$params.outdir/circRNA_Report", mode:'copy'
	
	input:
		file(phenotype) from ch_phenotype
		tuple val(base), file(targetscan), file(miranda), file(bed), file(parent_gene), file(mature_len), file(circRNA), file(rnaseq) from ch_report.combine(ch_DESeq2_dirs)

	output: 
		file("chr*") into circRNA_plots
		
	script:
	up_reg = "${circRNA}/*up_regulated_differential_expression.txt"
	down_reg = "${circRNA}/*down_regulated_differential_expression.txt"
	circ_counts = "${circRNA}/DESeq2_normalized_counts.txt"
	gene_counts = "${rnaseq}/DESeq2_normalized_counts.txt"
	"""
	# create file for circos plot
	bash "$projectDir"/bin/prep_circos.sh $bed
	
	# merge upreg, downreg info 
	cat $up_reg $down_reg > de_circ.txt
	
	# remove 6mers from TargetScan 
	grep -v "6mer" $targetscan > targetscan_filt.txt
	
	# Make plots and generate circRNA info
	Rscript "$projectDir"/bin/circ_report.R de_circ.txt $circ_counts $gene_counts $parent_gene $bed $miranda targetscan_filt.txt $mature_len $phenotype circlize_exons.txt 
	"""
}

// collect all from previous process
//master_ch = circRNA_plots.collect()
(test, test1) = circRNA_plots.into(2)
test.view()
// delete text files in process script, left with only dirs. 


process master_report{
	publishDir "$params.outdir/circRNA_Report", mode:'copy'
	
	input:
		file(reports) from test1.collect()
		
	output:
		file("DE_circRNA_Report.txt") into master_report
		
	script:
	"""
	## extract reports
	for dir in '*/'; do cp \$dir/*_Report.txt .; done
	
	# remove header, add manually
	cat *.txt > merged.txt
	grep -v "Log2FC" merged.txt > no_headers.txt
	echo "circRNA_ID Parent_Gene Mature_Length Strand Log2FC pvalue Adjusted_pvalue" | tr ' ' '\t' > headers.txt
	cat headers.txt no_headers.txt > merged_reports.txt
	
	Rscript "$projectDir"/bin/annotate_report.R
	"""
}

// Check parameter existence
def checkParameterExistence(it, list) {
    if (!list.contains(it)) {
        log.warn "Unknown parameter: ${it}"
        return false
    }
    return true
}

// Compare each parameter with a list of parameters
def checkParameterList(list, realList) {
    return list.every{ checkParameterExistence(it, realList) }
}


// Define list of available tools
def defineToolList() {
    return [
        'ciriquant',
        'circexplorer2',
        'find_circ',
        'circrna_finder',
        'dcc',
        'mapsplice',
        'uroborus',
	'combine'
        ]
}
