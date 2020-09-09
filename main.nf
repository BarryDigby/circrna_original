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

            print_yellow('    Mandatory arguments:\n') +
            print_cyan('      --inputdir <path>            ') + print_green('Path to input data\n') +
            print_cyan('      --input_type <str>           ') + print_green('Input data type. Supported: \'fastq\', \'bam\'\n') +
            print_cyan('      --fastq_glob <str>         ') + print_green('Glob pattern of fastq files e.g: \'_R{1,2}.fastq.gz\'\n') +
            print_cyan('      --bam_glob <str>           ') + print_green('Glob pattern of bam files expected: \'*.bam\'\n') +
            print_cyan('      --tool <str>              ') + print_green('circRNA tool to use for analysis. Supported: \'CIRCexplorer2\', \'CIRIquant\', \'find_circ\', \'UROBORUS\', \'mapsplice\'\n') +
            '\n' +
            print_yellow('    Input Files:            if left empty will be generated\n') +
            print_cyan('      --fasta <path>               ') + print_green('Path to genome fasta file\n') +
            print_cyan('      --fasta_fai <path>           ') + print_green('Path to genome fasta fai file\n') +
            print_cyan('      --gencode_gtf <path>         ') + print_green('Path to genocde gtf file\n') + 
            print_cyan('      --gene_annotation <path>     ') + print_green('Path to gene annotation file \n') + 
            print_cyan('      --star_index <str>         ') + print_green('Path to STAR index\n') +
            print_cyan('      --bwa_index <str>    ') + print_green('Path to BWA index\n') +
            print_cyan('      --bowtie_index <str>    ') + print_green('Path to Bowtie index (must include glob for files)\n') +
            print_cyan('      --bowtie2_index <str>    ') + print_green('Path to Bowtie2 index (must include glob for files)\n') +
            print_cyan('      --hisat2_index <str>    ') + print_green('Path to Hisat2 index\n') +
            print_cyan('      --ciriquant_yml <str>    ') + print_green('Path to CIRIquant yml configuration file\n') +
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


params.outdir ='.'
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
params.fasta_chr = ''
params.ciriquant_yml = ''
params.inputdir = '/data/bdigby/circTCGA/fastq/'
params.input_type = 'fastq'
params.fastq_glob = '*_R{1,2}.fastq.gz'
params.bam_glob = '*.bam'
params.adapters = '/data/bdigby/grch38/adapters.fa'

toolList = defineToolList()
tool = params.tool ? params.tool.split(',').collect{it.trim().toLowerCase()} : []
if (!checkParameterList(tool, toolList)) exit 1, 'Unknown tool, see --help for more information'


/*
================================================================================
                          Download Reference Files
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


ch_fasta.view()
ch_gene_annotation.view()
ch_gencode_gtf.view()

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
            val("$launchDir/index/hisat2") into hisat2_path
            
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
            file ("${fasta.baseName}.*") into bowtie_built
            val("$launchDir/index/bowtie") into bowtie_path
            
        when: !(params.bowtie_index) && ('mapsplice' in tool || 'uroborus' in tool)

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
            
        when: !(params.bowtie2_index) && ('find_circ' in tool || 'uroborus' in tool)

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
             
        when: ('mapsplice' in tool || 'find_circ' in tool)
        
        shell:
        '''
        awk '/^>/ {F=substr($0, 2, length($0))".fa"; print >F;next;} {print >> F;}' < !{fasta}
        rm !{fasta}
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

        when: !(params.ciriquant_yml) && 'ciriquant' in tool

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

(circexplorer2_reads, find_circ_reads, ciriquant_reads, mapsplice_reads, uroborus_reads) = trim_reads_built.into(5)

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
         
        when: 'circexplorer2' in tool
        
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
            tuple val(base), file("${base}.txt") into circexplorer2_results
        
        when: 'circexplorer2' in tool
        
        script:
        """
        CIRCexplorer2 parse -t STAR $chimeric_reads -b ${base}.STAR.junction.bed
        CIRCexplorer2 annotate -r $gene_annotation -g $fasta -b ${base}.STAR.junction.bed -o ${base}.txt
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
            
        when: 'find_circ' in tool
        
        script:
        """
        bowtie2 -p 8 --very-sensitive --mm -D 20 --score-min=C,-15,0 \
        -x ${fasta.baseName} -q -1 ${fastq[0]} -2 ${fastq[1]} \
        | samtools view -hbuS - | samtools sort --threads 8 -m 2G - > ${base}.bam

        samtools view -hf 4 ${base}.bam | samtools view -Sb - > ${base}_unmapped.bam

        python ${baseDir}/bin/unmapped2anchors.py ${base}_unmapped.bam | gzip > ${base}_anchors.qfa.gz
        """
}


process find_circ{

        publishDir "$params.outdir/circrna_dicovery/find_circ", mode:'copy'
        
        input:
            tuple val(base), file(anchors) from ch_anchors
            file(bowtie2_index) from ch_bowtie2_index.collect()
            file(fasta) from ch_fasta
            val(fasta_chr_path) from ch_fasta_chr
        
        output:
            tuple val(base), file("${base}.bed") into find_circ_results
         
        when: 'find_circ' in tool
        
        script:
        """
        bowtie2 -p 8 --reorder --mm -D 20 --score-min=C,-15,0 -q -x ${fasta.baseName} \
        -U $anchors | python ${baseDir}/bin/find_circ.py -G $fasta_chr_path -p ${base} -s ${base}.sites.log > ${base}.sites.bed 2> ${base}.sites.reads

        echo "#chrom:start:end:name:n_reads:strand:n_uniq:best_qual_A:best_qual_B:spliced_at_begin:spliced_at_end:tissues:tiss_counts:edits:anchor_overlap:breakpoints" > tmp.txt

        cat tmp.txt | tr ':' '\t' > ${base}.bed

        grep circ ${base}.sites.bed | grep -v chrM | python ${baseDir}/bin/sum.py -2,3 | python ${baseDir}/bin/scorethresh.py -16 1 | python ${baseDir}/bin/scorethresh.py -15 2 | python ${baseDir}/bin/scorethresh.py -14 2 | python ${baseDir}/bin/scorethresh.py 7 2 | python ${baseDir}/bin/scorethresh.py 8,9 35 | python ${baseDir}/bin/scorethresh.py -17 100000 >> ${base}.bed
        """
}


// CIRIquant


process ciriquant{

        publishDir "$params.outdir/circrna_discovery/ciriquant", mode:'copy'
        
        input:
            tuple val(base), file(fastq) from ciriquant_reads
            file(ciriquant_yml) from ch_ciriquant_yml

        output:
            tuple val(base), file("${base}.gtf") into ciriquant_results
            
        when: 'ciriquant' in tool
        
        script:
        """
        CIRIquant -t 8 \
        -1 ${fastq[0]} \
        -2 ${fastq[1]} \
        --config $ciriquant_yml \
        --no-gene \
        -o ${base} \
        -p ${base}
        
        mv ${base}/${base}.gtf ./
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

        when: 'mapsplice' in tool

        script:
        prefix = gtf.toString() - ~/.gtf/
        """
        gzip -d --force ${fastq[0]}
        gzip -d --force ${fastq[1]}
        
        mapsplice.py \
        -c $mapsplice_ref \
        -x $prefix \
        -1 ${base}_1.fastq \
        -2 ${base}_2.fastq \
        -p 8 \
        --bam \
        --seglen 20 \
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
            tuple val(base), file("${base}.txt") into mapsplice_results
        
        when: 'mapsplice' in tool
        
        script:
        """
        CIRCexplorer2 parse -t MapSplice $raw_fusion -b ${base}.mapsplice.junction.bed

        CIRCexplorer2 annotate -r $gene_annotation -g $fasta -b ${base}.mapsplice.junction.bed -o ${base}.txt
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
        mv ${base}/unmapped_bam ./ 
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
            file(uroborus_ref) from ch_fasta_chr
            file(fasta) from ch_fasta
            
        output:
            file("${base}.txt") into uroborus_results
            
        when: 'uroborus' in tool
        
        script:
        """
        samtools view $unmapped_bam > unmapped.sam
        
        perl ${baseDir}/bin/UROBORUS.pl \
        -index ${fasta.baseName} \
        -gtf $gtf \
        -fasta $uroborus_ref \
        unmapped.sam $accepted_hits &> uroborus_logs.txt
        
        mv circRNA_list.txt ${base}.txt
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
        'mapsplice',
        'uroborus'
        ]
}
