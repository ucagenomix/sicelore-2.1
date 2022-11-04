#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

workflow {

	STEP1_readscan()
	STEP1_validbarcodes(STEP1_readscan.out.scancsv)
	STEP2_mapping(STEP1_readscan.out.fastqgz)
	STEP3_umis(STEP2_mapping.out.mappingbam)
	
	STEP4a_matrix(STEP1_validbarcodes.out.csv, STEP3_umis.out.parsedbam)
	
	STEP4b_addsequence(STEP3_umis.out.parsedbam, STEP1_readscan.out.fastqgz)
	STEP4b_getchrs(STEP4b_addsequence.out.parsedbamseq) | splitText | map{it -> it.trim()} | STEP4b_splitbam | STEP4b_consensus | STEP4b_concatenate | collectFile | STEP4b_deduplicate | STEP4b_mapping | STEP4b_addtags | STEP4b_addgenes
	STEP4b_matrix(STEP1_validbarcodes.out.csv, STEP4b_addgenes.out.bam)
}

process STEP1_readscan {
    
    output:
    path './passed/ReadScanner.html'	, emit: scanhtml
    path './passed/BarcodesAssigned.tsv'	, emit: scancsv
    path 'fastq_pass.fastq.gz'	, emit: fastqgz
    
    publishDir "${params.outdir}/${params.scandir}", mode: 'copy'
    
    """
    mkdir ./passed
    $params.java -jar $params.javaXmx $params.nanopore scanfastq -d $params.fastqdir -o ./passed --bcEditDistance 1 --compress
	find ./passed/passed/ -type f -name '*' | xargs pigz -dc |  pigz > fastq_pass.fastq.gz
	"""
}

process STEP1_validbarcodes {
   
    input:
	path(barcodeassigned)

    output:
    path 'BarcodesValidated.csv'	, emit: csv
    
    publishDir "${params.outdir}/${params.scandir}", mode: 'copy'
    
    """
	$params.java -jar $params.javaXmx $params.sicelore SelectValidCellBarcode I=$barcodeassigned O=BarcodesValidated.csv MINUMI=$params.MINUMI ED0ED1RATIO=$params.ED0ED1RATIO
	"""
}

process STEP2_mapping {
   
    input:
	path(fastqgz)

    output:
    path 'passed.bam'	, emit: mappingbam
    
    publishDir "${params.outdir}/${params.mappingdir}", mode: 'copy'
    
    """
	$params.minimap2 -ax splice -uf --sam-hit-only -t $params.max_cpus --junc-bed $params.juncbed $params.minimapfasta $fastqgz | $params.samtools view -bS -@ $params.max_cpus - | $params.samtools sort -m 2G -@ $params.max_cpus -o passed.bam -&& $params.samtools index passed.bam
	"""
}

process STEP3_umis {
    
    input:
	path(mappingbam)
 
    output:
    path 'passedParsed.bam'	, emit: parsedbam
    path 'passedParsed.bai'	, emit: parsedbai
    
    publishDir "${params.outdir}/${params.umisdir}", mode: 'copy'
    
    """
 	$params.java -jar $params.javaXmx $params.nanopore assignumis --inFileNanopore $mappingbam -o passedParsed.bam --annotationFile $params.refflat
	"""
}

process STEP4a_matrix {

    input:
 	path(csv)
 	path(bam)
 	
 	output:
 	path("*.txt")			, emit: isotxt
 	path("*_isobam.bam")	, emit: isobam
 	path("*_isobam.bam.bai")	, emit: isobai
 	
 	publishDir "${params.outdir}/${params.matrixdir}", mode: 'copy'
 	
 	
    """
	$params.java -jar $params.javaXmx $params.sicelore IsoformMatrix I=$bam REFFLAT=$params.refflat CSV=$csv OUTDIR=. PREFIX=$params.PREFIX CELLTAG=$params.CELLTAG UMITAG=$params.UMITAG GENETAG=$params.GENETAG TSOENDTAG=$params.TSOENDTAG POLYASTARTTAG=$params.POLYASTARTTAG CDNATAG=$params.CDNATAG USTAG=$params.USTAG RNTAG=$params.RNTAG MAPQV0=$params.MAPQV0 DELTA=$params.DELTA METHOD=$params.METHOD ISOBAM=$params.ISOBAM AMBIGUOUS_ASSIGN=$params.AMBIGUOUS_ASSIGN VALIDATION_STRINGENCY=SILENT
    $params.samtools index -@ $params.max_cpus ${params.PREFIX}_isobam.bam
    """
}

process STEP4b_addsequence {

    input:
	path(bam)
	path(fastqgz)
 
    output:
    path 'parsedbamseq.bam'	, emit: parsedbamseq
    
    publishDir "${params.outdir}/${params.matrixconsdir}", mode: 'copy'
    
    """
    $params.java -jar $params.javaXmx $params.nanopore tagbamwithread --inFastq $fastqgz --inBam $bam --outBam parsedbamseq.bam --readTag US --qvTag QS
	"""
}

process STEP4b_getchrs {
    input:
	path(bam)
    
    output:
    path("chromo.csv")
    
    """
	$params.samtools view -H $bam | grep SQ | awk '{ print \$2 }' | sed 's/SN://' | grep -v 'ERCC\\|SIRV\\|phiX174' > chromo.csv
	"""
}

process STEP4b_splitbam {
    input:
 	val(chromo)
 	
 	output:
 	path 'chromosome.bam'	, emit: bam
 	
    """
    $params.samtools view -Sb $params.matrixconsdirfull/parsedbamseq.bam $chromo -o chromosome.bam
	$params.samtools index -@ $params.max_cpus chromosome.bam
    """
}

process STEP4b_consensus {  
    input:
 	path(bam)
 	
 	output:
 	path 'chr.fq'	, emit: fq
 	
    """
    $params.java -jar $params.javaXmx $params.sicelore ComputeConsensus T=$params.max_cpus I=$bam O=chr.fq CELLTAG=$params.CELLTAG UMITAG=$params.UMITAG GENETAG=$params.GENETAG TSOENDTAG=$params.TSOENDTAG POLYASTARTTAG=$params.POLYASTARTTAG CDNATAG=$params.CDNATAG USTAG=$params.USTAG RNTAG=$params.RNTAG MAPQV0=$params.MAPQV0 TMPDIR=$params.tmpdir VALIDATION_STRINGENCY=SILENT MAXREADS=$params.MAXREADS MINPS=$params.MINPS MAXPS=$params.MAXPS DEBUG=$params.DEBUG
    """
}

process STEP4b_concatenate {
    input:
  	path x
  
  	output:
  	path 'consensus.fq'	, emit: cons
  
  	script:
  	"""
  	< $x cat > consensus.fq
  	"""
}

process STEP4b_deduplicate {
    
    input:
 	path(fq)
 	
 	output:
 	path 'molecules.fastq'	, emit: dedup
 	
    publishDir "${params.outdir}/${params.matrixconsdir}", mode: 'copy'

    """
 	$params.java -jar $params.javaXmx $params.sicelore DeduplicateMolecule I=$fq O=molecules.fastq SELECT=true VALIDATION_STRINGENCY=SILENT
    """
}

process STEP4b_mapping {
    
    input:
 	path(dedup)
 	
 	output:
 	path 'molecules.bam'		, emit: bam
 	path 'molecules.bam.bai'	, emit: bai

    publishDir "${params.outdir}/${params.matrixconsdir}", mode: 'copy'
 	
    """
	$params.minimap2 -ax splice -uf --sam-hit-only -t $params.max_cpus --junc-bed $params.juncbed $params.minimapfasta $dedup > molecules.sam
	$params.samtools view -Sb -@ $params.max_cpus molecules.sam -o deleted.bam
	$params.samtools sort -@ $params.max_cpus deleted.bam -o molecules.bam
	$params.samtools index -@ $params.max_cpus molecules.bam
    """
}

process STEP4b_addtags {
   
    input:
 	path(bam)
 	path(bai)
 	
 	output:
 	path 'molecules.tags.bam'		, emit: bam
 	path 'molecules.tags.bam.bai'	, emit: bai
 	
    publishDir "${params.outdir}/${params.matrixconsdir}", mode: 'copy'

    """
	$params.java -jar $params.javaXmx $params.sicelore AddBamMoleculeTags I=$bam O=molecules.tags.bam CELLTAG=$params.CELLTAG UMITAG=$params.UMITAG RNTAG=$params.RNTAG
	$params.samtools index -@ $params.max_cpus molecules.tags.bam
    """
}

process STEP4b_addgenes {
    
    input:
 	path(bam)
 	path(bai)
 	
 	output:
 	path 'molecules.bam'		, emit: bam
 	path 'molecules.bam.bai'	, emit: bai
 	
    publishDir "${params.outdir}/${params.matrixconsdir}", mode: 'copy'

    """
	$params.java -jar $params.javaXmx $params.sicelore AddGeneNameTag I=$bam O=molecules.tags.GE.bam REFFLAT=$params.refflat GENETAG=$params.GENETAG ALLOW_MULTI_GENE_READS=$params.ALLOW_MULTI_GENE_READS USE_STRAND_INFO=$params.USE_STRAND_INFO VALIDATION_STRINGENCY=SILENT
    $params.samtools index -@ $params.max_cpus molecules.bam
    """
}

process STEP4b_matrix {
    
    input:
	path(csv)
 	path(bam)
 	
 	output:
 	path("*.txt")			, emit: isotxt
 	path("*_isobam.bam")	, emit: isobam
 	path("*_isobam.bam.bai")	, emit: isobai
 	
    publishDir "${params.outdir}/${params.matrixconsdir}", mode: 'copy'

    """
	$params.java -jar $params.javaXmx $params.sicelore IsoformMatrix I=$bam REFFLAT=$params.refflat CSV=$csv OUTDIR=./ PREFIX=$params.PREFIX CELLTAG=$params.CELLTAG UMITAG=$params.UMITAG GENETAG=$params.GENETAG TSOENDTAG=$params.TSOENDTAG POLYASTARTTAG=$params.POLYASTARTTAG CDNATAG=$params.CDNATAG USTAG=$params.USTAG RNTAG=$params.RNTAG MAPQV0=$params.MAPQV0 DELTA=$params.DELTA METHOD=$params.METHOD ISOBAM=$params.ISOBAM AMBIGUOUS_ASSIGN=$params.AMBIGUOUS_ASSIGN VALIDATION_STRINGENCY=SILENT
    $params.samtools index -@ $params.max_cpus ${params.PREFIX}_isobam.bam
    """
}





