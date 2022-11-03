#!/bin/bash

start=$SECONDS
readonly analysisdir=${PWD}/outputdir_4b/
readonly readscandir=${analysisdir}readscan/
readonly mappingdir=${analysisdir}mapping/
readonly umidir=${analysisdir}umis/
readonly fastqdir=Data/fastq_pass/
readonly tmpdir=${analysisdir}tmp/
readonly siceloredir=${analysisdir}sicelore/

mkdir $analysisdir
mkdir $readscandir
mkdir $umidir
mkdir $mappingdir
mkdir $siceloredir
mkdir $tmpdir

# need at least Java >= 12
java=`which java`
spoa=`which spoa`
minimap2=`which minimap2`
samtools=`which samtools`

if [ -z "$java" ] || [ -z "$spoa" ] || [ -z "$samtools" ] || [ -z "$minimap2" ]
then
    echo -e "\nMissing path to required softwares:"
    echo -e "\tjava=$java"
    echo -e "\tspoa=$spoa"
    echo -e "\tsamtools=$samtools"
    echo -e "\tminimap2=$minimap2"
    echo -e "\nPlease update your \$PATH and rerun.\n\n"
    exit
fi

### Step 1 - Adapter search, stranding, barcode assignment to reads ###
$java -jar -Xmx4G Jar/NanoporeBC_UMI_finder-2.1.jar scanfastq -d $fastqdir -o $readscandir --bcEditDistance 1 --compress 
find ${readscandir}/passed/ -type f -name '*' | xargs pigz -dc |  pigz > ${analysisdir}fastq_pass.fastq.gz

### Step 2 - Mapping ###
$minimap2 -ax splice -uf --sam-hit-only -t 4 --junc-bed Data/gencode.v38.chr12.bed Data/chr12.fa.gz ${analysisdir}fastq_pass.fastq.gz | $samtools view -bS -@ 60 - | $samtools sort -m 2G -@ 4 -o ${mappingdir}passed.bam -&& $samtools index ${mappingdir}passed.bam

### Step 3 - UMI assignment to Nanopore SAM records ###
$java -jar  -Xmx4G Jar/NanoporeBC_UMI_finder-2.1.jar assignumis --inFileNanopore ${mappingdir}passed.bam -o ${umidir}passedParsed.bam --annotationFile Data/gencode.v38.chr12.refFlat

### Step 4 - option b ###
$java -jar -Xmx4g Jar/Sicelore-2.1.jar SelectValidCellBarcode I=${readscandir}BarcodesAssigned.tsv O=${siceloredir}barcodes.csv MINUMI=1 ED0ED1RATIO=1
$java -jar -Xmx4g Jar/NanoporeBC_UMI_finder.jar tagbamwithread --inFastq ${analysisdir}fastq_pass.fastq.gz --inBam ${umidir}passedParsed.bam --outBam ${siceloredir}passedParsedWithSequences.bam --readTag US --qvTag QS
$java -jar -Xmx4g Jar/Sicelore-2.1.jar ComputeConsensus I=${siceloredir}passedParsedWithSequences.bam O=${siceloredir}molecules.fastq T=4 TMPDIR=${tmpdir}

$minimap2 -ax splice -uf --sam-hit-only -t 4 --junc-bed Data/gencode.v38.chr12.bed Data/chr12.fa.gz ${siceloredir}molecules.fastq > ${siceloredir}molecules.sam
$samtools view -Sb ${siceloredir}molecules.sam -o ${siceloredir}unsorted.bam
$samtools sort ${siceloredir}unsorted.bam -o ${siceloredir}molecules.bam
$samtools index ${siceloredir}molecules.bam

$java -jar -Xmx4g Jar/Sicelore-2.1.jar AddGeneNameTag I=${siceloredir}molecules.bam O=${siceloredir}molecules.GE.bam REFFLAT=Data/gencode.v38.chr12.refFlat TAG=GE ALLOW_MULTI_GENE_READS=true USE_STRAND_INFO=true VALIDATION_STRINGENCY=SILENT
$samtools index ${siceloredir}molecules.GE.bam

$java -jar -Xmx44g Jar/Sicelore-2.1.jar AddBamMoleculeTags I=${siceloredir}molecules.GE.bam O=${siceloredir}molecules.GE.tags.bam
$samtools index ${siceloredir}molecules.GE.tags.bam

$java -jar -Xmx4g Jar/Sicelore-2.1.jar IsoformMatrix I=${siceloredir}molecules.GE.tags.bam REFFLAT=Data/gencode.v38.chr12.refFlat CSV=${readscandir}BarcodesAssigned.tsv OUTDIR=$siceloredir PREFIX=mols ISOBAM=true VALIDATION_STRINGENCY=SILENT

end=$SECONDS
duration=$(( end - start ))
echo "took $duration seconds to complete"
