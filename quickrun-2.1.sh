#!/bin/bash

start=$SECONDS
readonly analysisdir=${PWD}/outputdir_a/
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
minimap2=`which minimap2`
samtools=`which samtools`

if [ -z "$java" ] || [ -z "$samtools" ] || [ -z "$minimap2" ]
then
    echo -e "\nMissing path to required softwares:"
    echo -e "\tjava=$java"
    echo -e "\tsamtools=$samtools"
    echo -e "\tminimap2=$minimap2"
    echo -e "\nPlease update your \$PATH and rerun.\n\n"
    exit
fi

### Step 1 - Adapter search, stranding, barcode assignment to reads ###
$java -jar -Xmx4G Jar/NanoporeBC_UMI_finder-2.1.jar scanfastq -d $fastqdir -o $readscandir --bcEditDistance 1 --compress 
find ${readscandir}/passed/ -type f -name '*' | xargs pigz -dc |  pigz > ${analysisdir}fastq_pass.fastq.gz

### Step 2 - Mapping ###
$minimap2 -ax splice -uf --sam-hit-only -t 4 --junc-bed Data/gencode.v38.chr12.bed Data/chr12.fa.gz ${analysisdir}fastq_pass.fastq.gz | $samtools view -bS -@ 4 - | $samtools sort -m 2G -@ 4 -o ${mappingdir}passed.bam -&& $samtools index ${mappingdir}passed.bam

### Step 3 - UMI assignment to Nanopore SAM records ###
$java -jar  -Xmx4G Jar/NanoporeBC_UMI_finder-2.1.jar assignumis --inFileNanopore ${mappingdir}passed.bam -o ${umidir}passedParsed.bam --annotationFile Data/gencode.v38.chr12.refFlat

### Step 4 - option a ###
$java -jar -Xmx4g Jar/Sicelore-2.1.jar SelectValidCellBarcode -I ${readscandir}BarcodesAssigned.tsv -O ${siceloredir}barcodes.csv -MINUMI 1 -ED0ED1RATIO 1
$java -jar -Xmx4g Jar/Sicelore-2.1.jar IsoformMatrix -I ${umidir}passedParsed.bam -REFFLAT Data/gencode.v38.chr12.refFlat -CSV ${siceloredir}barcodes.csv -OUTDIR $siceloredir -PREFIX reads -ISOBAM true -VALIDATION_STRINGENCY SILENT
$samtools index ${siceloredir}reads_isobam.bam

end=$SECONDS
duration=$(( end - start ))
echo "took $duration seconds to complete"
