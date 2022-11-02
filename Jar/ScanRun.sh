#!/bin/bash
start=$SECONDS
readonly analysisdir=./analysis/
readonly readscandir=${analysisdir}readscan/
readonly mappingdir=${analysisdir}mapping/
readonly umidir=${analysisdir}umis/
readonly fastqdir=./fastq_pass/
mkdir $analysisdir
mkdir $readscandir
mkdir $umidir
mkdir $mappingdir

nice -n 10 java -jar -Xmx250G ~/apps/new/NanoporeBC_UMI_finder-2.1.jar scanfastq  -d $fastqdir,/data/20221011_D534_6_BiopINT_CTRL/20221011_D534_6_BiopINT_CTRL/20221011_1052_2G_PAM23499_8fff69c0/fastq_pass,/data/20221011_D534_6_BiopINT_CTRL/20221011_D534_6_BiopINT_CTRL/20221011_1429_2G_PAM23499_c7bab6e0/fastq_pass  -o $readscandir  --bcEditDistance 1 --compress 
rm -r ${readscandir}/failed
count=$(find ./analysis/readscan/passed/ -maxdepth 1 -type f|wc -l)
echo $count fastq files found
find ${readscandir}/passed/ -type f -name '*' | xargs pigz -dc |  pigz > ${analysisdir}fastq_pass.fastq.gz
#cat ./analysis/readscan/passed/* > ./analysis/fastq_pass.fastq.gz
#pigz -dc ./analysis/readscan/passed/* |  pigz > ./analysis/fastq_pass.fastq.gz
minimap2 -ax splice -uf --sam-hit-only -t 90 --junc-bed /data/REFERENCES/minimap/hg38/gencode.v38.annotation.chr.bed \
 /data/REFERENCES/minimap/hg38covid/hg38covid.mmi ${analysisdir}fastq_pass.fastq.gz | \
   samtools view -bS -@ 60 - | samtools sort -m 2G -@ 90 -o ${mappingdir}passed.bam -&& samtools index ${mappingdir}passed.bam

rm ${analysisdir}fastq_pass.fastq.gz

nice -n 10 java -jar  -Xmx250G ~/apps/new/NanoporeBC_UMI_finder-2.1.jar assignumis  --inFileNanopore ${mappingdir}passed.bam -o ${umidir}passedParsed.bam   --annotationFile /data/REFERENCES/minimap/hg38covid/cellranger-GRCh38-2020-A-covid.refFlat

end=$SECONDS
duration=$(( end - start ))
echo "took $duration seconds to complete"
