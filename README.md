
<img src="sicelore.png">


SiCeLoRe (Single Cell Long Read) is a suite of tools dedicated to cell barcode / UMI (unique molecular identifier) assignment and bioinformatics analysis of highly multiplexed single cell Nanopore or PacBIo long read sequencing data.

# Installation

just copy files.

requires: 

* Java >= 12

* [spoa](https://github.com/rvaser/spoa)

* [fastp](https://github.com/OpenGene/fastp)

* samtools

# Workflow

[1) Adapter search, stranding, barcode assignment to reads](#nanopore-scan)

 

              Mapping against reference genome

 

[2) UMI assignment to Nanopore SAM records](#umi-assignment)


 
[6) Consensus sequence calculation for RNA molecules (UMIs)](#ComputeConsensus)

Generates consensus sequence for transcripts from multiple reads for UMI.

[7) Mapping of molecules consensus sequences to the reference genome with minimap2](#molecule-consensus-mapping)

Consensus sequences are mapped to the reference genome

[8) Tag molecule SAM records with gene names, cell barcodes and UMI sequence](#molecule-tagging)

Adds gene names, cell barcode and UMI sequence. Required for [cell x genes/isoforms/junctions] matrices generation

[9) Transcript isoform expression quantification](#IsoformMatrix)

Identifies matching Gencode transcript isoforms and generates [cell x genes/isoforms/junctions] matrices.

[10) SNP calling](#snp-calling)

Calling nucleotide polymorphism cell by cell

[11) Fusions gene calling](#fusion-calling)

Detecting fusion transcripts cell by cell

[12) Novel transcripts isoforms detection](#collapse-model)

Indentifying novel transcripts isoforms 

## Authors

* Kevin Lebrigand <[lebrigand@ipmc.cnrs.fr](mailto:lebrigand@ipmc.cnrs.fr)>

[![Twitter Follow](https://img.shields.io/twitter/follow/kevinlebrigand.svg?style=social&logo=twitter)](https://twitter.com/kevinlebrigand)

* Rainer Waldmann <[waldmann@ipmc.cnrs.fr](mailto:waldmann@ipmc.cnrs.fr)>

<a id="parsing-illumina-data"></a>

## Quick run analysis

We provide test data as a subsampling of reads for the Mus musculus Clta locus for the 190 cells dataset.
It requires java 1.8 (JAVA_HOME), minimap2, samtools in your PATH as well as racon and poa (including blosum80.mat in same folder) for consensus calling part.
This test script should takes under 5mn to run, output files are located in ./output_dir directory.

```

git clone https://github.com/ucagenomix/sicelore.git
cd sicelore
chmod +x quickrun.v2.2.sh
dos2unix quickrun.v2.2.sh
export JAVA_HOME=<path to Java 1.8>
export PATH=$PATH:<minimap2path>:<samtoolspath>:<spoapath>
./quickrun.v2.2.sh

```


<a id="nanopore-scan"></a>

 

# 1) Scan Nanopore reads - assign cell barcodes.

 

1) Scans for chimeric cDNA generated during library preparation and splits chimeric reads . Reads that contain two adjacent internal (>200 nt from end)  polyA/T sequences flanked by an adapter (pAad), one  internal pAad with adjacent TSO or two adjacent TSOs are split into separate reads.

 

2) Scans the Nanopore fastq reads for poly(A/T) and adapter sequence and generates stranded (forward) reads for reads with found polyA and adapter.

Scans by default for >= 15 nt. polyA (or T) with >= 75% As within 100 nt from both ends of the read. If poly(A) was found, Searches for a 10xGenomics adapter sequence "CTTCCGATCT" downstream of the poly(A).

Scan for  polyA is optional and can be skipped (--noPolyARequired option) for samples that don't have polyA tails (e.g. RACE PCR combined with 5' barcoding).  Splitting of chimeric transcripts is skipped when no polyA is searched.

3) Assigns cell barcodes: First, a list of used barcodes is generated. The software only uses high quality reads during the first pass. A list of all possible 10xGenomics barcodes (e.g. 3M-february-2018.txt.gz from Cellranger) should be supplied to guide the definition of used barcodes. If Illumina data are available the file with the list of used barcodes can be supplied (--cellRangerBCs < path to file >) -> the first pass will be skipped and the supplied barcode list will be used.

In a second pass the software assigns the cell barcodes. With the current read qualities a maximal edit distance of one is largely sufficient (--bcEditDistance 1).  Barcode assignment accuracy is typically > 99.9% for edit distance 1, edit distance 2 yields about 5 - 10% more barcode assigned reads but the tradeoff is a lower assignment accuracy.

 

Adds info required for the Illumina-guided BC UMI assignment (e.g. sequence upstream of polyA) to the read name. Info on polyA, adapter position and info on found cellBC are added.

 

When strandedness of read can be determined (adapter and polyA[if opted for] were found) the barcode assigned read is written stranded (forward) into a "pass" folder. Failed reads are written unstranded into a "failed" folder.

 

Scan of 100M reads takes about 80 min on the 96 core Promethion. In our experience running it does not interfere with ongoing sequencing runs. When started in parallel with sequencing runs we use nice -n +10 to decrease priority and allow a max of 200G of RAM ( -Xmx200G)

 

**Use the reads in the "pass" folder to proceed**

 

Generates a html file with stats.

 

## Required files

 

* NanoporeBC_UMI_finder.jar

* Libraries in the ./lib folder

* Config file: config.xml (Most default settings can be changed there). Structure of the config file needs revision - will be updated in the near future.

* If no config.xml file is found in the current path (working directory), the software takes the default config file from the directory where the applications (jars) are installed.

*

* If a 10xGenomics system was used, a file with a list of all possible barcodes (e.g.  3M-february-2018.txt.gz or 737K-august-2016.txt from Cellranger (gz compressed is o.k.) should be provided in the the same folder as the jar file. The path of the file must be in the config.xml  e.g. \<fileWithAllPossibleTenXbarcodes\>3M-february-2018.txt.gz\</fileWithAllPossibleTenXbarcodes\> will use the 3M-february-2018.txt.gz file. 

 

 

 

## Usage

 

```bash

 

java -jar -Xmx300G <path>/NanoporeBC_UMI_finder.jar scanfastq -d <directory to start recursive search for fastq files> -o outPutDirectory --bcEditDistance 1

 

-Xmx : allow the RAM you have available -Xmx300G is an example for running it on the Promethion

```

## Comments

Input: fastq files or gzipped fastq files.  Use directly the fastq files generated by the Nanopore sequencer, **don't merge them into a single fastq file**. The software reads the fastqs highly parallelized - multiple fastqs are processed much faster than a single fastq. If compressed files are used they need to be individually compressed - no tar gz archive. Multiple input directories must be seperated by ","

Typically just provide the path to the fastq_pass folder of the sequencing data

 

## Parameters

 

## **required:**

 

**-b,--bcEditDistance**

edit distance for barcode assignment.  --bcEditDistance 1 is a good tradeoff between efficiency and specifity (typically >99.9%). Specifity can be estimated by re-running the program with the -e option where the read bc sequence gets replaced by a random sequence and comparing the number of barcode assigned reads.

 

 

**-d,--inDir**

<,> separated list of directories to start file search. starting at this directory, takes recursively (if  --nonrecursive is not set) fastq files with names that match the RegEx pattern provided in <-–pattern> option, defaults to "".{1,}\\.fastq(\\.gz)?"

Gzipped fastq files are also o.k.

 

**-o,--outDir**

path to output directory

 

### **Optional:** 

**-g,--cellRangerBCs**

If Illumina data are available a tsv file with cell barcodes (e.g. from Cellranger) can be provided.

When provided the software won't define the used cell barcodes, it  assigns to one of the barcodes in the provided list.

If this list of barcodes is not provided, the software uses the **10xgenomics  list of all possible barcodes** (e.g.  3M-february-2018.txt.gz from Cellranger (gz compressed is o.k.) and defines the used barcodes. The name of the file with all possible must be in the config.xml :\<fileWithAllPossibleTenXbarcodes\>3M-february-2018.txt.gz\</fileWithAllPossibleTenXbarcodes\>.  The file has to be in the application root directory.

If neither file is provided (e.g. systems other than 10xGenomics), the software assigns barcodes without any guidance (still experimental).

 

**-h,--fivePbc** (no arguments)

Set this flag if **5' barcoding** was used.

Important: Set the right \<fileWithAllPossibleTenXbarcodes> in the config.xml

Current 3' barcoding kits use 3M barcode set:

\<fileWithAllPossibleTenXbarcodes>3M-february-2018.txt.gz\</fileWithAllPossibleTenXbarcodes>

5' barcoding might still use the 737K set

\<fileWithAllPossibleTenXbarcodes>737K-august-2016.txt\</fileWithAllPossibleTenXbarcodes>

If the efficiency of barcode assignment is low, check whether a  wrong whitelist was used.

 

**-y,--noPolyARequired**

No arguments. Set this flag only for targeted sequencing by RACE PCR when 5' barcoding was used. RACE PCR products won't have a polyA if 5' barcoding was used. Don't set it for 3' barcoding.

 

**-c,--compress**

compress fastq output files.

 

**-v,--pattern**

RegEx pattern for fastq input files , Defaults to: ".{1,}\\.fastq(\\.gz)?"

 

**-e,--randomBarcode** (no arguments)

Performs a random Barcode Simulation. For accuracy testing. When set barcode sequences are replaced by random sequences during barcode assignment. No argument.

 

**-w,--windowAT**

window to search for  AT, defaults to 150. Default can be modified:  <windowSearchForPolyA> of config.xml

 

**-p,--polyAlength**

min length of polyA.  Defaults to the value provided in \<polyATlength\> of the config.xml  

 

**-f,--fractionAT**

Min fraction of A or Ts in polyA or polyT tail respectively. Defaults to \<fractionATInPolyAT\> in config.xml

 

**-k,--skipNfastqs**

skips the given number of fastq files. Can be used in combination with "-z" to define batches if jobs are distributed on a cluster. Barcode assignment is rather fast (< 120 min on the Promethion). Splitting into batches is normally not required.

 

**-z,--onlyNfastqs**

only processes the given number of fastq files. Can be used in combination with "-k" to define batches if jobs are distributed on a cluster.


**-l,--logFile**

path to log file, defaults to terminal.

 

**-n,--nonrecursive**

if set searches for fastq files only in specified folder(s), ignores subdirectories. No argument.

 

## Output:

- stranded barcode assigned reads in a "passed" folder in the output directory

- unassigned reads (not stranded) in "failed" folder

- html file with stats and a knee plot (reads per cell). The lower curve of the knee plot are just the reads used for barcode definition. The upper curve is for all reads.

- BarcodeList.tsv: This file is generated during the first pass (definition of used barcodes). First column is the barcode sequence, the second column are the number of high quality reads that match this barcode without mismatches. The 10x Genomics barcode whitelists contain barcodes that are just one Levenshtein distance apart. Barcodes with a sequencing error can be falsely assigned to such a barcode.  To avoid this, the following strategy is used. If, for two barcodes with a Levenshtein distance of just one, one barcode has at least 30x more reads than the other barcode the other barcode is not considered as valid barcode. The third column shows barcodes in the whitelist that have a Levenshtein distance of just one and in parenthesis the number of full matching barcode sequences.

- BarcodesAssigned.tsv: Barcode assignment data for all reads. First column , the barcode sequence; 2nd column, the number of reads. The following columns show for each barcode the number of reads assigned at different Levenshtein distances.

Examples of read name extensions:

_FWD_PS=566_PE=590_AE=619_bc=TCCGATCGTGCCAAGA_ed=0_ed_sec=2147483647_bcStart=618_bcEnd=603_rk=2987_X=AAAAAAAAAAAATGGCGTGTATTGTCTTGGCACGATCGGAAGA_Q=27.1

FWD: read was forward

PS=566: polyA start in stranded read

PE=590: polyA end in stranded read

AE=619: Adapter end

bc=TCCGATCGTGCCAAGA assigned barcode sequence

ed=0: Levenshtein distance of barcode

ed_sec=2147483647:  Levenshtein distance of secondary barcode match, INTMAX(2147483647) means none found

bcStart=618 :barcode start

bcEnd=603 :barcode end

rk=2987 ;rank of this barcode . Rank 1 is barcode with highest number of reads

X=AAAAAAAAAAAATGGCGTGTATTGTCTTGGCACGATCGGAAGA  contains 3 bases of adapter add 3' + upstream sequences

Q=27.1 : mean quality of the sequence in X=

**sp2**_REV_PS=1257_PE=1305_AE=1327_T=40_bc=GAGTGAGGTTGGGTAG_ed=1_ed_sec=2147483647_bcStart=1326_bcEnd=1311_rk=3883_X=AAAAAAAAAAACAAACCAAGTAACCAACCCAACCTCACTCAGA_Q=15.9

sp2 indicates this is the second read obtained after splitting the initial chimeric read

 

## Options in config.xml

 

```xml

              <--        reads below this length will be discarded-->

         <minReadLength>200</minReadLength>

               <!--        Levenshtein distance (ED) with which to merge barcodes when creating list of used barcodes

         If BC a has an ED from BC b equal or less this value and minCountFold x less reads it will be merge to Barcode b and not added to the list of used barcodes

         if not provided or null will be set to barcode edit distance providid in the command line

         should be null or at least as big as the ED for barcode search-->

        <mergeBCsED>null</mergeBCsED>

        <!--         minamal fold more counts in BC a than in BC b to merge b into a during generation of list of used barcodes-->

        <minCountFold>20</minCountFold>

        <!--        Minimal mean BC qv to consider barcode for barcode list-->

        <minMeanBCqv>10</minMeanBCqv>

        <!--        Minimal mean Read qv to consider barcode for barcode list-->

        <minMeanReadqv>10</minMeanReadqv>

        <!--        Minimal consecutive 3p adapter matches to consider barcode for barcode list-->

        <minAdapter3pMatches>8</minAdapter3pMatches>

       

        <!--         cell barcodes with read counts this fold below the max reads for a cell won't be assigned-->

         <cellsWithReadsnFoldBelowMaxToKeep>500</cellsWithReadsnFoldBelowMaxToKeep>

        <!--        search for barcode at position predicted from adapter finding +/- this value-->

        <testPlusMinusPos>2</testPlusMinusPos>

        <!--File with all possible 10x barcodes, one BC per line, can contain -1 after the BC sequence , can be gz compressed , SHOULD BE IN APPLICATION ROOT DIRECTORY-->

       <!--   <fileWithAllPossibleTenXbarcodes>3M-february-2018.txt.gz</fileWithAllPossibleTenXbarcodes>-->

        <fileWithAllPossibleTenXbarcodes>737K-august-2016.txt</fileWithAllPossibleTenXbarcodes>

                             Prefixes in Read names:

        <!--Prefix before polyA position, first A after cDNA -->

        <pa_start_prefix>PS=</pa_start_prefix>

        <!--Prefix before polyA position, last A after cDNA -->

        <pa_end_prefix>PE=</pa_end_prefix>

        <!--Prefix before Adapter position, last adapter base before cellBC -->

        <adapter_pos_prefix>AE=</adapter_pos_prefix>

        <!--Prefix before TSO position in read, last TSO base before cDNA-->

        <tso_pos_prefix>T=</tso_pos_prefix>

        <!--        Prefix before seq is between last adapter base and 42 bases further  downstream

        sequence is forward on stranded read-->

        <seq_prefix>X=</seq_prefix>

        <!--        Prefix before mean qv -->

        <qv_prefix>Q=</qv_prefix>

        <!--number of bases of adapter (the one b4 the BC) to include into the read sequence (adapter end, bc, umi, start of polyA) in readname-->

        <nbasesOfAdapterSeqInReadname>3</nbasesOfAdapterSeqInReadname>

```

 

 

 

 

# 3) Mapping

Map the reads from the pass folder generated in the barcode assignment step against the reference genome.

We typically use minimap2.

UMI assignment requires sorted bam files. When you split files to start multiple UMI assignment jobs keep genomic regions together (e.g. split by chromosome)

 

 


<a id="umi-assignment"></a>

# 4) UMI assignment

 

- Searches for  UMIs in  Nanopore BAM files.

-Gets sequence info and adapter, polyA position from the info added to the read name during [Nanopore read scan](#nanopore-scan).

Groups and assigns UMIs by clustering.

## Strategy

 

Uses a clustering approach that does not require Illumina data. UMI read sequences for each cell that match to the same genomic region (default within 500 nt, <max_GenomeDistance_forGrouping> in config.xml) are clustered with an ED limit of two. 

For each cluster the assigned UMI is: 

a) if  two reads were found for the UMI, the UMI sequence with the highest read quality defines the UMI sequence for the cluster.

b) if > 2 reads were found for the UMI the center of the cluster (least square sum distance from others) is assigned as the UMI

c) The reads for which the UMIs don't cluster with any other UMI for the same genomic region at ed <=2 are considered as UMIs backed by just one read. -> the read sequence is assigned as UMI.

 

 

# Considerations

The input bam file must be sorted. **if you use separate jobs, don't split records that correspond to the same genomic region over several batches ** For clustering, all reads corresponding to the same genomic region must be treated in the same batch. If you split the reads for mapping, merge the Bam files, sort them and then split them again (e.g. one bam per chromosome). 

 

For UMI clustering, records corresponding to the same genomic region must be kept in RAM.  The number of records treated in each batch can be defined in the config.xml (\<sam_records_chunk_size\>).  A rather big chunk size is important for targeted sequencing where a big fraction of the records match the same genomic region and should be kept together.

 

I'll add a BAM splitter in a later release, that splits a Bam and keeps data for each cell barcode identified during the read scan together.

 

RAM requirement:  one GB of RAM per 10,000 - 20,000 records (\<sam_records_chunk_size\>) should be o.k. I did not test the RAM requirement extensively since I was running it mostly directly on the machine that pilots the Promethion (384 Gb RAM, 96 cores) in just one batch. I use a batch size of 1 - 2 million.  Just allow the RAM you have available. **Adjust the JVM -Xmx option accordingly.**

 

## Bam output file

 

 

Cell barcodes are in the BC tag, UMI sequence is in U8 tag in the output BAM file specified  in the -o option,

 The tags can be customized in the config.xml.

 Info on barcode and UMI edit distances and other info is written to various tags. See "config.xml" for more details and customization.


# Required files

 

* Nanopore_BC_UMI_finder.jar

 

* Libraries in the ./lib folder

 

* Config.xml: Config file to set various parameters.

If no config file is found in the current path (working directory) , the software takes the config files from the directory where the application (jar) is. Path to the config file can be provided on the command line (-c,--config).

 

* bcMaxEditDistances.xml and umiMaxEditDistances.xml: contain max edit distances for various barcode and UMI false assignment percentages respectively.

 

* refFlat file for gene assignment. If supplied will generate a cell/gene umi count table and add the GE tag to SAM record.  Uses the  “TagReadWithGeneExonFunction” from the DropSeq package. Not required if the GE tag was added to BAM records with the Sicelore package. 

 

 

## Usage

 

```bash

 

java -jar -Xmx30g NanoporeBC_UMI_finder.jar assignumis --inFileNanopore <Nanopore Bam>

--outfile <cell bc and UMI assigned output bam file>

```

Not all of the options in the config.xml are currently used. Some options are for Illlumina guided assignment.

Structure of the config file is currently being revised.

 

## Command line arguments (some override defaults in config xml).

 

### Required

 

* **-i,--inFileNanopore**

 

Nanopore input BAM file. The genome aligned Nanopore BAM input file with the read sequence in “US” tag and the gene name (if any) in the “GE” tag. All tags can be changed in “config.xml”.

 

* **-o,--outfile**

 

Output bam file.

 

### Optional (some defaults can be modified in config xml)

 

* **-p,--fivePbc** (no arguments)

 

Set this flag when **5' barcoding** was used. Software default is 3' barcoding.

 

 

* **-a,--annotationFile**  <name.refFlat> <a id="annotationFile"></a>

 

path to refflat or GTF file.  Name has to end with ".refFlat" , ".refflat" or ".gtf".  Can be compressed, ".bz2" or ".gz" file name extension.  If supplied will generate a cell/gene umi count table and add the GE tag to SAM record.  Uses the  “TagReadWithGeneExonFunction” from the DropSeq package (https://github.com/broadinstitute/Drop-seq ).

 

* **-c,--config**

 

path to config.xml.  By default tests current path (working directory) , if not found, the software takes the config files from the directory where the application (jar) is.

 

* **-l,--logFile**

 

If not provided will create a default log with same name as outfile with “.log” appended in folder where output bam goes.

 

* **-e,--randomBarcode; -f,--randomUMI** (no arguments)

 

For benchmarking. performs random Barcode or UMI Simulation. If a previously scanned and BC and UMI annotated file (the file with all entries not just the entries with found BC/UMIs) is rescanned with this option the max edit distance is the ED of the previously found match if any. Provide the barcode and UMI assigned BAM file (the file with all sam entries not just the BC/UMI assigned) as input BAM for the random barcode or umi scan (retrieves edit sdistance info for found BCs and UMIs there).

 

* **-g,--ONTgene**

 

2 char SAM attribute for gene name in Nanopore SAM. Default: GE (from config.xml)

 

* **-h,--help**

 

displays list of option.

 

 

### SAM tags in output file

 

SAM tags in output bam can be modified in the config xml. If entry for tag is commented out in XML it won’t be printed.

By default cell barcode is in BC tag (\<CELL_BC\>) in config.xml. UMI is in U8  (\<UMI_sequence\>) tag.

 

Other info is added to the following tags by default:

### Infos from read scan:

positions are 1-based, in <> the tag in config.xml is shown where it can be customized.

**RE** Read was reverse complemented (read scan) <ReadReversed>

**PS**  polyA start  <POLYAT_START>

**PE**  polyA end  <POLYAT_END>

**AE** Adapter end facing cell BC  <ADAPTER_END>

**TE** 3' end of TSO <TSO_END>

**BU** BC sequence from readscan <BC_SEQ_READSCAN>

**BV**  start of cellBC <BC_SEQ_READSCAN_BEGIN>

**BE**  end of cellBC <BC_SEQ_READSCAN_END>

**BW**  BC edit distance read scan <BC_SEQ_READSCAN_ED>

**BX**  ED of second best readscan BC match <BC_SEQ_READSCAN_ED_SECOND>

**BH** Rank of cell barcode (higher rank-> lower read number) <CELL_BC_SEQ_FROM_READSCAN_RANK>

 

### Umi Finding

**UC**   if set indicates that UMI is from clustering                   

**UZ**   UMI is just the post BC read sequence               

**U1** UMI edit distance (mutation cycles)

**U2** UMI edit distance (mutation cycles) of second best match (if any)

**UB** UMI start (1 based)

**UE** UMI end (one based)

**U8** UMI sequence

**U7** UMI read sequence


<a id="computeConsensus"></a>

## 6) Generate consensus sequences

**use ComputeConsensus (sicelore.v2.jar)**

The pipeline allows to compute the consensus sequence for molecule detected  in input .bam file. First step is loading the molecules, the cDNA sequence is defined as [tsoEnd(TE tag) ... umiEnd(UE tag)] for consensus sequence computation.

Briefly, each molecule is processed as follows depending the number of reads the molecule has: (i) just one read per molecule (UMI), the consensus sequence is set to the read sequence; (ii) 2 reads per molecule, the consensus sequence is set as the cDNA sequence of the best mapping read according to the "de" minimap2 SAMrecord tag value; (iii) More than two reads per molecule, a consensus sequence is comppute using [spoa](https://github.com/rvaser/spoa) multiple alignment using a set of reads for the molecule. 

By default the top MAXREADS (default=20) reads having the lowest divergence to the reference ("de" minimap2 SAMrecord tag value) are taken for consensus calling.



The speed of consensus sequence computation is dependent of the sequencing depth which induce a low/high number of multi-reads molecules for which a consensus sequence is computed. In a standard whole transcriptome experiment having roughly 50% of the molecules having more than 2 reads the speed is about 600k UMIs/hour on a 20 core compute node using 20 threads. For time calculation optimization, this step could be parrallelized, for instance on a per chromosome basis, and dispense on a calcul cluster.

SPOA executable needs to be in your PATH variable.


**splitting bam by chromosomes**

```

@jobs = ('1','2','X','MT','3','4',..all chromosomes..);
for($i=0; $i<@jobs; $i++){
    samtools view -Sb GEUS10xAttributes.umifound.bam ".$jobs[$i]." -o GEUS10xAttributes.umifound.chr".$jobs[$i].".bam
    samtools index GEUS10xAttributes.umifound.chr".$jobs[$i].".bam
}

```

**ComputeConsensus pipeline**

**INPUT=** (required)

GEUS10xAttributes.umifound.chr1.bam file containing cellBC (BC tag) and UMI (U8 tag) assigned Nanopore SAM entries with the gene mame in the IG tag. SAMrecords lacking any of those 3 required fields are discarded from further analysis.

**OUTPUT=** (required)

Consensus sequence in fastq format for all molecules detected. Name of each molecules set to CellBC(BC), UMIs(U8) and read number RN) ">BC-UMI-RN"

**THREADS=,T=** (required)

Number of threads for multi-threading (typically number of cores of compute node)

**MAXREADS=,** (required)

Maximum number of reads per UMI to use for consensus sequence calling (default=20)

**TMPDIR=** (required)

Full path to temporary directory

**example below is for chromosome 1, repeat for all chromosomes**

```

java -jar -Xmx44g sicelore.v2.jar ComputeConsensus I=GEUS10xAttributes.umifound.chr1.bam O=molecules.chr1.fq T=20 TMPDIR=/scracth/tmp/

```

<a id="molecule-consensus-mapping"></a>

## 7) Mapping of molecules consensus sequences to the reference genome with minimap2

If ComputeConsensus step has been split per chromosome, we need to deduplicate molecules from genes having multiple copies in the genome.

First we need to concatenate all chromosomes fastq files then use ***DeduplicateMolecule*** pipeline (sicelore.v2.jar)

```

cat */molecules.chr*.fa > molecules.fa
java -jar -Xmx22g sicelore.v2.jar DeduplicateMolecule I=molecules.fq O=deduplicate.fq

```

Molecule consensus sequences can then be mapped to the reference genome to generate a molecules.bam file for further analysis.

```

minimap2 -ax splice -uf --sam-hit-only -t 20 --junc-bed junctions.bed $BUILD.mmi molecules.fq > molecules.sam
samtools view -Sb molecules.sam -o unsorted.bam
samtools sort unsorted.bam -o molecules.bam
samtools index molecules.bam

```

<a id="molecule-tagging"></a>

## 8) Tag molecule SAM records with gene names, cell barcodes, UMI sequence and reads number

**uses AddGeneNameTag (sicelore.v2.jar)** 

Add gene name tag (GE) to molecules SAM records

```

java -jar -Xmx22g sicelore.v2.jar AddGeneNameTag I=molecules.bam O=molecules.GE.bam REFFLAT=refFlat.txt TAG=GE ALLOW_MULTI_GENE_READS=true USE_STRAND_INFO=true VALIDATION_STRINGENCY=SILENT
samtools index molecules.GE.bam

```

**uses AddBamMoleculeTags (sicelore.v2.jar)**

Add CellBC (BC), UMIs (U8) and read number (RN) tags from molecule read name to molecules SAM records

``` 

java -jar -Xmx22g sicelore.v2.jar AddBamMoleculeTags I=molecules.GE.bam O=molecules.GE.tags.bam
samtools index molecules.GE.tags.bam

```

<a id="isoformmatrix"></a>

## 9) Transcript isoform expression quantification

**uses IsoformMatrix (sicelore.v2.jar)**

This step can be done at the reads level but prefer using it at the molecules level after consensus calling as describes in the protocole. In case of processing at the read level here is the algotithm used.

SAM records matching known genes are grouped by UMI and analyzed for matching Gencode v18 transcripts. SAM records with extensive non-matching sequences at either end, hard- or soft clipping are discarded (MAXCLIP parameter, defaults to > 150 nt ). To assign a UMI to a Gencode transcript a STRICT and a SCORE strategy is available(METHOD parameter). Using ***STRICT*** strategy (default mode), a read is assigned to a given transcript when it recapitulates the full exon-exon junction layout authorizing a DELTA (default 2) bases margin of added or lacking sequences at exon boundaries, to allow for indels at exon junctions and imprecise mapping by Minimap2. For each UMI, all its reads are analyzed and the UMI is assigned to the Gencode transcript supported by the majority of the reads.

If an equal number of reads supports two different Gencode transcript, the UMI is considered as ambiguous and randomly assigned to one of the top scoring isoforms (AMBIGUOUS_ASSIGN=true) or not assigned to an isoform (default mode, AMBIGUOUS_ASSIGN=false).

The ***SCORE*** strategy allows to assign an isoform to an UMI that does not recapitulate the full exon-exon structure of the isoform. Some single cell reads are incomplete due to limitations of the 10xGenomics library preparation. Those incomplete cDNAs are likely mainly due to partial degradation of RNA before reverse transcription (5' truncation) or internal priming of the oligo(dT) RT primer on A rich sequences within the transcript since reverse transcription starts at room temperature (3' truncation). This mode allows to assign such truncated reads if the available structure allows to assign it unanimously to one Genecode isoform. In SCORE mode, for every UMI each genome mapped read is examined for exon junctions matching the junctions present in the Gencode transcripts for the corresponding gene, authorizing a DELTA (2 default) bases margin of added or lacking sequences at exon boundaries, to allow for indels at exon junctions and imprecise mapping by Minimap2. For each found exon-exon junction, a score of 1 is added to each Gencode isoform carrying the junction. Finally, if just one isoform has the highest score it gets selected as the isoform for the UMI. If more than one Gencode transcript obtain the highest score, the UMI is considered as ambiguous and assigned (AMBIGUOUS_ASSIGN=true) or not assigned (AMBIGUOUS_ASSIGN=false) to one of the best matching isoforms.

**INPUT=** (required)

Barcode and UMI assigned GEUS10xAttributes.umifound.bam (or <b>molecules.GE.tags.bam</b> file) file containing the cellBC (BC tag) and UMI (U8 tag) with the gene mame in the IG tag. SAMrecords lacking any of those 3 required fields are not analyzed.

**CSV=** (required)

.csv file listing, one per line, the cell barcodes that need to be quantified (e.g. brain951.barcodes.tsv in github /barcodes/ directory)

**REFFLAT=** (required)

Can be generated base on Gencode GTF file for genome build used for mapping with ***gtfToGenePred*** from [UCSC](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/)

```

gtfToGenePred -genePredExt -geneNameAsName2 gencode.vM18.primary_assembly.annotation.gtf gencode.vM18.primary_assembly.annotation.refflat.txt
paste <(cut -f 12 gencode.vM18.primary_assembly.annotation.refflat.txt) <(cut -f 1-10 gencode.vM18.primary_assembly.annotation.refflat.txt) > refFlat.txt

```

**OUTDIR**

Output directory where output files are created

**PREFIX=** (required)

prefix used for output file names

**DELTA=** (required)

the number of extra or lacking bases allowed at exon-exon junctions (default = 2)

**MAXCLIP=**

Maximum number of extensive non-matching sequences at either end, hard- or soft clipping to call the read as chimeric and discards it (default = 150)

**GENETAG=**

Gene name tag (default = GE)

**UMITAG=**

UMI sequence tag (default = U8)

**CELLTAG=**

Cell tag (default = BC)

**METHOD=** (required)

STRICT full exon-exon structure required for assignation (default mode)

SCORE more lenient way of isoform assignation

**AMBIGUOUS_ASSIGN=**

Only active for SCORE mode, whether or not to assign an UMI that has 2 or more matching transcript model (default=false)

**ISOBAM=**

Wether or not to output a BAM file containg a flag (IT) with the transcriptID called during this process (default=false)

**output files**

**PREFIX**_cellmetrics.txt

cell by cell metrics (cellBC, nbReads, nbUmis, nbIsoformSet, nbIsoformNotSet).

**PREFIX**_genemetrics.txt

gene by gene metrics (geneId, nbUmis, nbIsoformSet, nbIsoformNotSet)

**PREFIX**_isometrics.txt

isoform by isoform metrics (geneId, transcriptId, nbExons, nbUmis)

**PREFIX**_juncmetrics.txt

exon-exon junction by junction metrics (junctionId, nbUmis)

**PREFIX**_genematrix.txt

gene level [geneId x cellBC] count matrix 

**PREFIX**_isomatrix.txt

isoform level [transcriptId x cellBC] count matrix 

**PREFIX**_juncmatrix.txt

junction level [junctionId x cellBC] count matrix 

**PREFIX**_molinfos.txt

molecule per molecules information (cellBC, UMI, nbReads, nbSupportingReads, geneId, transcriptId)

```

java -jar -Xmx44g sicelore.v2.jar IsoformMatrix I=molecules.tags.GE.bam GENETAG=GE UMITAG=U8 CELLTAG=BC REFFLAT=refFlat.txt CSV=barcodes.csv DELTA=2 MAXCLIP=150 METHOD=STRICT AMBIGUOUS_ASSIGN=false OUTDIR=. PREFIX=sicelore

```

<a id="snp-calling"></a>

## 10) Calling Single Nucleotide Polymorphisms cell by cell

**use SNPMatrix (sicelore.v2.jar)**

Consensus sequence show higher sequence accuracy than raw nanopore reads and we can now call SNPs using the generated molecules bam file (molecules.tags.GE.bam)

**INPUT=** (required)

Molecules bam file (molecules.tags.GE.bam) including cell barcode and UMI tags

**CSV=** (required)

.csv file listing, one per line, the cell barcodes in quantification (e.g. brain951.barcodes.tsv in /barcodes/ dir)

**SNP=** (required)

SNPs descriptor comma-separated .csv file (snp_description.csv), 1-position or x-positions per line as follow:

```
chromosome,position,strand,name
3,80692286,-,Gria2_RG                   // SNP call at one position
3,80706912,-,Gria2_QR                   // SNP call at one position
3,80692286|80706912,-,Gria2_RGQR        // SNP association call at 2-positions ("|" separator, no limit on number of positions)
...
```

**MINRN** (required)

Minimum read number to keep the molecule for SNP calling (default=0, means all)

**MINQV** (required)

Minimum QV score at position to keep the molecule for SNP calling (default=0, means all)

**OUTPUT=** (required)

Output directory.

**PREFIX**

Prefix for _matrix.txt/_metrics.txt/_molinfos.txt tab-delimited output text files

```

java -jar -Xmx22g sicelore.v2.jar SNPMatrix I=molecule.tags.bam MINRN=0 MINQV=0 CSV=barcodes.csv SNP=snp_description.csv O=/path/to/output/directory/ PREFIX=snp

```

**output**

**PREFIX_matrix.txt**

[Cell] x [SNP position - base] matrix molecule count.

**PREFIX_metrics.txt**

Total number of molecules per SNP position - base observed in dataset

**PREFIX**_molinfos.txt

Per molecule per SNP position information

<a id="fusion-calling"></a>

## 11) Detecting fusion transcripts cell by cell

**use ExportClippedReads and FusionDetector (sicelore.v2.jar)**

First step is exporting hard and soft clipped reads from dataset. A fusion transcript read map the genome in more than one unique location. It gives two differents minimap2 SAM records each including a portion of starting or ending clipping. Those reads might also comes from chimeric molecules produced during amplification step (i.e. PCR artefacts). Those PCR artefacts are arising mainly between transcripts having a high sequence similarity  and might preferentially happens wih highly expressed genes.

### ExportClippedReads pipeline

**INPUT=** (required)

GEUS10xAttributes.umifound.bam file containing cellBC (BC tag) and UMI (U8 tag)

**OUTPUT=,O=** (required)

Output fastq file of clipped reads.

**MINCLIP**

Minimum length of start/end hard/soft clipping needed far calling a read as clipped

```

java -jar -Xmx22g sicelore.v2.jar ExportClippedReads I=GEUS10xAttributes.umifound.bam O=clipped_reads.fastq

```

### Preparing files for Fusion detection 

```

# Minimap2 clipped reads mapping 
minimap2 -ax splice <b>-k13</b> -t 20 --junc-bed junctions.bed $BUILD.mmi clipped_reads.fastq > clipped_reads.sam
samtools view -Sb clipped_reads.sam -o unsorted.bam
samtools sort unsorted.bam -o clipped_reads.bam
samtools index clipped_reads.bam

# add gene names to Nanopore SAM records
java -jar -Xmx12g sicelore.v2.jar AddGeneNameTag I=clipped_reads.bam O=clipped_reads.GE.bam REFFLAT=refFlat.txt TAG=GE ALLOW_MULTI_GENE_READS=true USE_STRAND_INFO=true VALIDATION_STRINGENCY=SILENT
samtools index clipped_reads.GE.bam

# add read sequence and QV values to Nanopore SAM records
java -jar -Xmx12g sicelore.v2.jar AddBamReadSequenceTag I=clipped_reads.GE.bam O=clipped_reads.GEUS.bam FASTQ=clipped_reads.fastq
samtools index clipped_reads.GEUS.bam

# Add Nanopore Barcode/UMI
java -jar -Xmx22000m NanoporeBC_UMI_finder.jar -i clipped_reads.GEUS.bam -o clipped_reads.GEUS10xAttributes.umifound.bam -k parsedForNanopore_v0.2.obj --maxUMIfalseMatchPercent 1 --maxBCfalseMatchPercent 5 --logFile out.log

```

### FusionDetector pipeline

**INPUT=** (required)

clipped_reads.GEUS10xAttributes.umifound.bam file containing cellBC (BC tag) and UMI (U8 tag) for clipped reads only

**CSV=** (required)

.csv file listing, one per line, the cell barcodes in quantification (e.g. brain951.barcodes.tsv in /barcodes/ dir)

**OUTPUT=,O=** (required)

Output directory
 
**PREFIX**

Prefix for _matrix.txt, _metrics.txt and _molinfos.txt tab-delimited output text files

```

java -jar -Xmx22g sicelore.v2.jar FusionDetector I=clipped_reads.GEUS10xAttributes.umifound.bam O=/path/to/output/directory PREFIX=fus CSV=barcodes.csv

```

**output**

**PREFIX_matrix.txt**

Cell Barcode / fusion transcript matrix molecule counts

**PREFIX_metrics.txt**

Total number of molecules per fusion transcripts detected in dataset

**PREFIX**_molinfos.txt

per molecules (UMI/BC) fusion information



<a id="collapse-model"></a>

## 12) Detecting novel transcripts isoforms

**use CollaspeModel (sicelore.v2.jar)**

SAM records are grouped per cellBC and UMI by transcript isoform based on exon makeup considering only UMI supported by RNMIN (default=1) reads. 
Transcripts isoforms showing an exon structure supported by less than MINEVIDENCE (default=5) UMIs are filtered out. 
Transcripts isoforms detected as potential 3p and /or 5p degradated sequence of Gencode reference transcripts or longer potential novels transcripts isoforms are filtered out.
The set of novel isoforms are then validated using **CAGE** / **SHORT** / **POLYA** files if provided. Novel isoform is classify as valid if: 
(i) all exon-exon junctions either described in Gencode or confirmed in an external short read data given as STAR aligned bam file by at least **juncCo** reads; 
(ii) a 5’ start within **cageCo** nucleotides of a transcription start site identified by CAGE peaks given as .bed file; 
(iii) a 3’ end within **polyaCo** nucleotides of a polyadenylation site given as .bed file.


**INPUT=** (required)

Isoforms bam file (isobam.bam) produced by IsoformMatrix pipeline

**CSV=** (required)

.csv file listing, one per line, the cell barcodes in quantification (e.g. brain951.barcodes.tsv in /barcodes/ dir)

**REFFLAT=** (required)

Can be generated base on Gencode GTF file for genome build used for mapping with ***gtfToGenePred*** from [UCSC](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/)

**OUTDIR**

Output directory where output files are created

**PREFIX=** (required)

prefix used for output file names (default=CollapseModel)

**DELTA=** (required)

Number of extra or lacking bases allowed at exon-exon junctions (default = 2)

**MINEVIDENCE=** (required)

Minimum number of UMIs supporting the potential novel isoform to keep it in set (default = 5)

**RNMIN=** (required)

Minimum number of reads supporting the UMI to consider the UMI in novel isoform detection (default = 1)

**CELLTAG=**

Cell tag (default = BC)

**UMITAG=**

UMI sequence tag (default = U8)

**GENETAG=**

Gene name tag (default = GE)

**ISOFORMTAG=**

Gene name tag (default = IT)

**RNTAG=**

Read number tag (default=RN)

**SHORT=**

The short read SAM or BAM file (STAR aligned external short reads dataset)

**CAGE=**

CAGE peaks file (.bed) (Ressources/mm10.liftover.Fantom5.cage_peaks.bed)

**POLYA=**

POLYA sites file (.bed) (Ressources/mm10.gencode.vM24.polyAs.bed)

**cageCo=**

CAGE validation cutoff (default = 50 bases)

**polyaCo=**

PolyA validation cutoff (default = 50 bases)

**juncCo=**

Junction validation cutoff (default = 1 read)

```

java -jar -Xmx40g sicelore.v2.jar CollapseModel I=isobam.bam CSV=barcodes.csv REFFLAT=refFlat.txt O=. PREFIX=CollapseModel MINEVIDENCE=5 DELTA=2 RNMIN=1 SHORT=SRA.E18brain.star.bam CAGE=Ressources/mm10.liftover.Fantom5.cage_peaks.bed POLYA=Ressources/mm10.gencode.vM24.polyAs.bed

```

**output**

**PREFIX.d'DELTA'.rn'RNMIN'.e'MINEVIDENCE'.txt**

Gencode and novel isoforms detected classification .txt file

**PREFIX.d'DELTA'.rn'RNMIN'.e'MINEVIDENCE'.gff**

Gencode and all novel isoforms .gff file

**PREFIX.d'DELTA'.rn'RNMIN'.e'MINEVIDENCE'.final.gff**

Gencode and all validated novel isoforms .gff file

**PREFIX.d'DELTA'.rn'RNMIN'.e'MINEVIDENCE'.refflat.txt**

Gencode and all novel isoforms .refflat file

**PREFIX.d'DELTA'.rn'RNMIN'.e'MINEVIDENCE'.final.refflat.txt**

Gencode and all validated novel isoforms .refflat file (can be used for IsoformMatrix Sicelore pipeline quantification)

**PREFIX.d'DELTA'.rn'RNMIN'.e'MINEVIDENCE'.fas**

Representative sequence fasta file, poa/racon consensus sequence using top 20 best qualities UMIs (based on "de" minimap2 tag value)


