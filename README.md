## **ChIP-seq and ATAC-seq Processing and Peak calling Software (CAPPS)**

Python based software for the processing and peak calling of ATAC-seq and ChIP-seq Datasets. The program relies on the following key python modules:
```
tkinter
datetime
pandas
argparse
subprocess
turtle
```
CAPPS will attempts to install the necessary python packages on first launch. At first, it will check for `pip ` installation. If **pip** is installed, then it will attempt installation of the above mentioned programs, if they are not imported. 

## Table of Contents
1. [Genome Alignment](#genome-alignment)
2. [Requirement for Bioinformatic packages](#requirement-for-bioinformatic-packages)

&nbsp;

### **Genome Alignment
----------------------------

CAPPS currently does not offer `fastq` file alignment for SAM file generation or SAM to BAM file conversion. Use the code below for alignment of fastq files, constituting raw sequence data, with mouse or human genome and SAM to BAM conversion:

&nbsp;
&nbsp;

1. Install bowtie2 and samtools
> sudo apt-get install bowtie2
> sudo apt-get install samtools

2. Mouse genome preparation:
  1. Open terminal and make directory
  
  > mkdir genome_dir/
  
  2. Download genome. Provided is a code for downloading mouse genome:
  
  > wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
  
  3. Index mouse genome 
  
  > bowtie2-build mm39.fa mm39bowtie

3. Align the fastq file with the appropriate genome:

> bowtie2 --very-sensitive –dovetail --quiet --maxins 1000  --met-file ./file_metrics.txt -x ./path/to/genome/indexes -1 ./path/to/R1_reads -2 ./path/to/R2/reads  -S output.sam

4. Convert the resulting `.SAM` file to `.BAM` file:

>Samtools view -bS 58_DKO.sam > 58_DKO.bam


&nbsp;
&ensp;
&nbsp;
&ensp;

### Requirement for Bioinformatic packages
------------------------------------------
CAPPS automatically integrates shell and Java scripting into the python code, so the user does not have to be familiar either of the programming languages. The following packages need to be installed for CAPPS to work:

1. for BAM file filtering
```
Samtools 
Bamtools 
Bedtools 

```