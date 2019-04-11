<h1 align="center">Preparing High Quality Data for Assembly Construction<a id="top"></a></h1>

<p>Once we have our data in our raw_reads directory we are ready to evaluate the quality of that data, <a href="#tirm">trim contaminants</a> and prepare the reads for assembly construction.  In the project directory make another directory.  Call it something like "fastqc_reports".  This is where we will put our fastqc files generated from running fastqc.</p>

<h2 align="center">Running Fastqc on Marconi</h2>

<p>Use the following script to run fastqc on the shrimp data set.  The script can be easily modified to run with the lobster and yeast data sets.</p>

```
#!/bin/bash -l
#PBS -q haswell
#PBS -N fastqc
#PBS -l nodes=1:ppn=2
#PBS -l walltime=00:30:00
#PBS -o out.txt
#PBS -e err.txt

cd #PBS_O_WORKDIR

module load java
module load fastqc

fastqc -f fastq -o /home/nbumpus/shrimp/fastqc_reports/ /home/nbumpus/shrimp/raw_reads/SRR4341164_1.fastq
fastqc -f fastq -o /home/nbumpus/shrimp/fastqc_reports/ /home/nbumpus/shrimp/raw_reads/SRR4341164_2.fastq
fastqc -f fastq -o /home/nbumpus/shrimp/fastqc_reports/ /home/nbumpus/shrimp/raw_reads/SRR4341163_1.fastq
fastqc -f fastq -o /home/nbumpus/shrimp/fastqc_reports/ /home/nbumpus/shrimp/raw_reads/SRR4341163_2.fastq
fastqc -f fastq -o /home/nbumpus/shrimp/fastqc_reports/ /home/nbumpus/shrimp/raw_reads/SRR4341162_1.fastq
fastqc -f fastq -o /home/nbumpus/shrimp/fastqc_reports/ /home/nbumpus/shrimp/raw_reads/SRR4341162_2.fastq
fastqc -f fastq -o /home/nbumpus/shrimp/fastqc_reports/ /home/nbumpus/shrimp/raw_reads/SRR4341161_1.fastq
fastqc -f fastq -o /home/nbumpus/shrimp/fastqc_reports/ /home/nbumpus/shrimp/raw_reads/SRR4341161_2.fastq

```

<p>This script along with most of the others that follow it takes advantage of #PBS_O_WORKDIR variable which takes the value of the directory from which the script is called.</p>

<p>When the script finishes, sftp the html files in the fastqc_reports directory to your hard drive.  Once these files have been reviewed create another directory in the project directory.  Call it something like trimmed_reads.  We are now ready to trim the raw RNA seq reads of contamination such as adaptors to obtain high quality data from which to build an assembly. <a href="#top">back to top </a><a href="#contents">table of contents</a></p>

<h2 align="center">Running Trimmomatic on Marconi<a id="trim"></a></h2>

<p>To run Trimmomatic on Marconi use the following script for each set of the raw RNA seq reads</p>

```
#!/bin/bash -l
#PBS -q haswell
#PBS -N Trim1
#PBS -l nodes=1:ppn=2
#PBS -l walltime=1:00:00
#PBS -o out.txt
#PBS -e err.txt

cd #PBS_O_WORKDIR

module load java
module load trimmomatic

java -jar /opt/modules/universal/trimmomatic/0.36/trimmomatic-0.36.jar PE -phred33 \
/home/nbumpus/shrimp/raw_reads/SRR4341161_1.fastq \
/home/nbumpus/shrimp/raw_reads/SRR4341161_2.fastq \
/home/nbumpus/shrimp/trimmed_reads/SRR4341161_1.trim.paired.fastq \
/home/nbumpus/shrimp/trimmed_reads/SRR4341161_1.trim.unpaired.fastq \
/home/nbumpus/shrimp/trimmed_reads/SRR4341161_2.trim.paired.fastq \
/home/nbumpus/shrimp/trimmed_reads/SRR4341161_2.trim.unpaired.fastq \
ILLUMINACLIP:/opt/modules/universal/trimmomatic/0.36/adapters/TruSeq2-PE.fa:2:30:10 HEADCROP:12 LEADING:25 TRAILING:25 MINLEN:40

```

<p> Trimmomatic is a java program so we need to point our script at the jar file.  We then provide the read pairs that we want to trim and then direct their output to the trimmed_reads directory.  The ILLUNINACLIP line is where we specify how we want to trim the reads.  This step will perform the arguments in the order given so trim the adapters first.  The above script tells trimmomatic to remove the Illumina adapters, then trim the first 12 bases from reads, trim bases with a quality score less than 25 from the leading end and trailing ends of the reads and, lastly, drop reads that are less than 40 bases long.  Repeat this with lobster data as well. Then go back and run fastqc again. This time using the trimmed reads</p>

<p>Any over represented sequences that do not have any data available should be blasted against the NCBI database</p>

<p>If we now perform a head command one of the fastq files in the trimmed_reads directory we should see something like this.</p>

```
@D1317JN1:268:C5A93ACXX:4:1101:1719:2086_forward/1
AAGCATTGGAATCCGCCTCTTGGGAATCACGAGAATATGGACTGGCGCTTGTGGCGAAACATCCCGGAAAGCTAAGCACTGGTCATCTT
+SRR4341161.163 D1317JN1:268:C5A93ACXX:4:1101:1719:2086 length=101
CGF@GEHGIHF:C?CCFGBHAGDGD?F3?FHFH@FGGIII<FGG9@GE/AA@@?;939;=ACCCBBBBCBB9>:4:3ACACBCAC(:AA
@D1317JN1:268:C5A93ACXX:4:1101:1698:2209_forward/1
CAAACTAAACTACTTGTCATTCCTCAAAAAGAGATGTAATGTGAACCCCAAGCGTGGTCCTTTCCATCACCGATCTCCTGCCAAGATTT
+SRR4341161.164 D1317JN1:268:C5A93ACXX:4:1101:1698:2209 length=101
FIFIIEIEFAH<C@EG?FGFFFIFFIIEEFFCCDEFFCCF<BDCFFFF;CEFAECE<;7@@AAA;@BABABB/3=BBBB:A?A:??ABA
@D1317JN1:268:C5A93ACXX:4:1101:1916:2050_forward/1
ATTAACCTAAGTAAAAATATATTTGTTCATCATAATGGATAGCTACGTAATAGATTTAAATTTACTAGTTTTTTTTACTAAGAAATCTA

```

<p>This looks good except for one thing.  At the end of the first line of each entry there is a _forward/1.  The /1 tells trinity that the read is from the forward or left strand.  Trinity will not recognize the _forward so we need to remove it so that the only thing left is /1.  First make a new directory in the project folder called trinity_reads and then run the following script.</p>

```
#!/bin/bash -l
#PBS -q haswell
#PBS -N dirty-data
#PBS -l nodes=1:ppn=2
#PBS -l walltime=01:00:00
#PBS -o out.txt
#PBS -e err.txt

cd #PBS_O_WORKDIR

module load perl

cat /home/nbumpus/shrimp/trimmed_reads/SRR4341161_1.trim.paired.fastq | perl -lane 's/_forward//; print;' > /home/nbumpus/shrimp/trinity_reads/SRR4341161_1.trim.paired.adj.fastq
cat /home/nbumpus/shrimp/trimmed_reads/SRR4341162_1.trim.paired.fastq | perl -lane 's/_forward//; print;' > /home/nbumpus/shrimp/trinity_reads/SRR4341162_1.trim.paired.adj.fastq
cat /home/nbumpus/shrimp/trimmed_reads/SRR4341163_1.trim.paired.fastq | perl -lane 's/_forward//; print;' > /home/nbumpus/shrimp/trinity_reads/SRR4341163_1.trim.paired.adj.fastq
cat /home/nbumpus/shrimp/trimmed_reads/SRR4341164_1.trim.paired.fastq | perl -lane 's/_forward//; print;' > /home/nbumpus/shrimp/trinity_reads/SRR4341164_1.trim.paired.adj.fastq
cat /home/nbumpus/shrimp/trimmed_reads/SRR4341161_2.trim.paired.fastq | perl -lane 's/_reverse//; print;' > /home/nbumpus/shrimp/trinity_reads/SRR4341161_2.trim.paired.adj.fastq
cat /home/nbumpus/shrimp/trimmed_reads/SRR4341162_2.trim.paired.fastq | perl -lane 's/_reverse//; print;' > /home/nbumpus/shrimp/trinity_reads/SRR4341162_2.trim.paired.adj.fastq
cat /home/nbumpus/shrimp/trimmed_reads/SRR4341163_2.trim.paired.fastq | perl -lane 's/_reverse//; print;' > /home/nbumpus/shrimp/trinity_reads/SRR4341163_2.trim.paired.adj.fastq
cat /home/nbumpus/shrimp/trimmed_reads/SRR4341164_2.trim.paired.fastq | perl -lane 's/_reverse//; print;' > /home/nbumpus/shrimp/trinity_reads/SRR4341164_2.trim.paired.adj.fastq

```
<p>If we go to the trinity_reads directory and use the head command on one of the files the _forward and _reverse should be gone.  Notice we do not have to do this for the lobster data.  We are now ready to build an assembly using trinity. <a href="#top">back to top></a></p>


<h2 align="center">Table of Contents<a id="contents"></a></h2>
* [Home](README.md)
* [Obtaining Data](data.md)
* [Data Quality](dataqc.md)
* [Building an Assembly](assembly.md)
* [Abundance Estimation](abundance.md)
* [Assembly Quality](assemblyqc.md)
* [Differential Expression](DE.md)


