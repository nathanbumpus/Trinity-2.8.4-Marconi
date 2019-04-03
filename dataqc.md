<h1 align="center">Preparing High Quality Data for Assembly Construction</h1>

<p>Once we have our data in our raw_reads directory we are ready to evaluate the quality of that data, trim contaminants and prepare the reads for assembly construction.  In the project directory make another directory.  Call it something like "fastqc_reports".  This is where we will put our fastqc files generated from running fastqc.</p>

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

<p>When the script finishes, sftp the html files in the fastqc_reports directory to your hard drive.  Once these files have been reviewed create another directory in the project directory.  Call it something like trimmed_reads.  We are now ready to trim the raw RNA seq reads of contamination such as adaptors to obtain high quality data from which to build an assembly.</p>

<h2 align="center">Running Trimmomatic on Marconi</h2>
