<h1 align="center">Scripts Used for Lobster Assembly</h1>

<h2 align="center">Dowload and Prepare the Raw Reads from NCBI</h2>

<p>First run the following script to download each set of raw reads using the SRA toolkit.</p>

```
#!/bin/bash -l
#PBS -q haswell
#PBS -N data
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00
#PBS -o out.txt
#PBS -e err.txt

cd /home/nbumpus/sratoolkit.2.9.4-1-ubuntu64/bin/

module load java

./prefetch.2.9.4 -v SRR7156183
./prefetch.2.9.4 -v SRR7156182
./prefetch.2.9.4 -v SRR7156181
./prefetch.2.9.4 -v SRR7156180
./prefetch.2.9.4 -v SRR7156179
./prefetch.2.9.4 -v SRR7156178
./prefetch.2.9.4 -v SRR7156177
./prefetch.2.9.4 -v SRR7156176
./prefetch.2.9.4 -v SRR7156175
./prefetch.2.9.4 -v SRR7156174
./prefetch.2.9.4 -v SRR7156173
./prefetch.2.9.4 -v SRR7156172
```
<p>The downloaded reads are in SRA format.  To separate the reads into left and right reads use the following scripts for each group of reads. Your file paths will vary.  The important thing is that the first --outdir argument specifies where to put the resulting left and right raw reads and the --split-files argument defines the path to the downloaded SRA files.  The --defline-seq argument will always be the same when downloading SRA files from NCBI. First for the cardiac ganglion group.</p>

```
#!/bin/bash -l
#PBS -q haswell
#PBS -N data-prep-cg
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00
#PBS -o out.txt
#PBS -e err.txt

cd /home/nbumpus/sratoolkit.2.9.4-1-ubuntu64/bin/

module load java

./fastq-dump.2.9.4 --outdir /home/nbumpus/lobster/raw_reads/cg/ --defline-seq '@$sn[_$rn]/$ri' --split-files /home/nbumpus/ncbi/public/sra/SRR7156177.sra
./fastq-dump.2.9.4 --outdir /home/nbumpus/lobster/raw_reads/cg/ --defline-seq '@$sn[_$rn]/$ri' --split-files /home/nbumpus/ncbi/public/sra/SRR7156176.sra
./fastq-dump.2.9.4 --outdir /home/nbumpus/lobster/raw_reads/cg/ --defline-seq '@$sn[_$rn]/$ri' --split-files /home/nbumpus/ncbi/public/sra/SRR7156175.sra
./fastq-dump.2.9.4 --outdir /home/nbumpus/lobster/raw_reads/cg/ --defline-seq '@$sn[_$rn]/$ri' --split-files /home/nbumpus/ncbi/public/sra/SRR7156174.sra
```
<p>Next for the Premotor Neuron group</p>

```
#!/bin/bash -l
#PBS -q haswell
#PBS -N data-prep-pmn
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00
#PBS -o out.txt
#PBS -e err.txt

cd /home/nbumpus/sratoolkit.2.9.4-1-ubuntu64/bin/

module load java

./fastq-dump.2.9.4 --outdir /home/nbumpus/lobster/raw_reads/pmn/ --defline-seq '@$sn[_$rn]/$ri' --split-files /home/nbumpus/ncbi/public/sra/SRR7156183.sra
./fastq-dump.2.9.4 --outdir /home/nbumpus/lobster/raw_reads/pmn/ --defline-seq '@$sn[_$rn]/$ri' --split-files /home/nbumpus/ncbi/public/sra/SRR7156182.sra
./fastq-dump.2.9.4 --outdir /home/nbumpus/lobster/raw_reads/pmn/ --defline-seq '@$sn[_$rn]/$ri' --split-files /home/nbumpus/ncbi/public/sra/SRR7156173.sra
./fastq-dump.2.9.4 --outdir /home/nbumpus/lobster/raw_reads/pmn/ --defline-seq '@$sn[_$rn]/$ri' --split-files /home/nbumpus/ncbi/public/sra/SRR7156172.sra
```

<p>And lastly for the the Motor Neuron group</p>

```
#!/bin/bash -l
#PBS -q haswell
#PBS -N data-prep-mn
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00
#PBS -o out.txt
#PBS -e err.txt

cd /home/nbumpus/sratoolkit.2.9.4-1-ubuntu64/bin/

module load java

./fastq-dump.2.9.4 --outdir /home/nbumpus/lobster/raw_reads/mn/ --defline-seq '@$sn[_$rn]/$ri' --split-files /home/nbumpus/ncbi/public/sra/SRR7156181.sra
./fastq-dump.2.9.4 --outdir /home/nbumpus/lobster/raw_reads/mn/ --defline-seq '@$sn[_$rn]/$ri' --split-files /home/nbumpus/ncbi/public/sra/SRR7156180.sra
./fastq-dump.2.9.4 --outdir /home/nbumpus/lobster/raw_reads/mn/ --defline-seq '@$sn[_$rn]/$ri' --split-files /home/nbumpus/ncbi/public/sra/SRR7156179.sra
./fastq-dump.2.9.4 --outdir /home/nbumpus/lobster/raw_reads/mn/ --defline-seq '@$sn[_$rn]/$ri' --split-files /home/nbumpus/ncbi/public/sra/SRR7156178.sra
```

<h2 align="center">FastQC</h2>

<p>We should now have sets of left and right reads for each sample.  We can now use the following scripts in parallel to run fastqc for each set of reads.</p>

```
#!/bin/bash -l
#PBS -q haswell
#PBS -N fastqc-cg
#PBS -l nodes=1:ppn=1
#PBS -l walltime=02:00:00
#PBS -o out.txt
#PBS -e err.txt

cd #PBS_O_WORKDIR

module load java
module load fastqc

fastqc -f fastq -o /home/nbumpus/lobster/fastqc_reports/cg/ /home/nbumpus/lobster/raw_reads/cg/SRR7156177_1.fastq
fastqc -f fastq -o /home/nbumpus/lobster/fastqc_reports/cg/ /home/nbumpus/lobster/raw_reads/cg/SRR7156177_2.fastq
fastqc -f fastq -o /home/nbumpus/lobster/fastqc_reports/cg/ /home/nbumpus/lobster/raw_reads/cg/SRR7156176_1.fastq
fastqc -f fastq -o /home/nbumpus/lobster/fastqc_reports/cg/ /home/nbumpus/lobster/raw_reads/cg/SRR7156176_2.fastq
fastqc -f fastq -o /home/nbumpus/lobster/fastqc_reports/cg/ /home/nbumpus/lobster/raw_reads/cg/SRR7156175_1.fastq
fastqc -f fastq -o /home/nbumpus/lobster/fastqc_reports/cg/ /home/nbumpus/lobster/raw_reads/cg/SRR7156175_2.fastq
fastqc -f fastq -o /home/nbumpus/lobster/fastqc_reports/cg/ /home/nbumpus/lobster/raw_reads/cg/SRR7156174_1.fastq
fastqc -f fastq -o /home/nbumpus/lobster/fastqc_reports/cg/ /home/nbumpus/lobster/raw_reads/cg/SRR7156174_2.fastq
```
```
#!/bin/bash -l
#PBS -q haswell
#PBS -N fastqc-pmn
#PBS -l nodes=1:ppn=1
#PBS -l walltime=02:00:00
#PBS -o out.txt
#PBS -e err.txt

cd #PBS_O_WORKDIR

module load java
module load fastqc

fastqc -f fastq -o /home/nbumpus/lobster/fastqc_reports/pmn/ /home/nbumpus/lobster/raw_reads/pmn/SRR7156183_1.fastq
fastqc -f fastq -o /home/nbumpus/lobster/fastqc_reports/pmn/ /home/nbumpus/lobster/raw_reads/pmn/SRR7156183_2.fastq
fastqc -f fastq -o /home/nbumpus/lobster/fastqc_reports/pmn/ /home/nbumpus/lobster/raw_reads/pmn/SRR7156182_1.fastq
fastqc -f fastq -o /home/nbumpus/lobster/fastqc_reports/pmn/ /home/nbumpus/lobster/raw_reads/pmn/SRR7156182_2.fastq
fastqc -f fastq -o /home/nbumpus/lobster/fastqc_reports/pmn/ /home/nbumpus/lobster/raw_reads/pmn/SRR7156173_1.fastq
fastqc -f fastq -o /home/nbumpus/lobster/fastqc_reports/pmn/ /home/nbumpus/lobster/raw_reads/pmn/SRR7156173_2.fastq
fastqc -f fastq -o /home/nbumpus/lobster/fastqc_reports/pmn/ /home/nbumpus/lobster/raw_reads/pmn/SRR7156172_1.fastq
fastqc -f fastq -o /home/nbumpus/lobster/fastqc_reports/pmn/ /home/nbumpus/lobster/raw_reads/pmn/SRR7156172_2.fastq
```
```
#!/bin/bash -l
#PBS -q haswell
#PBS -N fastqc-mn
#PBS -l nodes=1:ppn=1
#PBS -l walltime=02:00:00
#PBS -o out.txt
#PBS -e err.txt

cd #PBS_O_WORKDIR

module load java
module load fastqc

fastqc -f fastq -o /home/nbumpus/lobster/fastqc_reports/mn/ /home/nbumpus/lobster/raw_reads/mn/SRR7156181_1.fastq
fastqc -f fastq -o /home/nbumpus/lobster/fastqc_reports/mn/ /home/nbumpus/lobster/raw_reads/mn/SRR7156181_2.fastq
fastqc -f fastq -o /home/nbumpus/lobster/fastqc_reports/mn/ /home/nbumpus/lobster/raw_reads/mn/SRR7156180_1.fastq
fastqc -f fastq -o /home/nbumpus/lobster/fastqc_reports/mn/ /home/nbumpus/lobster/raw_reads/mn/SRR7156180_2.fastq
fastqc -f fastq -o /home/nbumpus/lobster/fastqc_reports/mn/ /home/nbumpus/lobster/raw_reads/mn/SRR7156179_1.fastq
fastqc -f fastq -o /home/nbumpus/lobster/fastqc_reports/mn/ /home/nbumpus/lobster/raw_reads/mn/SRR7156179_2.fastq
fastqc -f fastq -o /home/nbumpus/lobster/fastqc_reports/mn/ /home/nbumpus/lobster/raw_reads/mn/SRR7156178_1.fastq
fastqc -f fastq -o /home/nbumpus/lobster/fastqc_reports/mn/ /home/nbumpus/lobster/raw_reads/mn/SRR7156178_2.fastq
```
<p>Note that we can split the sets of reads into as many scripts as we want or run them all in one script.  I found it helpful to keep everything organized by condition.  After reviewing the fastqc reports we can use trimmomatic to clean the reads of contaminants.</p>

<h2 align="center">Trimmomatic</h2>



<h1>Work Being Done Here</h1>


* [Home](README.md)
* [Obtaining Data](data.md)
* [Data Quality](dataqc.md)
* [Building an Assembly](assembly.md)
* [Abundance Estimation](abundance.md)
* [Assembly Quality](assemblyqc.md)
* [Differential Expression](DE.md)
* [Annotations](annotations.md)
* [Results](results.md)
* [References](references.md)

