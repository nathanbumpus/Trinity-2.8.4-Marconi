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
<p>The downloaded reads are in SRA format.  To separate the reads into left and right reads use the following scripts for each group of reads. Your file paths will vary.  The important thing is that first --outdir specifies where to put the resulting left and right raw reads and the --split-files argument is the path to the downloaded SRA files.  The --defline-seq argunment will stay the same. First for the cardiac ganglion group.</p>

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

