<h1 align="center">How do I download data from NCBI or ENA<a id="top"></a></h1>

<p>If you do not already have data from an experiment NCBI and ENA have raw RNA seq data publicly available for download.</p>

<h2 align="center">Downloading data from the NCBI SRA database with the SRA Toolkit</h2>

<p>The easiest way to obtain data from the NCBI SRA database is to use the SRA Toolkit.  The Toolkit can be downloaded from <a href="https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software" target="_blank">SRA Toolkit Download</a> and then sftp to Marconi.</p>

<p>Once the packaged is unzipped we can use prefetch to obtain experimental data sets in SRA format.  I will use the shrimp data as an example</p>

```
#!/bin/bash -l
#PBS -q haswell
#PBS -N data
#PBS -l nodes=1:ppn=1
#PBS -l walltime=01:00:00
#PBS -o out.txt
#PBS -e err.txt

cd /home/nbumpus/sratoolkit.2.9.4-1-ubuntu64/bin/

module load java

./prefetch.2.9.4 -v SRR4341164
./prefetch.2.9.4 -v SRR4341163
./prefetch.2.9.4 -v SRR4341162
./prefetch.2.9.4 -v SRR4341161

```
<p>Since this is not a high memory job we will run this on haswell cluster.  The bio cluster should be reserved for jobs requiring high memory such as building the assembly.  This will output four SRA files into the ~/ncbi/public/sra/ directory. Next make a project directory called "shrimp" and inside that directory make a directory "raw_reads" using the mkdir command.  We can then use the following script to split the sra file into forward and reverse strands and output the resulting fastq files designated with a _1 for the left reads and a _2 for the right reads in the raw_reads directory.</p>
```
#!/bin/bash -l
#PBS -q haswell
#PBS -N data-prep
#PBS -l nodes=1:ppn=4
#PBS -l walltime=02:00:00
#PBS -o out.txt
#PBS -e err.txt

cd /home/nbumpus/sratoolkit.2.9.4-1-ubuntu64/bin/

module load java

./fastq-dump.2.9.4 --outdir /home/nbumpus/shrimp/raw_reads/ --defline-seq '@$sn[_$rn]/$ri' --split-files /home/nbumpus/ncbi/public/sra/SRR4341164.sra
./fastq-dump.2.9.4 --outdir /home/nbumpus/shrimp/raw_reads/ --defline-seq '@$sn[_$rn]/$ri' --split-files /home/nbumpus/ncbi/public/sra/SRR4341163.sra
./fastq-dump.2.9.4 --outdir /home/nbumpus/shrimp/raw_reads/ --defline-seq '@$sn[_$rn]/$ri' --split-files /home/nbumpus/ncbi/public/sra/SRR4341162.sra
./fastq-dump.2.9.4 --outdir /home/nbumpus/shrimp/raw_reads/ --defline-seq '@$sn[_$rn]/$ri' --split-files /home/nbumpus/ncbi/public/sra/SRR4341161.sra
```
<p>Repeat these steps for the SRA lobster data using SRR5428110, SRR5428111, SRR5428112, SRR5428113 or any other SRA data set you want to work with from NCBI.  If instead you want to download data from ENA you can do so directly from their website <a href="https://www.ebi.ac.uk/ena/browse/download" target="_blank">ENA website</a>.  The yeast data can be download via the command line using wget as described in the paper here <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3875132/" target="_blank">yeast data</a> Once we have our data we can move on to cleaning up the data before building the assembly. <a href="#top">back to top</a></p>

<h2 align="center">Table of Contents</h2>
* [Home](README.md)
* [Obtaining Data](data.md)
* [Data Quality](dataqc.md)
* [Building an Assembly](assembly.md)
* [Abundance Estimation](abundance.md)
* [Assembly Quality](assemblyqc.md)
* [Differential Expression](DE.md)
* [Annotations](annotations.md)
* [References](references.md)
