<h1 align="center"><a id="top"></a>Scripts Used for Lobster Assembly</h1>

<p>This page contains the scripts to use when working with the lobster data.  Scripts contained on this page may be of use in <a href="#download">downloading data</a> from NCBI, <a href="#separating">separating reads</a> in SRA format into left and right reads, performing <a href="#fastqc">fastqc </a>on the reads, <a href="#trim">trimming reads</a> with trimmomatic, <a href="#assembly">building the assembly</a>, <a href="#abundances">abundance estimation</a>, assesing the <a href="#quality">quality </a>of the assembly, <a href="#map">mapping </a>the reads back to the assembly, viewing the mappings in <a href="#igv">IGV </a> and determining which isoforms are <a href="#DE">differentially expressed</a>.

<h2 align="center"><a id="download"></a>Download and Prepare the Raw Reads from NCBI</h2>

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
<p><a id="separating"></a>The downloaded reads are in SRA format.  To separate the reads into left and right reads use the following scripts for each group of reads. Your file paths will vary.  The important thing is that the first --outdir argument specifies where to put the resulting left and right raw reads and the --split-files argument defines the path to the downloaded SRA files.  The --defline-seq argument will always be the same when downloading SRA files from NCBI. First for the cardiac ganglion group.</p> 
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
<a href="#top">back to top</a> <a href="#table">Table of Contents</a>

<h2 align="center"><a id="fastqc"></a>FastQC</h2>

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

<p>Note that we can split the sets of reads into as many scripts as we want or run them all in one script.  I found it helpful to keep everything organized by condition.  After reviewing the fastqc reports we can use trimmomatic to clean the reads of contaminants.</p> <a href="#top">back to top</a> <a href="#table">Table of Contents</a>

<h2 align="center"><a id="trim"></a>Trimmomatic</h2>

<p>Below is one of the trimmoatic scripts used to clean one of the sets of reads.  We can adjust this script for each set of reads and run all of the scripts in parallel.  Here is the script used for one set of the Premotor Neuron reads</p>

```
#!/bin/bash -l
#PBS -q haswell
#PBS -N trim72
#PBS -l nodes=1:ppn=1
#PBS -l walltime=1:00:00
#PBS -o out.txt
#PBS -e err.txt

cd #PBS_O_WORKDIR

module load java
module load trimmomatic

java -jar /opt/modules/universal/trimmomatic/0.36/trimmomatic-0.36.jar PE -phred33 \
/home/nbumpus/lobster/raw_reads/pmn/SRR7156172_1.fastq \
/home/nbumpus/lobster/raw_reads/pmn/SRR7156172_2.fastq \
/home/nbumpus/lobster/trimmed_reads/pmn/SRR7156172_1.trim.paired.fastq \
/home/nbumpus/lobster/trimmed_reads/pmn/SRR7156172_1.trim.unpaired.fastq \
/home/nbumpus/lobster/trimmed_reads/pmn/SRR7156172_2.trim.paired.fastq \
/home/nbumpus/lobster/trimmed_reads/pmn/SRR7156172_2.trim.unpaired.fastq \
ILLUMINACLIP:/opt/modules/universal/trimmomatic/0.36/adapters/TruSeq2-PE.fa:2:30:10 HEADCROP:9 LEADING:20 TRAILING:20 MINLEN:50
```
<p>Next perform fastqc again only this time using the resulting trimmed reads as inputs.  Once we have determined that reads have been sufficiently trimmed we can use these reads to build the assembly.</p> <a href="#top">back to top</a> <a href="#table">Table of Contents</a>

<h2 align="center"><a id="assembly"></a>Lobster Assembly</h2>

<p>Instead of catonating all of the left reads together and all of the right reads together we can create a tab deliminited assembly_samples_file.txt file describing each pair of left and right reads like so.</p>

```
PMN     PMN2    /home/nbumpus/lobster/trinity_reads/left/SRR7156183_1.trim.paired.fastq /home/nbumpus/lobster/trinity_reads/right/SRR7156183_2.trim.paired.fastq
PMN     PMN1    /home/nbumpus/lobster/trinity_reads/left/SRR7156182_1.trim.paired.fastq /home/nbumpus/lobster/trinity_reads/right/SRR7156182_2.trim.paired.fastq
PMN     PMN4    /home/nbumpus/lobster/trinity_reads/left/SRR7156173_1.trim.paired.fastq /home/nbumpus/lobster/trinity_reads/right/SRR7156173_2.trim.paired.fastq
PMN     PMN3    /home/nbumpus/lobster/trinity_reads/left/SRR7156172_1.trim.paired.fastq /home/nbumpus/lobster/trinity_reads/right/SRR7156172_2.trim.paired.fastq
MN	MN2     /home/nbumpus/lobster/trinity_reads/left/SRR7156181_1.trim.paired.fastq /home/nbumpus/lobster/trinity_reads/right/SRR7156181_2.trim.paired.fastq
MN	MN1     /home/nbumpus/lobster/trinity_reads/left/SRR7156180_1.trim.paired.fastq /home/nbumpus/lobster/trinity_reads/right/SRR7156180_2.trim.paired.fastq
MN	MN4     /home/nbumpus/lobster/trinity_reads/left/SRR7156179_1.trim.paired.fastq /home/nbumpus/lobster/trinity_reads/right/SRR7156179_2.trim.paired.fastq
MN	MN3     /home/nbumpus/lobster/trinity_reads/left/SRR7156178_1.trim.paired.fastq /home/nbumpus/lobster/trinity_reads/right/SRR7156178_2.trim.paired.fastq
CG	CG2     /home/nbumpus/lobster/trinity_reads/left/SRR7156177_1.trim.paired.fastq /home/nbumpus/lobster/trinity_reads/right/SRR7156177_2.trim.paired.fastq
CG	CG1     /home/nbumpus/lobster/trinity_reads/left/SRR7156176_1.trim.paired.fastq /home/nbumpus/lobster/trinity_reads/right/SRR7156176_2.trim.paired.fastq
CG	CG4     /home/nbumpus/lobster/trinity_reads/left/SRR7156175_1.trim.paired.fastq /home/nbumpus/lobster/trinity_reads/right/SRR7156175_2.trim.paired.fastq
CG	CG3     /home/nbumpus/lobster/trinity_reads/left/SRR7156174_1.trim.paired.fastq /home/nbumpus/lobster/trinity_reads/right/SRR7156174_2.trim.paired.fastq
```
<p>Each line in the file has the group name followed by the replicate name followed by the left reads followed by the right reads.  This file can then be fed into the trinity script to construct the assembly as shown below.</p>

```
#!/bin/bash -l
#PBS -q bio
#PBS -N stuckonchrysallis
#PBS -l nodes=fast1:ppn=12
#PBS -l walltime=240:00:00
#PBS -o out4.txt
#PBS -e err4.txt
cd #PBS_O_WORKDIR

module load python
module load trinity/2.8.4

Trinity --seqType fq --SS_lib_type RF --normalize_max_read_cov 50 --min_contig_length 324 \
--samples_file /home/nbumpus/lobster/assembly_samples_file.txt \
--output /home/nbumpus/lobster/trinity_out_dir3 \
--CPU 12 \
--max_memory 250G
```
<p>This script needs to be run on one of the high memory bio nodes.  By specifying that this job should run on the fast1 node we take all of the available CPU's and ensure that we have access to all of the memory on that node.  This also leaves the larger bio node free for others to use.  This job should take between three and four days to run.</p> <a href="#top">back to top</a> <a href="#table">Table of Contents</a>

<h2 align="center"><a id="abundances"></a>Abundance Estimates</h2>

<p>The abundance estimations can be calculated by running the following series of scripts after setting up an abundances directory.  We can use salmon for this task but RSEM will also work.  The first step when using salmon is to prepare a reference.  We can accomplish this by using the following script.</p>

```
#!/bin/bash -l
#PBS -q haswell
#PBS -N referenceprep
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00
#PBS -o out.txt
#PBS -e err.txt

cd #PBS_O_WORKDIR

module load trinity/2.8.4

$TRINITY_HOME/util/align_and_estimate_abundance.pl \
--transcripts /home/nbumpus/lobster/trinity_out_dir3/Trinity.fasta \
--est_method salmon \
--gene_trans_map /home/nbumpus/lobster/trinity_out_dir3/Trinity.fasta.gene_trans_map \
--SS_lib_type RF \
--prep_reference
```
<p>Create a tab delimiited text file describing the sample to replicate relationships along with file paths to the reads.  The file should look similar to this</p>

```
SRR7156183	/home/nbumpus/lobster/abundances/PMN2   /home/nbumpus/lobster/trimmed_reads/pmn/SRR7156183_1.trim.paired.fastq  /home/nbumpus/lobster/trimmed_reads/pmn/SRR7156183_2.trim.paired.fastq
SRR7156182	/home/nbumpus/lobster/abundances/PMN1   /home/nbumpus/lobster/trimmed_reads/pmn/SRR7156182_1.trim.paired.fastq  /home/nbumpus/lobster/trimmed_reads/pmn/SRR7156182_2.trim.paired.fastq
SRR7156173	/home/nbumpus/lobster/abundances/PMN4   /home/nbumpus/lobster/trimmed_reads/pmn/SRR7156173_1.trim.paired.fastq  /home/nbumpus/lobster/trimmed_reads/pmn/SRR7156173_2.trim.paired.fastq
SRR7156172	/home/nbumpus/lobster/abundances/PMN3   /home/nbumpus/lobster/trimmed_reads/pmn/SRR7156172_1.trim.paired.fastq  /home/nbumpus/lobster/trimmed_reads/pmn/SRR7156172_2.trim.paired.fastq
SRR7156181	/home/nbumpus/lobster/abundances/MN2   /home/nbumpus/lobster/trimmed_reads/mn/SRR7156181_1.trim.paired.fastq  /home/nbumpus/lobster/trimmed_reads/mn/SRR7156181_2.trim.paired.fastq
SRR7156180	/home/nbumpus/lobster/abundances/MN1   /home/nbumpus/lobster/trimmed_reads/mn/SRR7156180_1.trim.paired.fastq  /home/nbumpus/lobster/trimmed_reads/mn/SRR7156180_2.trim.paired.fastq
SRR7156179	/home/nbumpus/lobster/abundances/MN4   /home/nbumpus/lobster/trimmed_reads/mn/SRR7156179_1.trim.paired.fastq  /home/nbumpus/lobster/trimmed_reads/mn/SRR7156179_2.trim.paired.fastq
SRR7156178	/home/nbumpus/lobster/abundances/MN3   /home/nbumpus/lobster/trimmed_reads/mn/SRR7156178_1.trim.paired.fastq  /home/nbumpus/lobster/trimmed_reads/mn/SRR7156178_2.trim.paired.fastq
SRR7156177	/home/nbumpus/lobster/abundances/CG2   /home/nbumpus/lobster/trimmed_reads/cg/SRR7156177_1.trim.paired.fastq  /home/nbumpus/lobster/trimmed_reads/cg/SRR7156177_2.trim.paired.fastq
SRR7156176	/home/nbumpus/lobster/abundances/CG1   /home/nbumpus/lobster/trimmed_reads/cg/SRR7156176_1.trim.paired.fastq  /home/nbumpus/lobster/trimmed_reads/cg/SRR7156176_2.trim.paired.fastq
SRR7156175	/home/nbumpus/lobster/abundances/CG4   /home/nbumpus/lobster/trimmed_reads/cg/SRR7156175_1.trim.paired.fastq  /home/nbumpus/lobster/trimmed_reads/cg/SRR7156175_2.trim.paired.fastq
SRR7156174	/home/nbumpus/lobster/abundances/CG3   /home/nbumpus/lobster/trimmed_reads/cg/SRR7156174_1.trim.paired.fastq  /home/nbumpus/lobster/trimmed_reads/cg/SRR7156174_2.trim.paired.fastq
```
<p>Feed this file into the the following script as the --samples_file argument like so.</p>

```
#!/bin/bash -l
#PBS -q haswell
#PBS -N abundances
#PBS -l nodes=1:ppn=12
#PBS -l walltime=10:00:00
#PBS -o abundout.txt
#PBS -e abunderr.txt

cd #PBS_O_WORKDIR

module load trinity/2.8.4

$TRINITY_HOME/util/align_and_estimate_abundance.pl \
--transcripts /home/nbumpus/lobster/trinity_out_dir3/Trinity.fasta \
--seqType fq \
--est_method salmon \
--samples_file /home/nbumpus/lobster/samples_file.txt \
--salmon_add_opts --validateMappings \
--gene_trans_map /home/nbumpus/lobster/trinity_out_dir3/Trinity.fasta.gene_trans_map \
--output_dir /home/nbumpus/lobster/abundances/ \
--thread_count 12
```
<p>This will generate the quant.sf files to be used to generate the matrices.  First make a matrix directory inside the abundances directory and then make a text file listing the file paths to each of the quant.sf files.  The file should look similar to this</p>

```
/home/nbumpus/lobster/abundances/PMN2/quant.sf
/home/nbumpus/lobster/abundances/PMN1/quant.sf
/home/nbumpus/lobster/abundances/PMN4/quant.sf
/home/nbumpus/lobster/abundances/PMN3/quant.sf
/home/nbumpus/lobster/abundances/MN2/quant.sf
/home/nbumpus/lobster/abundances/MN1/quant.sf
/home/nbumpus/lobster/abundances/MN4/quant.sf
/home/nbumpus/lobster/abundances/MN3/quant.sf
/home/nbumpus/lobster/abundances/CG2/quant.sf
/home/nbumpus/lobster/abundances/CG1/quant.sf
/home/nbumpus/lobster/abundances/CG4/quant.sf
/home/nbumpus/lobster/abundances/CG3/quant.sf
```
<p>Feed this file into the --quant_files argument in the following script to generate the matrices.</p>

```
#!/bin/bash -l
#PBS -q haswell
#PBS -N matrix
#PBS -l nodes=1:ppn=1
#PBS -l walltime=05:00:00
#PBS -o matrixout.txt
#PBS -e matrixerr.txt

cd #PBS_O_WORKDIR

module load trinity/2.8.4
module load R/3.5.2

$TRINITY_HOME/util/abundance_estimates_to_matrix.pl \
--est_method salmon \
--gene_trans_map /home/nbumpus/lobster/trinity_out_dir3/Trinity.fasta.gene_trans_map2 \
--name_sample_by_basedir \
--out_prefix /home/nbumpus/lobster/abundances/matrix/lobster \
--quant_files /home/nbumpus/lobster/quant_files.txt
```
<p>You may notice that Trinity.fasta.gene_trans_map file has been modified.  This is due to an artifact that is sometimes generated when Trinity works with Salmon and the auto removal fails.  We can grep the Trinity id of the offending abundance and remove the entry from the Trinity.fasta.gene_trans_map file</p>

<p>Once this script finishes we should have counts matrices and EXPR matrices for both the isoforms and the genes.  We can now evaluate the quality of our assembly</p> <a href="#top">back to top</a> <a href="#table">Table of Contents</a>

<h2 align="center"><a id="quality"></a>Assembly Quality</h2>

<p>The following script can be used to generate the basic assembly stats including the N50 statistic.</p>

```
#!/bin/bash -l
#PBS -q haswell
#PBS -N AssemblyStats
#PBS -l nodes=1:ppn=1
#PBS -l walltime=02:00:00
#PBS -o /home/nbumpus/lobster/assembly_quality/assembly_stats.txt
#PBS -e err.txt

cd #PBS_O_WORKDIR

module load trinity/2.8.4

$TRINITY_HOME/util/TrinityStats.pl \
/home/nbumpus/lobster/trinity_out_dir3/Trinity.fasta
```
<p><a id="map"></a>To produce alignmentment statistics using bowtie2 first use the following script to build an index</p>

```
#!/bin/bash -l
#PBS -q haswell
#PBS -N Bowtie2-build
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -o out.txt
#PBS -e err.txt

cd #PBS_O_WORKDIR

module load Bowtie2/2.3.4.3

bowtie2-build /home/nbumpus/lobster/trinity_out_dir3/Trinity.fasta \
/home/nbumpus/lobster/assembly_quality/lobster_ref
```
<p>Next align the reads to the assembly by modifying the following script for each individual set of reads.  We can run each of the scripts in parallel.  Here is an example using the PMN1 replicate.  This will generate the alignment statistics shown <a href="https://nathanbumpus.github.io/Trinity-2.8.4-Marconi/results.html#alignment" target="_blank">here.</a></p>

```
#!/bin/bash -l
#PBS -q haswell
#PBS -N bowtie2-PMN1-dove
#PBS -l nodes=1:ppn=16
#PBS -l walltime=02:00:00
#PBS -o doveout.txt
#PBS -e doveerr.txt

cd #PBS_O_WORKDIR

module load perl
module load Bowtie2/2.3.4.3
module load samtools

bowtie2 -p 16 -q \
--no-unal --dovetail -k 20 -x /home/nbumpus/lobster/assembly_quality/lobster_ref \
-1 /home/nbumpus/lobster/trinity_reads/left/SRR7156182_1.trim.paired.fastq \
-2 /home/nbumpus/lobster/trinity_reads/right/SRR7156182_2.trim.paired.fastq \
2>/home/nbumpus/lobster/assembly_quality/PMN1_dove_align_stats.txt| \
samtools view -@16 -Sb \
-o /home/nbumpus/lobster/assembly_quality/PMN1_dove.bowtie2.bam
```
<p><a id="igv"></a>If we want to view the alignments in IGV we can use the following script to create the necessary coordinate sorted bam files and indices.</p>

```
#!/bin/bash -l
#PBS -q haswell
#PBS -N igvprep
#PBS -l nodes=1:ppn=16
#PBS -l walltime=24:00:00
#PBS -o out.txt
#PBS -e err.txt

cd #PBS_O_WORKDIR

module load samtools

samtools sort -@16 /home/nbumpus/lobster/assembly_quality/PMN1_dove.bowtie2.bam \
-o /home/nbumpus/lobster/assembly_quality/PMN1_dove.bowtie2.coordsort.bam
samtools sort -@16 /home/nbumpus/lobster/assembly_quality/PMN2_dove.bowtie2.bam \
-o /home/nbumpus/lobster/assembly_quality/PMN2_dove.bowtie2.coordsort.bam
samtools sort -@16 /home/nbumpus/lobster/assembly_quality/PMN3_dove.bowtie2.bam \
-o /home/nbumpus/lobster/assembly_quality/PMN3_dove.bowtie2.coordsort.bam
samtools sort -@16 /home/nbumpus/lobster/assembly_quality/PMN4_dove.bowtie2.bam \
-o /home/nbumpus/lobster/assembly_quality/PMN4_dove.bowtie2.coordsort.bam
samtools sort -@16 /home/nbumpus/lobster/assembly_quality/MN1_dove.bowtie2.bam \
-o /home/nbumpus/lobster/assembly_quality/MN1_dove.bowtie2.coordsort.bam
samtools sort -@16 /home/nbumpus/lobster/assembly_quality/MN2_dove.bowtie2.bam \
-o /home/nbumpus/lobster/assembly_quality/MN2_dove.bowtie2.coordsort.bam
samtools sort -@16 /home/nbumpus/lobster/assembly_quality/MN3_dove.bowtie2.bam \
-o /home/nbumpus/lobster/assembly_quality/MN3_dove.bowtie2.coordsort.bam
samtools sort -@16 /home/nbumpus/lobster/assembly_quality/MN4_dove.bowtie2.bam \
-o /home/nbumpus/lobster/assembly_quality/MN4_dove.bowtie2.coordsort.bam
samtools sort -@16 /home/nbumpus/lobster/assembly_quality/CG1_dove.bowtie2.bam \
-o /home/nbumpus/lobster/assembly_quality/CG1_dove.bowtie2.coordsort.bam
samtools sort -@16 /home/nbumpus/lobster/assembly_quality/CG2_dove.bowtie2.bam \
-o /home/nbumpus/lobster/assembly_quality/CG2_dove.bowtie2.coordsort.bam
samtools sort -@16 /home/nbumpus/lobster/assembly_quality/CG3_dove.bowtie2.bam \
-o /home/nbumpus/lobster/assembly_quality/CG3_dove.bowtie2.coordsort.bam
samtools sort -@16 /home/nbumpus/lobster/assembly_quality/CG4_dove.bowtie2.bam \
-o /home/nbumpus/lobster/assembly_quality/CG4_dove.bowtie2.coordsort.bam
samtools index /home/nbumpus/lobster/assembly_quality/PMN1.bowtie2.coordsort.bam
samtools index /home/nbumpus/lobster/assembly_quality/PMN2.bowtie2.coordsort.bam
samtools index /home/nbumpus/lobster/assembly_quality/PMN3.bowtie2.coordsort.bam
samtools index /home/nbumpus/lobster/assembly_quality/PMN4.bowtie2.coordsort.bam
samtools index /home/nbumpus/lobster/assembly_quality/MN1.bowtie2.coordsort.bam
samtools index /home/nbumpus/lobster/assembly_quality/MN2.bowtie2.coordsort.bam
samtools index /home/nbumpus/lobster/assembly_quality/MN3.bowtie2.coordsort.bam
samtools index /home/nbumpus/lobster/assembly_quality/MN4.bowtie2.coordsort.bam
samtools index /home/nbumpus/lobster/assembly_quality/CG1.bowtie2.coordsort.bam
samtools index /home/nbumpus/lobster/assembly_quality/CG2.bowtie2.coordsort.bam
samtools index /home/nbumpus/lobster/assembly_quality/CG3.bowtie2.coordsort.bam
samtools index /home/nbumpus/lobster/assembly_quality/CG4.bowtie2.coordsort.bam
samtools faidx /home/nbumpus/lobster/trinity_out_dir/Trinity.fasta
```
<p>Then use the coordinate sorted bam files, the index files for the bam files, the Trinity.fasta file and the Trinity.fasta files index as inputs in IGV.</p> <a href="#top">back to top</a> <a href="#table">Table of Contents</a>

<h2 align="center"><a id="DE"></a>Differential Expression</h2>

<p>To determine which isoforms are differentially expressed we can use edgeR with the following two scripts.</p>

```
#!/bin/bash -l
#PBS -q haswell
#PBS -N lob-DE-iso
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00
#PBS -o out.txt
#PBS -e err.txt

cd #PBS_O_WORKDIR

module load trinity/2.8.4
module load R/3.5.2

$TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl \
--matrix /home/nbumpus/lobster/abundances/matrix/lobster.isoform.counts.matrix \
--method edgeR \
--samples_file /home/nbumpus/lobster/samples_described.txt \
--output /home/nbumpus/lobster/DifferentialExpression/isoform/
```
<p>Second, from with the directory containing the differential expression results (the same directory as the --output from the previous script) run the following script to obtain subsets of up regulated isoforms and heatmaps.</p>

```
#!/bin/bash -l
#PBS -q haswell
#PBS -N lob-DE-iso-anlys
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:10:00
#PBS -o out.txt
#PBS -e err.txt

cd /home/nbumpus/lobster/DifferentialExpression/isoform/

module load trinity/2.8.4
module load R/3.5.2

$TRINITY_HOME/Analysis/DifferentialExpression/analyze_diff_expr.pl \
--matrix /home/nbumpus/lobster/abundances/matrix/lobster.isoform.TMM.EXPR.matrix \
--samples /home/nbumpus/lobster/samples_described.txt \
--max_genes_clust 18000
```
<h2 align="center"><a id="table"></a>Table of Contents</h2>
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

<a href="#top">back to top</a>
