<h1 align="center">Estimating Transcript Abundance with Trinity</h1>

<p>When estinmating transcript abundance with Trinity on Marconi we need to choose whether to use an alignment based or a non-alignment based approach.  Trinity has a script align_and_estimate_abundance.pl which has direct support for the RSEM and Salmon modules on Marconi.  In this project I used Salmon as a non-alignment based method to estimate transcript abundance.  A comparison of methods including Salmon and RSEM can be found here.  <a href="https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-4002-1">RNA Seq Computational Tools</a>  Go into the shrimp project directory and make a new directory called abundances.  Inside the abundances directory make another directory called matrix.  We will use this matrix directory immediately after estimating our abundances.</p>

<h2 align="center">Non-Alignment Based Abundance Estimation with Salmon</h2>

<p>I will use the shrimp data set to demonstrate the use of Salmon to produce abundance estimates but this procedure can be used for any data set.  It is useful to create a tab delimited file countaining the file paths to all of the samples.  We can use nano to do this.  Open nano in the shrimp project directory and build the samples file.  It should look similar to the following.</p>

```
SRR4341161	/home/nbumpus/shrimp/abundances/larvae1 /home/nbumpus/shrimp/trinity_reads/SRR4341161_1.trim.paired.adj.fastq   /home/nbumpus/shrimp/trinity_reads/SRR4341161_2.trim.paired.adj.fastq
SRR4341162	/home/nbumpus/shrimp/abundances/larvae2 /home/nbumpus/shrimp/trinity_reads/SRR4341162_1.trim.paired.adj.fastq   /home/nbumpus/shrimp/trinity_reads/SRR4341162_2.trim.paired.adj.fastq
SRR4341163	/home/nbumpus/shrimp/abundances/adult3  /home/nbumpus/shrimp/trinity_reads/SRR4341163_1.trim.paired.adj.fastq   /home/nbumpus/shrimp/trinity_reads/SRR4341163_2.trim.paired.adj.fastq
SRR4341164	/home/nbumpus/shrimp/abundances/adult4  /home/nbumpus/shrimp/trinity_reads/SRR4341164_1.trim.paired.adj.fastq   /home/nbumpus/shrimp/trinity_reads/SRR4341164_2.trim.paired.adj.fast
```
<p>The first column contains the name of the sample, followed by the output directory for the abundance estimates and finally the left reads followed by the right reads for the sample.  Next prepare the reference like so.</p>

```
#!/bin/bash -l
#PBS -q bio
#PBS -N referenceprep
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:10:00
#PBS -o out.txt
#PBS -e err.txt

cd #PBS_O_WORKDIR

module load trinity/2.8.4

$TRINITY_HOME/util/align_and_estimate_abundance.pl \
--transcripts /home/nbumpus/shrimp/trinity_out_dir/Trinity.fasta \
--est_method salmon \
--gene_trans_map /home/nbumpus/shrimp/trinity_out_dir/Trinity.fasta.gene_trans_map \
--prep_reference
```

<p>If you have an --SS_lib_type, as is the case with the lobster and yeast data, it should be added to the arguments here.  Next run Salmon with the following script.</p>

```
#!/bin/bash -l
#PBS -q bio
#PBS -N abundances
#PBS -l nodes=1:ppn=8
#PBS -l walltime=02:00:00
#PBS -o out.txt
#PBS -e err.txt

cd #PBS_O_WORKDIR

module load trinity/2.8.4

$TRINITY_HOME/util/align_and_estimate_abundance.pl \
--transcripts /home/nbumpus/shrimp/trinity_out_dir/Trinity.fasta \
--seqType fq \
--est_method salmon \
--samples_file /home/nbumpus/shrimp/samples_file.txt \
--salmon_add_opts --validateMappings \
--gene_trans_map /home/nbumpus/shrimp/trinity_out_dir/Trinity.fasta.gene_trans_map \
--output_dir /home/nbumpus/shrimp/abundances/
--thread_count 8
```
<p>The current version of Salmon on Marconi is 0.12.0 which requires us to add the --validateMappings argument which we can do through the --salmon_add_opts argument.  --validateMappings is on by default versions of Salmon following 0.12.0.  We can dictate resources to be used by the --thread_count argument making sure that the --thread_count value is the same as the ppn in the Torque settings.</p>

<p>Once the abundances have been estimated we can generate counts and expression matrices by running the followin script.  First we should create a txt file containing the paths to all of the quant.sf files generated during abundance estimation.  This file does not need to be tab dilimited and is simply a list.  Use nano to create a quant_files.txt file in the shrimp project directory like so.</p>

```
/home/nbumpus/shrimp/abundances/larvae1/quant.sf
/home/nbumpus/shrimp/abundances/larvae2/quant.sf
/home/nbumpus/shrimp/abundances/adult3/quant.sf
/home/nbumpus/shrimp/abundances/adult4/quant.sf
```
<p>Before running the trinity abundance__estimates_to_matrix.pl script we should build temporary R libraries.  Once these libraries are installed they will be available whenrever you are logged into Marconi.  In the shrimp project library use the module load command to import R/3.5.2 and then enter the R environment. Install the following temporary libraries using the commands below.  Some of the packages are needed later on during differential expression but now is a good time to install everything we will need</p>

```
source("http://bioconductor.org/biocLite.R")
biocLite('edgeR')
biocLite('ctc')
biocLite('Biobase')
install.packages('gplots')
install.packages('ape')
```

<p>Exit R and run the following script to generate the matrices.</p>

```
#!/bin/bash -l
#PBS -q bio
#PBS -N matrix
#PBS -l nodes=1:ppn=1
#PBS -l walltime=01:00:00
#PBS -o out.txt
#PBS -e err.txt

cd #PBS_O_WORKDIR

module load trinity/2.8.4
module load R/3.5.2

$TRINITY_HOME/util/abundance_estimates_to_matrix.pl \
--est_method salmon \
--gene_trans_map /home/nbumpus/shrimp/trinity_out_dir/Trinity.fasta.gene_trans_map \
--name_sample_by_basedir \
--out_prefix /home/nbumpus/shrimp/abundances/matrix/shrimp \
--quant_files /home/nbumpus/shrimp/quant_files.txt
```


