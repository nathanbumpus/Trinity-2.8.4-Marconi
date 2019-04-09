<h1 align="center">Estimating Transcript Abundance with Trinity</h1>

<p>When estinmating transcript abundance with Trinity on Marconi we need to choose whether to use an alignment based or a non-alignment based approach.  Trinity has a script align_and_estimate_abundance.pl which has direct support for the RSEM and Salmon modules on marconi.  If an alignment based method is prefered RSEM should be used.  Salmon should be used for non-alignment based estimation.  Before deciding which method to use go into the shrimp project directory and make a new directory called abundances.  Inside the abundances directory make another directory called matrix.  We will use this matrix directory immediately after estimating our abundances.</p>

<h2 align="center">Non-Alignment Based Abundance Estimation with Salmon</h2>

<p>I will use the shrimp data set to demonstrate the use of Salmon to produce abundance estimates but this procedure can be used for any data set.  It is useful to create a tab delimited file countaining the file paths to all of the samples.  We can use nano to do this.  Open nano in the shrimp project directory and build the samples file.  It should look similar to the following.</p>

```
SRR4341161	/home/nbumpus/shrimp/abundances/larvae1 /home/nbumpus/shrimp/trinity_reads/SRR4341161_1.trim.paired.adj.fastq   /home/nbumpus/shrimp/trinity_reads/SRR4341161_2.trim.paired.adj.fastq
SRR4341162	/home/nbumpus/shrimp/abundances/larvae2 /home/nbumpus/shrimp/trinity_reads/SRR4341162_1.trim.paired.adj.fastq   /home/nbumpus/shrimp/trinity_reads/SRR4341162_2.trim.paired.adj.fastq
SRR4341163	/home/nbumpus/shrimp/abundances/adult3  /home/nbumpus/shrimp/trinity_reads/SRR4341163_1.trim.paired.adj.fastq   /home/nbumpus/shrimp/trinity_reads/SRR4341163_2.trim.paired.adj.fastq
SRR4341164	/home/nbumpus/shrimp/abundances/adult4  /home/nbumpus/shrimp/trinity_reads/SRR4341164_1.trim.paired.adj.fastq   /home/nbumpus/shrimp/trinity_reads/SRR4341164_2.trim.paired.adj.fast
```
<p>The first column contains the name of the sample, followed by the output directory for the abundance estimates and finally the left reads followed by the right reads for the sample.</p>

<p>To use Salmon run the following script.</p>

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
