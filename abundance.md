<h1 align="center">Estimating Transcript Abundance with Trinity</h>

<p>When estinmating transcript abundance with Trinity on Marconi we need to choose whether to use an alignment based or a non-alignment based approach.  Trinity has a script align_and_estimate_abundance.pl which has direct support for the RSEM and Salmon modules on marconi.  If an alignment based method is prefered RSEM should be used.  Salmon should be used for non-alignment based estimation.  Before deciding which method to use go into the shrimp project directory and make a new directory called abundances.  Inside the abundances directory make another directory called matrix.  We will use this matrix directory immediately after estimating our abundances.</p>

<h2 align="center">Non-Alignment Based Abundance Estimation with Salmon</h2>

<p>I will use the shrimp data set to demonstrate the use of Salmon to produce abundance estimates but this procedure can be used for any data set.  To use Salmon run the following script.</p>

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
