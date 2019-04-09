<h1 align="center"><a id="top"></a>Evaluation of Assembly Quality</h1>

<p>Once transcript abundance estimates have been obtained there are a number of tools we can use to determine if our assemblies are acceptable for downstream analysis.  These include the generation of the N50 statistic, <a href="#align">realigning the reads back to the assemblies</a>, <a href="#specificity">verifying strand specificity</a> (if not already done), building a blastable database and blasting the assemblies and evaluating the Ex90N50 statistic.  First go to the shrimp project directory and make a new directory called assembly_quality.</p>

<h2 align="center">Generation of the N50 statistic</h2>

<p>To generate the Nx50 statistic for the shrimp data run the following script.</p>

```
#!/bin/bash -l
#PBS -q bio
#PBS -N AssemblyStats
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:05:00
#PBS -o /home/nbumpus/shrimp/assembly_quality/assembly_stats.txt
#PBS -e err.txt

cd #PBS_O_WORKDIR

module load trinity/2.8.4

$TRINITY_HOME/util/TrinityStats.pl \
/home/nbumpus/shrimp/trinity_out_dir/Trinity.fasta
```

<p>This will generate the assembly_stats.txt file that looks something like this.</p>

```

################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':  124180
Total trinity transcripts:	186206
Percent GC: 39.17

########################################
Stats based on ALL transcript contigs:
########################################

        Contig N10: 9401
        Contig N20: 7201
        Contig N30: 5808
        Contig N40: 4663
        Contig N50: 3673

        Median contig length: 460
        Average contig: 1372.72
        Total assembled bases: 255608356


#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

        Contig N10: 7728
        Contig N20: 5416
        Contig N30: 3974
        Contig N40: 2908
        Contig N50: 2050

        Median contig length: 362
        Average contig: 876.75
        Total assembled bases: 108874612

```
<p>If we look at only the longest isoform per gene and line up all of our assembled contigs by decreasing length the N50 is the length of the contig at the at the point of coverage at 50 percent of the length of the total assembly.</p>

<h2 align="center">Realigning the Reads Back to the Assembly<a id="align"></a></h2>

<p>To map the reads back to the assembly we can use Bowtie2/2.3.4.3 on Marconi.  First build a Bowtie2 index.</p>

```
#!/bin/bash -l
#PBS -q bio
#PBS -N Bowtie2-build
#PBS -l nodes=1:ppn=1
#PBS -l walltime=02:00:00
#PBS -o out.txt
#PBS -e err.txt

cd #PBS_O_WORKDIR

module load Bowtie2/2.3.4.3

bowtie2-build /home/nbumpus/shrimp/trinity_out_dir/Trinity.fasta \
/home/nbumpus/shrimp/assembly_quality/shrimp_ref
```
<p>Then align each set of reads back to the assembly modifying the following script as necessary for each pair of left and right reads</p>

```
#!/bin/bash -l
#PBS -q bio
#PBS -N bowtie2-shrimp1
#PBS -l nodes=1:ppn=5
#PBS -l walltime=05:00:00
#PBS -o out.txt
#PBS -e err.txt

cd #PBS_O_WORKDIR

module load perl
module load Bowtie2/2.3.4.3
module load samtools

bowtie2 -p 4 -q \
--no-unal -k 20 -x /home/nbumpus/shrimp/assembly_quality/shrimp_ref \
-1 /home/nbumpus/shrimp/trinity_reads/SRR4341161_1.trim.paired.adj.fastq \
-2 /home/nbumpus/shrimp/trinity_reads/SRR4341161_2.trim.paired.adj.fastq \
2>/home/nbumpus/shrimp/assembly_quality/SRR4341161_align_stats.txt| \
samtools view -@4 -Sb \
-o /home/nbumpus/shrimp/assembly_quality/SRR4341161.bowtie2.bam
```
<p>This script tell Bowtie2 to use 4 threads and when aligning the fastq files (-q), then suppress SAM records for all unaligned reads and to report up to 20 alignments per read.  The results are written to an align_stats.txt file and piped to samtools to create a bam file.  The thread count from Bowtie2 is also piped to samtools and should be less than or equal to the ppn in the Torque settings.  We then specify that the output file should be .bam.</p>

<p>If we open the align_stats.txt file we will see something like this.</p>

```
45248026 reads; of these:
  45248026 (100.00%) were paired; of these:
    3135976 (6.93%) aligned concordantly 0 times
    15302515 (33.82%) aligned concordantly exactly 1 time
    26809535 (59.25%) aligned concordantly >1 times
    ----
    3135976 pairs aligned concordantly 0 times; of these:
      543470 (17.33%) aligned discordantly 1 time
    ----
    2592506 pairs aligned 0 times concordantly or discordantly; of these:
      5185012 mates make up the pairs; of these:
        2010855 (38.78%) aligned 0 times
        414129 (7.99%) aligned exactly 1 time
        2760028 (53.23%) aligned >1 times
97.78% overall alignment rate
```
<p>This file tells us that ~93% of the reads from this sample aligned concordantly at least one time.  Bowtie2 then took the ~7% of the reads that did not align concordantly and tried to align them discordantly.  Bowtie2 was able to align ~17 percent of these reads discordantly.  Bowtie2 then attempted to align the remaining unaligned reads to the assembly as singletons and was successful with about 61% of these reads.  This gives an overall alignment rate of ~98%.  The mapping shown is an example of a sucessful alignment showing greater than 70% of the reads aligning concordantly and a high overall alignment rate greater than 80%. <a href="#top">back to top </a><a href="#contents">table of contents</a></p>

<h2 align="center">Determining Strand Specificity<a id="specificity"></a></h2>



