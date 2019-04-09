<h1 align="center">Evaluation of Assembly Quality</h1>

<p>Once transcript abundance estimates have been obtained there are a number of tools we can use to determine if our assemblies are acceptable for downstream analysis.  These include the generation of the N50 statistic, realigning the reads back to the assemblies, verifying strand specificity (if not already done), building a blastable database and blasting the assemblies and evaluating the Ex90N50 statistic.  First go to the shrimp project directory and make a new directory called assembly_quality.</p>

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

<h2 align="center">Realigning the Reads Back to the Assembly</h2>
