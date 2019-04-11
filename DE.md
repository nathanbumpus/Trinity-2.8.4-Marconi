<h1 align="center">Differential Expression with Trinity on Marconi<a id="top"></a></h1>

<p>Once we have determined that our assemblies are of high quality we can look at which transcripts are differentially expressed.  Trinity supports the use of edgeR, DESeq2 and limma.  I chose to use edgeR because the yeast and lobster data sets do not contain biological replicates.  In this case edgeR is the only option as we can set a dispersion estimate.  The recommended dispersion estimate is between 0.1 and 0.4 depending on the organism. For the yeast I picked a dispersion value of 0.1 but more can read about dispersion in the edgeR documentation <a href="http://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf" target="_blank">here.</a> In the case of the shrimp data we have two larvae replicates and two adult replicates.  We need to create a txt file describing these sample to replicate relationships to feed into edgeR.  Use nano to create the following tab delimited samples_described.txt file.</p>

```
larvae  larvae1
larvae  larvae2
adult   adult3
adult   adult4
```

<p>Next make a new directory within the project directory called DifferentialExpression and use the following script to determine Differential Expression.</p>

```
#!/bin/bash -l
#PBS -q bio
#PBS -N shr-DE-iso
#PBS -l nodes=1:ppn=1
#PBS -l walltime=02:00:00
#PBS -o out.txt
#PBS -e err.txt

cd #PBS_O_WORKDIR

module load trinity/2.8.4
module load R/3.5.2

$TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl \
--matrix /home/nbumpus/shrimp/abundances/matrix/shrimp.isoform.counts.matrix \
--method edgeR \
--samples_file /home/nbumpus/shrimp/samples_described.txt \
--output /home/nbumpus/shrimp/DifferentialExpression/isoform/
```
<p>Note that this script calculates the differential expression at the isoform level.  You can also calculate differential expression at the gene level by using the shrimp.gene.counts.matrix and output to a new gene folder within the DifferentialExpression directory.  For the yeast data add the --dispersion 0.1 argument after the --method edgeR argument.  The output will give a set of differentially expressed isoforms, sets of up-regulated isoforms and MA and volcano plots describing differential expression between all pairwise comparisons.  Below are the MA and volcano plots for the shrimp data</p>

<h2 align="center">Table of Contents<a id="contents"></a></h2>
* [Home](README.md)
* [Obtaining Data](data.md)
* [Data Quality](dataqc.md)
* [Building an Assembly](assembly.md)
* [Abundance Estimation](abundance.md)
* [Assembly Quality](assemblyqc.md)
* [Differential Expression](DE.md)
