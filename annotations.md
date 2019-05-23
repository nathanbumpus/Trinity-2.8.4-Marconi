<h1 align="center"><a id="top"></a>Annotations</h1>

<P>Once we have the differentially expressed isoforms we can begin to functionally annotate them.  There are a few ways to do this. Here we will use blast+ to get functional annotations.  For the shrimp data we already have some functional annotations that we obtained from blasting our assembly against the uniprot database.  We can use the command line to extract the Trinity id's of the differentially expressed isoforms, add the annotations we already have to these id's and extract the id's that do not yet have an annotation to create a new list of id's to blast against another database.  I will demonstrate using the Trembl database.  First make a new directory in the shrimp directory called annotations and enter the directory.  Then use the awk command to extract a list of differentially expressed isoforms from the DE.subset file located in the DifferentialExpression/isoform/ directory.</p>

```
awk 'BEGIN{ FS=OFS="\t" }{ print $1 }' ../DifferentialExpression/isoform/shrimp.isoform.counts.matrix.adult_vs_larvae.edgeR.DE_results.P0.001_C2.DE.subset | grep TRINITY > isoDElist.txt
```
<p>Use awk again to extract the trinity ids that already have functional annotations into a toberemoved file.</p>

```
awk 'BEGIN{ FS=OFS="\t" }{ print $1 }' ../assembly_quality/pooled.blastx.outfmt6.w_pct_hit_length | grep TRINITY > toberemoved.txt
```
<p> awk extracts the first column of each of the files then pipes to grep to remove the header and write each set of id's to their respective files. Next use grep to make a newblastfile.txt file by removing the ids in the toberemoved.txt file from the isoDElist.txt file like so</p>

```
grep -vwF -f toberemoved.txt isoDElist.txt > newblastfile.txt
```

<p>We should now have a new list of trinity id's that did not get a hit during the initial blasting.  Next add the annotations that we already have to the toberemoved.txt file.</p>

```
grep -wF -f toberemoved.txt ../assembly_quality/pooled.blastx.outfmt6.w_pct_hit_length > DEisoAnnotations
```

<p>The DEisoAnnotations file contains the annotations we already have for the differentially expressed isoforms.  We can now extract sequences associated with the id's in the newblastfile.txt from the Trinity.fasta file to make a new fasta file for blasting with the following script</p>

```
#!/bin/bash -l
#PBS -q haswell
#PBS -N Fasta_Extract
#PBS -l nodes=1:ppn=1
#PBS -l walltime=01:00:00
#PBS -o /home/nbumpus/shrimp/annotations/
#PBS -e err.txt

cd #PBS_O_WORKDIR

module load trinity/2.8.4

$TRINITY_HOME/util/misc/acc_list_to_fasta_entries.pl \
/home/nbumpus/shrimp/annotations/newblastfile.txt \
/home/nbumpus/shrimp/trinity_out_dir/Trinity.fasta
```
<p>Next, after downloading and unzipping the Trembl database, build the blastable database.</p>

```
#!/bin/bash -l
#PBS -q haswell
#PBS -N builddb
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:05:00
#PBS -o out.txt
#PBS -e err.txt

cd #PBS_O_WORKDIR

module load blast+

makeblastdb -in /home/nbumpus/databases/uniprot_trembl.fasta \
-dbtype prot
```
<p>Then blast the new fasta file against the trembl database</p>

```
#!/bin/bash -l
#PBS -q haswell
#PBS -N blast-trembl
#PBS -l nodes=1:ppn=16
#PBS -l walltime=48:00:00
#PBS -o out.txt
#PBS -e err.txt

cd #PBS_O_WORKDIR

module load blast+

blastx \
-query /home/nbumpus/shrimp/annotations/<name of new fasta file> \
-db /home/nbumpus/databases/uniprot_trembl.fasta \
-out /home/nbumpus/shrimp/annotations/trembl.blastx.outfmt6 \
-evalue 1e-20 \
-num_threads 16 \
-max_target_seqs 1 \
-outfmt 6
```
<p>Next use the blastcoverage script below to add the annotations.</p>

```
#!/bin/bash -l
#PBS -q haswell
#PBS -N blastcoverage
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:05:00
#PBS -o out.txt
#PBS -e err.txt

cd #PBS_O_WORKDIR

module load trinity/2.8.4

$TRINITY_HOME/util/analyze_blastPlus_topHit_coverage.pl \
/home/nbumpus/shrimp/assembly_quality/trembl.blastx.outfmt6 \
/home/nbumpus/shrimp/annotations/<name of new fasta file> \
/home/nbumpus/databases/uniprot_trembl.fasta
```

<p>We can repeat these steps and blast any database we wish. Finally, use the cat command to add the new annotations to the DEisoAnnotations file.</p>

```
cat DEisoAnnotations trembl.blast.outfmt6 > totalDEannotations
```
<p>We can keep repeating this process to reduce the size of the unannotated Trinity id's.  Note that when blasting a large database like NR it is helpful to break the newblastfile.txt file into multiple files containing smaller amounts of unannotated Trinity id's.  This can be accomplished using the split command like so</p>

```
split -l 500 newblastfile.txt
```
<p>This splits the newblastfile.txt into a number of files each containing 500 ids.  The sequence extractions can then be performed for each of these new files and the blast scripts for each file can be run in parallel on the cluster to save time.  Then catonate all of the resulting outfmt6 files into one file and run the blastcoverage script to add the annotations.</p>

<h2 align="center">Table of Contents<a id="contents"></a></h2>
* [Home](README.md)
* [Obtaining Data](data.md)
* [Data Quality](dataqc.md)
* [Building an Assembly](assembly.md)
* [Abundance Estimation](abundance.md)
* [Assembly Quality](assemblyqc.md)
* [Differential Expression](DE.md)
* [Results](results.md)
* [Scripts for Lobster Assembly](Lobster_Scripts.md)
* [References](references.md)

<p><a href="#top">back to top</a></p>





