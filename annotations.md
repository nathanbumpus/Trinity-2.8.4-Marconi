<h1 align="center">Annotations</h1>

<P>Once we have the differentially expressed isoforms we can begin to functionally annotate them.  There are a few ways to do this.
Here we will use blast+ to get functional annotations.  For the shrimp data we already have some functional annotations that we 
obtained from blasting our assembly against the uniprot database.  We can use command line to remove these annotated trinity id's 
from the Trinity.fasta.gene_trans_map file and then extract the sequences from the Trinity.fasta file to create a new list for 
blasting.  Make a new directory in the shrimp directory called annotations and make a directory inside the annotations directory
called isoform.  Then use the awk command to extract the second column of the gene_trans_map file into a new file like so.</p>

```
awk 'BEGIN{ FS=OFS="\t" }{ print $2 }' ../../trinity_out_dir/Trinity.fasta.gene_trans_map > isolist.txt
```
<p>Use awk again to extract the trinity ids that already have functional annotations into a toberemoved file.</p>

```
awk 'BEGIN{ FS=OFS="\t" }{ print $1 }' ../../assembly_quality/pooled.blastx.outfmt6.w_pct_hit_length | grep TRINITY > toberemoved.txt
```
<p> awk takes the first column of the hit_length file.  The pipe to grep then removes the header so that all that is left to be written
are the trinity id's.  Then use grep to make a newblastfile.txt file by removing the ids in the toberemoved.txt file from the 
isolist.txt file like so</p>

```
grep -vwF -f toberemoved.txt isolist.txt > newblastfile.txt
```

<p>We should now have a new list of trinity id's that did not get a hit during the initial blasting.  We can now extract the sequences
from the Trinity.fasta file associated with our list of trinity ids with the following script</p>

```
#!/bin/bash -l
#PBS -q haswell
#PBS -N Fasta_Extract
#PBS -l nodes=1:ppn=1
#PBS -l walltime=01:00:00
#PBS -o /home/nbumpus/shrimp/annotations/isoform/
#PBS -e err.txt

cd #PBS_O_WORKDIR

module load trinity/2.8.4

$TRINITY_HOME/util/misc/acc_list_to_fasta_entries.pl \
/home/nbumpus/shrimp/annotations/isoform/newblastfile.txt \
/home/nbumpus/shrimp/trinity_out_dir/Trinity.fasta
```
<p>We can now blast any data base we wish by downloading and building using the steps outlined in the assembly quality
section. Once we get the new set of blast hits we can use the blastcoverage.sh script again to add annotations.  Finally use the 
cat command to add the new annotations to the orignal annotations file like so</p>

```
cat ../../assembly_quality/pooled.blastx.outfmt6.w_pct_hit_length shrimp.trembl.blast.outfmt6 > shrimp.annotations
```
<p>We can keep repeating this process to reduce the size of the unannotated trinity id's.  I used the Trembl database.  The next 
database I used was the NR database.  When blasting a large database like NR it is helpful to break the newblastfile.txt file into 
multiple files containing smaller amounts of unannotated trinity id's.  This can be accomplished using the split command like so</p>

```
split -l 10000 newblastfile.txt
```
<p>This splits the newblastfile.txt into a number of files each containing 10000 ids.  The sequence extractions can then be performed
for each of these new files and the blast scripts for each file can be run in parallel on the cluster.  Again the blast script is 
outlined in the assembly quality section.</p>






