<h1 align="center">Building a Trinity Assembly on Marconi<a id="top"></a></h1>

<p>Now that we have quality sets of data for the yeast, lobster and shrimp we can build a trinity assembly for each. First we need to know how the cDNA library was constructed.  Often the literature that supplied the data will tell us this, as was the case for the yeast and lobster, but sometimes, in the case of the shrimp data this information is omitted.  When the latter occurs there are a few things we can do.  First, we can take a sub set of our reads and align them back to a reference genome of a similar organism.  Then view the alignments in IGV by viewing as pairs and coloring alignments by first-of-pair-strand.  Secondly, we could build the assembly, align the reads back to the assembly and then view in IGV.  I will demonstrate this in the Assembly Quality section.  But for now let's build an assembly.</p>

<p>First go into the trinity_reads directory inside the shrimp project directory.  We can group all of left reads and right reads into two separate files using the cat command.</p>

```
cat *_1.trim.paired.adj.fastq > reads.all.left.fastq
```

<p>Then</p>

```
cat *_2.trim.paired.adj.fastq > reads.all.right.fastq
```
Now go back one directory to the shrimp project directory and run the following script.

```
#!/bin/bash -l
#PBS -q bio
#PBS -N pooled-unstranded
#PBS -l nodes=1:ppn=12
#PBS -l walltime=24:00:00
#PBS -o out.txt
#PBS -e err.txt

cd #PBS_O_WORKDIR

module load python
module load trinity/2.8.4

Trinity --seqType fq \
--left /home/nbumpus/shrimp/trimmed_reads/all.left.adj.fastq \
--right /home/nbumpus/shrimp/trimmed_reads/all.right.adj.fastq \
--output /home/nbumpus/shrimp/trinity_out_dir \
--CPU 12 \
--max_memory 240G
```
<p>This version of trinity requires the python module to run one of its scripts.  We then give trinity the file types 
(fastq = --seqType fq).  If this were strand specific data such as the case with the lobster and yeast data and used dUTP 
in the cDNA library construction we add the --SS_lib_type RF argument after the --seqType fq.  We then present trinity with all of our left reads.  If we had not previously combined all of the left reads into a single file we could provide each file individually separated by a comma with no spaces.  Do the same for the right reads.  Make sure that the --CPU argument matches the #PBS ppn parameter to not overload the node.  A rule of thumb for the --max_memory setting is 1G for every 1 million reads.  This run will take anywhere between  12 and 20 hours to complete.  --CPU can be set anywhere between 8 and 16 depending on time constraints but will only be utilized during the chrysallis phase.  If the run crashes during jellyfish due to memory issues you can try reducing the --max_memory parameter.  Make sure the script is running on the bio cluster to avoid this pitfall.  Note that by default trinity 2.8.4 will normalize your reads.  2.2.0 will not.</p>

<h2 align="center">Monitoring the Assembly Construction</h2>

<p>The only way that I have found to monitor the assembly construction is to go into the trinity_out_dir and compare the length of the recursive_trinity.cmds file to the recursive_trinity.cmds.completed file.  If all the commands are completed you can figure on about one or two hours for completion of the assembly Next we need to check how good our assembly is and if it is appropriate for downstream analysis. <a href="#top">back to top</a></p>


<h2 align="center">Table of Contents</h2>
* [Home](README.md)
* [Obtaining Data](data.md)
* [Data Quality](dataqc.md)
* [Building an Assembly](assembly.md)
* [Abundance Estimation](abundance.md)
* [Assembly Quality](assemblyqc.md)
* [Differential Expression](DE.md)
