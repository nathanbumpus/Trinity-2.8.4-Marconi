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

<p>To check the status of the assembly construction we can us the top command.  First obtain the job id by using the qstat command.  Then use the checkjob command and provide the job id.  The terminal output will contain an "Allocated Nodes" section telling us where the job is running.  Next ssh gromit which is the job scheduler.  From gromit we can ssh the node the job is running on.  Then use the top command.  We can use the -u argument and provide our username to show only the jobs we are running on the node.  Here we can see the %CPU which divided by 100 will tell us how many CPU's we are using.  RES shows the memory usage for the the job.  If we type "1" we can toggle and see which CPU's on the node are actively running our job.  If we type "c" the command column will expand to show exactly which script is running and which files are being used.  Below is the output for a top command for a job running on the n151 node.</p>

```
top - 11:01:44 up 86 days, 15:34,  1 user,  load average: 1.64, 2.08, 2.32
Tasks: 279 total,   2 running, 277 sleeping,   0 stopped,   0 zombie
%Cpu0  :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu1  :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu2  : 86.4 us, 13.6 sy,  0.0 ni,  0.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu3  :  0.0 us,  0.3 sy,  0.0 ni, 99.7 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu4  :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu5  :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu6  :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu7  :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu8  :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu9  :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu10 :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu11 :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu12 :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu13 :  0.3 us,  0.3 sy,  0.0 ni, 99.3 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu14 :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu15 :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu16 :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu17 :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu18 :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu19 :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu20 :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu21 :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu22 :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu23 :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
KiB Mem : 26404700+total,  5243064 free,  5749108 used, 25305483+buff/cache
KiB Swap:        0 total,        0 free,        0 used. 25601696+avail Mem 

  PID USER      PR  NI    VIRT    RES    SHR S  %CPU %MEM     TIME+ COMMAND                                                                                                                                        
 9619 nbumpus   20   0 4252636 3.824g   2296 R 100.0  1.5  54:46.03 perl /opt/modules/centos7/trinity/2.8.4/util/insilico_read_normalization.pl --seqType fq --JM 240G --max_cov 200 --min_cov 1 --CPU 16 --outpu+ 
10769 nbumpus   20   0  157868   2428   1548 R   0.3  0.0   0:00.15 top -u nbumpus                                                                                                                                 
 9454 nbumpus   20   0   11628   1572   1316 S   0.0  0.0   0:00.01 -bash                                                                                                                                          
 9477 nbumpus   20   0    9504   1496   1228 S   0.0  0.0   0:00.00 /bin/bash -l /var/spool/torque/mom_priv/jobs/102250.gromit.SC                                                                                  
 9502 nbumpus   20   0   36440   8640   2320 S   0.0  0.0   0:00.11 perl /opt/modules/centos7/trinity/2.8.4/Trinity --seqType fq --SS_lib_type RF --left /home/nbumpus/lobster/trimmed_reads/all/all.left.fastq -+ 
10739 nbumpus   20   0  136716   2200    984 S   0.0  0.0   0:00.00 sshd: nbumpus@pts/0                                                                                                                            
10740 nbumpus   20   0  115516   2096   1616 S   0.0  0.0   0:00.00 -bash          
```
  
 <a href="#top">back to top</a></p>


<h2 align="center">Table of Contents</h2>
* [Home](README.md)
* [Obtaining Data](data.md)
* [Data Quality](dataqc.md)
* [Building an Assembly](assembly.md)
* [Abundance Estimation](abundance.md)
* [Assembly Quality](assemblyqc.md)
* [Differential Expression](DE.md)
