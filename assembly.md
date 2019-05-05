<h1 align="center">Building a Trinity Assembly on Marconi<a id="top"></a></h1>

<p>Now that we have quality sets of data for the yeast, lobster and shrimp we can build a Trinity assembly for each organism. First we need to know how the cDNA library was constructed.  Often the literature that supplied the data will tell us this, as is the case for the yeast and lobster.  However, in the case of the shrimp, this information is omitted.  When this occurs there are a couple things we can do.  First, we can take a sub set of our reads and align them back to a reference genome of a similar organism.  Then view the alignments in IGV by viewing as pairs and coloring alignments by first-of-pair-strand.  Secondly, we could build the assembly, align the reads back to the assembly and then view in IGV and construct violin plots predicting strand specificity.  I will demonstrate this in the <a href="https://nathanbumpus.github.io/Trinity-2.8.4-Marconi/assemblyqc.html">Assembly Quality</a> section, but for now let's build an assembly.</p>

<p>First go into the trinity_reads directory inside the shrimp project directory.  We can group all of the left reads together and all of the right reads together using the cat command.</p>

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
<p>This version of Trinity requires the python module to run one of its scripts.  We then give Trinity the file types 
(fastq = --seqType fq).  If this were strand specific data such as the case with the lobster and yeast data and used dUTP 
in the cDNA library construction we add the --SS_lib_type RF argument after the --seqType fq.  We then present Trinity with all of our left reads.  If we had not previously combined all of the left reads into a single file we could provide each file individually separated by a comma with no spaces.  Do the same for the right reads.  Make sure that the --CPU argument matches the #PBS ppn parameter to not overload the node.  A rule of thumb for the --max_memory setting is 1G for every 1 million reads.  This run will take anywhere between  12 and 20 hours to complete. Note the lobster assembly takes between three and four days.  --CPU can be set anywhere between 8 and 16 depending on time constraints however, all of the allocated resources will only be used during certain phases of the assembly construction.  The --max_memory parameter should be set to below the available memory of the node.  If you run the assembly on fast1 you can take advantage of all of the ~260 GB of memory by setting the --CPU to 12.  This will leave the other bio node free for others to use.  Note that by default trinity 2.8.4 will normalize your reads.  2.2.0 will not.</p>

<h2 align="center">Monitoring the Assembly Construction</h2>

<p>To check the status of the assembly construction we can us the top command.  First obtain the job id by using the qstat command.  Then use the checkjob command and provide the job id.  The terminal output will contain an "Allocated Nodes" section telling us where the job is running.  Next ssh gromit which is the job scheduler.  From gromit we can ssh the node the job is running on.  Then use the top command.  We can use the -u argument and provide our username to show only the jobs we are running on the node.  Here we can see the %CPU which divided by 100 will tell us how many CPU's we are using.  RES shows the memory usage for the the job.  If we type "1" we can toggle and see which CPU's on the node are actively running our job.  If we type "c" the command column will expand to show exactly which script is running and which files are being used.  Below is the output for a top command for a job running on the n151 node.</p>

```
[nbumpus@n151 ~]$ top -u nbumpus

top - 11:46:02 up 86 days, 16:18,  1 user,  load average: 5.37, 4.38, 5.08
Tasks: 289 total,   2 running, 287 sleeping,   0 stopped,   0 zombie
%Cpu0  : 95.0 us,  4.6 sy,  0.0 ni,  0.3 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu1  :  2.0 us,  0.3 sy,  0.0 ni, 97.7 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu2  : 94.7 us,  5.3 sy,  0.0 ni,  0.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu3  :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu4  :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu5  :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu6  :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu7  :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu8  : 95.7 us,  4.0 sy,  0.0 ni,  0.3 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu9  : 94.7 us,  5.0 sy,  0.0 ni,  0.3 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu10 :  0.0 us,  1.3 sy,  0.0 ni, 98.0 id,  0.0 wa,  0.0 hi,  0.7 si,  0.0 st
%Cpu11 : 95.4 us,  4.3 sy,  0.0 ni,  0.3 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu12 : 89.0 us,  8.3 sy,  0.0 ni,  2.7 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu13 :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu14 :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu15 :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu16 :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu17 :  0.0 us,  0.3 sy,  0.0 ni, 99.7 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu18 :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu19 :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu20 :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu21 :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu22 :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
%Cpu23 :  0.0 us,  0.0 sy,  0.0 ni,100.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
KiB Mem : 26404700+total, 13341980+free, 13475548 used, 11715164+buff/cache
KiB Swap:        0 total,        0 free,        0 used. 24829128+avail Mem 

  PID USER      PR  NI    VIRT    RES    SHR S  %CPU %MEM     TIME+ COMMAND                                                                                                                                        
11086 nbumpus   20   0 11.351g 0.011t   1616 R 597.4  4.4  37:45.93 /opt/modules/centos7/trinity/2.8.4/Inchworm/bin//inchworm --kmers jellyfish.kmers.fa --run_inchworm -K 25 -L 25 --monitor 1 --num_threads 6 -+ 
 9454 nbumpus   20   0   11628   1572   1316 S   0.0  0.0   0:00.01 -bash                                                                                                                                          
 9477 nbumpus   20   0    9504   1496   1228 S   0.0  0.0   0:00.00 /bin/bash -l /var/spool/torque/mom_priv/jobs/102250.gromit.SC                                                                                  
 9502 nbumpus   20   0  188024  12180   1684 S   0.0  0.0   0:06.02 perl /opt/modules/centos7/trinity/2.8.4/Trinity --seqType fq --SS_lib_type RF --left /home/nbumpus/lobster/trimmed_reads/all/all.left.fastq -+ 
11049 nbumpus   20   0  136716   2204    992 S   0.0  0.0   0:00.07 sshd: nbumpus@pts/0                                                                                                                            
11050 nbumpus   20   0  115516   2144   1656 S   0.0  0.0   0:00.02 -bash                                                                                                                                          
11085 nbumpus   20   0    9512   1156    964 S   0.0  0.0   0:00.00 sh -c /opt/modules/centos7/trinity/2.8.4/Inchworm/bin//inchworm --kmers jellyfish.kmers.fa --run_inchworm -K 25 -L 25 --monitor 1   --num_thr+ 
11152 nbumpus   20   0  157972   2548   1636 R   0.0  0.0   0:00.04 top -u nbumpus                     
```
<p>We can also get some more information about the job by using the cat command with the PID value.  First use ctrl-c to exit top then use the cat command for the for the job like so.</p>
```
[nbumpus@n151 ~]$ cat /proc/11086/status
Name:	inchworm
State:	R (running)
Tgid:	11086
Ngid:	11088
Pid:	11086
PPid:	11085
TracerPid:	0
Uid:	1647	1647	1647	1647
Gid:	1557	1557	1557	1557
FDSize:	64
Groups:	1557 
VmPeak:	12258676 kB
VmSize:	12198288 kB
VmLck:	       0 kB
VmPin:	       0 kB
VmHWM:	11968540 kB
VmRSS:	11968540 kB
RssAnon:	11966924 kB
RssFile:	     948 kB
RssShmem:	     668 kB
VmData:	12176828 kB
VmStk:	     140 kB
VmExe:	     168 kB
VmLib:	    4652 kB
VmPTE:	   23468 kB
VmSwap:	       0 kB
Threads:	6
SigQ:	0/1031307
SigPnd:	0000000000000000
ShdPnd:	0000000000000000
SigBlk:	0000000000000000
SigIgn:	0000000000001000
SigCgt:	0000000180000000
CapInh:	0000000000000000
CapPrm:	0000000000000000
CapEff:	0000000000000000
CapBnd:	0000001fffffffff
CapAmb:	0000000000000000
Seccomp:	0
Cpus_allowed:	ffffff
Cpus_allowed_list:	0-23
Mems_allowed:	00000000,00000000,00000000,00000000,00000000,00000000,00000000,00000000,00000000,00000000,00000000,00000000,00000000,00000000,00000000,00000000,00000000,00000000,00000000,00000000,00000000,00000000,00000000,00000000,00000000,00000000,00000000,00000000,00000000,00000000,00000000,00000003
Mems_allowed_list:	0-1
voluntary_ctxt_switches:	31111
nonvoluntary_ctxt_switches:	1969
```
<p><a href="#top">back to top</a></p>


<h2 align="center">Table of Contents</h2>
* [Home](README.md)
* [Obtaining Data](data.md)
* [Data Quality](dataqc.md)
* [Abundance Estimation](abundance.md)
* [Assembly Quality](assemblyqc.md)
* [Differential Expression](DE.md)
* [Annotations](annotations.md)
* [Results](results.md)
* [References](references.md)
