<h1 align="center">Results for the Yeast, Shrimp and Lobster Assemblies</h1>

<p>Below are the results generated from the trinity pipeline on marcroni for each of the three assemblies.</p>

<h2 align="center">Yeast Assembly</h2>

<p>A small set of of RNA seq data used to construct a trinity assembly from samples from four conditions with no replicates.</p>

<h3 align="center">Assembly Quality</h3>
<p>Below are the stats basic stats for the yeast trinity assembly including the N50 statistic</p>
```
################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':  8698
Total trinity transcripts:	9245
Percent GC: 38.09

########################################
Stats based on ALL transcript contigs:
########################################

        Contig N10: 3196
        Contig N20: 2601
        Contig N30: 2173
        Contig N40: 1861
        Contig N50: 1603

        Median contig length: 740
        Average contig: 1022.39
        Total assembled bases: 9452030


#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

        Contig N10: 3194
        Contig N20: 2599
        Contig N30: 2174
        Contig N40: 1869
        Contig N50: 1603

        Median contig length: 717
        Average contig: 1011.08
        Total assembled bases: 8794382
```
<p>Next the graph depicting the Ex90N50 statistic suggesting that deep enough sequencing was achieved to construct the assembly</p>

<p align="center">
  <img src="spEx90N50.25.jpg" alt="Ex90N50"> 
</p>

<h3 align="center">Alignment Statistics</h3>

<p>Below are the alignment rates of the reads from each of the samples to the trinity assembly</p>

```
ds reads
1000000 reads; of these:
  1000000 (100.00%) were paired; of these:
    74023 (7.40%) aligned concordantly 0 times
    811523 (81.15%) aligned concordantly exactly 1 time
    114454 (11.45%) aligned concordantly >1 times
    ----
    74023 pairs aligned concordantly 0 times; of these:
      9992 (13.50%) aligned discordantly 1 time
    ----
    64031 pairs aligned 0 times concordantly or discordantly; of these:
      128062 mates make up the pairs; of these:
        99862 (77.98%) aligned 0 times
        21441 (16.74%) aligned exactly 1 time
        6759 (5.28%) aligned >1 times
95.01% overall alignment rate

hs reads
1000000 reads; of these:
  1000000 (100.00%) were paired; of these:
    72134 (7.21%) aligned concordantly 0 times
    790420 (79.04%) aligned concordantly exactly 1 time
    137446 (13.74%) aligned concordantly >1 times
    ----
    72134 pairs aligned concordantly 0 times; of these:
      8381 (11.62%) aligned discordantly 1 time
    ----
    63753 pairs aligned 0 times concordantly or discordantly; of these:
      127506 mates make up the pairs; of these:
        98520 (77.27%) aligned 0 times
        21950 (17.21%) aligned exactly 1 time
        7036 (5.52%) aligned >1 times
95.07% overall alignment rate

log reads
1000000 reads; of these:
  1000000 (100.00%) were paired; of these:
    78228 (7.82%) aligned concordantly 0 times
    783620 (78.36%) aligned concordantly exactly 1 time
    138152 (13.82%) aligned concordantly >1 times
    ----
    78228 pairs aligned concordantly 0 times; of these:
      10453 (13.36%) aligned discordantly 1 time
    ----
    67775 pairs aligned 0 times concordantly or discordantly; of these:
      135550 mates make up the pairs; of these:
        105156 (77.58%) aligned 0 times
        23562 (17.38%) aligned exactly 1 time
        6832 (5.04%) aligned >1 times
94.74% overall alignment rate

plat reads
1000000 reads; of these:
  1000000 (100.00%) were paired; of these:
    76864 (7.69%) aligned concordantly 0 times
    822698 (82.27%) aligned concordantly exactly 1 time
    100438 (10.04%) aligned concordantly >1 times
    ----
    76864 pairs aligned concordantly 0 times; of these:
      10743 (13.98%) aligned discordantly 1 time
    ----
    66121 pairs aligned 0 times concordantly or discordantly; of these:
      132242 mates make up the pairs; of these:
        100536 (76.02%) aligned 0 times
        25258 (19.10%) aligned exactly 1 time
        6448 (4.88%) aligned >1 times
94.97% overall alignment rate
```
<h3 align="center">Strand Specificity</h3>

<p>Violin plot suggesting that library construction was performed in a strand specific fashion</p>

<p align="center">
  <img src="spviolinplot.25.jpg" alt="violinplot">
</p>

<h3 align="center">Differential Expression and Annotations</h3>

<p>Below are the MA and volcano plots for visualizing the differentially expressed isoforms between each of the pairwise comparisons and annotation statistics for the obtained differentially expressed isoforms.  The uniprot and trembl databases were used for the annotations.</p>

<h4 align="center">ds versus hs</h4>

<p align="center">
  <img src="dsVShsMA50.jpg" alt="MAplot">
  <img src="dsVShsVolcano50.jpg" alt="Volcanoplot">
</p>

<p>Blastp functional annotations were obtained for 86 of the 94 differentially expressed isoforms in this group</p>

<h4 align="center">ds versus log</h4>

<p align="center">
  <img src="dsVSlogMA50.jpg" alt="MAplot">
  <img src="dsVSlogVolcano50.jpg" alt="Volcanoplot">
</p>

<p>Blastp functional annotations were obtained for 99 of the top 100 differentially expressed isoforms in this group</p>

<h4 align="center">ds versus plat</h4>

<p align="center">
  <img src="dsVSplatMA50.jpg" alt="MAplot">
  <img src="dsVSplatVolcano50.jpg" alt="Volcanoplot">
</p>

<p>Blastp functional annotations were obtained for 68 of the 77 differentially expressed isoforms in this group</p>

<h4 align="center">hs versus log</h4>

<p align="center">
  <img src="hsVSlogMA50.jpg" alt="MAplot">
  <img src="hsVSlogVolcano50.jpg" alt="Volcanoplot">
</p>

<p>Blastp functional annotations were obtained for 94 of the top 100 differentially expressed isoforms in this group</p>

<h4 align="center">hs versus plat</h4>

<p align="center">
  <img src="hsVSplatMA50.jpg" alt="MAplot">
  <img src="hsVSplatVolcano50.jpg" alt="Volcanoplot">
</p>

<p>Blastp functional annotations were obtained for 98 of the top 100 differentially expressed isoforms in this group</p>

<h4 align="center">log versus plat</h4>

<p align="center">
  <img src="logVSplatMA50.jpg" alt="MAplot">
  <img src="logVSplatVolcano50.jpg" alt="Volcanoplot">
</p>

<p>Blastp functional annotations were obtained for 97 of the top 100 differentially expressed isoforms in this group</p>

<h4 align="center">Heatmaps for Determining Cluster Extractions</h4>

<p>Below is the heatmap to help visualize isoform clusters of interest.</p>

<p align="center">
  <img src="spHeatmap50.jpg" alt="Heatmap">
</p>

<h2 align="center">Shrimp Assembly</h2>

<p>Since most of the results for the shrimp data are shown in the building of the assembly I will not repost them here.  However, a third database was used to try to recover more functional annotations for the differentially expressed isoforms.  The non redundant NCBI data was used for this purpose. Using the uniprot, trembl and nr databases functional annotations were obtained for 9411 of the 16836 differenially expressed isoforms.</p>

<h3 align="center">Lobster Assembly</h3>





