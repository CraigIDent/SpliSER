<br>

<p align="center">
  <img src="Images/SpliSER.png" width="200">
</p>

<p align="center">
  <b>Splice-site Strength Estimation using RNA-seq</b>
</p>

<br>
<br>

**version 1.0 (24th April 2024)**

This "version 1" release comes with performance improvements and several quality-of-life updates:

* ~5x speedup (thanks to Pysam and @chenkenbio)

* No more regtools: measure directly from BAM files (via Pysam)

* Duplicate introns on both strands (template switching errors), are now merged by default.

* A new pre-combine command (details below). 

* Replaced some output columns with new info to help downstream processing.


If you are looking for the SpliSER version mentioned in the published article (v0.1.8), you'll find it in the *archive* directory.
<br>

<hr>

SpliSER quantifies the utilisation of splice sites across the genome. Generating a Splice-site Strength Estimate (SSE) for each individual site.

SpliSER has three commands which are used in succession across an analysis: *process*, *combine*, and *output*.
There is an additonal, optional,  *preCombineIntron* command which may save time, especially for big genomes/BAMs.
There is also a second version of *combine* called *combineShallow* with an IO workaround for experiments with many samples.

**Some tips to get started:**

* If you have stranded data, use the --isStranded flag consistently for all of your commands. the default --strandedType "rf" means that that (first-in-pair) read strands are the opposite of the gene that they map to. 

* Always check your outputs for a few sites against the BAM file itself (using an alignment viewer like IGV) to see if the SpliSER output makes sense.  Count how many uses of a splice site you see: does it match the alpha counts that SpliSER gave? How many reads map directly across the splice site without a gap, does it match the beta1 counts?

You can catch a lot of issues with strandedness and annotation this way.

<br>

### Installation
---
*Required python modules: pysam, HTSeq*

Clone the SpliSER repository:

```bash
git clone https://github.com/YOURNAME/SpliSER.git
cd SpliSER
pip install .
```

## process

The *process* command generates a file containing the SSE values (an associated information) of each splice site observed in an RNA-seq alignment. You will need to run this command once on each individual sample; it will then produce a .SpliSER.bed file which can be taken as input for the *combine* command. 

**Which of my samples do I need to process?**<br>
All of them. If you have 3 RNA-seq replicates across 2 conditions, you will need to call *process* once for each of your 6 samples.

**What will I need for this step?**<br>
1. You will need your RNA-seq alignment files in BAM format
2. You may also want an annotation file for your organism (from the same reference genome used for the alignment step) in GFF or GTF format.
<br>

### parameters

The *process* command requires the following two input parameters:

| Required Parameter      | Description |
| ----------- | ----------- |
| -B &nbsp;    \--BAMFile      | The path to an RNA-seq alignment file, in BAM format       |
| -o &nbsp;    \--outputPath  | The path to a directory (including sample prefix) where the resulting .SpliSER.tsv file will be written      |

* The **BAM file** can the the output of any RNA-seq alignment program

* The **outputPath** needs to end with the sample prefix, so if you are processing sample1.bam your output path might read '-o /path/to/directory/sample1'; this will produce a file sample1.SpliSER.tsv in the folder /path/to/directory.
<br>
<br>

There are several optional parameters, which add gene annotations, flag strandedness in your data,  and allow the analysis to be limited to a particular region:

| Optional Parameter      | Description |
| ----------- | ----------- |
| -A &nbsp;    \--annotationFile | The path to a GFF or GTF file matching the genome to which your RNA-seq data was aligned|
| -t &nbsp;   \--annotationType | The type of feature to be extracted from the annotation file. (Default: gene). |
| --isStranded | Include this flag if your RNA-seq data is stranded, prevents opposite-strand reads from contributing to a site's SSE|
| -s &nbsp; \--strandedType | REQUIRED IF USING --isStranded. Strand specificity of RNA library preparation, where \"rf\" is first-strand/RF and \"fr\" is second-strand/FR.|
| -c &nbsp;    \--chromosome | Limit the analysis to one chromosome/scaffold, given by name matching the annotation file *eg.* '-c Chr1'. **required if using -g** |
| -g &nbsp; \--gene | Limit the analysis to one locus, given by name matching the annotation file *eg.* '-g ENSMUSG00000024949'. (If using this parameter you must also specify the --chromosome and --maxIntronSize) |
| -m &nbsp; \--maxIntronSize | **only required if using -g** This is the maximum intron size used in your alignment (If you're unsure, take a maximum intron size for your species *eg.* '-m 6000' for *A.thaliana* or '-m 500000' for *M.musculus*).  |
| -I &nbsp; \--intronFilePath |  Absolute path to a .introns.tsv file output by the preCombineIntrons command.  |

* Add an **annotationFILE** so that you can see which genes your splice sites belong to. SpliSER is annotation-independent by design - when SpliSER reads in an annotation file, all it is really doing is identifying the 'left' and 'right' boundaries of each gene, so it can determine whether a splice-site falls within it or not. This works best with stranded data.

* SpliSER uses the HTSeq package to interpret the GFF/GTF files, thus the first item in the attributes column is taken as the Gene Name.

* The **annotationType** is what will be searched for in the GFF/GTF file. This is 'gene' by default, but these files vary between species - check to make sure you don't need to change it to 'Gene'.

* The **chromosome** parameter allows you to restrict your analysis to a single genomic region. You'll need your input it to match however it appears in the first column of the GFF/GTF annotation file.

* The **gene** parameter allows you to only assess splicing of your favourite locus, this will save you a lot of time compared to the genome-wide approach. Sometimes there are splicing events spanning across annotated gene boundaries, so you'll also need to provide a **maxIntronSize** to ensure that all splice-site strength scores for sites inside the locus are correctly calculated.

<br>

### result

A SpliSER.tsv file containing information about all of the splice sites measured in this sample
| Column  |  Name     | Description |
| ----------- | ----------- | ----------- |
| 1 | Region | The chromosome or scaffold where the splice site occurs (ie Chr1) |
| 2 | Site | The position of the splice site. The position immediately left of the splice junction in the reference genome, for both donors and acceptors. |
| 3 | Strand| The strand on which the site occurred, if known, '+', '-', or '?' |
| 4 | Gene | The gene in which the splice site falls. Otherwise "NA" |
| 5 | SSE | The Splice-site Strength Estimate for this site |
| 6 | alpha_count | How many times did we see this site being used in the sample |
| 7 | beta1_count | How many reads (on the same strand, if applicable) map directly across the splice site without a gap, and don't also show splicing of a known competitor site |
| 8 | beta2_count | How many reads show usage of a competitor site which means that this site probably wasn't used. Because the intron spans over this site, or because the flanking parts of the read cover this site in an ungapped (beta1) manner.  |
| 9 | MultiGeneFlag | Flags 'True' if the SSE for this site uses sites from other genes. Otherwise 'False' (or 'NA' if running for a single gene). |
| 10 | Other | Positions of sites involved in the SSE calculation of this site not covered by the next two columns, in a list format eg. [7228, 8479]  |
| 11 | Partners | Positions of sites which form introns with this site, and their associated counts in this sample, in a dictionary format eg. {7863: 8, 7869: 2 } |
| 12 | Competitors | Positions of sites which form introns with the Partners of this site, in a list format eg. [7710, 7642] |

* A *Partner* site, is any site that this site forms an intron with.

* A *Competitor* site, is any site which also forms an intron with a Partner of this site.

*Other sites, are any sites involved in the calculation of SSE which aren't a Partner or a Competitor. This can happen in a mutually exclusive splicing scenario, where the two sets of splices site never actually form junctions with one another. You'll also see it a lot when a spuriously mapped intron is spanning 10 genes.

<br>
Help for this command can also be viewed in terminal using:
```bash
spliser process -h
```

<br>
<br>
<br>

## combine


The *combine* command takes a set of *processed* files and combines them into a single experiment. *Combine* interpolates all of the splice-site data across your samples. If a splice-site was observed in one sample but not another, SpliSER will get the associated read counts in that second sample for that site. The *combine* command produces a .combined.tsv file which holds a complete set of information for all splice-sites observed across all of your samples. This is usually the longest step.

**Note:**given how long this step can take (it is expensive to lookup individual sites across multiple BAM files one by one). I've added a new command **preCombineIntrons** (described further down) which might speed up the pipeline. 
<br>

**What will I need for this step?**<br>
1. You will need your *processed* files from the previous step.
2. You will need to make your own samples file (detailed below). You can make this in a text editor and save it under whatever name you like. I prefer *SamplesFile.tsv*.
3. You will need to have a quick look at the first few rows of the annotation file (if you used one), and record what the first genomic region is usually 'chr1' or '1'.


### parameters

The *combine* command requires the following three input parameters:

| Required Parameter      | Description |
| ----------- | ----------- |
| -S &nbsp;    \--samplesFile      | The path to a user generated file with three tab-separated columns providing information for each of the processed samples to be combined (see example below)|
| -o &nbsp;    \--outputPath  | The path to a directory (including sample prefix) where the resulting .combined.tsv file will be written|

The **samplesFile** needs to have three columns and no header row. Each row stores the information for one of the samples you wish to combine.



| Column | content |
| --- | --- |
| 1 | **Sample name** - this can be whatever you like, so long as it uniquely identifies each sample |
| 2 | **Path to processed file** - this is the path to the *processed* SpliSER.tsv file for this sample  |
| 3 | **Path to BAM** - this is the path to the BAM file for this sample (SPliSER will need to access it if there was a splice-site observed in another sample, which will then need to be assessed in this sample)|

Here is an example for combining four RNAseq samples:
```
Sample1 /path/to/Sample1.SpliSER.tsv  /path/to/bams/Sample1.bam
Sample2 /path/to/Sample2.SpliSER.tsv  /path/to/bams/Sample2.bam
Sample3 /path/to/Sample3.SpliSER.tsv  /path/to/bams/Sample3.bam
Sample4 /path/to/Sample4.SpliSER.tsv  /path/to/bams/Sample4.bam

```

* The **outputPath** needs to end with the sample prefix, so if you are processing samples for a WT vs mutant analysis, your output path might read '-o /path/to/directory/WTvsMut'; this will produce a file WTvsMut.combined.tsv in the folder /path/to/directory.

| Optional Parameter      | Description |
| ----------- | ----------- |
| --isStranded | Include this flag if your RNA-seq data is stranded, prevents opposite-strand reads from contributing to a site's SSE|
| -s &nbsp; \--strandedType | REQUIRED IF USING --isStranded. Strand specificity of RNA library preparation, where \"rf\" is first-strand/RF and \"fr\" is second-strand/FR.|
| -g &nbsp; \--gene | Limit the analysis to one locus *eg.* '-g ENSMUSG00000024949'(only use this if you also applied the --gene parameter in the previous process step) (Default: All) |

* The -1 / \--firstChrom parameter is redundant as of v0.1.3. The combine command now uses a topological sort to infer the order of genomic regions present in the input files.

<br>

### result

An interleaved combined.tsv file containing information about all of the splice sites measured in all samples for this experiment. 
Very similar to the output of the process command (above) , just with an extra column at the start to denote which sample the measurement came from
| Column  |  Name     | Description |
| ----------- | ----------- | ----------- |
| 1 | Sample | The alias for the sample in which the site was measured, taken from the your samplesFile |
| 2 | Region | The chromosome or scaffold where the splice site occurs (ie Chr1) |
| 3 | Site | The position of the splice site. The position immediately left of the splice junction in the reference genome, for both donors and acceptors. |
| 4 | Strand| The strand on which the site occurred, if known, '+', '-', or '?' |
| 5 | Gene | The gene in which the splice site falls. Otherwise "NA" |
| 6 | SSE | The Splice-site Strength Estimate for this site |
| 7 | alpha_count | How many times did we see this site being used in the sample |
| 8 | beta1_count | How many reads (on the same strand, if applicable) map directly across the splice site without a gap, and don't also show splicing of a known competitor site |
| 9 | beta2_count | How many reads show usage of a competitor site which means that this site probably wasn't used. Because the intron spans over this site, or because the flanking parts of the read cover this site in an ungapped (beta1) manner.  |
| 10 | MultiGeneFlag | Flags 'True' if the SSE for this site uses sites from other genes. Otherwise 'False' (or 'NA' if running for a single gene). |
| 11 | Other | Positions of sites involved in the SSE calculation of this site not covered by the next two columns, in a list format eg. [7228, 8479]  |
| 12 | Partners | Positions of sites which form introns with this site, and their associated counts in this sample, in a dictionary format eg. {7863: 8, 7869: 2 } |
| 13 | Competitors | Positions of sites which form introns with the Partners of this site, in a list format eg. [7710, 7642] |

<br>

Help for this command can also be viewed in terminal using:
```bash
spliser combine -h
```
<br>
<br>
<br>

## preCombineIntrons

The *preCombineIntrons* command generates a file containing the introns seen across all BAM files in an experiment. This means that the **process** command will already know which sites to measure, and the **combine** command won't spend so much time filling in missing sites in each sample. This is more efficient and saves time. If you want to use this, it should be applied *before* running the **process** command. You will need to run this command once per experiment; it will then produce a .introns.tsv file which can be taken as input for the *process* command. 

**What will I need for this step?**<br>
1. A list of all the paths to your BAM files in a comma-separated list (eg. /path/to/control1.bam,/path/to/control2.bam,/path/to/test1.bam,/path/to/test2.bam)
<br>

### parameters

The *process* command requires the following two input parameters:

| Required Parameter      | Description |
| ----------- | ----------- |
| -L &nbsp;    \--ListOfBAMFiles      | A comma-separated (no spaces) list of BAM files in this experiment (which will be later combined)     |
| -o &nbsp;    \--outputPath  | The path to a directory (including sample prefix) where the resulting .introns.tsv file will be written      |

* The **outputPath** needs to end with a filename  prefix, so if your output path reads '-o /path/to/directory/WTvsMUT'; this will produce a file WTvsMUT.introns.tsv in the folder /path/to/directory.

| Optional Parameter      | Description |
| ----------- | ----------- |
| --isStranded | Include this flag if your RNA-seq data is stranded, prevents opposite-strand reads from contributing to a site's SSE|
| -s &nbsp; \--strandedType | REQUIRED IF USING --isStranded. Strand specificity of RNA library preparation, where \"rf\" is first-strand/RF and \"fr\" is second-strand/FR.|
| -c &nbsp;    \--chromosome | Limit the analysis to one chromosome/scaffold, given by name matching the annotation file *eg.* '-c Chr1'. **required if using -g** |
<br>

### result

A 4-column tsv file summarising introns seen in any sample for this experiment.
| Column  | Description |
| ----------- | ----------- |
| 1 | The chromosome or scaffold where the intron occurs (ie Chr1) |
| 2 | The leftmost position of the intron. The position immediately left of the left site of the splice junction, in the reference genome. |
| 3 | The rightmost position of the intron. The position immediately left of the right site of the splice junction, in the reference genome. |
| 4 |  The strand on which the intron occurred, if known, '+', '-', or '?' |

<br>
<br>

### combineShallow

For GWAS-type analyses, you may wish to combine hundreds (if not thousands) of RNA-seq samples into a single experiment. 
The standard combine command will keep each file open as it combines them, but large numbers of files may exceed the file handle limit of your operating system. 
In UNIX systems you can check this limit with:
```
ulimit -n
```
In these circumstances you can use the *combineShallow* command, this works the same as the *combine* but works around the file handle limit. The downside being this is more memory intensive.

To improve performance you can run this command with optional parameters, which will filter out some low-coverage sites and save time. It is up to you to decide which thresholds are appropriate for your downstream analyses.

| Optional Parameter      | Description |
| ----------- | ----------- |
| -m &nbsp;    \--minSamples  | For any given splice site, the minimum number of samples passing the --minReads filter in order for a site to be kept in the analysis - default: 0 |
| -r &nbsp;    \--minReads  | The minimum number of reads giving evidence for a splice site needed for downstream analyses - default: 10 |
| -e &nbsp;    \--minSSE  | The minimum SSE for a splice site to count towards the --minSamples filter - default: 0.00 |

For example: If you don't plan to analyse splice sites which are not supported by 10+ reads in at least 50 samples. You could select "-m 50 -r 10" to skip over these sites during the combineShallow run. If your sample number is in the 1000s, this will considerably speed up the command. 
Further, if you don't care to analyse sites which vary from 0.00 to 0.002, you can select "-e 0.05" to ignore those sites which never get an SSE above that threshold.

<br>

## output

The *output* command takes a *combined* file and outputs the data in one of two formats for downstream analysis.
<br>
**What will I need for this step?**<br>
By this step you should already have everything you need
<br>

### parameters

| Required Parameter      | Description |
| ----------- | ----------- |
| -S &nbsp;    \--samplesFile      | The path to the samples file you used in the previous step, used here to access the sample names|
| -C &nbsp;    \--combinedFile  | The path to the combined.tsv file that you generated in the previous step|
| -t &nbsp;    \--outputType  | 'GWAS' for a SpliSE-QTL analysis, 'DiffSpliSER' for traditional comparison between groups|
| -o &nbsp;    \--outputPath  | The path to a directory (including sample prefix) where the resulting output file will be written|

* The *outputType* command tells SpliSER to format results in one of two ways. 
  * **GWAS** Will produce an individual file for each splice site. This is a two coloumn file containing only the sample name, and the associated SSE value. These files can then be further processed by the user as input for their preferred GWAS pipeline.
  * **DiffSPliSER** Will produce a single file with one row for each splice site, this file can be directly taken as input for the DiffSpliSER.Rmd pipeline script for differential splicing analysis between groups.
  * The samplesFile parameter changed from -s to -S in v0.1.3 onwards
<br>

| Optional Parameter      | Description |
| ----------- | ----------- |
| -r &nbsp;    \--minReads  | The minimum number of reads giving evidence for a splice site in a given sample, below which SpliSER will report NA. The default is 10 reads.|
| -g &nbsp;    \--gene  | Output a the data only from a particular locus *eg.* '-g ENSMUSG00000024949' |
| -m &nbsp;    \--minSamples  | **Only when using --outputType GWAS** - the minimum number of samples passing the read filter for a splice site file to be written |

* **minReads** specifies how few reads you are willing to trust for a SSE value. By default this is set at 10 reads giving evidence for usage or non-usage of the site. We've seen biologically meaningful results lower than this, but the data will tend to get noisier. If you are interested in a particular gene, but don't have sufficient read coverage on your first pass, try 5 reads here.

* **minSamples** This is only for GWAS-type analyses. You may not be interested in  generating an output file for a splice site if there are less than 100 samples passing the read-filter for it, in this case you can specify '-m 100' and those files will not be produced.

<br>  
  
Help for this command can also  be viewed in terminal using:
```bash
spliser output -h
```

## diffSpliSER

We provide an R markdown script for identifying statistically significant differences in splice site usage between samples aligned to the same reference genome. 

### required libraries
- edgeR
- dplyr
- tidyr
- ggplot2

<br>
  
You'll need to edit some of the code to work for your particular sample names/numbers of replicates, just follow the instrcutions in the comments of the R markdown script. 

You'll also need to make tab-separated 'target' file with some additional details about the samples that you are analyzing. We've included a template file with dummy values which you can replace with your own.

Here is the basic layout of the 'target' file:

| Column | content |
| --- | --- |
| 1 | **Sample** - this can be whatever you like, so long as it uniquely identifies each sample. These should be in the same order as the samples appear in your SpliSER output. |
| 2 | **Group** - Each sample should belong to one of two groups, you can call the groups whatver you like.  |
| 3 | **Library_size** - This is the number of mapped reads in your bam files for these samples. This is used for read count normalisation. You can get this by calling *samtools flagstat <your.bam>*|
| 4 | **Description** - A column where you can make any notes you want, you can also leave these entries blank with empty quotation marks ("").


<br>  
<br>  

## Citation

Please cite the original article describing SpliSER:

> Craig I. Dent, Shilpi Singh, Sourav Mukherjee, Shikhar Mishra, Rucha D. Sarwade, Nawar Shamaya, Kok Ping Loo, Paul Harrison, Sridevi Sureshkumar, David Powell, Sureshkumar Balasubramanian.  
> *Quantifying splice-site usage: a simple yet powerful approach to analyze splicing.*  
> **NAR Genomics and Bioinformatics**, Volume 3, Issue 2, June 2021, lqab041.  
> [https://doi.org/10.1093/nargab/lqab041](https://doi.org/10.1093/nargab/lqab041)


## Further Information

Please don't hesitate to get in contact with me at cdent @ mpipz . mpg . de
