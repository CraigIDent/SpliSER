<img src="Images/SpliSER.png" width="200">
Splice-site Strength Estimation using RNA-seq

<br>

SpliSER quantifies the utilisation of splice sites across the genome. Generating a Splice-site Strength Estimate (SSE) for each individual site.

SpliSER has three commands which are used in succession across an analysis: *process*, *combine*, and *output*

## process

The *process* command generates a file containing the SSE values (an associated information) of each splice site observed in an RNA-seq alignment. You will need to run this command once on each individual sample; it will then produce a .SpliSER.bed file which can be taken as input for the *combine* command. 

**Which of my samples do I need to process?**<br>
All of them. If you have 3 RNA-seq replicates across 2 conditions, you will need to call *process* once for each of your 6 samples.

**What will I need for this step?**<br>
1. You will need your RNA-seq alignment files in BAM format
2. An accompanying BED file for each BAM. 
3. You may also want an annotation file for your organism (from the same reference genome used for the alignment step) in GFF or GTF format.
<br>

### parameters

The *process* command requires the following three input parameters:

| Required Parameter      | Description |
| ----------- | ----------- |
| -B &nbsp;    \--BAMFile      | The path to an RNA-seq alignment file, in BAM format       |
| -b &nbsp;    \--BEDFile  | The path to a corresponding splice junctions file, in BED12 format        |
| -o &nbsp;    \--outputPath  | The path to a directory (including sample prefix) where the resulting .SpliSER.bed file will be written      |

* The **BAM file** can the the output of any RNA-seq alignment program

* The **BED file** is a TopHat-style summary of splice junctions detected in the alignment. This can be generated from any BAM file using the RegTools (https://regtools.readthedocs.io/en/latest/) 'junction extract' command.

* The **outputPath** needs to end with the sample prefix, so if you are processing sample1.bam your output path might read '-o /path/to/directory/sample1'; this will produce a file sample1.SpliSER.tsv in the folder /path/to/directory. (In version 0.1.1 this was called a .SpliSER.bed file).
<br>
<br>

There are several optional parameters, which add gene annotations to the output and allow the analysis to be limited to a particular region:

| Optional Parameter      | Description |
| ----------- | ----------- |
| -A &nbsp;    \--annotationFile | The path to a GFF or GTF file mathcing the genome to which your RNA-seq data was aligned|
| -t &nbsp;   \--annotationType | The type of feature to be extracted from the annotation file. 'gene' by default. |
| -c &nbsp;    \--chromosome | Limit the analysis to one chromosome/scaffold, given by name matching the annotation file *eg.* '-c Chr1'. **required if using -g** |
| -g &nbsp; \--gene | Limit the analysis to one locus, given by name matching the annotation file *eg.* '-g ENSMUSG00000024949'. (If using this parameter you must also specify the --chromosome and --maxIntronSize) |
| -m &nbsp; \--maxIntronSize | **only required if using -g** This is the maximum intron size used in your alignment (If you're unsure, take a maximum intron size for your species *eg.* '-m 6000' for *A.thaliana* or '-m 500000' for *M.musculus*).  |

* Add an **annotationFILE** so that you can see which genes your splice sites belong to. SpliSER is annotation-independent by design - when SpliSER reads in an annotation file, all it is really doing is identifying the 'left' and 'right' boundaries of each gene, so it can determine whether a splice-site falls within it or not. In version 0.1.1 a custom annotation file format was required (see the attached TAIR10_genes.tsv file as an example).

* SpliSER uses the HTSeq package to interpret the GFF/GTF files, thus the first item in the attributes column is taken as the Gene Name.

* The **annotationType** is what will be searched for in the GFF/GTF file. This is 'gene' by default, but these files vary between species - check to make sure you don't need to change it to 'Gene'.

* The **chromosome** parameter allows you to restrict your analysis to a single genomic region. You'll need to input it to match however it appears in the first coloumn of the GFF/GTF annotation file.

* The **gene** parameter allows you to only assess splicing of your favourite locus, this will save you a lot of time compared to the genome-wide approach. Sometimes there are splicing events spanning across annotated gene boundaries, so you'll also need to provide a **maxIntronSize** to ensure that all splice-site strength scores for sites inside the locus are correctly calculated.

<br>

Help for this command can also be viewed in terminal using:
```
python SpliSER_v0.1.2.py process -h
```
<br>
<br>

## combine


The *combine* command takes a set of *processed* files and combines them into a single experiment. *Combine* interpolates all of the splice-site data across your samples. If a splice-site was observed in one sample but not another, SpliSER will get the associated read counts in that second sample for that site. The *combine* command produces a .combined.tsv file which holds a complete set of information for all splice-sites observed across all of your samples. This is usually the longest step.
<br>

**What will I need for this step?**<br>
1. You will need your *processed* files from the previous step.
2. You will need to make your own samples file (detailed below). You can make this in a text editor and save it under whatever name you like. I prefer *SamplesFile.tsv*.
3. You will need to have a quick look at the first few rows of the annotation file (if you used one), and record what the first genomic region is usually 'chr1' or '1'.


### parameters

The *combine* command requires the following three input parameters:

| Required Parameter      | Description |
| ----------- | ----------- |
| -s &nbsp;    \--samplesFile      | The path to a user generated file with three tab-separated columns providing information for each of the processed samples to be combined (see example below)|
| -1 &nbsp;    \--firstChrom  | The name of the chromosome/scaffold that appears first in the original annotation file (usually '1' or 'chr1')|
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
* **firstChrom** gives the first genomic region (usually a chromosome) that should appear in your data. This is provided in case some samples did not record splicing in the first region, SpliSER still knows where to start combining.

* The **outputPath** needs to end with the sample prefix, so if you are processing samples for a WT vs mutant analysis, your output path might read '-o /path/to/directory/WTvsMut'; this will produce a file WTvsMut.combined.tsv in the folder /path/to/directory.

<br>

Help for this command can also be viewed in terminal using:
```
python SpliSER_v0.1.2.py combine -h
```
<br>

### combineShallow

For GWAS-type analyses, you may wish to combine hundreds (if not thousands) of RNA-seq samples into a single experiment. 
The standard combine command will keep each file open as it combines them, but large numbers of files may exceed the file handle limit of your operating system. 
In UNIX systems you can check this limit with:
```
ulimit -n
```
In these circumstances you can use the *combineShallow* command, this works the same as the *combine* but works around the file handle limit. The downside being this runs considerably slower.

<br>
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
| -s &nbsp;    \--samplesFile      | The path to the samples file you used in the previous step, used here to access the sample names|
| -C &nbsp;    \--combinedFile  | The path to the combined.tsv file that you generated in the previous step|
| -t &nbsp;    \--outputType  | 'GWAS' for a SpliSE-QTL analysis, 'diffSpliSE' for traditional comparison between groups|
| -o &nbsp;    \--outputPath  | The path to a directory (including sample prefix) where the resulting output file will be written|

* The *outputType* command tells SpliSER to format results in one of two ways. 
  * **GWAS** Will produce an individual file for each splice site. This is a two coloumn file containing only the sample name, and the associated SSE value. These files can then be further processed by the user as input for their preferred GWAS pipeline.
  * **DiffSpliSE** Will produce a single file with one row for each splice site, this file can be directly taken as input for the DiffSpliSE.Rmd pipeline script for differential splicing analysis between groups.
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
```
python SpliSER_v0.1.2.py output -h
```

## Further Information

Please don't hesitate to get in contact with me at Craig. Dent @ Monash. edu