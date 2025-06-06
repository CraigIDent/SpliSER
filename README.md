<img src="Images/SpliSER.png" width="200">
Splice-site Strength Estimation using RNA-seq

<br>

version 0.1.8 (16th November 2022)

<br>

**An updated (faster) version is now available on the branch 'speedups' - https://github.com/CraigIDent/SpliSER/tree/speedups**

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
**Note**: This step seems to be a common place to trip up, so I suggest checking the files at each step. I've included some trouble shooting tips at the bottom of this section


* The **outputPath** needs to end with the sample prefix, so if you are processing sample1.bam your output path might read '-o /path/to/directory/sample1'; this will produce a file sample1.SpliSER.tsv in the folder /path/to/directory. (In version 0.1.1 this was called a .SpliSER.bed file).
<br>
<br>

There are several optional parameters, which add gene annotations to the output and allow the analysis to be limited to a particular region:

| Optional Parameter      | Description |
| ----------- | ----------- |
| -A &nbsp;    \--annotationFile | The path to a GFF or GTF file mathcing the genome to which your RNA-seq data was aligned|
| -t &nbsp;   \--annotationType | The type of feature to be extracted from the annotation file. (Default: gene). |
| --isStranded | Include this flag if your RNA-seq data is stranded, prevents opposite-strand reads from contributing to a site's SSE|
| -s &nbsp; \--strandedType | REQUIRED IF USING --isStranded. Strand specificity of RNA library preparation, where \"rf\" is first-strand/RF and \"fr\" is second-strand/FR.|
| --beta2Cryptic | Calculate SSE of sites taking into account the weighted utilisation of competing splice sites as indirect evidence of site non-utilisation (Legacy).|
| -c &nbsp;    \--chromosome | Limit the analysis to one chromosome/scaffold, given by name matching the annotation file *eg.* '-c Chr1'. **required if using -g** |
| -g &nbsp; \--gene | Limit the analysis to one locus, given by name matching the annotation file *eg.* '-g ENSMUSG00000024949'. (If using this parameter you must also specify the --chromosome and --maxIntronSize) |
| -m &nbsp; \--maxIntronSize | **only required if using -g** This is the maximum intron size used in your alignment (If you're unsure, take a maximum intron size for your species *eg.* '-m 6000' for *A.thaliana* or '-m 500000' for *M.musculus*).  |

* Add an **annotationFILE** so that you can see which genes your splice sites belong to. SpliSER is annotation-independent by design - when SpliSER reads in an annotation file, all it is really doing is identifying the 'left' and 'right' boundaries of each gene, so it can determine whether a splice-site falls within it or not. In version 0.1.1 a custom annotation file format was required (see the attached TAIR10_genes.tsv file as an example).

* SpliSER uses the HTSeq package to interpret the GFF/GTF files, thus the first item in the attributes column is taken as the Gene Name.

* The **annotationType** is what will be searched for in the GFF/GTF file. This is 'gene' by default, but these files vary between species - check to make sure you don't need to change it to 'Gene'.

* The **chromosome** parameter allows you to restrict your analysis to a single genomic region. You'll need to input it to match however it appears in the first coloumn of the GFF/GTF annotation file.

* The **gene** parameter allows you to only assess splicing of your favourite locus, this will save you a lot of time compared to the genome-wide approach. Sometimes there are splicing events spanning across annotated gene boundaries, so you'll also need to provide a **maxIntronSize** to ensure that all splice-site strength scores for sites inside the locus are correctly calculated.

<br>
**Troubleshooting guide for splice juction .bed files**

* After running Regtools, open up the .bed files in a text editor and have a look at the 6th column. If you chose unstranded, these should all be question marks (?). If you chose a stranded analysis, then you should see some pluses (+) and minuses (-). 

* These pluses and minuses should mostly occur in groups of consecutive pluses or minuses. Splice sites from the same genes should be in the same orientation as eachother. If your data is unstranded, but you incorrectly tell Regtools it is stranded, then you will still get pluses and minuses but you will see an apparently random alternation of + and - in the 6th column of the bed file (because Regtools will base the orientation on the first read it encounters at each junction).

* Regtools sometimes doesn't like stranded bam file and will only correctly annotate one strand, and the other strand will all be question marks. If about 50% of your junctions are '+', and the other 50% are '?', then you may need to manually overwrite the '?' symbols as '-' (conversely if you see 50% '-' and  '?'). 

* The rf, fr options can be confusing, but there is a rule of thumb to help: Open up your Bam in IGV, colour the reads by first-in-pair, and find a splice junction. Find the corresponding junction in the bed file. If the orientation of the junction +/- is the same orientation as the reads +/-, then you are all set. From there, if the Reads are in the same orientation as the gene annotation, then use SpliSER's 'fr' option, if they are opposite orientation then use SpliSER's 'rf option.

*Issues with the above steps can look like the following:
    *Half of the splice sites didn't get assigned to a gene
    *There are no Beta1 read counts for any splice site

Any issues please don't hesitate to contact me at cdent @ mpipz . mpg . de

<br>
Help for this command can also be viewed in terminal using:
```
python SpliSER_v0.1.8.py process -h
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
| --beta2Cryptic | Calculate SSE of sites taking into account the weighted utilisation of competing splice sites as indirect evidence of site non-utilisation (Legacy).|

* The -1 / \--firstChrom parameter is redundant as of v0.1.3. The combine command now uses a topological sort to infer the order of genomic regions present in the input files.

<br>

Help for this command can also be viewed in terminal using:
```
python SpliSER_v0.1.8.py combine -h
```
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
```
python SpliSER_v0.1.8.py output -h
```

## diffSpliSER

We provide an R markdown script for identifying statistically significant differences in splice site usage between samples aligned to the same reference genome. 

You'll need to edit some of the code to work for your particular sample names/numbers of replicates, just follow the instrcutions in the comments of the R markdown script. 

You'll also need to make tab-separated 'target' file with some additional details about the samples that you are analyzing. We've included a template file with dummy values which you can replace with your own.

Here is the basic layout of the 'target' file:

| Column | content |
| --- | --- |
| 1 | **Sample** - this can be whatever you like, so long as it uniquely identifies each sample. These should be in the same order as the samples appear in your SpliSER output. |
| 2 | **Group** - Each sample should belong to one of two groups, you can call the groups whatver you like.  |
| 3 | **Library_size** - This is the number of mapped reads in your bam files for these samples. This is used for read count normalisation. You can get this by calling *samtools flagstat <your.bam>*|
| 4 | **Description** - A column where you can make any notes you want, you can also leave these entries blank with empty quotation marks ("").

## Further Information

Please don't hesitate to get in contact with me at cdent @ mpipz . mpg . de
