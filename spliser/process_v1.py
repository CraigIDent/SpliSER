from spliser.Gene_Site_Iter_Graph_v1 import Gene, Site, Iter, Graph
from spliser.site_ops_v1 import calculateSSE
from spliser.gene_creator_v1 import createGenes
from spliser.bam_parser_v1 import checkBam_pysam, findAlphaCounts_pysam

def outputBedFile(outputPath,chrom_index, site2D_array,qGene):
	outBed = open(outputPath+".SpliSER.tsv","w+")
	#Write the header line
	outBed.write("Region\tSite\tStrand\tGene\tSSE\talpha_count\tbeta1_count\tbeta2_count\tMultiGeneFlag\tOthers\tPartners\tCompetitors\n")


	for c_index, c in enumerate(chrom_index):
		for site in site2D_array[c_index]:
			outBed.write(str(site.getChromosome())+"\t")
			outBed.write(str(site.getPos())+"\t")
			outBed.write(str(site.getStrand())+"\t")
			outBed.write(str(site.getGeneName())+"\t")
			outBed.write("{0:.3f}\t".format(site.getSSE(0)))
			outBed.write(str(site.getAlphaCount(0))+"\t")
			outBed.write(str(site.getBeta1Count(0))+"\t")
			outBed.write(str(site.getBeta2SimpleCount(0))+"\t")
			if qGene !='All':
				outBed.write("NA\t")
			else:
				outBed.write(str(site.getMultiGeneFlag())+"\t")
			outBed.write(str(site.getMutuallyExclusivePos())+"\t")
			#outBed.write(str(site.getBeta2WeightedCount(0))+"\t")
			outBed.write(str(site.getPartnerCount(0))+"\t")
			outBed.write(str(site.getCompetitorPos())+"\n")
	outBed.close()

def buildSiteWhitelist_fromSite2DArray(site2D_array,chrom_index):
	whitelist=set()
	#Read in the samples file and get the SpliSER.tsv file paths
	for c in chrom_index:
		for idx, site in enumerate(site2D_array[chrom_index.index(c)]):
			whitelist.add(str(c)+"_"+str(site.getPos())+"_"+str(site.getStrand()))
	return whitelist


def processSites(inBAM, qChrom, isStranded, strandedType, chrom_index, gene2D_array,site2D_array, sample = 0, numsamples = 1):
	capCounts=False
	print('Processing sample '+ str(int(sample)+1)+' out of '+str(numsamples))
	Whiteset=buildSiteWhitelist_fromSite2DArray(site2D_array,chrom_index)
	for c in chrom_index:
		if qChrom == c or qChrom =="All":
			print("Processing region "+str(c))
			for idx, site in enumerate(site2D_array[chrom_index.index(c)]):
				#Go assign Beta 1 type reads from BAM file
				checkBam_pysam(inBAM, site, sample, isStranded, strandedType, capCounts, Whiteset)
			#Once this is done for all sites, we can calculate SSE
			for idx, site in enumerate(site2D_array[chrom_index.index(c)]):
				#findBeta2Counts(site, numsamples)
				calculateSSE(site)
	return chrom_index, gene2D_array,site2D_array

def flagSiteGenes(chrom_index,site2D_array):
	for c_index, c in enumerate(chrom_index):
		SiteGeneDex={}

		for site in site2D_array[c_index]:
			SiteGeneDex[str(site.getPos())+"_"+site.getStrand()]=site.getGeneName()

		for site in site2D_array[c_index]:
			thisGene = site.getGeneName()
			AssociatedGenes=set()
			AssociatedGenes.add(thisGene)
			#check all partner sites
			for pSite in site.getPartners():
				if str(pSite.getPos())+"_"+site.getStrand() in SiteGeneDex: #In rare cases might not be true where template switching co-occurs with anti-sense introns and they appear to be in competition. In which case silently skip.
					AssociatedGenes.add(pSite.getGeneName())
				else:
					print("warning, skipping ",str(pSite.getPos())+"_"+site.getStrand()," for multiGeneFlagging")
			#check all competitor sites
			for cpos in site.getCompetitorPos():
				if str(cpos)+"_"+site.getStrand() in SiteGeneDex:
					AssociatedGenes.add(SiteGeneDex[str(cpos)+"_"+site.getStrand()])
				else:
					print("warning skipping ",str(cpos)+"_"+site.getStrand()," for multiGeneFlagging") 

			for mpos in site.getMutuallyExclusivePos():
				if str(mpos)+"_"+site.getStrand() in SiteGeneDex:
					AssociatedGenes.add(SiteGeneDex[str(mpos)+"_"+site.getStrand()])
				else:
					print("warning skipping ",str(mpos)+"_"+site.getStrand()," for multiGeneFlagging")

			#check now if we have multiple genes involved.
			AssociatedGenes.discard('NA')
			if len(AssociatedGenes) >1:#if there are multiple genes here
				#print(site.getPos(),AssociatedGenes)
				site.setMultiGeneFlag(True)

def findCompetitorPos(site2D_array, chrom_index):
	for c_index, c in enumerate(chrom_index):
		for site in site2D_array[c_index]:
			sPos = site.getPos()
			for p in site.getPartners():
				for c in p.getPartners():
					cPos = c.getPos()
					if cPos != sPos:
						site.addCompetitorPos(cPos)

def process(inBAM, outputPath, qGene, qChrom, maxIntronSize, annotationFile,aType, isStranded, strandedType, site2D_array=[], intronFilePath = ''):
	capCounts=False
	print('Processing')
	if isStranded:
		print('Stranded Analysis {}'.format(strandedType))
	else:
		print('Unstranded Analysis')

	if annotationFile is not None:
		print('\n\nStep 0: Creating Genes from Annotation...')
		chrom_index, gene2D_array, QUERY_gene, NA_gene = createGenes(annotationFile, aType, qGene)

	#TODO: Add a checker here for qGene not being found (ie mispelled)

	print(('\n\nPreparing Splice Site Arrays'))
	for x in chrom_index:
		site2D_array.append([])

	print('\n\nStep 1A: Finding Splice Sites / Counting Alpha reads...')
	chrom_index, gene2D_array, site2D_array, = findAlphaCounts_pysam(inBAM,qChrom, qGene, int(maxIntronSize), isStranded, strandedType, NA_gene, QUERY_gene=QUERY_gene,chrom_index=chrom_index, gene2D_array=gene2D_array, site2D_array=site2D_array, intronFilePath=intronFilePath, annotationFile=annotationFile) #We apply theqGene filter here, where the splice site objects are made
	print('\n\nStep 2: Finding Beta reads')
	#if capCounts:
	#	print("\t(Capping beta read counts at 2000, or when SSE >0.000)")
	findCompetitorPos(site2D_array,chrom_index) #Add competitors that have been observed in data so far.‚	
	chrom_index, gene2D_array,site2D_array=processSites(inBAM,qChrom, isStranded, strandedType, chrom_index, gene2D_array,site2D_array)

	
	if annotationFile is not None and qGene == 'All':
		print('\n\nFlagging sites involved in multi-gene splicing')
		flagSiteGenes(chrom_index,site2D_array)

	print('\nOutputting .tsv file')
	outputBedFile(outputPath, chrom_index,site2D_array,qGene)
