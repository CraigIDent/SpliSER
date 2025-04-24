from ast import literal_eval
from src.Gene_Site_Iter_Graph_v1 import Gene, Site, Iter, Graph
from src.site_ops_v1 import makeSingleSpliceSite, calculateSSE
from src.bam_parser_v1 import checkBam_pysam

def deduceChromOrder(bedPaths):

	#First iterate through all files and make a list of all regions
	#make a graph of all relationships in the list
	allChroms = []
	beforeList = []
	afterList = []
	
	for b in bedPaths: # for each processed file
		before = "-1" # set an arbitary initial region
		for idx, line in enumerate(open(b,'r')):
			if idx >0: # skip headers
				chrom = line.split("\t")[0] # get the genomic region
				if chrom != before: # if this is a new region
					beforeList.append(before)
					afterList.append(chrom)
					before = chrom
					if chrom not in allChroms:
						allChroms.insert(0,chrom)
	if not allChroms:
		print("No genomic regions found - EXITING")
		sys.exit()

	allChroms.insert(0,"-1")

	g = Graph(allChroms)
	for idx, b in enumerate(beforeList):
		g.addEdge(b, afterList[idx])

	chromsInOrder = g.topologicalSort()[1:] #droppping off the arbitrary first region again
	
	return chromsInOrder

def outputCombinedLines(outTSV, site, gene,qGene,allTitles):
	for idx, t in enumerate(allTitles): # for each sample
		outTSV.write(str(t)+"\t")
		outTSV.write(str(site.getChromosome())+"\t")
		outTSV.write(str(site.getPos())+"\t")
		outTSV.write(str(site.getStrand())+"\t")
		outTSV.write(gene+"\t")
		outTSV.write("{0:.3f}\t".format(site.getSSE(idx)))
		outTSV.write(str(site.getAlphaCount(idx))+"\t")
		outTSV.write(str(site.getBeta1Count(idx))+"\t")
		outTSV.write(str(site.getBeta2SimpleCount(idx))+"\t")
		if qGene !='All':
			outTSV.write("NA\t")
		else:
			outTSV.write(str(site.getMultiGeneFlag())+"\t")
		outTSV.write(str(site.getMutuallyExclusivePos())+"\t")
		outTSV.write(str(site.getPartnerCount(idx))+"\t")
		outTSV.write(str(site.getCompetitorPos())+"\n")



def combine(samplesFile, outputPath,qGene, isStranded, strandedType):
	print('Combining samples...')
	outTSV = open(outputPath+".combined.tsv", 'w+')
	outTSV.write("Sample\tRegion\tSite\tStrand\tGene\tSSE\talpha_count\tbeta1_count\tbeta2_count\tMultiGeneFlag\tOthers\tPartners\tCompetitors\n")
	#Process the input paths file
	samples = 0
	allTitles=[]
	bedPaths = [] # stores the absolute path to each SpliSER.bed file
	BAMPaths = [] # stores the absolute path to each orginal bam file
	for line in open(samplesFile,'r'):
		values = line.split("\t")
		if len(values) ==3:
			samples +=1
			allTitles.append(values[0]) # record the sample moniker
			bedPaths.append(values[1]) # record the bed file paths
			BAMPaths.append(values[2].rstrip()) # record BAM file paths
		elif len(values) >= 0:
			print(str(allTitles), str(bedPaths), str(BAMPaths))
			raise Exception('Samples File contains lines that do not have exactly 3 tab-separated columns')


	print('Establishing order of genomic regions...')
	chromsInOrder = deduceChromOrder(bedPaths)
	print("order of genomic regions deduced: {}".format(chromsInOrder))

	#Iterate bed files
	print('Iterating through files in parallel, to interleave lines and fill gaps.')
	iters = []
	for file in bedPaths: # make a bed file line generator for each bed file
		iters.append((open(file, 'r')))

	for iter in iters: #skip header line
		temp = next(iter, None)


	#initialise data structures needed for iterating along bed files
	currentVals = [['']*11]*samples #array to hold the current line for each bed file - 11 elements
	currentChromPos=0
	maxChromPos = len(chromsInOrder)-1
	currentChrom = chromsInOrder[currentChromPos] # The chromosome currently being assessed
	lowestPos = -1
	lowestPosStrand="?"
	assocGene = "" #use this to store the gene associated with the lowestPos
	chroms = ['']*samples
	iterGo = [True]*samples #initialise iterGo for each sample
	iterDone = [False]*samples


	#load up intial values
	count = 0
	filledCount = 0
	print('updated currentChrom to {}'.format(currentChrom))
	print("Combining data for site# "+str(count)+"...")
	while not(len(set(iterDone)) <=1 and iterDone[0] ==True): #Stop if all files are done iterating
		#refresh values
		sSite = None
		assocGene = ""
		partners = []
		competitors = []
		others = []
		filledGap = False

		for idx, iter in enumerate(iters):
			#	print('processing Site: ',count,'file: ',idx)
			if iterDone[idx] ==False: #unless we have run out of lines for this file, get the values for this file
				#Get the next value for appropriate iters
				if iterGo[idx] == True:
					nextLine = next(iter,None)
					if nextLine is not None:
						currentVals[idx] = nextLine.rstrip().split("\t")
						#Update chroms seen this round
						chroms[idx] = currentVals[idx][0]
						iterGo[idx] = False #pause iterator until we know this value was used
					else: #if the iterator is done
						iterDone[idx] = True #recognise it is done
						chroms[idx] = None #set Chrom to None so we can ignore this file in our currentChrom equation

				#check position(if on curent chrom) - record if it's the lowest we've seen so far
				if chroms[idx] == currentChrom:
					pos = int(currentVals[idx][1])
					strand = currentVals[idx][2]
					if pos < lowestPos or lowestPos == -1 or (isStranded==True and pos==lowestPos and strand=="+"): #check if this is a new lowest position
						lowestPos = int(pos)
						lowestPosStrand = strand
						assocGene = currentVals[idx][3]

		#if we didn't exhaust all iterators at the start of this loop
		if not(len(set(iterDone)) <=1 and iterDone[0] ==True):
			#Create a Splice Site for lowestPos
			sSite = makeSingleSpliceSite(currentChrom, lowestPos, samples, '', isStranded)

			chromsCheck = [i for i in chroms if i == currentChrom]
			if len(chromsCheck)==0: # if we have exhausted the current region
				currentChromPos = currentChromPos + 1
				if currentChromPos > maxChromPos:
					currentChrom = None
				else:
					currentChrom = chromsInOrder[currentChromPos]
					count = 0
					print('updated currentChrom to {}'.format(currentChrom))
					print("Combining data for site# "+str(count)+"...")
			else: #if we are still adding new values
				#add values into SpliceSite Object
				for idx, vals in enumerate(currentVals):
						if vals[0] == currentChrom and int(vals[1]) == lowestPos and iterDone[idx] ==False and (isStranded==False or vals[2] == lowestPosStrand): #if this sample has values for the splice site (and if it's a stranded analysis we are looking at the same strand)
							iterGo[idx] = True #if we take values from this file, then we want to get a new line next time
							#add details in for splice sites
							sSite.setStrand(str(vals[2])) # set the strand for this site
							#sSite.setSSE(float(vals[4]),idx) # set SSE for this site
							sSite.addAlphaCount(int(vals[5]), idx) # add alpha Counts
							sSite.addBeta1Count(int(vals[6]), idx) # add beta1 Counts
							sSite.addBeta2SimpleCount(int(vals[7]), idx) # add beta2Simple Counts
							sSite.setMultiGeneFlag(str(vals[8]))
							#get other (mutually exclusive sites)
							mPosList = literal_eval(str(vals[9]))
							for m in mPosList:
								others.append(m)
								sSite.addMutuallyExclusivePos(m)
							#recalculate SSE, since we might have used the isbeta2Cryptic flag differently in this step
							try:
								calculateSSE(sSite)
							except:
								print("Could not recalculate SSE.")
							#read partner counts as a dictionary and update the splice
							pCounts = literal_eval(str(vals[10]))
							for key, val in pCounts.items():
								partners.append(key)
								sSite.addPartnerCount(key, val ,idx)
							#read competitor positions as a list and add to the
							cPosList = literal_eval(str(vals[11]))
							for c in cPosList:
								competitors.append(c)
								sSite.addCompetitorPos(c)

						else: #if this sample doesn't have values for the spliceSite
							#find beta1 and beta2Simple counts for the site, using partners and competitors
							if qGene == 'All' or qGene == assocGene:
								filledGap = True
								checkBam_pysam(BAMPaths[idx], sSite, idx, isStranded, strandedType)
								sSite.setSSE(0.000,idx)
				#output lines for this splice site
				if qGene == 'All' or qGene == assocGene:
					outputCombinedLines(outTSV,sSite, assocGene,qGene,allTitles)
			#reset the lowestPos COUNTER
			lowestPos = -1
		if filledGap:
			filledCount += 1

		count += 1
		if count%10000 == 0:
			print("Combining data for site# "+str(count))
	#loop back
	print('Filled in Beta read counts for {} Sites not detected in some samples'.format(filledCount))


def combineShallow(samplesFile, outputPath, qGene, isStranded, minSamples, minReads, minSSE, strandedType):
	print('Combining samples...')
	outTSV = open(outputPath+".combined.tsv", 'w+')
	outTSV.write("Sample\tRegion\tSite\tStrand\tGene\tSSE\talpha_count\tbeta1_count\tbeta2_count\tMultiGeneFlag\tOthers\tPartners\tCompetitors\n")
	#Process the input paths file
	samples = 0
	allTitles=[]
	bedPaths = [] # stores the absolute path to each SpliSER.bed file
	BAMPaths = [] # stores the absolute path to each orginal bam file
	print('Reading in Samples File')
	for line in open(samplesFile,'r'):
		values = line.split("\t")
		if len(values) ==3:
			samples +=1
			allTitles.append(values[0]) # record the sample moniker
			bedPaths.append(values[1]) # record the bed file paths
			BAMPaths.append(values[2].rstrip()) # record BAM file paths

	print('Reading SpliSER processed files into memory')
	#instead of an iterator, make a list of a list of lines
	bedLines = []
	line_counter = []
	max_lines = []
	for file in bedPaths:
		line_list = []
		maxi = 0
		for idx,line in enumerate(open(file,"r")):
			if idx >0: # Skip headers

				if qGene != "All": #if we are prefiltering by qGene
					g = line.split("\t")[3]
					if g == qGene:
						maxi += 1
						line_list.append(line)
				else: # otherwise add as normal
					maxi += 1
					line_list.append(line)
		bedLines.append(line_list)
		line_counter.append(0)
		max_lines.append(maxi)

	print('Establishing order of genomic regions...')
	chromsInOrder = deduceChromOrder(bedPaths)
	print("order of genomic regions deduced: {}".format(chromsInOrder))

	#initialise data structures needed for iterating along bed files
	currentVals = [['']*11]*samples #array to hold the current line for each bed file - 11 columns
	currentChromPos=0
	maxChromPos = len(chromsInOrder)-1
	currentChrom = chromsInOrder[currentChromPos] # The chromosome currently being assessed
	lowestPos = -1
	lowestPosStrand="?"
	assocGene = "" #use this to store the gene associated with the lowestPos
	chroms = ['']*samples
	iterGo = [True]*samples #initialise iterGo for each sample
	iterDone = [False]*samples

	#load up intial values
	count = 0
	filledCount = 0
	print('updated currentChrom to {}'.format(currentChrom))
	print("Combining data for site# "+str(count)+"...")
	while not(len(set(iterDone)) <=1 and iterDone[0] ==True): #Stop if all files are done iterating
		#refresh values
		sSite = None
		assocGene = ""
		partners = []
		competitors = []
		filledGap = False
		posCounter = 0

		for idx, bed in enumerate(bedLines):
			#print('processing Site: ',count,'file: ',idx)
			if iterDone[idx] ==False: #unless we have run out of lines for this file, get the values for this file
				#Get the next linr from this bed file
				if iterGo[idx] == True:
					if line_counter[idx] < max_lines[idx]:
						nextLine = bed[line_counter[idx]]
						line_counter[idx] += 1
					else:
						nextLine = None

					if nextLine is not None:
						#get the values of the next line
						if qGene == 'All': #If we are assessing all sites - proceed as normal
							currentVals[idx] = nextLine.rstrip().split("\t")
							#Update chroms seen this round
							chroms[idx] = currentVals[idx][0]
							#pause iterator until we know this value was used
							iterGo[idx] = False
						else: #If we are looking for sites from a particular gene
							qGeneFound = False
							#check the gene associated with this site
							while qGeneFound == False and iterDone[idx] == False: #If the site we are looking at isn't from the qGene,  and the iterator isn't yet exhausted
								currentVals[idx] = nextLine.rstrip().split("\t") # get values for this site
								#print('{} {} {}'.format(currentVals[idx][0],currentVals[idx][1],currentVals[idx][3]))
								if currentVals[idx][3] == qGene: #If this is a site from qGene - proceed as normal
									qGeneFound = True
									chroms[idx] = currentVals[idx][0]
									iterGo[idx] = False
								elif iterGo[idx] == True: #else, if the iterator is still going
									if line_counter[idx] < max_lines[idx]:
										nextLine = bed[line_counter[idx]]
										line_counter[idx] += 1
									else:
										nextLine = None
									if nextLine is None: #if that is empty
										iterDone[idx] = True #recognise it is done
										chroms[idx] = None #set Chrom to None so we can ignore this file in our currentChrom equation
									#Otherwise the while loop will start again until we find something from qGene..


					else: #if the iterator is done
						iterDone[idx] = True #recognise it is done
						chroms[idx] = None #set Chrom to None so we can ignore this file in our currentChrom equation

				#check position(if on curent chrom) - record if it's the lowest we've seen so far

			if chroms[idx] == currentChrom:
				pos = int(currentVals[idx][1])
				strand = currentVals[idx][2]
				if pos < lowestPos or lowestPos == -1 or (isStranded==True and pos==lowestPos and strand != lowestPosStrand and strand=="+"):

					lowestPos = int(pos)
					lowestPosStrand = strand
					assocGene = currentVals[idx][3]
					#If our first encounter of a new lowestPos has more than 'minReads' reads, record that we've seen it
					reads = int(currentVals[idx][5])+ int(currentVals[idx][6])+int(currentVals[idx][7])
					sse = float(currentVals[idx][4])
					if reads >= minReads and sse >= minSSE:
						posCounter =  1
					else:
						posCounter = 0
				#elseagain, count up how many times it is seen with more than 'minReads' reads
				elif pos == lowestPos:
					reads = int(currentVals[idx][5])+ int(currentVals[idx][6])+int(currentVals[idx][7])
					sse = float(currentVals[idx][4]) #fixed 15Mar2022
					if reads >= minReads and sse >= minSSE:
						#If there are at least 10 reads, record that we've seen it
						posCounter = posCounter + 1

		#print("{} {} {} ".format(str(currentVals[0][1]),str(currentVals[1][1]),str(currentVals[2][1])))
		#print("{} {} {} ".format(str(iterGo[0]),str(iterGo[1]),str(iterGo[2])))
		#if we didn't exhaust all iterators at the start of this loop
		if not(len(set(iterDone)) <=1 and iterDone[0] ==True):



				#Create a Splice Site for lowestPos
				sSite = makeSingleSpliceSite(currentChrom, lowestPos, samples, '',isStranded)

				chromsCheck = [i for i in chroms if i == currentChrom]
				if len(chromsCheck)==0: # if we have exhausted the current region
					currentChromPos = currentChromPos + 1
					if currentChromPos > maxChromPos:
						currentChrom = None
					else:
						currentChrom = chromsInOrder[currentChromPos]
						count = 0
						print('updated currentChrom to {}'.format(currentChrom))
						print("Combining data for site# "+str(count)+"...")
				else: #if we are still adding new values
					#add values into SpliceSite Object
					if posCounter >= minSamples:  #If we have seen this site enough times for it to be worth processing.
						for idx, vals in enumerate(currentVals):
							if vals[0] == currentChrom and int(vals[1]) == lowestPos and iterDone[idx] ==False and (isStranded==False or vals[2] == lowestPosStrand): #if this sample has values for the splice site (and if it's a stranded analysis we are looking at the same strand)
								iterGo[idx] = True #if we take values from this file, then we want to get a new line next time
								#add details in for splice sites
								sSite.setStrand(str(vals[2])) # set the strand for this site
								#sSite.setSSE(float(vals[4]),idx) # set SSE for this site
								sSite.addAlphaCount(int(vals[5]), idx) # add alpha Counts
								sSite.addBeta1Count(int(vals[6]), idx) # add beta1 Counts
								sSite.addBeta2SimpleCount(int(vals[7]), idx)  # add beta2Simple Counts
								sSite.setMultiGeneFlag(str(vals[8])) 
								#get other (mutually exclusive sites)
								mPosList = literal_eval(str(vals[9]))
								for m in mPosList:
									others.append(m)
									sSite.addMutuallyExclusivePos(m)
								#recalculate SSE, redundant
								try:
									calculateSSE(sSite)
								except:
									print("Could not recalculate SSE.")
								#read partner counts as a dictionary and update the splice
								pCounts = literal_eval(str(vals[10]))
								#print(pCounts, type(pCounts))
								for key, val in pCounts.items():
									partners.append(key)
									sSite.addPartnerCount(key, val ,idx)
								#read competitor positions as a list and add to the
								cPosList = literal_eval(str(vals[11]))
								for c in cPosList:
									competitors.append(c)
									sSite.addCompetitorPos(c)

							else: #if this sample doesn't have values for the spliceSite
								#find beta1 and beta2Simple counts for the site, using partners and competitors
								if qGene == 'All' or qGene == assocGene:
									filledGap = True
									checkBam_pysam(BAMPaths[idx], sSite, idx, isStranded, strandedType)
									sSite.setSSE(0.000,idx)
								#output lines for this splice site
						if qGene == 'All' or qGene == assocGene:
							outputCombinedLines(outTSV,sSite, assocGene,qGene,allTitles)	
					else:	# If there were not enough samples recording the splice site to pass minSamples
						print("Skipped site {} for insufficient evidence, only {} samples with Site using minimum reads".format(lowestPos,posCounter))
						#print(lowestPos)
						#print("{} {} {} ".format(str(iterGo[0]),str(iterGo[1]),str(iterGo[2])))
						for idx, vals in enumerate(currentVals):
							if iterDone[idx] == False and int(currentVals[idx][1]) == lowestPos:
								iterGo[idx] = True #we want to get a new line next time
				#reset the lowestPos COUNTER
				lowestPos = -1

		if filledGap:
			filledCount += 1

		count += 1
		if count%10000 == 0:
			print("Combining data for site# "+str(count)+"...")
	#loop back
	print('Filled in Beta read counts for {} Sites not detected in some samples'.format(filledCount))

