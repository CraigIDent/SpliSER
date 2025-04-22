"""
SpliSER- Splice-site Strength Estimation from RNA-seq
"""
#Version 1.0 - 18th April 2025
version = "v1.0"
import sys
import timeit
import time
import subprocess
import argparse
import HTSeq
import re
from operator import truediv
from Gene_Site_Iter_Graph_v0_1_8 import Gene, Site, Iter, Graph
import numpy
from operator import add, truediv, mul, sub
import bisect
from ast import literal_eval

from site_ops_v1_dev import findAlphaCounts_pysam, findBeta2Counts, calculateSSE
from gene_creator_v1_dev import createGenes
from bam_parser_v1_dev import checkBam
digit_pattern = re.compile(r'\D')

#chrom_index = []
#gene2D_array = []
#site2D_array = []
#Sites_2Darray = []
Genes = []
allChroms = []
allCounts = [] #stores the total number of evidencing reads for a given splice site, for a given sample.
allNames =[]
allSEs = [] #stores all the effiency values
allTitles = [] #stores all the names of files
geneCounter = 0


digit_pattern = re.compile(r'\D') # pattern, non-digit
char_pattern = re.compile(r'\d') # pattern, digit
sSite = None
#QUERY_gene = None

def outputBedFile(outputPath,isbeta2Cryptic,chrom_index, site2D_array):
	outBed = open(outputPath+".SpliSER.tsv","w+")
	#Write the header line
	outBed.write("Region\tSite\tStrand\tGene\tSSE\talpha_count\tbeta1_count\tbeta2Simple_count\tbeta2Cryptic_count\tbeta2Cryptic_weighted\tPartners\tCompetitors\n")
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
			if isbeta2Cryptic:
				outBed.write(str(site.getBeta2CrypticCount(0))+"\t")
				outBed.write("{0:.5f}\t".format(site.getBeta2WeightedCount(0)))
			else:
				outBed.write("NA\t")
				outBed.write("NA\t")
			#outBed.write(str(site.getBeta2WeightedCount(0))+"\t")
			outBed.write(str(site.getPartnerCount(0))+"\t")
			outBed.write(str(site.getCompetitorPos())+"\n")
	outBed.close()

#Check if all elements of the list are equivalent.
def checkEqual2(list):
   return len(set(list)) <= 1

def makeSingleSpliceSite(iChrom, iPos, numSamples, iStrand, isStranded):
	return(Site(
				chromosome = iChrom,
				pos = iPos,
				samples = numSamples,
				strand = iStrand,
				source = '',
				isStranded=isStranded
				)
			)

def processSites(inBAM, qChrom, isStranded, strandedType,isbeta2Cryptic, chrom_index, gene2D_array,site2D_array, sample = 0, numsamples = 1):
	print('Processing sample '+ str(int(sample)+1)+' out of '+str(numsamples))
	for c in chrom_index:
		print("Processing region "+str(c))
		if qChrom == c or qChrom =="All":
			for idx, site in enumerate(site2D_array[chrom_index.index(c)]):
				#Go assign Beta 1 type reads from BAM file
				checkBam(inBAM, site, sample, isStranded, strandedType)
			#Once this is done for all sites, we can calculate SSE
			for idx, site in enumerate(site2D_array[chrom_index.index(c)]):
				findBeta2Counts(site, numsamples)
				calculateSSE(site,isbeta2Cryptic)
	return chrom_index, gene2D_array,site2D_array


def process(inBAM, inBed, outputPath, qGene, qChrom, maxIntronSize, annotationFile,aType, isStranded, strandedType, isbeta2Cryptic, site2D_array=[]):
	print('Processing')
	if isStranded:
		print('Stranded Analysis {}'.format(strandedType))
	else:
		print('Unstranded Analysis')

	if annotationFile is not None:
		print('\n\nStep 0: Creating Genes from Annotation...')
		chrom_index, gene2D_array, QUERY_gene, NA_gene = createGenes(annotationFile, aType, qGene)

	#TODO: Add a checker here for qGene not being found

	print(('\n\nPreparing Splice Site Arrays'))
	for x in chrom_index:
		site2D_array.append([])

	print('\n\nStep 1A: Finding Splice Sites / Counting Alpha reads...')
	chrom_index, gene2D_array, site2D_array = findAlphaCounts_pysam(inBAM,qChrom, qGene, int(maxIntronSize), isStranded, strandedType, NA_gene, QUERY_gene=QUERY_gene,chrom_index=chrom_index, gene2D_array=gene2D_array, site2D_array=site2D_array) #We apply theqGene filter here, where the splice site objects are made
	print('\n\nStep 2: Finding Beta reads')
	chrom_index, gene2D_array,site2D_array=processSites(inBAM,qChrom, isStranded, strandedType, isbeta2Cryptic,chrom_index, gene2D_array,site2D_array)

	print('\nOutputting .tsv file')
	outputBedFile(outputPath,isbeta2Cryptic,chrom_index,site2D_array)

def outputCombinedLines(outTSV, site, gene,isbeta2Cryptic):
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
		if isbeta2Cryptic:
			outTSV.write(str(site.getBeta2CrypticCount(idx))+"\t")
			outTSV.write(str(site.getBeta2WeightedCount(idx))+"\t")
		else:
			outTSV.write("NA\t")
			outTSV.write("NA\t")
		outTSV.write(str(site.getPartnerCount(idx))+"\t")
		outTSV.write(str(site.getCompetitorPos())+"\n")

def combine(samplesFile, outputPath,qGene, isStranded, strandedType, isbeta2Cryptic):
	print('Combining samples...')
	outTSV = open(outputPath+".combined.tsv", 'w+')
	outTSV.write("Sample\tRegion\tSite\tStrand\tGene\tSSE\talpha_count\tbeta1_count\tbeta2Simple_count\tbeta2Cryptic_count\tbeta2_weighted\tPartners\tCompetitors\n")
	#Process the input paths file
	samples = 0
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

	#First iterate through all files and make a list of all regions
	#make a graph of all relationships in the list
	allchroms = []
	beforeList = []
	afterList = []

	print('Establishing order of genomic regions.')
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
							if vals[8]!= "NA":
								sSite.addBeta2CrypticCount(int(vals[8]), idx) # add beta2Cryptic Counts
								sSite.addBeta2Weighted(float(vals[9]), idx) # add beta2WeightedCounts
							#else do nothing, we don't need these values

							#recalculate SSE, since we might have used the isbeta2Cryptic flag differently in this step
							try:
								calculateSSE(sSite, isbeta2Cryptic)
							except:
								print("Could not recalculate SSE. You might be trying to use --beta2Cryptic flag without using it in the process step")
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
								checkBam(BAMPaths[idx], sSite, idx, isStranded, strandedType)
								sSite.setSSE(0.000,idx)
				#output lines for this splice site
				if qGene == 'All' or qGene == assocGene:
					outputCombinedLines(outTSV,sSite, assocGene,isbeta2Cryptic)
			#reset the lowestPos COUNTER
			lowestPos = -1
		if filledGap:
			filledCount += 1

		count += 1
		if count%10000 == 0:
			print("Combining data for site# "+str(count))
	#loop back
	print('Filled in Beta read counts for {} Sites not detected in some samples'.format(filledCount))


def combineShallow(samplesFile, outputPath, qGene, isStranded, minSamples, minReads, minSSE, strandedType, isbeta2Cryptic):
	print('Combining samples...')
	outTSV = open(outputPath+".combined.tsv", 'w+')
	outTSV.write("Sample\tRegion\tSite\tStrand\tGene\tSSE\talpha_count\tbeta1_count\tbeta2Simple_count\tbeta2Cryptic_count\tbeta2_weighted\tPartners\tCompetitors\n")
	#Process the input paths file
	samples = 0
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

	#First iterate through all files and make a list of all regions
	#make a graph of all relationships in the list
	allchroms = []
	beforeList = []
	afterList = []

	print('Establishing order of genomic regions.')
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
	#Build the graph of genomic regions
	g = Graph(allChroms)
	for idx, b in enumerate(beforeList):
		g.addEdge(b, afterList[idx])
	#Make a sorted list of genomic regions
	chromsInOrder = g.topologicalSort()[1:] #droppping off the arbitrary first region again
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
								sSite.addBeta2SimpleCount(int(vals[7]), idx) # add beta2Simple Counts
								if vals[8]!= "NA":
									sSite.addBeta2CrypticCount(int(vals[8]), idx) # add beta2Cryptic Counts
									sSite.addBeta2Weighted(float(vals[9]), idx) # add beta2WeightedCounts
								#Else do nothing, we don't need these values unless they've already been specified

								#recalculate SSE, since we might have used the isbeta2Cryptic flag differently in this step
								try:
									calculateSSE(sSite, isbeta2Cryptic)
								except:
									print("Could not recalculate SSE. You might be trying to use --beta2Cryptic flag without using it in the process step")

								#read partner counts as a dictionary and update the splice site
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
									checkBam(BAMPaths[idx], sSite, idx, isStranded, strandedType)
									sSite.setSSE(0.000,idx)
								#output lines for this splice site
						if qGene == 'All' or qGene == assocGene:
							outputCombinedLines(outTSV,sSite, assocGene,isbeta2Cryptic)
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

def DiffSpliSER_output(samplesFile,combinedFile, outputPath, minReads, qGene):

	outDiff = open(outputPath+str(qGene)+".DiffSpliSER.tsv", "w+")
	correct = True
	fileFinished = False
	#record the samples we're assessing
	samples = 0
	for line in open(samplesFile,'r'):
		values = line.split("\t")
		if len(values) ==3:
			samples +=1
			allTitles.append(values[0]) # record the sample moniker

	#print the headerLine
	outDiff.write("Region\tSite\tStrand\tGene")
	for t in allTitles:
		outDiff.write("\t"+str(t)+"_alpha")
		outDiff.write("\t"+str(t)+"_beta")
		outDiff.write("\t"+str(t)+"_SSE")
	outDiff.write("\n")

	#make an iterator for the combined file
	comboFile = open(combinedFile, 'r')
	#skip the header
	#comboFile.next()
	nextLine = next(comboFile,None)

	while not fileFinished: #go until the file is exhausted
		currentVals = []
		for idx, t in enumerate(allTitles):
			#Skip the header
			nextLine = next(comboFile,None)
			if nextLine is not None:
				currentVals.append(nextLine.rstrip().split("\t")) #store an array of values for each sample
				if currentVals[idx] == t: #check that each title matches up with the line we've just extracted
					correct = True
			else:
				fileFinished = True
		#if we didnt hit the end of the file, output some valuuuues
		if not fileFinished:
			if not correct:
				print("Samples aren\'t match up, please check there are no missing lines in your combined file, or missing lines in your samples file" )
				break

			outDiff.write(str(currentVals[0][1])+"\t"+str(currentVals[0][2])+"\t"+str(currentVals[0][3])+"\t"+str(currentVals[0][4])) #write the region, splice site, and gene
			for idx, t in enumerate(allTitles):
				t_alpha = int(currentVals[idx][6])# get alpha Values
				t_beta = float(currentVals[idx][7])+float(currentVals[idx][8]) # add beta1 and beta2Simple
				if currentVals[idx][10] != "NA":
					t_WeightedCrypticBeta = float(currentVals[idx][10])
				else:
					t_WeightedCrypticBeta =0

				t_SSE = float(currentVals[idx][5])
				if t_alpha+t_beta >= minReads: # if this sample passes the minimum read count for this site

					t_beta = t_beta+t_WeightedCrypticBeta
					outDiff.write("\t"+str(t_alpha)+"\t"+"{0:.2f}".format(t_beta)+"\t"+"{0:.2f}".format(t_SSE))
				else:
					outDiff.write("\tNA\tNA\tNA")
			outDiff.write("\n")

def GWAS_output(samplesFile,combinedFile, outputPath, minReads, qGene, minSamples):
	print(qGene)
	fileFinished = False
	correct = True
	#record the samples we're assessing
	samples = 0
	for line in open(samplesFile,'r'):
		samples +=1
		values = line.split("\t")
		allTitles.append(values[0]) # record the sample moniker

	#open the datafile
	comboFile = open(combinedFile, 'r')
	#skip the header
	#comboFile.next()
	nextLine = next(comboFile,None)
	while not fileFinished: #go until the file is exhausted
		currentVals = []
		for idx, t in enumerate(allTitles):
			nextLine = next(comboFile,None)
			if nextLine is not None:
				currentVals.append(nextLine.rstrip().split("\t")) #store an array of values for each sample
				if currentVals[idx] == t: #check that each title matches up with the line wes've just extracted
					correct = True
			else:
				fileFinished = True
		#if we didnt hit the end of the file, output some valuuuues
		if not fileFinished:
			if not correct:
				print("Samples aren\'t match up, please check there are no missing lines in your combined file, or missing lines in your samples file" )
				break

			currentGene = str(currentVals[0][4])
			currentSite = str(currentVals[0][2])
			bufferString = ''
			samplesPassing = 0
			if qGene == currentGene or qGene == 'All':
				filtered = open(str(outputPath+currentGene+"_"+currentSite+"_filtered.log"),'w+') # otherwise write into a filter log file specific for this gene
				for idx, t in enumerate(allTitles):
					t_alpha = int(currentVals[idx][6])# get alpha Values
					t_beta = float(currentVals[idx][7])+float(currentVals[idx][8]) # add beta1 and beta2Simple values
					if t_alpha+t_beta >= minReads: # if this sample passes the minimum read count for this site
						samplesPassing +=1
						bufferString = bufferString+str(currentVals[idx][0])+"\t"+str(currentVals[idx][5])+"\n" #store sample name and SSE in buffer
					else:
						filtered.write(str(t)+" did not pass minReads for "+currentSite+"\n")
				if samplesPassing >= minSamples:
					out = open(outputPath+currentGene+"_"+currentSite+".tsv","w+") #open a file named after the splice site
					#out.write("Sample\tSSE\n")#Write in a header
					out.write(bufferString) # write stored info
					out.close()
				else:
					filtered.write("Site: "+currentSite+" minSamples not met - venting buffer\n")
					filtered.write(bufferString+"\n")

				filtered.close()

def output(outputType, samplesFile,combinedFile, outputPath, minReads, qGene, minSamples):
	if outputType == 'DiffSpliSER':
		DiffSpliSER_output(samplesFile,combinedFile, outputPath, minReads, qGene)
	if outputType == 'GWAS':
		GWAS_output(samplesFile,combinedFile, outputPath, minReads, qGene, minSamples)

#Main Method
if __name__ == "__main__":
	print("\nSpliSER "+version+" SKB LAB\n")
	# this won't be run when imported
	start = timeit.default_timer()

	#Parse arguments form the command line
	parser = argparse.ArgumentParser(description="SpiSER - Splice Site Strength Estimates from RNA-seq")
	subparsers = parser.add_subparsers(dest="command")
	#Parser for arguments when user calls command 'process'
	parser_process = subparsers.add_parser('process')
	parser_process.add_argument('-B', '--BAMFile', dest='inBAM', help="The mapped RNA-seq file in BAM format", required=True)
	parser_process.add_argument('-b', '--bedFile', dest='inBed', help="The Tophat-style splice junction bed file", required=True)
	parser_process.add_argument('-o', '--outputPath', dest='outputPath', help="Absolute path, including file prefix where SpliSER will write the .SpliSER.tsv file", required=True)
	parser_process.add_argument('-A', '--annotationFile', dest='annotationFile', help="optional: gff3 or gtf file matching the reference genome used for alignment", required=False)
	parser_process.add_argument('-t', '--annotationType', dest='aType', nargs='?', default='gene', type=str, help="optional: The feture to be extracted from the annotation file - default: gene", required=False)
	parser_process.add_argument('-c', '--chromosome', dest='qChrom', nargs='?', default='All', type=str,help="optional: limit SpliSER to one chromosome/scaffold -  default: All", required=False)
	parser_process.add_argument('-g', '--gene', dest='qGene', nargs='?', default='All', type=str, help="optional:Limit SpliSER to splice sites falling in a single locus (requires --chromosome and --annotationFile and --maxIntronSize)", required=False)
	parser_process.add_argument('-m', '--maxIntronSize', dest='maxIntronSize', nargs='?', default=0, type=int, required=False, help="optional: required if passing the --gene parameter, the max intron size used in aligning the bam file")
	#parser_process.add_argument('-S', '--isStranded', dest='isStranded', nargs='?', default=False, type=bool, required=False, help="optional: True/False Perform a strand-aware analysis - default: False")
	parser_process.add_argument('--isStranded', dest='isStranded', default=False, action='store_true')
	parser_process.add_argument('-s', '--strandedType', dest='strandedType', nargs='?', type=str, required=False, help="optional: Strand specificity of RNA library preparation, where \"rf\" is first-strand/RF and \"fr\" is second-strand/FR - default : fr")
	parser_process.add_argument('--beta2Cryptic', dest='isbeta2Cryptic', default=False, action='store_true', help="optional: Calculate SSE of sites taking into account the weighted utilisation of competing splice sites as indirect evidence of site non-utilisation (Legacy).")
	#Parser for arguments when user calls command 'combine'
	parser_combine = subparsers.add_parser('combine')
	parser_combine.add_argument('-S', '--samplesFile', dest='samplesFile', required=True, help="A three-column .tsv file, each line containing a sample name, the absolute path to a processed .SpliSER.tsv file input, and the absolute path to the original bam file")
	parser_combine.add_argument('-o', '--outputPath', dest='outputPath', required=True, help="Absolute path to a folder, including file_prefix where SpliSER will write the output .tsv file")
	parser_combine.add_argument('-g', '--gene', dest='qGene', nargs='?', default='All', type=str, help="optional:Limit SpliSER to splice sites falling in a single locus - default All", required=False)
	parser_combine.add_argument('--isStranded', dest='isStranded', default=False, action='store_true')
	parser_combine.add_argument('-s', '--strandedType', dest='strandedType', nargs='?', default="fr", type=str, required=False, help="optional: Strand specificity of RNA library preparation, where \"rf\" is first-strand/RF and \"fr\" is second-strand/FR - default : fr")
	parser_combine.add_argument('--beta2Cryptic', dest='isbeta2Cryptic', default=False, action='store_true', help="optional: Calculate SSE of sites taking into account the weighted utilisation of competing splice sites as indirect evidence of site non-utilisation (Legacy).")

	parser_combineShallow = subparsers.add_parser('combineShallow')
	parser_combineShallow.add_argument('-S', '--samplesFile', dest='samplesFile', required=True, help="A three-column .tsv file, each line containing a sample name, the absolute path to a processed .SpliSER.tsv file input, and the absolute path to the original bam file")
	parser_combineShallow.add_argument('-g', '--gene', dest='qGene', nargs='?', default='All', type=str, help="optional:Limit SpliSER to splice sites falling in a single locus - default All", required=False)
	parser_combineShallow.add_argument('-o', '--outputPath', dest='outputPath', required=True, help="Absolute path to a folder, including file_prefix where SpliSER will write the output .tsv file")
	parser_combineShallow.add_argument('--isStranded', dest='isStranded', default=False, action='store_true')
	parser_combineShallow.add_argument('-m', '--minSamples', dest='minSamples', required=False, nargs='?', default=0, type=int, help="For optional filtering: For any given splice site, the minimum number of samples passing the --minReads filter in order for a site to be kept in the analysis (Raising this value increases speed of combineShallow step)")
	parser_combineShallow.add_argument('-r', '--minReads', dest='minReads',required=False, nargs='?', default=10, type=int, help="For optional filtering: The minimum number of reads giving evidence for a splice site needed for downstream analyses - default: 10")
	parser_combineShallow.add_argument('-e','--minSSE', dest='minSSE',required=False, nargs='?', default=0.00, type=float, help="For optional filtering: The minimum SSE of a site for a given sample, for it to be considered in the --minSamples filter - default: 0.00")
	parser_combineShallow.add_argument('-s', '--strandedType', dest='strandedType', nargs='?', type=str, required=False, help="optional: Strand specificity of RNA library preparation, where \"rf\" is first-strand/RF and \"fr\" is second-strand/FR - default : fr")
	parser_combineShallow.add_argument('--beta2Cryptic', dest='isbeta2Cryptic', default=False, action='store_true', help="optional: Calculate SSE of sites taking into account the weighted utilisation of competing splice sites as indirect evidence of site non-utilisation (Legacy).")

	parser_output = subparsers.add_parser('output')
	parser_output.add_argument('-S', '--samplesFile', dest='samplesFile', required=True, help="the three-column .tsv file you used to combine the samples in the previous step")
	parser_output.add_argument('-C', '--combinedFile', dest='combinedFile', required=True, help="a SpliSER .combined.tsv file containing the splice site information for each sample")
	parser_output.add_argument('-t', '--outputType', dest='outputType', required = True, help="Type of output file: -t DiffSpliSER will output a file ready for differential splicing analysis. -t GWAS will output an SSE phenotype file for each Splice Site (writes to outputPath folder, ignoring file prefix)")
	parser_output.add_argument('-o', '--outputPath', dest='outputPath', required=True, help="Absolute path to an output folder(ie ending in a slash), iff using -t DiffSpliSER also provide a file_prefix where SpliSER will write the output .tsv file")
	parser_output.add_argument('-r', '--minReads', dest='minReads',required=False, nargs='?', default=10, type=int, help="The minimum number of reads giving evidence for a splice site in a given sample, below which SpliSER will report NA - default: 10")
	parser_output.add_argument('-g', '--gene', dest='qGene', required=False, nargs='?', default='All', type=str, help="optional:Limit SpliSER to splice sites falling in a single locus")
	parser_output.add_argument('-m', '--minSamples', dest='minSamples', required=False, nargs='?', default=50, type=int, help="optional: when using --outputType GWAS: the minimum number of samples passing the read filter for a splice site file to be written")
	#Parse arguments
	kwargs = vars(parser.parse_args())
	command = kwargs.pop('command')

	#in process command, --g requires --A
	if command == 'process' and (kwargs.get('qGene') != "All") and ((kwargs.get('annotationFile') is None) or (kwargs.get('maxIntronSize') is None)):
		print(kwargs.get('qGene'))
		print(kwargs.get('annotationFile'))
		parser.error("--gene requires --annotationFile and --maxIntronSize")
	elif (command == 'process' or command == 'combine' or command =='combineShallow') and (kwargs.get('isStranded') is True) and (kwargs.get('strandedType') is None):
		parser.error("--isStranded requires parameter --strandedType/-s as fr or rf")
	else: #otherwise
		# call appropriate functions
		globals()[command](**kwargs)

	stop = timeit.default_timer()
	print("Total runtime (s): \t"+ str(stop - start))
#EOF
