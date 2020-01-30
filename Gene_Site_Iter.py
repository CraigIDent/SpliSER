#Gene.py
#defines a 'Gene' object. Containing annotation information and all splice sites detected through RNA-seq analysis.
from operator import add, truediv, mul, sub
import bisect
import numpy as np
class Gene:

	def __init__(self, chromosome, name, leftPos, rightPos, readNums, samples, strand, source):
		self.chromosome = str(chromosome)
		self.name = str(name)
		if readNums == None:
			self.reads = [None]*int(samples)
		else:
			self.reads = readNums[:samples]

		self.leftPos = int(leftPos)
		self.rightPos = int(rightPos)
		self.source = str(source)
		self.sites = []
		self.strand = str(strand)

    #Operators for genes
    	def __lt__(self, other):
		return self.leftPos < other.leftPos
	def __le__(self, other):
		return self.leftPos <= other.leftPos
	def __eq__(self, other):
		return self.leftPos == other.leftPos
	def __ne__(self, other):
		return self.leftPos != other.leftPos
	def __gt__(self, other):
		return self.leftPos > other.leftPos
	def __ge__(self, other):
		return self.leftPos >= other.leftPos

	def getChromosome(self):
		return self.chromosome

	def getName(self):
		return self.name

	def getReads(self):
			return self.reads

	def getLeftPos(self):
		return self.leftPos

	def getRightPos(self):
		return self.rightPos

	def getSource(self):
		return self.source

	def getStrand(self):
		return self.strand

	def addSite(self, l):
		#if l not in self.sites:
		self.sites.append(l)

	def getSites(self):
		return self.sites

	def getReadNum(self, sample):
		return self.reads[sample]

	def getReadNums(self):
		return self.reads

	def addReadNum(self, newReads, sample):
		self.reads[sample] = newReads

	def popSite(self):
		return self.sites.pop()


#	def addJunction(self, j):
#		uniqueJunction = True # binary true/false
#		for junction in self.junctions:
#			if j.getLeftPos() == junction.getLeftPos() and j.getRightPos() == junction.getRightPos():
#				uniqueJunction = False
#		if uniqueJunction == True:
#			self.junctions.append(j)

#	def addCountToJunction(self, l, r, sample, strand, count):
#		for j in self.junctions:
#			if int(j.getLeftPos()) == int(l) and int(j.getRightPos()) == int(r) and str(j.getStrand()) == str(strand):
#				j.addCount(count, sample)


#END OF GENE CLASS
#--------------------------------------------------------------------------------------------------------------

class Site:
	#Slot attributes so __dict__ is not required for each
	__slots__ = ('chromosome','samples','pos','source','alphaCounts', 'beta1Counts','beta2Counts','beta2Weighted','Partners','PartnerCounts','SEs','strand','beta2weights','Gene')

	def __init__(self, chromosome, pos, samples, strand, source):
		self.chromosome = chromosome
		self.samples = samples
		self.pos = int(pos)
		self.source = source
		self.alphaCounts = [0]*int(samples)
		#self.alphaCounts = np.empty(int(samples), np.int32)
		self.beta1Counts = [0]*int(samples)
		#self.beta1Counts = np.empty(int(samples), np.int32)
		self.beta2SoftCounts = [0]*int(samples)
		self.beta2HardCounts = [0]*int(samples) #counts where we observe non-usage of the target site, along with usage of a partner and a competitor.
		#self.beta2Counts = np.empty(int(samples), np.int32)
		self.beta2Weighted = [0.000]*int(samples)
		self.Partners = [] # array of Site objects
		self.CompetitorPos = [] #array of Site Positions
		self.PartnerCounts = {} #Dictionary, where key is Partner Position and value is an array of Counts across samples.
		self.PartnerBeta2DoubleCounts = {}
		self.SSEs= [0.000]*int(samples)
		self.strand = strand
		self.beta2weights = {} # Dictionary where key is Partner Position and value is a weight between 0 and 1
		#self.id = str(chromosome)+"_"+str(pos)
		self.Gene = None
#Operators for Sites
	def __lt__(self, other):
		return self.pos < other.pos
	def __le__(self, other):
		return self.pos <= other.pos
	def __eq__(self, other):
		return self.pos == other.pos
	def __ne__(self, other):
		return self.pos != other.pos
	def __gt__(self, other):
		return self.pos > other.pos
	def __ge__(self, other):
		return self.pos >= other.pos

#SETTERS
	def setGene(self, g):
		self.Gene = g

	def setSSEs(self, values):
		self.SSEs = values

	def setSSE(self, value, sample):
		self.SSEs[sample] = value

	def updateBeta2Weighted(self, values):
		self.beta2Weighted = values

#BETTERS
	def addAlphaCount(self, count, sample):
		self.alphaCounts[sample] += count

	def addBeta1Count(self, count, sample):
		self.beta1Counts[sample] += count

	def addBeta2SoftCount(self, count, sample):
		self.beta2SoftCounts[sample] += count

	def addBeta2HardCount(self, count, sample):
		self.beta2HardCounts[sample] += count

	def addBeta2SoftCounts(self, values):
		self.beta2SoftCounts = map(add, self.beta2SoftCounts, values)

	def addBeta2Weighted(self, count, sample):
		self.beta2Weighted[sample] += count

	def addPartnerCount(self, sitePos, count, sample):
		if sitePos not in self.PartnerCounts:
			self.PartnerCounts[sitePos] = [0]*int(self.samples)
		self.PartnerCounts[sitePos][sample] = self.PartnerCounts[sitePos][sample] + count

	def addPartnerBeta2DoubleCount(self, sitePos, count, sample):
		if sitePos not in self.PartnerBeta2DoubleCounts:
			self.PartnerBeta2DoubleCounts[sitePos] = [0]*int(self.samples)
		self.PartnerBeta2DoubleCounts[sitePos][sample] = self.PartnerBeta2DoubleCounts[sitePos][sample] + count

	def addPartner(self, site):
		if site not in self.Partners:
			self.Partners.append(site)

	def addCompetitorPos(self, pos):
		if pos not in self.CompetitorPos:
			bisect.insort(self.CompetitorPos,int(pos))
#GETTERS

	def getAlphaCount(self, sample):
		return self.alphaCounts[sample]

	def getAlphaCounts(self):
		return map(int, self.alphaCounts)

	def getBeta1Count(self, sample):
		return self.beta1Counts[sample]

	def getBeta1Counts(self):
		return map(int, self.beta1Counts)

	def getBeta2SoftCount(self, sample):
		return self.beta2SoftCounts[sample]

	def getBeta2SoftCounts(self):
		return self.beta2SoftCounts

	def getBeta2HardCount(self, sample):
		return self.beta2HardCounts[sample]

	def getBeta2HardCounts(self):
		return self.beta2HardCounts

	def getBeta2WeightedCount(self, sample):
		return self.beta2Weighted[sample]

	def getBeta2WeightedCounts(self):
		return self.beta2Weighted

#	def getAllCount(self, sample):
#		return [self.alphaCounts[sample], self.beta1Counts[sample], self.beta2Weighted[sample]]

	def getGeneName(self):
		return self.Gene.getName()

	def getPos(self):
		return self.pos

	def getChromosome(self):
		return self.chromosome

	def getStrand(self):
		return self.strand

	def getPartners(self):
		return self.Partners

	def getPartners(self):
		return self.Partners

	def getPartnerCounts(self):
		return self.PartnerCounts

#	def getPartnerBeta2HardCounts(self):
#		return self.PartnerBeta2HardCounts

	def getPartnerBeta2DoubleCounts(self):
		return self.PartnerBeta2DoubleCounts

	def getPartnerCount(self, sample):
		val = {}
		for p in self.PartnerCounts:
			val[p] = self.PartnerCounts[p][sample]
		return val

	def getCompetitorPos(self):
		return self.CompetitorPos

	def getSSE(self, sample):
		return self.SSEs[sample]

	def getSSEs(self):
		return self.SSEs

	def getSource(self):
		return self.source

#END OF JUNCTION CLASS
class Iter:

	def __init__(self, File):
		self.File = str(File)
		self.pos = 0
		done = False

	def getNext(self):
			self.pos+=1
			try:
				return [x for i, x in enumerate(open(self.File,'r')) if i == (self.pos-1)][0]
			except IndexError:
				return None
