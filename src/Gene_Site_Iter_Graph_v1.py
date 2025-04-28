
from operator import add, truediv, mul, sub
import bisect
from collections import defaultdict


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
	__slots__ = ('chromosome','samples','pos','source','alphaCounts', 'beta1Counts','beta2SimpleCounts','Partners','CompetitorPos','MutuallyExclusivePos','PartnerCounts','PartnerBeta2DoubleCounts','SSEs','strand','Gene','isStranded','multiGeneFlag')

	def __init__(self, chromosome, pos, samples, strand, source, isStranded):
		self.chromosome = chromosome
		self.samples = samples
		self.pos = int(pos)
		self.source = source
		self.alphaCounts = [0]*int(samples)
		self.beta1Counts = [0]*int(samples)
		self.beta2SimpleCounts = [0]*int(samples) #counts where we observe non-usage of the target site, along with usage of a partner and a competitor.
		self.Partners = [] # array of Site objects
		self.CompetitorPos = [] #array of Site Positions
		self.MutuallyExclusivePos = [] #array of Site Positions
		self.PartnerCounts = {} #Dictionary, where key is Partner Position and value is an array of Counts across samples.
		self.PartnerBeta2DoubleCounts = {}
		self.SSEs= [0.000]*int(samples)
		self.strand = strand
		self.Gene = None
		self.isStranded = isStranded #marks whether or not this site is being used in a stranded analysis
		self.multiGeneFlag = False
#Operators for Sites
	def __lt__(self, other):
		if self.isStranded: #if a stranded analysis
			if self.pos==other.pos: #if sites have equal position
				if self.strand == other.strand: #if strands are same, then equal
					return False
				elif self.strand =="+" and other.strand =="-": ## + < -
					return True
				elif self.strand =="-" and other.strand =="+": ## - not < +
					return False
			else:
				return self.pos < other.pos
					
		else: #if not a stranded analysis
			return self.pos < other.pos

	def __le__(self, other):
		if self.isStranded:
			if self.pos==other.pos: #if sites have equal position
				if self.strand == other.strand: #if strands are same, then equal
					return True
				elif self.strand =="+" and other.strand =="-": ## + <= -
					return True
				elif self.strand =="-" and other.strand =="+": ## - not <= +
					return False
			else:
				return self.pos <= other.pos	
		else:
			return self.pos <= other.pos
	def __eq__(self, other):
		if self.isStranded:
			if self.pos==other.pos: #if sites have equal position
				if self.strand == other.strand: #if strands are same, then equal
					return True
				elif self.strand =="+" and other.strand =="-": ## + not == -
					return False
				elif self.strand =="-" and other.strand =="+": ## - not == +
					return False	
			else:
				return self.pos == other.pos					
		else:
			return self.pos == other.pos
	def __ne__(self, other):
		if self.isStranded:
			if self.pos==other.pos: #if sites have equal position
				if self.strand == other.strand: #if strands are same, then equal
					return False
				elif self.strand =="+" and other.strand =="-": ## + not == -
					return True
				elif self.strand =="-" and other.strand =="+": ## - not == +
					return True	
			else:
				return self.pos != other.pos			
		else:
			return self.pos != other.pos
	def __gt__(self, other):
		if self.isStranded:
			if self.pos==other.pos: #if sites have equal position
				if self.strand == other.strand: #if strands are same, then equal
					return False
				elif self.strand =="+" and other.strand =="-": ## + not >= -
					return False
				elif self.strand =="-" and other.strand =="+": ## - <= +
					return True
			else:
				return self.pos > other.pos
		else:
			return self.pos > other.pos
	def __ge__(self, other):
		if self.isStranded:
			if self.pos==other.pos: #if sites have equal position
				if self.strand == other.strand: #if strands are same, then equal
					return True
				elif self.strand =="+" and other.strand =="-": ## + not >= -
					return False
				elif self.strand =="-" and other.strand =="+": ## - <= +
					return True
			else:
				return self.pos >= other.pos
		else:
			return self.pos >= other.pos

#SETTERS
	def setStrand(self, s ):
		self.strand = s

	def setGene(self, g):
		self.Gene = g

	def setSSEs(self, values):
		self.SSEs = values

	def setSSE(self, value, sample):
		self.SSEs[sample] = value

	def setMultiGeneFlag(self, value):
		self.multiGeneFlag = value

#BETTERS
	def addAlphaCount(self, count, sample):
		self.alphaCounts[sample] += count

	def addBeta1Count(self, count, sample):
		self.beta1Counts[sample] += count

	def addBeta2SimpleCount(self, count, sample):
		self.beta2SimpleCounts[sample] += count

	def addBeta2SimpleCounts(self, values):
		self.beta2SimpleCounts =[x + y for x, y in zip(self.beta2SimpleCounts, values)]

	def addPartnerCount(self, sitePos, count, sample):
		if sitePos not in self.PartnerCounts:
			self.PartnerCounts[sitePos] = [0]*int(self.samples)
		self.PartnerCounts[sitePos][sample] = self.PartnerCounts[sitePos][sample] + count

	def addPartner(self, site):
		if site not in self.Partners:
			self.Partners.append(site)

	def addCompetitorPos(self, pos):
		if pos not in self.CompetitorPos:
			bisect.insort(self.CompetitorPos,int(pos))

	def addMutuallyExclusivePos(self, pos):
		if pos not in self.MutuallyExclusivePos:
			bisect.insort(self.MutuallyExclusivePos,int(pos))
#GETTERS

	def getAlphaCount(self, sample):
		return self.alphaCounts[sample]

	def getAlphaCounts(self):
		return [int(a) for a in self.alphaCounts]

	def getBeta1Count(self, sample):
		return self.beta1Counts[sample]

	def getBeta1Counts(self):
		return [int(b) for b in self.beta1Counts]

	def getBeta2SimpleCount(self, sample):
		return self.beta2SimpleCounts[sample]

	def getBeta2SimpleCounts(self):
		return self.beta2SimpleCounts

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

	def getPartnerCount(self, sample):
		val = {}
		for p in self.PartnerCounts:
			val[p] = self.PartnerCounts[p][sample]
		return val

	def getCompetitorPos(self):
		return self.CompetitorPos

	def getMutuallyExclusivePos(self):
		return self.MutuallyExclusivePos

	def getSSE(self, sample):
		return self.SSEs[sample]

	def getSSEs(self):
		return self.SSEs

	def getSource(self):
		return self.source

	def getMultiGeneFlag(self):
		return self.multiGeneFlag

#END OF SITE CLASS

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

#GRAPH CLASS ADAPTED FROM CODE CONTRIBUTED TO 'GEEKS FOR GEEKS' BY Neelam Yadav, for topolgical sorting
#I've altered the code to handle and store strings and indices in parallel
class Graph:
	#Class to represent a graph
	def __init__(self,chromList):
		self.graph = defaultdict(list) #dictionary containing adjacency List
		self.R = len(chromList) #No. of genomic regions
		self.chromList = chromList

    # function to add an edge to graph
	def addEdge(self,before,after):
		if after not in self.graph[before]:
			self.graph[before].append(after)
		else:
			pass

    # A recursive function used by topologicalSort
	def topologicalSortUtil(self,region,visited,stack):
		idx = self.chromList.index(region)
		visited[idx] = True

        # Recur for all the vertices adjacent to this vertex
		for i in self.graph[region]:
			if visited[self.chromList.index(i)] == False:
				self.topologicalSortUtil(i,visited,stack)

        # Push current vertex to stack which stores result
		stack.insert(0,region)
    # The function to do Topological Sort. It uses recursive
    # topologicalSortUtil()
	def topologicalSort(self):
        # Mark all the vertices as not visited
		visited = [False]*self.R
		stack =[]
        # Call the recursive helper function to store Topological
        # Sort starting from all vertices one by one
		for idx, i in enumerate(self.chromList):
			if visited[idx] == False:
				self.topologicalSortUtil(i,visited,stack)

		return(stack)

#print(Site.__module__)