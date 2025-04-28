# bam_parser.py

import sys
import pysam
import bisect
from collections import Counter
from spliser.Gene_Site_Iter_Graph_v1 import Site
from spliser.binary_searches_v1 import binary_gene_search, binary_site_search

def check_strand(strandedType, SAMflag, siteStrand):
    '''
    Takes the flag element of a SAM and determines if that read originated from the same strand as the splice site
    '''
    if strandedType == "fr":
        if (SAMflag & 64) or not (SAMflag & 1): #If first in pair, or not paired read
            if (SAMflag & 16):                      #If read reverse strand
                readStrand = "-"
            else:
                readStrand = "+"
        else:                                   #else (if paired and second in pair)
            if (SAMflag & 16):                      #If read is reverse strand
                readStrand = "+"    
            else:
                readStrand = "-"

    if strandedType == "rf":                  #Flip everything if it is "rf" library prep
        if (SAMflag & 64) or not (SAMflag & 1):
            if (SAMflag & 16):
                readStrand = "+"
            else:
                readStrand = "-"
        else:
            if (SAMflag & 16):
                readStrand = "-"
            else:
                readStrand = "+"

    return (readStrand == siteStrand)


def collapse_duplicate_introns(introns_plus: Counter, introns_minus: Counter):
    # Create a copy to avoid modifying inputs directly
    introns_plus = introns_plus.copy()
    introns_minus = introns_minus.copy()

    # Check for overlapping introns by (start, end) tuple
    shared_introns = set(introns_plus.keys()) & set(introns_minus.keys())

    for intron in shared_introns:
        plus_count = introns_plus[intron]
        minus_count = introns_minus[intron]

        if plus_count >= minus_count:
            introns_plus[intron] += minus_count
            del introns_minus[intron]
        else:
            introns_minus[intron] += plus_count
            del introns_plus[intron]

    return introns_plus, introns_minus


def preCombineIntrons(BAMPathList,outputPath,qChrom,isStranded,strandedType):
    dontCollapse=False
    BAMPaths = BAMPathList.split(",") # stores the absolute path to each orginal bam file

    print("Finding all introns across:", len(BAMPaths),"samples")
    print("Standed analysis: ", isStranded)
    intronSet = set()
    #Find all introns for this experiment
    for bdx,bamFile in enumerate(BAMPaths):
        bam = pysam.AlignmentFile(bamFile, "rb")

        chroms = bam.references

        for chrom in chroms:
            if qChrom == chrom or qChrom == "All":
                print("sample",bdx+1, "chromosome",chrom, "processing observed introns..")
                #If this is a stranded analysis, then make a generator for the forward and reverse strand reads (accounting for paired-end reads)
                #and package them up in tuples with the appropriate strand designation
                if isStranded:
                    bamgen_plus = (read for read in bam.fetch(chrom) if check_strand(strandedType,read.flag,"+") and not read.is_unmapped and not read.is_secondary and not read.is_supplementary)
                    bamgen_minus = (read for read in bam.fetch(chrom) if check_strand(strandedType,read.flag,"-") and not read.is_unmapped and not read.is_secondary and not read.is_supplementary)
                    introns_plus=bam.find_introns(bamgen_plus)
                    introns_minus=bam.find_introns(bamgen_minus)
                    if not dontCollapse:
                        introns_plus, introns_minus = collapse_duplicate_introns(introns_plus,introns_minus)
                    Intron_info =[(introns_plus,"+"),(introns_minus,"-")]
                else:
                    bamgen=(read for read in bam.fetch(chrom) if not read.is_unmapped and not read.is_secondary and not read.is_supplementary)
                    introns=bam.find_introns(b)
                    Intron_info = [(introns,"?")]
            
                for introns,strand in Intron_info:
                    for i in introns:
                        #print(i,strand)
                        intronSet.add((chrom, i[0],i[1],strand))
    print(len(intronSet))
    intronOut = open(outputPath+".introns.tsv","w+")
    for line in sorted(intronSet):
        #print(line)
        intronOut.write("\t".join([str(line[0]),str(line[1]),str(line[2]),str(line[3])]))
        intronOut.write("\n")



def findAlphaCounts_pysam(bamFile, qChrom, qGene, maxIntronSize, isStranded,strandedType, NA_gene, sample=0, numsamples=1, QUERY_gene=None, chrom_index=None, gene2D_array=None, site2D_array=None, intronFilePath=''):
    """
    Find junctions in bam file, record the number of reads showing usage (alpha1 reads) of each splice site forming the junction.

    Parameters
    ----------
    bAMFile : str
        The absolute path to a .BAM file for a given alignment.
    qChrom : str
        Chromosome to restrict analysis to, or 'All' for full genome.
    qGene : str
        Gene of interest, or 'All'.
    maxIntronSize : int
        Max intron size used for filtering junctions near gene boundaries.
    isStranded : bool
        Whether to consider strand in site-gene mapping.
    NA_gene : Gene
        Placeholder gene object used when no gene is found.
    sample : int, default=0
        The index of this sample in the total sample set.
    numsamples : int, default=1
        Total number of samples being processed.
    QUERY_gene : Gene or None
        If a specific gene is targeted, this should be its object.
    chrom_index : list
        List of chromosomes already seen (can be empty).
    gene2D_array : list
        2D list of genes indexed by chromosome.
    site2D_array : list
        2D list of splice sites indexed by chromosome.

    Returns
    -------
    tuple
        (chrom_index, gene2D_array)
    """
    print('Processing sample '+ str(int(sample)+1)+' out of '+str(numsamples))
    #Initialisations
    lineCounter = 0
    left_alphaCount = 0
    right_alphaCount = 0
    foundCounter = 0
    newCounter = 0
    duplicateCounter = 0
    assessedCounter = 0
    left_site_new = True
    right_site_new = True
    left_site_idx = -1
    right_site_idx = -1
    dontCollapse = False
    if qGene == "All":
        print("Query gene:","All")
    else:
        print("Query gene:",QUERY_gene.getName())
        
    #Load up the bam file
    bam = pysam.AlignmentFile(bamFile, "rb")
    #Get names of chromosomes
    chroms = bam.references

    #Get addtional introns from file if user has specificed
    additionalIntrons_plus = {} #dictionary to store introns seen in each chromosome across the whoel experiment (if user gave --intronFilePath)
    additionalIntrons_minus = {} 
    additionalIntrons_unstranded = {}
    counter =0
    if intronFilePath != '':
        print("Loading introns from pre-combined intron file..")
        for c in chroms:
            additionalIntrons_plus[c]=[]
            additionalIntrons_minus[c]=[]
            additionalIntrons_unstranded[c]=[]
        for line in open(intronFilePath,"r"):
            chrom, left, right, strand = line.rstrip().split("\t")
            if isStranded:
                if strand == "+":
                    counter +=1
                    additionalIntrons_plus[chrom].append((int(left),int(right)))
                if strand == "-": 
                    counter+=1
                    additionalIntrons_minus[chrom].append((int(left),int(right)))
            else:
                counter+=1
                additionalIntrons_unstranded[chrom].append((int(left),int(right)))
        #print('additionalIntronss:',additionalIntrons_unstranded )
        print('picked up list of ',counter, 'introns from preCombined intron file')
            
    #Now time to process the bams
    for chrom in chroms:
        if chrom not in chrom_index: #Add chromsome to index if not already there (sometimes there are no annotated genes in a scaffold)
                chrom_index.append(chrom)
                site2D_array.append([])
                gene2D_array.append([])
        if qChrom == chrom or qChrom == "All":
            chrom_idx = chrom_index.index(str(chrom))

            #If this is a stranded analysis, then make a generator for the forward and reverse strand reads (accounting for paired-end reads)
            #and package them up in tuples with the appropriate strand designation
            if isStranded:
                bamgen_plus = (read for read in bam.fetch(chrom) if check_strand(strandedType,read.flag,"+") and not read.is_unmapped and not read.is_secondary and not read.is_supplementary)
                bamgen_minus = (read for read in bam.fetch(chrom) if check_strand(strandedType,read.flag,"-") and not read.is_unmapped and not read.is_secondary and not read.is_supplementary)
                
                introns_plus=bam.find_introns(bamgen_plus)
                introns_minus=bam.find_introns(bamgen_minus)
                if not dontCollapse:
                    print(chrom, "Collapsing duplicate introns to the majority strand..")
                    introns_plus, introns_minus = collapse_duplicate_introns(introns_plus,introns_minus)
                Intron_info =[(introns_plus,"+"),(introns_minus,"-")]
                if intronFilePath != '': # if we are adding some additional introns from file
                    #print(additionalIntrons_plus[chrom])
                    #print({name: 0 for name in additionalIntrons_plus[chrom] if name not in introns_plus})
                    introns_supp_plus = {name: 0 for name in additionalIntrons_plus[chrom] if name not in introns_plus}
                    introns_supp_minus = {name: 0 for name in additionalIntrons_minus[chrom] if name not in introns_minus}
                    print('added ', len(introns_supp_plus)+len(introns_supp_minus),' new introns from list')
                    Intron_info.append((introns_supp_plus,"+"))
                    Intron_info.append((introns_supp_minus,"-"))
            else:
                bamgen=(read for read in bam.fetch(chrom) if not read.is_unmapped and not read.is_secondary and not read.is_supplementary)
                introns_unstranded=bam.find_introns(bamgen)
                Intron_info = [(introns_unstranded,"?")]
                if intronFilePath != '': # if we are adding some introns from file
                    introns_supp_unstranded = {name: 0 for name in additionalIntrons_unstranded[chrom] if name not in introns_unstranded}
                    print('added ', len(introns_supp_unstranded),' new introns from list')
                    Intron_info.append((introns_supp_unstranded,"?"))
    
            #If user has provided an intron file path, then add these to the introns list with a zero count. 
            if intronFilePath != '':
                print("Loading introns from pre-combined intron file..")

                pass

            #Now iterate over introns to fill in splice sites, where b is the bamgen tuple that we made before (linking the bam gen to the strand value)
            for introns,strand in Intron_info:
                print("processing introns in region:",chrom, " strand:", strand)
                for i in introns:
                    LinQGene = False
                    RinQGene = False
                    leftpos = i[0]
                    rightpos = i[1]
                    alpha = introns[i]
                    strand = strand

                    if qGene != 'All': # if we are chasing a particular gene
                        lPosAdj = leftpos+maxIntronSize
                        rPosAdj = rightpos-maxIntronSize
                        #check if either site is within an intron length of the gene boundaries
                        if lPosAdj >= QUERY_gene.getLeftPos() and leftpos <= QUERY_gene.getRightPos():
                            LinQGene = True
                        if rPosAdj <= QUERY_gene.getRightPos() and rightpos >= QUERY_gene.getLeftPos():
                            RinQGene = True

                    #if we are taking any gene, or if the sites are query-gene-associated
                    if qGene =="All" or LinQGene or RinQGene:
    
                        if len(site2D_array[chrom_idx])>0:
                            left_site_idx = binary_site_search(site2D_array[chrom_idx],leftpos,strand,isStranded)
                            right_site_idx = binary_site_search(site2D_array[chrom_idx],rightpos,strand,isStranded)
    
                        else:
                            left_site_idx = -1
                            right_site_idx = -1
    
                        assessedCounter += 1
                        #FOR EACH SITE
                        site_left = None
                        site_right = None
                        for num, site_info in enumerate([[leftpos,left_site_idx],[rightpos,right_site_idx]]):
                            ss = None
                            # if this is a new site
                            if site_info[1] == -1:
    
                                if num == 0:
                                    left_site_new = True
                                if num == 1:
                                    right_site_new = True
    
                                newCounter = newCounter + 1
                                gene_idx = binary_gene_search(gene2D_array[chrom_idx], site_info[0], strand,isStranded)
                                #Check if we found a gene for site
                                if gene_idx >=0:
                                    gene = gene2D_array[chrom_idx][gene_idx]
                                    foundCounter +=1
                                else:
                                    gene = NA_gene
                                #Make the site
                                ss = Site(
                                        chromosome = str(chrom),
                                        pos = site_info[0],
                                        samples = numsamples,
                                        strand = strand,
                                        source = '',
                                        isStranded=isStranded
                                        )
                                ss.setGene(gene)
                                #add the site in it's ordered position
                                #bisect.insort(site2D_array[chrom_idx],ss)
                            else: #If the site has been recorded already
                                if num == 0:
                                    left_site_new = False
                                    duplicateCounter += 1
                                if num == 1:
                                    right_site_new = False
                                    duplicateCounter += 1
                                ss = site2D_array[chrom_idx][site_info[1]]
                            #record the alpha read counts
                            ss.addAlphaCount(int(alpha), sample)
                            if num == 0: # if the left site
                                site_left = ss
                            else: #if the right site
                                site_right = ss
    
                        if left_site_new:
                            bisect.insort(site2D_array[chrom_idx],site_left)
                        if right_site_new:
                            bisect.insort(site2D_array[chrom_idx],site_right)
                            #add sites as partners for eachother, and record the shared alpha
                        site_left.addPartner(site_right)
                        site_left.addPartnerCount(site_right.getPos(), alpha, sample)
                        site_right.addPartner(site_left)
                        site_right.addPartnerCount(site_left.getPos(), alpha, sample)
    bam.close()
    print("Done")



    print("Introns assessed:\t"+str(assessedCounter))
    print("Sites found:\t\t\t"+str(newCounter))
    print("Sites assigned to a Gene:\t"+str(foundCounter))
    num =0
    for s in site2D_array:
        num = num + len(s)
    #print("Sites:\t\t\t"+str(num))

    #Find competitorPositions of each site
    def findCompetitorPos():
        for c_index, c in enumerate(chrom_index):
            for site in site2D_array[c_index]:
                sPos = site.getPos()
                for p in site.getPartners():
                    for c in p.getPartners():
                        cPos = c.getPos()
                        if cPos != sPos:
                            site.addCompetitorPos(cPos)

    print('\n\nStep 1B: Identifying Competitors of each splice site...')
    findCompetitorPos()
    print('Competitors updated')
    return chrom_index, gene2D_array, site2D_array


def parse_cigar(cigarstring):
    """
    Parse a CIGAR string into a list of (operation, length) tuples.

    Mimics the behavior of pysam's AlignedSegment.cigartuples but avoids relying
    on internal integer opcodes. This provides a future-proof alternative that 
    directly returns the character-based operations (e.g. 'M', 'I', 'D').

    Parameters
    ----------
    cigarstring : str
        A standard CIGAR string (e.g., '76M2I3D20M').

    Returns
    -------
    List[Tuple[str, int]]
        A list of (operation, length) tuples, e.g. [('M', 76), ('I', 2), ...].
    """
    digits = ''
    ops = []

    for char in cigarstring:
        if char.isdigit():
            digits += char
        else:
            ops.append((char,(int(digits))))
            digits = ''
    return ops


def checkBam_pysam(bamFile, sSite, sample, isStranded, strandedType):
	targetPos = sSite.getPos()
	competitors = sSite.getCompetitorPos()
	siteStrand = sSite.getStrand()
	siteChrom = sSite.getChromosome()
	partners = []
	for partner, counts in sSite.getPartnerCounts().items():
		partners.append(partner)

	mode = sys.argv[1]

	bam = pysam.AlignmentFile(bamFile, "rb")
	bamgen = (read for read in bam.fetch(siteChrom, targetPos, targetPos + 2) if not read.is_unmapped and not read.is_secondary and not read.is_supplementary)

	for read in bamgen:
		cPos = -1 #stores position of competitor site, if found
		leftBound = read.reference_start
		rightBound = read.reference_end

		if leftBound <= targetPos and rightBound >= targetPos:
			flag = read.flag
			cigar = read.cigarstring
			cigar_tups = parse_cigar(cigar)



			#spliceSites = []
			partnerUsed = ""
			compSplicing = False
			alpha_read = False
			beta1_read = False
			compSplicing_read = False
			SimpleBeta2_flanking_read = False
			SimpleBeta2_beta1type_read = False
			SimpleBeta2_mutuallyExclusive_read = False
			
			currentPos = int(leftBound)
			#if sSite.getPos()==27281:
			#	print(read.query_name, leftBound,rightBound, cigar, cigar_tups, currentPos)
			for idx, t in enumerate(cigar_tups):
				case = t[0]
				d = t[1]
				if case in ['M', 'X', '=']:
					mappedRegion = True
					progression = True
				if case in ['N', 'D']:
					mappedRegion = False
					progression = True
				if case in ['I', 'S', 'H', 'P']:
					progression = False

				if progression:
					currentPos += int(d)
					#if sSite.getPos()==27281:
					#	print("\t",currentPos)
					#if read.query_name == "SRR3462015.1716895":
					#	print("\t",case, d, mappedRegion, progression, currentPos, (currentPos - int(d)),targetPos)

					#If we have mapped across the splice site, then record a beta1 read
					if mappedRegion and targetPos > (currentPos - int(d)) and currentPos > targetPos and currentPos >= targetPos + 1:
						if not isStranded or check_strand(strandedType,read.flag,siteStrand):
							beta1_read = True
						#if read.query_name == "SRR3462015.1716895":
						#	print("\t","\t",targetPos, currentPos, (currentPos - int(d)))
						#	print("\t","\t","beta1 read: ",beta1_read)

					if case == 'N': #Otherwise only care if we have split/gapped reads 
						lSite = currentPos - int(d)
						rSite = currentPos
						#spliceSites.append(lSite)
						#spliceSites.append(rSite)
						#if read.query_name == "SRR3462008.769256":
						#	print("\t","\t",lSite,rSite)

						if lSite == targetPos:
							partnerUsed = rSite
							alpha_read = True
						if rSite == targetPos:
							partnerUsed = lSite
							alpha_read = True

						if rSite in competitors and lSite in partners:
							compSplicing = True
							cPos = rSite
						if lSite in competitors and rSite in partners:
							compSplicing = True
							cPos = lSite

						if compSplicing and targetPos > lSite and targetPos < rSite: #strandedness protected by usage of partner site in compSplicing check
							SimpleBeta2_flanking_read = True

						if not alpha_read and not compSplicing and targetPos > lSite and targetPos < rSite: #If these are two seemingly unrelated sites with a gap spanning the targetPos
							if not isStranded or check_strand(strandedType,read.flag,siteStrand):
								sSite.addMutuallyExclusivePos(lSite)
								sSite.addMutuallyExclusivePos(rSite)
								SimpleBeta2_mutuallyExclusive_read = True

			if beta1_read and compSplicing: #strandedness protected by earlier beta1_read check
				SimpleBeta2_beta1type_read = True
			#if sSite.getPos()==27281:
			#	print(alpha_read, beta1_read, SimpleBeta2_flanking_read,SimpleBeta2_mutuallyExclusive_read,SimpleBeta2_beta1type_read)


			if SimpleBeta2_flanking_read: 
				#if sSite.getPos()==36685:
				#	print("SimpleBeta2_flanking_read",read.query_name)
				#if mode == 'combine' or mode == 'combineShallow':
				sSite.addBeta2SimpleCount(1, sample)
				if compSplicing:
					sSite.addCompetitorPos(cPos)

			elif SimpleBeta2_mutuallyExclusive_read:
				#if sSite.getPos()==36685:
				#	print("SimpleBeta2_mutuallyExclusive_read",read.query_name)
				sSite.addBeta2SimpleCount(1, sample)

			elif SimpleBeta2_beta1type_read:
				#if sSite.getPos()==36685:
				#	print("SimpleBeta2_beta1type_read",read.query_name)
				sSite.addBeta2SimpleCount(1, sample)
				if compSplicing:
					sSite.addCompetitorPos(cPos)

			elif beta1_read and not SimpleBeta2_beta1type_read:
				sSite.addBeta1Count(1, sample)
					

	#if sSite.getPos()==27281:
	#	print(sSite.getAlphaCount(0), sSite.getBeta1Count(0), sSite.getBeta2SimpleCount(0))
	bam.close()
