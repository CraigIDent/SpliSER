#site_ops.py
import pysam
import bisect
from operator import add, truediv, mul, sub
from Gene_Site_Iter_Graph_v0_1_8 import Site
from binary_searches_v1_dev import binary_gene_search, binary_site_search

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

def findAlphaCounts(bedFile, qChrom, qGene, maxIntronSize, isStranded, NA_gene, sample=0, numsamples=1, QUERY_gene=None, chrom_index=None, gene2D_array=None, site2D_array=None):
    """
    Take each junction from the bed file, If it falls within a recorded Gene object,
    record the number of reads evidencing usage (alpha1 reads) of each splice site forming the junction.

    Parameters
    ----------
    bedFile : str
        The absolute path to a junctions.bed file for a given alignment.
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
    print(QUERY_gene)
    if QUERY_gene.getLeftPos != None:
        print("Query gene:",QUERY_gene)
    else:
        print("Query gene:","All")
    #check each junction in this bed file
    for line_idx, line in enumerate(open(bedFile, 'r')):
        #check for a header line
        values = str(line).split("\t")

        if len(values) == 12: # if not a header line
            lineCounter = lineCounter+1
            #values = str(line).split("\t")
            chrom = str(values[0])
            LinGene = False
            RinGene = False
            if chrom not in chrom_index: #Add chromsome to index if not already there (sometimes there are no annotated genes in a scaffold)
                chrom_index.append(chrom)
                site2D_array.append([])
                gene2D_array.append([])
            if qChrom == chrom or qChrom == "All":
                # find the index, which we'll use to specify the array we'll search.
                chrom_idx = chrom_index.index(str(chrom))
                strand = str(values[5])
                #identify the positions of each splice site
                flank = values[10].split(",")
                leftpos = int(values[1])+int(flank[0])
                rightpos = int(values[2])-int(flank[1])
                alpha = int(values[4])

                if qGene != 'All': # if we are chasing a particular gene
                    lPosAdj = leftpos+maxIntronSize
                    rPosAdj = rightpos-maxIntronSize
                    #check if either site is within an intron length of the gene boundaries
                    if lPosAdj >= QUERY_gene.getLeftPos() and leftpos <= QUERY_gene.getRightPos():
                        LinGene = True
                    if rPosAdj <= QUERY_gene.getRightPos() and rightpos >= QUERY_gene.getLeftPos():
                        RinGene = True
                #if we are taking any gene, or if the sites are query-gene-associated
                if qGene =="All" or LinGene or RinGene:

                    if len(site2D_array[chrom_idx])>0:
                        left_site_idx = binary_site_search(site2D_array[chrom_idx],leftpos,strand,isStranded)
                        right_site_idx = binary_site_search(site2D_array[chrom_idx],rightpos,strand,isStranded)

                    else:
                        left_site_idx = -1
                        right_site_idx = -1

                    assessedCounter += 2
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
                                    chromosome = str(values[0]),
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
                        ss.addAlphaCount(int(values[4]), sample)
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
    print("Sites assessed:\t"+str(assessedCounter))
    print("Sites found:\t\t\t"+str(newCounter))
    print("Sites assigned to a Gene:\t"+str(foundCounter))
    num =0
    for s in site2D_array:
        num = num + len(s)
    print("Sites:\t\t\t"+str(num))

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


def findAlphaCounts_pysam(bamFile, qChrom, qGene, maxIntronSize, isStranded,strandedType, NA_gene, sample=0, numsamples=1, QUERY_gene=None, chrom_index=None, gene2D_array=None, site2D_array=None):
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
    print(QUERY_gene)
    if QUERY_gene.getLeftPos != None:
        print("Query gene:",QUERY_gene)
    else:
        print("Query gene:","All")
    #Load up the bam file
    bam = pysam.AlignmentFile(bamFile, "rb")
    #Get names of chromosomes
    chroms = bam.references

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
                Bamgens =[(bamgen_plus,"+"),(bamgen_minus,"-")]
            else:
                bamgen=(read for read in bam.fetch(chrom) if not read.is_unmapped and not read.is_secondary and not read.is_supplementary)
                Bamgens = [(bamgen,"?")]
    
            #Now iterate over Bamgens to fill in splice sites, where b is the bamgen tuple that we made before (linking the bam gen to the strand value)
            for b,strand in Bamgens:
                print("processing region:",chrom," strand:",strand)
                introns=bam.find_introns(b)
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
    
                        assessedCounter += 2
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
    print("Sites assessed:\t"+str(assessedCounter))
    print("Sites found:\t\t\t"+str(newCounter))
    print("Sites assigned to a Gene:\t"+str(foundCounter))
    num =0
    for s in site2D_array:
        num = num + len(s)
    print("Sites:\t\t\t"+str(num))

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

def trueDivCatchZero(array1, array2):
    """
    Takes an Array (array1) and divides it by another array (array2), elementwise.
    Where an element of array2 is equal to zero, the resulting value will be Zero (rather than undefined).
    """

    array3 = [0.0]*int(len(array1))
    for i, a2 in enumerate(array2):
        if a2 > 0.0:
            array3[i] = truediv(array1[i],a2)
    return array3

def subIntNoNeg(element1, element2):
    ans = int(element1) - int(element2)
    if ans >= 0:
        return ans
    else:
        return 0

def findBeta2Counts(site, numSamps):
    #find the beta2 counts for this site
    Partners= site.getPartners() # array of Site objects the current site is known to partner with
    PartnerCounts = site.getPartnerCounts() #dictionary of known partners and alpha counts shared with the target site
    beta2CrypticCounts = [0]*numSamps
    beta2CrypticWeighted = [0.00]*numSamps
    #total alpha reads of the target site
    TotalAlphas = site.getAlphaCounts()

    for i, pSite in enumerate(Partners):

        for competitorPos, counts in pSite.getPartnerCounts().items(): # getting this from partner counts instead of sSite.getPartners().. .getPos() so the funciton is compatible with process and combine commands
            #Get all beta2 simple counts where partner and competitor flank the target site.
            if pSite.getPos() > site.getPos() and int(competitorPos) < site.getPos():
                site.addBeta2SimpleCounts(counts)
                site.addPartnerBeta2DoubleCounts(pSite.getPos(),counts)
            elif pSite.getPos() < site.getPos() and int(competitorPos) > site.getPos():
                site.addBeta2SimpleCounts(counts)
                site.addPartnerBeta2DoubleCounts(pSite.getPos(),counts)

        #get the total alpha of that partner, for all samples
        pAlphas = pSite.getAlphaCounts()
        #get the alpha reads shared between this site and the partner
        pCounts = PartnerCounts[pSite.getPos()]
        #beta 2 reads are all alpha reads of that partner, minus those shared between partner site and the target site
        b2 = [x - y for x, y in zip(pAlphas, pCounts)]
        #See if the b2 reads for this partner contain any double counts (that we observed in checkBam)
        if pSite.getPos() in site.getPartnerBeta2DoubleCounts():
            doubleCounts = site.getPartnerBeta2DoubleCounts()[pSite.getPos()]
            #adjust for double-counter reads
            b2 = [subIntNoNeg(x,y) for x, y in zip(b2, doubleCounts)]
        #tally Cryptic beta 2 counts
        beta2CrypticCounts = [x + y for x, y in zip(beta2CrypticCounts, b2)]
        #calculate weights for the given partner
        pWeights = trueDivCatchZero(pCounts, TotalAlphas)
        #tally weighted beta2 counts
        for i in range(0,numSamps):
            b2[i] = b2[i]*pWeights[i]
        beta2CrypticWeighted = [x + y for x, y in zip(beta2CrypticWeighted, b2)]

    #update values for this site
    site.addBeta2CrypticCounts(beta2CrypticCounts)
    site.updateBeta2Weighted(beta2CrypticWeighted)


def calculateSSE(site, isbeta2Cryptic):

    alpha = list(site.getAlphaCounts())
    beta1 = site.getBeta1Counts()
    beta2Simple = site.getBeta2SimpleCounts()
    betas = [x + y for x, y in zip(beta1, beta2Simple)]

    if isbeta2Cryptic:
        beta2w = site.getBeta2WeightedCounts()
        betas = [x + y for x, y in zip(betas, beta2w)]

    denominator = [x + y for x, y in zip(alpha, betas)]

    site.setSSEs(trueDivCatchZero(list(alpha), list(denominator)))