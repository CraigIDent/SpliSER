# bam_parser.py

import sys
import subprocess
import re
import pysam
digit_pattern = re.compile(r'\D')
char_pattern = re.compile(r'\d')
from site_ops_v1_dev import check_strand

def checkBam(bamFile, sSite, sample, isStranded, strandedType):
	##LOOK OUT FOR SECONDARY AND SUPPLEMENTARY ALIGNMENTS
	targetPos = sSite.getPos()
	competitors = sSite.getCompetitorPos()
	siteStrand = sSite.getStrand()

	partners = []
	for partner, counts in sSite.getPartnerCounts().items():
		partners.append(partner)

	bamview = subprocess.Popen(
		['samtools', 'view', str(bamFile), str(sSite.getChromosome()) + ':' + str(targetPos) + '-' + str((int(targetPos)+1))],
		stdout=subprocess.PIPE
	)
	bamstream = bamview.stdout
	mode=sys.argv[1]
	for line in bamstream:
		dline = line.decode('ascii')
		values = str(dline).split('\t')
		cPos = -1

		if len(values) > 3:
			leftBound = int(values[3])
			if leftBound <= int(targetPos):
				flag = int(values[1])
				cigar = str(values[5])
				digits = list(filter(None, digit_pattern.split(cigar)))
				chars = list(filter(None, char_pattern.split(cigar)))

				spliceSites = []
				partnerUsed = ""
				compSplicing = False

				alpha_read = False
				beta1_read = False
				compSplicing_read = False
				SimpleBeta2_flanking_read = False
				SimpleBeta2_beta1type_read = False
				SimpleBeta2_mutuallyExclusive_read = False
				currentPos = int(leftBound)

				for idx, d in enumerate(digits):
					case = chars[idx]
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

						if int(targetPos) >= (currentPos - int(d)) and currentPos > int(targetPos) and currentPos > int(targetPos) + 1:
							if mappedRegion:
								if isStranded:
									if check_strand(strandedType, flag, siteStrand):
										beta1_read = True
								else:
									beta1_read = True

						if case in ['N']:
							lSite = currentPos - int(d) - 1
							rSite = currentPos - 1
							spliceSites.append(lSite)
							spliceSites.append(rSite)

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

							if compSplicing and targetPos > lSite and targetPos < rSite:
								SimpleBeta2_flanking_read = True
							if not alpha_read and not compSplicing and targetPos > lSite and targetPos < rSite:
								if isStranded:
									if check_strand(strandedType, flag, siteStrand):
										SimpleBeta2_mutuallyExclusive_read = True
								else:
									SimpleBeta2_mutuallyExclusive_read = True

				if beta1_read and compSplicing:
					SimpleBeta2_beta1type_read = True

				if alpha_read and compSplicing:
					sP = set(partners)
					sS = set(spliceSites)
					sI = sP.intersection(sS)
					for p in sI:
						if p != partnerUsed:
							sSite.addPartnerBeta2DoubleCount(p, 1, sample)

				elif SimpleBeta2_flanking_read:
					if mode == 'combine' or mode == 'combineShallow':
						sSite.addBeta2SimpleCount(1, sample)
						if compSplicing:
							sSite.addCompetitorPos(cPos)

				elif SimpleBeta2_mutuallyExclusive_read:
					sSite.addBeta2SimpleCount(1, sample)

				elif SimpleBeta2_beta1type_read:
					sP = set(partners)
					sS = set(spliceSites)
					sI = sP.intersection(sS)
					for p in sI:
						sSite.addPartnerBeta2DoubleCount(p, 1, sample)
					sSite.addBeta2SimpleCount(1, sample)
					if compSplicing:
						sSite.addCompetitorPos(cPos)

				elif beta1_read and not SimpleBeta2_beta1type_read:
					sSite.addBeta1Count(1, sample)


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
	if isStranded:
		bamgen = (read for read in bam.fetch(siteChrom, targetPos, targetPos + 1) if check_strand(strandedType,read.flag,siteStrand) and not read.is_unmapped and not read.is_secondary and not read.is_supplementary)
	else:
		bamgen = (read for read in bam.fetch(siteChrom, targetPos, targetPos + 1) if not read.is_unmapped and not read.is_secondary and not read.is_supplementary)

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

						if compSplicing and targetPos > lSite and targetPos < rSite:
							SimpleBeta2_flanking_read = True
						if not alpha_read and not compSplicing and targetPos > lSite and targetPos < rSite:
							SimpleBeta2_mutuallyExclusive_read = True

			if beta1_read and compSplicing:
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
