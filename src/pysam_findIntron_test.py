#pysamtest.py
import pysam


def check_strand(strandedType, SAMflag, siteStrand):
	'''
	Takes the flag element of a SAM and determines if that read originated from the same strand as the splice site
	'''
	if strandedType == "fr":
		if (SAMflag & 64) or not (SAMflag & 1): #If first in pair, or not paired read
			if (SAMflag & 16):				#If read reverse strand
				readStrand = "-"
			else:
				readStrand = "+"
		else:								#else (if paired and second in pair)
			if (SAMflag & 16):					#If read is reverse strand
				readStrand = "+"    
			else:
				readStrand = "-"
	if strandedType == "rf":			#Flip everything if it is "rf" library prep
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



bamFile="/Users/craig/Documents/SpliSERv1/SpliSER/tests/test1_process_basic/9520_tophat2.bam"
bam = pysam.AlignmentFile(bamFile, "rb")
chroms = bam.references
print(chroms)
chrom = chroms[0]
#bamgen=(read for read in bam.fetch(chrom) if not read.is_unmapped and not read.is_secondary and not read.is_supplementary)

bamgen_plus = (read for read in bam.fetch(chrom) if check_strand("rf",read.flag,"+") and not read.is_unmapped and not read.is_secondary and not read.is_supplementary)

bamgen_plus_raw = (read for read in bam.fetch(chrom) if check_strand("rf",read.flag,"-"))
for rdx,read in enumerate(bamgen_plus_raw):
	#print(read.query_name,read.is_unmapped,read.is_secondary,read.is_supplementary)
	if read.query_name=="SRR3462010.3537636":
		print(read.query_name,read.is_unmapped,read.is_secondary,read.is_supplementary)

Bamgens =[(bamgen_plus,"+")]
for b,strand in Bamgens:
	print(b,strand)
I=bam.find_introns(bamgen_plus)
print(len(I))

for idx,intron in enumerate(I):
	if idx<10:
		print(intron)

for b,strand in Bamgens:
	introns=bam.find_introns(b)
	for idx, i in enumerate(introns):
		if idx < 10:
			LinGene = False
			RinGene = False
			left_pos = i[0]
			right_pos = i[1]
			alpha = introns[i]
			print(left_pos,right_pos,alpha,strand)























