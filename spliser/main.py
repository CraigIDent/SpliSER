"""
SpliSER- Splice-site Strength Estimation from RNA-seq
"""
#Version 1.0.0 - 24th April 2025
version = "v1.0.0"

import sys
import timeit
import time
import argparse
from ast import literal_eval

from spliser.Gene_Site_Iter_Graph_v1 import Gene, Site, Iter
from spliser.process_v1 import process
from spliser.site_ops_v1 import calculateSSE
from spliser.gene_creator_v1 import createGenes
from spliser.bam_parser_v1 import preCombineIntrons
from spliser.combine_v1 import combine, combineShallow

def DiffSpliSER_output(samplesFile,combinedFile, outputPath, minReads, qGene):

	outDiff = open(outputPath+str(qGene)+".DiffSpliSER.tsv", "w+")
	correct = True
	fileFinished = False
	#record the samples we're assessing
	samples = 0
	allTitles=[]
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
				t_beta = int(currentVals[idx][7])+int(currentVals[idx][8]) # add beta1 and beta2Simple
				t_SSE = float(currentVals[idx][5])
				if t_alpha+t_beta >= minReads: # if this sample passes the minimum read count for this site
					outDiff.write("\t"+str(t_alpha)+"\t"+str(t_beta)+"\t"+"{0:.3f}".format(t_SSE))
				else:
					outDiff.write("\tNA\tNA\tNA")
			outDiff.write("\n")

def GWAS_output(samplesFile,combinedFile, outputPath, minReads, qGene, minSamples):
	print(qGene)
	fileFinished = False
	correct = True
	#record the samples we're assessing
	samples = 0
	allTitles=[]
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
					t_beta = int(currentVals[idx][7])+int(currentVals[idx][8]) # add beta1 and beta2Simple values
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
def main():
	print("\nSpliSER "+version+" SKB LAB\n")
	# this won't be run when imported
	start = timeit.default_timer()

	#Parse arguments form the command line
	parser = argparse.ArgumentParser(description="SpiSER - Splice Site Strength Estimates from RNA-seq")
	subparsers = parser.add_subparsers(dest="command")

	#Parser for arguments when user calls command 'findIntrons'
	parser_introns = subparsers.add_parser('preCombineIntrons')
	parser_introns.add_argument('-L', '--ListOfBAMFiles', dest='BAMPathList', help="A comma-separated (no spaces) list of BAM files in this experiment (which will be later combined) ", required=True)
	parser_introns.add_argument('-o', '--outputPath', dest='outputPath', help="Absolute path, including file prefix where SpliSER will write the .introns.tsv file", required=True)
	parser_introns.add_argument('-c', '--chromosome', dest='qChrom', nargs='?', default='All', type=str,help="optional: limit SpliSER to one chromosome/scaffold -  default: All", required=False)
	parser_introns.add_argument('--isStranded', dest='isStranded', default=False, action='store_true')
	parser_introns.add_argument('-s', '--strandedType', dest='strandedType', nargs='?', type=str, required=False, help="optional: Strand specificity of RNA library preparation, where \"rf\" is first-strand/RF and \"fr\" is second-strand/FR - default : fr")
	parser_introns.add_argument('-A', '--annotationFile', dest='annotationFile', help="optional: gff3 or gtf file for finding strand of introns where there is conflicting strands in reads", required=False)

	#Parser for arguments when user calls command 'process'
	parser_process = subparsers.add_parser('process')
	parser_process.add_argument('-B', '--BAMFile', dest='inBAM', help="The mapped RNA-seq file in BAM format", required=True)
	parser_process.add_argument('-o', '--outputPath', dest='outputPath', help="Absolute path, including file prefix where SpliSER will write the .SpliSER.tsv file", required=True)
	parser_process.add_argument('-A', '--annotationFile', dest='annotationFile', help="optional: gff3 or gtf file matching the reference genome used for alignment", required=False)
	parser_process.add_argument('-t', '--annotationType', dest='aType', nargs='?', default='gene', type=str, help="optional: The feture to be extracted from the annotation file - default: gene", required=False)
	parser_process.add_argument('-c', '--chromosome', dest='qChrom', nargs='?', default='All', type=str,help="optional: limit SpliSER to one chromosome/scaffold -  default: All", required=False)
	parser_process.add_argument('-g', '--gene', dest='qGene', nargs='?', default='All', type=str, help="optional:Limit SpliSER to splice sites falling in a single locus (requires --chromosome and --annotationFile and --maxIntronSize)", required=False)
	parser_process.add_argument('-m', '--maxIntronSize', dest='maxIntronSize', nargs='?', default=0, type=int, required=False, help="optional: required if passing the --gene parameter, the max intron size used in aligning the bam file")
	parser_process.add_argument('--isStranded', dest='isStranded', default=False, action='store_true')
	parser_process.add_argument('-s', '--strandedType', dest='strandedType', nargs='?', type=str, required=False, help="optional: Strand specificity of RNA library preparation, where \"rf\" is first-strand/RF and \"fr\" is second-strand/FR - default : fr")
	#parser_process.add_argument('--capCounts', dest='capCounts', default=False, action='store_true',help="optional flag: Early stop beta-read counting (proportional to alpha reads). ie stop of read counting for rare sites in high-coverage regions once SSE reaches 0.000. Stop counting beta reads when they exceed 2000 (if alpha = 0) or otherwise 2000*alpha. SSE will still be accurate to 3 decimal places, beta1/beta2 ratios (internal values) may be biased.")
	parser_process.add_argument('-I','--intronFilePath', dest='intronFilePath',nargs='?',default='',type=str,required=False,help="optional: Path to the output of preCombineIntrons command for this experiment, containing all introns seen in this experiment")
	#Parser for arguments when user calls command 'combine'
	parser_combine = subparsers.add_parser('combine')
	parser_combine.add_argument('-S', '--samplesFile', dest='samplesFile', required=True, help="A three-column .tsv file, each line containing a sample name, the absolute path to a processed .SpliSER.tsv file input, and the absolute path to the original bam file")
	parser_combine.add_argument('-o', '--outputPath', dest='outputPath', required=True, help="Absolute path to a folder, including file_prefix where SpliSER will write the output .tsv file")
	parser_combine.add_argument('-g', '--gene', dest='qGene', nargs='?', default='All', type=str, help="optional:Limit SpliSER to splice sites falling in a single locus - default All", required=False)
	parser_combine.add_argument('--isStranded', dest='isStranded', default=False, action='store_true')
	parser_combine.add_argument('-s', '--strandedType', dest='strandedType', nargs='?', default="fr", type=str, required=False, help="optional: Strand specificity of RNA library preparation, where \"rf\" is first-strand/RF and \"fr\" is second-strand/FR - default : fr")
	#parser_combine.add_argument('--capCounts', dest='capCounts', default=False, action='store_true',help="optional flag: Early stop beta-read counting (proportional to alpha reads). ie stop of read counting for rare sites in high-coverage regions once SSE reaches 0.000. Stop counting beta reads when they exceed 2000 (if alpha = 0) or otherwise 2000*alpha. SSE will still be accurate to 3 decimal places, beta1/beta2 ratios (internal values) may be biased.")

	parser_combineShallow = subparsers.add_parser('combineShallow')
	parser_combineShallow.add_argument('-S', '--samplesFile', dest='samplesFile', required=True, help="A three-column .tsv file, each line containing a sample name, the absolute path to a processed .SpliSER.tsv file input, and the absolute path to the original bam file")
	parser_combineShallow.add_argument('-g', '--gene', dest='qGene', nargs='?', default='All', type=str, help="optional:Limit SpliSER to splice sites falling in a single locus - default All", required=False)
	parser_combineShallow.add_argument('-o', '--outputPath', dest='outputPath', required=True, help="Absolute path to a folder, including file_prefix where SpliSER will write the output .tsv file")
	parser_combineShallow.add_argument('--isStranded', dest='isStranded', default=False, action='store_true')
	parser_combineShallow.add_argument('-m', '--minSamples', dest='minSamples', required=False, nargs='?', default=0, type=int, help="For optional filtering: For any given splice site, the minimum number of samples passing the --minReads filter in order for a site to be kept in the analysis (Raising this value increases speed of combineShallow step)")
	parser_combineShallow.add_argument('-r', '--minReads', dest='minReads',required=False, nargs='?', default=10, type=int, help="For optional filtering: The minimum number of reads giving evidence for a splice site needed for downstream analyses - default: 10")
	parser_combineShallow.add_argument('-e','--minSSE', dest='minSSE',required=False, nargs='?', default=0.00, type=float, help="For optional filtering: The minimum SSE of a site for a given sample, for it to be considered in the --minSamples filter - default: 0.00")
	parser_combineShallow.add_argument('-s', '--strandedType', dest='strandedType', nargs='?', type=str, required=False, help="optional: Strand specificity of RNA library preparation, where \"rf\" is first-strand/RF and \"fr\" is second-strand/FR - default : fr")
	#parser_combineShallow.add_argument('--capCounts', dest='capCounts', default=False, action='store_true',help="optional flag: Early stop beta-read counting (proportional to alpha reads). ie stop of read counting for rare sites in high-coverage regions once SSE reaches 0.000. Stop counting beta reads when they exceed 2000 (if alpha = 0) or otherwise 2000*alpha. SSE will still be accurate to 3 decimal places, beta1/beta2 ratios (internal values) may be biased.")

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
	elif (command == 'process' or command == 'combine' or command =='combineShallow' or command == 'preCombineIntrons') and (kwargs.get('isStranded') is True) and (kwargs.get('strandedType') is None):
		parser.error("--isStranded requires parameter --strandedType/-s as fr or rf")
	else: #otherwise
		# Print the full command-line input for transparency
		print("Command line input:", ' '.join(sys.argv))
		print("")
		# call appropriate functions
		globals()[command](**kwargs)

	stop = timeit.default_timer()
	print("Total runtime (s): \t{:.3f}".format(stop - start))
#EOF

if __name__ == "__main__":
	main()