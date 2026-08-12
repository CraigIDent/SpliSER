# binary_searches.py

# Binary search for an insertion point into the array of genes, followed by (bounded) interval scan
def binary_gene_search(array, pos, strand, isStranded):

	pos = int(pos)

	if len(array) == 0:
		return -1

	# Find rightmost gene whose start is <= pos
	low = 0
	high = len(array) - 1
	idx = -1

	while low <= high:
		mid = (low + high) // 2

		if int(array[mid].getLeftPos()) <= pos:
			idx = mid
			low = mid + 1
		else:
			high = mid - 1

	if idx == -1:
		return -1

	# Check every earlier interval that could potentially contain pos
	for i in range(idx, -1, -1):
	
		gene = array[i]
	
		# Nothing at i or earlier can reach pos
		if gene.getMaxRightBefore() < pos:
			break
	
		strand_match = (
			not isStranded
			or strand not in ("+", "-")
			or strand == gene.getStrand()
		)
	
		if (int(gene.getLeftPos()) <= pos <= int(gene.getRightPos()) and strand_match):
			return i

	return -1


def binary_site_search(array, pos, strand, isStranded):
	'''
	Take an array and search for a Site within it whose position matches the query genomic position.
	Return its position within the Array. Otherwise return -1.
	'''
	length = int(len(array))
	idx = length // 2
	past_max = length
	past_min = 0
	last_idx = -1
	new_idx = idx
	stuck = False
	found = False

	while stuck == False and found == False:
		if int(pos) == int(array[idx].getPos()):
			if strand == array[idx].getStrand() or isStranded == False or (strand != '+' and strand != '-'):
				found = True
				break
			else:
				for a in [idx-1, idx+1]:
					if a >= 0 and a < len(array) and int(pos) == int(array[a].getPos()) and (strand == array[a].getStrand() or isStranded == False or (strand != '+' and strand != '-')):
						idx = a
						found = True
				break

		elif int(pos) >= int(array[idx].getPos()):
			new_idx = idx + ((past_max - idx) // 2)
			past_min = idx
		elif int(pos) <= int(array[idx].getPos()):
			new_idx = idx - ((idx - past_min) // 2)
			past_max = idx
			if idx == 1:
				new_idx = 0

		if idx != last_idx:
			last_idx = idx
			idx = new_idx
		else:
			stuck = True

	if found == True:
		return idx
	else:
		return -1
