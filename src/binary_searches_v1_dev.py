# binary_searches.py

def binary_gene_search(array, pos, strand, isStranded):
	'''
	Take an array and search for a Gene within it whose left and right bounds contain the given genomic position.
	Return its position within the Array. Otherwise return -1.

	If the search fails (e.g. overlapping genes not perfectly sorted), also checks three elements either side.
	'''
	length = int(len(array))
	if length == 0:
		return int(-1)

	idx = length // 2
	past_max = length
	past_min = 0
	last_idx = -1
	new_idx = idx
	stuck = False
	found = False
	while stuck == False and found == False:
		if int(pos) >= int(array[idx].getLeftPos()) and int(pos) <= int(array[idx].getRightPos()) and (strand == array[idx].getStrand() or isStranded==False or (strand != '+' and strand != '-')):
			found = True
			break
		elif int(pos) >= int(array[idx].getRightPos()):
			new_idx = idx + ((past_max-idx)//2)
			past_min = idx
		elif int(pos) <= int(array[idx].getLeftPos()):
			new_idx = idx - ((idx-past_min)//2)
			past_max = idx
			if idx == 1:
				new_idx = 0

		if idx != last_idx:
			last_idx = idx
			idx = new_idx
		else:
			stuck = True

	if found == False and stuck == True:
		for i in range(-3, 3):
			if idx+i >= 0 and idx+i < len(array)-1 and int(pos) >= int(array[idx+i].getLeftPos()) and int(pos) <= int(array[idx+i].getRightPos()):
				if strand == array[idx+i].getStrand() or isStranded==False or (strand != '+' and strand != '-'):
					found = True
					stuck = False
					idx = idx+i

	if found == True and stuck == False:
		return int(idx)
	else:
		return int(-1)


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
