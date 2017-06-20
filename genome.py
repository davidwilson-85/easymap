
def read_fasta(fp):
	name, seq = None, []
	for line in fp:
		line = line.rstrip()
		if line.startswith('>'):
			if name: yield (name, ''.join(seq))
			name, seq = line, []
		else:
			seq.append(line)
	if name: yield (name, ''.join(seq))

mut_list = [5925329,5864789,5807566,5769044,5637309,5622262,5497045,5405459,5085978,5067121]
mut_list = [5637309]

mut_list = [2137705,6306495,12385219,17669884,6751663,2336221,5710945,14171828,12310030,4246170,4916889,16553990,5163558,7100534,10003051,6231302,17546556,2185860,2384375,1788833,3575452,11864850,14295161,16875834,9133823,13774804,13501110,672934,1987346,2979912]


with open("user_data/4.at4.fa","r") as anaconda:
	n = 0
	for i in read_fasta(anaconda):
		y = i[1]
		for p in mut_list:
			p = p-1
			position = y[p:p+10]
			n += 1
			print n,position