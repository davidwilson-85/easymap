

contig_source = 'test.fa'


# Function to parse fasta file (based on one of the Biopython IOs)
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



# Read contig fasta file
with open(contig_source) as fp:
	fastalist = list()
	for name_contig, seq_contig in read_fasta(fp):
		innerlist = list()
		innerlist.append(name_contig)
		innerlist.append(len(seq_contig))
		fastalist.append(innerlist)
try:
	max_list = list()
	for c in fastalist:
		max_list.append(int(c[1]))
	max_length = max(max_list)
except:
	max_length = fastalist[0][1]

print max_length