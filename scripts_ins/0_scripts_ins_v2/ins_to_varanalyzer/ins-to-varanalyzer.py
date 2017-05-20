
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-a', action="store", dest='input')
parser.add_argument('-b', action="store", dest='output')

args = parser.parse_args()


#output
f2 = open(args.output, 'w')
f2.write('#data' + '\t' + 'contig' + '\t' + 'pos' + '\t' + 'ref' + '\t' + 'alt' + '\n')

#First I create a list of the insertions
insertion_list = list()
with open(args.input) as f1:
	for line in f1:
		if not line.startswith('@'):
			sp = line.split()
			if sp[0] == 'LOCAL':
				insertion = sp[2].strip()
				if insertion not in insertion_list:
					insertion_list.append(insertion)
				

#Then I loop through the insertion list, and take the highest RD value:
for ins in insertion_list:
	max_rd = 0
	with open(args.input) as f1:
		for line in f1: 
			if not line.startswith('@'):
				sp = line.split()
				if sp[0] == 'LOCAL' and sp[2].strip() == ins and max_rd < int(sp[4]):
					max_rd = int(sp[4])
					chrom = sp[1].strip()
					pos = sp[3].strip()
		f2.write('lim' + '\t' + chrom + '\t' + pos + '\t' + '-' + '\t' + '-' + '\n')

f2.close()