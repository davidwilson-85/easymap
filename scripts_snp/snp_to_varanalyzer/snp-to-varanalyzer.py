
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-a', action="store", dest='input')
parser.add_argument('-b', action="store", dest='output')

args = parser.parse_args()


#input
f1 = open(args.input, 'r')
lines = f1.readlines()

#output
f2 = open(args.output, 'w')
f2.write('#data' + '\t' + 'contig' + '\t' + 'pos' + '\t' + 'ref' + '\t' + 'alt' + '\n')


for i, line in enumerate(lines):
	sp = line.split()
	f2.write('snp' + '\t' + sp[0].strip() + '\t' + sp[1].strip() + '\t' + sp[2].strip() + '\t' + sp[3].strip() + '\n')
