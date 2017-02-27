#
#
# Hace falta ordenar por contigo y nt el archivo output-analysis.txt antes 
#
#

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-a', action="store", dest = 'input')
parser.add_argument('-b', action="store", dest = 'output')
args = parser.parse_args()

#Input file 
#input = str(args.input)
input = 'output_ordered.csv'
f1 = open(input, 'r')
lines = f1.readlines()	

#Output file
f2 = open('output_sorted.txt', 'w')
f2.write('@' + 'DATA\t' + 'Contig' + '\tins'+ '\t' + 'NT' + '\t' + '  RD' + '\t' + 'Direction' + '\n')

insertion_id = 1

for i, line in enumerate(lines):
	if not line.startswith('@'): 
		sp = line.split(',')
		if str(sp[0]).strip() == 'PAIRED' and  str(sp[4]).strip() == 'TOTAL':
			p = int(sp[2])
			contig = sp[1].strip('\t')
			try:
				d = abs(int(p2) - int(p))
				if d > 100 or contig != contig2:
					insertion_id = insertion_id + 1
					f2.write(sp[0] + '\t' + sp[1] + '\t' + str(insertion_id) + '\t' + sp[2] + '\t' + sp[3] + '\t' + sp[4])
					p2 = p
					contig2 = contig

				else:
					f2.write(sp[0] + '\t' + sp[1] + '\t' + str(insertion_id) + '\t' + sp[2] + '\t' + sp[3] + '\t' + sp[4] )
					p2 = p
					contig2 = contig
			except:
				f2.write(sp[0] + '\t' + sp[1] + '\t' + str(insertion_id) + '\t' + sp[2] + '\t' + sp[3] + '\t' + sp[4] )
				p2 = p
				contig2 = contig
		
		elif sp[0].strip('\t') == 'LOCAL':
			f2.write(sp[0] + '\t' + sp[1] + '\t' + ' - ' + '\t' + sp[2] + '\t' + sp[3] + '\t' + sp[4] )
		
		else:
			f2.write(sp[0] + '\t' + sp[1] + '\t' + str(insertion_id) + '\t' + sp[2] + '\t' + sp[3] + '\t' + sp[4] )

