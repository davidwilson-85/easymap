########################################################
#First we sort the data in the output_analysis.txt file#
########################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-a', action="store", dest = 'input')#Input 'output_analysis.txt'
parser.add_argument('-b', action="store", dest = 'finput')#Input '34k_genome_2c.fa'
parser.add_argument('-c', action="store", dest = 'output1')#Output1 'output_ordered.csv'
parser.add_argument('-d', action="store", dest = 'output2')#Output2 'sorted_insertions.txt'
parser.add_argument('-m', action="store", dest = 'mode', default = 'P')#Output2 'sorted_insertions.txt'
args = parser.parse_args()

#Input 'output_analysis.txt
input = str(args.input)
f1 = open(input, 'r')
lines = f1.readlines()	

#Input '34k_genome_2c.fa'
finput = str(args.finput)
f2 = open(finput, 'r') 
fasta_lines = f2.readlines()

#Output 'output_ordered.csv'
output1 = str(args.output1)
f3 = open(output1, 'w')
f3.write('@' + 'DATA\t' + 'Contig' + '\tins'+ '\t' + 'NT' + '\t' + '  RD' + '\t' + 'Direction' + '\n')

#Lists
contigs = []
data = []

#Create a list with all the genome contigs
for i, line in enumerate(fasta_lines):
	if line.startswith('>'): #fasta sequences start with '>'
		sp = line.split(' ')  #because some names have whitespaces and extra info that is not written to sam file
		cont = sp[0].strip()  #strip() is to remove the '\r\n' hidden chars
		cont = cont[1:]       #to remove the first char of string (>)
		if cont not in contigs:
			contigs.append(cont)

#Create a list from the input file
for i, line in enumerate(lines):
	if not line.startswith('@'):
		data.append(line.split())
		
#Sort list and write to file
import operator
sorted_data = sorted(data, key=lambda e: (e[1], e[2]))

import csv
with open(output1, 'wb') as f3:
    writer = csv.writer(f3)
    writer.writerows(sorted_data)
   
f1.close()
f2.close()
f3.close()

############################
#Now we sort the insertions#
############################

insertion_id = 1

#Input file 
input = str(args.output1)
f1 = open(input, 'r')
lines = f1.readlines()	

#Output file
output2 = str(args.output2)
f2 = open(output2, 'w')
f2.write('@' + 'DATA\t' + 'Contig' + '\tins'+ '\t' + 'NT' + '\t' + '  RD' + '\t' + 'Direction' + '\n')

insertion_id = 1

if args.mode == 'P': 
	#Insertion sorting
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
				f2.write(sp[0] + '\t' + sp[1] + '\t' + str(insertion_id) + '\t' + sp[2] + '\t' + sp[3] + '\t' + sp[4] )
		
			else:
				f2.write(sp[0] + '\t' + sp[1] + '\t' + str(insertion_id) + '\t' + sp[2] + '\t' + sp[3] + '\t' + sp[4] )

	f1.close()
	f2.close()



	#Now we create the candidate regions for each insertion
	f2 = open(args.output2, 'r')
	lines = f2.readlines()
	candidate_regions = list() #This list will have the format: list(list(d1, d2, e))

	for e in range(1, (insertion_id + 1)):
		d1 = float('inf')
		d2 = 0
		for i, line in enumerate(lines):
			if not line.startswith('@'):
				sp = line.split()
				if sp[0].strip() == 'PAIRED' and int(sp[2].strip()) == e: 

					#The following module creates a candidate region for the insertion
					#Delimiters:

					if sp[5].strip() == 'R': #reverse
						p = int(sp[3])
						if p < d1: 
							d1 = p

					if sp[5].strip() == 'F': #forward
						p2 = int(sp[3])
						if p2 > d2:
							d2 = p2
	
		cr = list()
		cr.append(d1)
		cr.append(d2)
		cr.append(e)
		
		candidate_regions.append(cr)

	f2 = open(args.output2, 'a')
	for i in candidate_regions: 
		f2.write('@#')
		f2.write(((str(i).strip('[')).strip(']')) + '\n')

elif args.mode == 'S': 
	#Insertion sorting
	for i, line in enumerate(lines):
		if not line.startswith('@'): 
			sp = line.split(',')
			if str(sp[0]).strip() == 'LOCAL' and  str(sp[4]).strip() == 'TOTAL':
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
		
#			elif sp[0].strip('\t') == 'PAIRED':
#				f2.write(sp[0] + '\t' + sp[1] + '\t' + str(insertion_id) + '\t' + sp[2] + '\t' + sp[3] + '\t' + sp[4] )
		
#			else:
#				f2.write(sp[0] + '\t' + sp[1] + '\t' + str(insertion_id) + '\t' + sp[2] + '\t' + sp[3] + '\t' + sp[4] )

	f1.close()
	f2.close()






















