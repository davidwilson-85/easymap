#Input file 
#input = str(args.input)
input = 'output_analysis.txt'
f1 = open(input, 'r')
lines = f1.readlines()	

finput = '34k_genome_2c.fa'
f2 = open(input, 'r') 
fasta_lines = f2.readlines()

#Output file
f3 = open('output_ordered.csv', 'w')
f3.write('@' + 'DATA\t' + 'Contig' + '\tins'+ '\t' + 'NT' + '\t' + '  RD' + '\t' + 'Direction' + '\n')

contigs = []

#create a list with all the genome contigs
for i, line in enumerate(fasta_lines):
	if line.startswith('>'): #fasta sequences start with '>'
		sp = line.split(' ')  #because some names have whitespaces and extra info that is not written to sam file
		cont = sp[0].strip()  #strip() is to remove the '\r\n' hidden chars
		cont = cont[1:]       #to remove the first char of string (>)
		if cont not in contigs:
			contigs.append(cont)


data = []

for i, line in enumerate(lines):
	if not line.startswith('@'):
		data.append(line.split())


import operator
sorted_data = sorted(data, key=lambda e: (e[1], e[2]))


import csv

with open('output_ordered.csv', 'wb') as f3:
    writer = csv.writer(f3)
    writer.writerows(sorted_data)

