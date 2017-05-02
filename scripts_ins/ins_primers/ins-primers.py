#python lin-primers.py -sam_in alignment4.sam -var_in variants.txt -sam_out 5_prime_end_reads


import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-sam_in', action="store", dest = 'input_sam')
parser.add_argument('-var_in', action="store", dest = 'input_var')
parser.add_argument('-sam_out', action="store", dest = 'output')

args = parser.parse_args()

import time
start_time = time.time()

#Input 
input1 = args.input_sam
f1 = open(input1, 'r')
sam_lines = f1.readlines()	

input2 = args.input_var
f2 = open(input2, 'r')
var_lines = f2.readlines()	


var_list = list()

for i, line in enumerate(var_lines):
	if not line.startswith('@'):
		sp = line.split('\t')
		var_list.append([sp[1].strip(), sp[2].strip()])

for insertion in var_list:
	
	#Output
	output = args.output
	f3 = open(output + '_' + str(insertion[0] + '-' + str(insertion[1])), 'w')

	ins_chromosome = insertion[0]
	ins_position = int(insertion[1])
	for i, line in enumerate(sam_lines):
		if line.startswith('@'):
			f3.write(line)

		if not line.startswith('@'):
			sp = line.split('\t')
			chromosome = (sp[2].strip()).lower()
			position = int(sp[3])
			cigar = sp[5]
			sequence = sp[9]

			if chromosome == ins_chromosome and position in range(ins_position - 200, ins_position + 200):

					x = ''
					x2 = ''

					for i in cigar: 
						if i == 'M' or i == 'D' or i == 'I' or i == 'N' or i == 'S' or i == 'H' or i == 'P' or i == 'X' : 
							x += str(i) + '\t'
						else:
							x += str(i)

					sp2 = x.split()
					for i in sp2:
						if 'M' in i:
							x2 += '1'
						if 'S' in i:
							x2 += '0'

					if x2.startswith('1'): 			 
						
						# CONTAR NUCLEOTIDOS LEIDOS
						l = 0
						l2 = 0
						for i in reversed(sp2): 				 							#Lee los elementos del CIGAR de atras a adelante
							if 'S' in i:
								num = i.replace('S', '')
								l = int(l) + int(num)
								break

						# SOBREESCRIBIR SECUENCIA CON NUCLEOTIDOS NO LEIDOS Y EL CIGAR
						sequence2 = sequence[len(sequence)-l: ]
						cigar2 = str(l) + 'M'												#Sobreescribimos los cigars con "matches" para que los siguientes programas no den errores
						
						# ESCRIBIR EN EL SAM DE SALIDA
						newline = line.replace(sequence, sequence2)
						newline2 = newline.replace(cigar, cigar2)
						f3.write(newline2)

					elif x2.endswith('0'):
						pass