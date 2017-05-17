#python lin-primers_v3.py -sam_in alignment4.sam -var_in variants.txt -sam_out out

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
input2 = args.input_var

var_list = list()

with open(input2) as f2:
	for line in f2:
		if not line.startswith('@'):
			sp = line.split('\t')
			var_list.append([sp[1].strip(), sp[2].strip()])

for insertion in var_list:
	
	#Output
	output = args.output
	f5 = open(output  + str(insertion[0] + '_' + str(insertion[1])) + '_5' + '.fq', 'w')
	f3 = open(output  + str(insertion[0] + '_' + str(insertion[1])) + '_3' + '.fq', 'w')

	ins_chromosome = insertion[0]
	ins_position = int(insertion[1])


	with open(input1) as f1:
		for line in f1:
			if not line.startswith('@'):
				sp = line.split('\t')
				chromosome = (sp[2].strip()).lower()
				position = int(sp[3])
				cigar = sp[5]
				sequence = sp[9]
				quality = sp[10]

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

						#Lecturas del extremo 5' de la insercion
						if x2.startswith('1') and x2.endswith('0'): 			 										#Lecturas cortadas en el 3', del extremo 5' de la insercion
							
							# CONTAR NUCLEOTIDOS LEIDOS
							for i in reversed(sp2): 				 													#Lee los elementos del CIGAR de atras a adelante
								if 'S' in i:
									num = i.replace('S', '')
									l = int(num)
									break

							# SOBREESCRIBIR SECUENCIA CON NUCLEOTIDOS NO LEIDOS Y EL CIGAR
							sequence2 = sequence[len(sequence)-l: ]
							quality2 = quality[len(quality)-l: ]
		
							# ESCRIBIR EN EL FQ DE SALIDA
							f5.write('@'+ sp[0] + '\n' + sequence2 + '\n' + '+' + '\n' + quality2 + '\n' ) 

						#Lecturas del extremo 3' de la insercion
						elif x2.startswith('0') and x2.endswith('1'):													#Lecturas cortadas en el 5', del extremo 3' de la insercion

							# CONTAR NUCLEOTIDOS LEIDOS
							for i in sp2: 				 													#Lee los elementos del CIGAR
								if 'S' in i:
									num = i.replace('S', '')
									l = int(num)
									break

							# SOBREESCRIBIR SECUENCIA CON NUCLEOTIDOS NO LEIDOS Y EL CIGAR
							sequence2 = sequence[ :l]
							quality2 = quality[ :l]
		
							# ESCRIBIR EN EL FQ DE SALIDA
							f3.write('@'+ sp[0] + '\n' + sequence2 + '\n' + '+' + '\n' + quality2 + '\n' ) 
