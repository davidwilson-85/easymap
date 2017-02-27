#This module will process the information in the .sam file to obtain the absolute frequency of aligments ending per nucleotide during local aligments. 

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-a', action="store", dest = 'input')
parser.add_argument('-b', action="store", dest = 'output')
parser.add_argument('-c', action="store", dest = 'finput')
parser.add_argument('-m', action="store", dest = 'mode', default='P')
args = parser.parse_args()

#Input file 
input = args.input
f1 = open(input, 'r')
lines = f1.readlines()	


#fasta input
fasta_input = str(args.finput)
fasta_f1 = open(fasta_input, 'r')
fasta_lines = fasta_f1.readlines()

#Output: If the paired-reads analysis is being performed we will oppen the output to append data to the file, else we create the output and open it in write mode
if args.mode == 'P':
	output = args.output
	f2 = open(output, 'a')
elif args.mode == 'S': 
	output = args.output
	f2 = open(output, 'w')


#create a list with all the genome contigs
contigs = []
for i, line in enumerate(fasta_lines):
	if line.startswith('>'):						 #fasta sequences start with '>'
		sp = line.split(' ') 						 #because some names have whitespaces and extra info that is not written to sam file
		cont = sp[0].strip() 						 #strip() is to remove the '\r\n' hidden chars
		cont = cont[1:]      						 #to remove the first character of a string (>)
		if cont not in contigs:
			contigs.append(cont)

#Analyze SAM file:									 We create three dictionaries in which we will compute the absolute frequency of how many times  
for c in contigs: 									 #a reads alignment finishes in each nucleotide, separating the information according to which  
	di_left = dict()								 #side of the read is aligned (left or right). d2_1 and d2_2 are lists of all the positions  
	di_right = dict()								 #in which a reads aligment finishes, repeated as many times as that happens, required for the creation of the dictionaries
	di_total = dict()
	d2_1 = []
	d2_2 = []
	for i, line in enumerate(lines):
		if not line.startswith('@'):
			sp = line.split('\t')
			cont = sp[2]
			if cont == c and cont != '*':				
				p = int(sp[3]) 						 #Read position
				cigar = sp[5].strip() 				 #Then we define the CIGAR parameter, from which we will extract the aligned nucleotides of each read 
				if cigar != '*':					 #and their relative position (left/right)
					x = ''							 #X is a string containing the cigar separating the M and S characters by a tabulation
					x2 = ''							 #X2 is a string containing ones and zeros that will map the position of the M and S characters 
					l = 0							 #to determine the part of the read that is aligned
					for i in cigar: 
						if i != 'M' and i != 'S' : 
							x += str(i) 
						else:
							x += str(i) + '\t'
					sp2 = x.split()
					for i in sp2:
						if 'M' in i:
							x2 += '1'
						if 'S' in i:
							x2 += '0'
					if x2.startswith('0'): 			 #In this case, the CIGAR starts with unaligned nucleotides so the first aligned nucleotide is 
						d2_2.append(str(p))			 #the one in which the aligment is cut 
					elif x2.endswith('0'): 			 #If the aligment ends with a zero, the unaligned nucleotides are at the end of the read so we
						for i in sp2: 				 #Need to calculate the aligned nucleotides in order to get the position of the last one aligned (pf)
							if 'M' in i:
								num = i.replace('M', '')
								l = int(l) + int(num)
							if 'D' in i:
								num = i.replace('D', '')
								l = int(l) + int(num)
							if 'I' in i:
								num = i.replace('I', '')
								l = int(l) - int(num)
						pf = int(p) + int(l) - 1
						d2_1.append(str(pf))				
					elif x2 == '1':
						pass

	#Key count										 #The "key count" is the transformation of the information in the lists (d2_1 and d2_2) in dictionaries
	#TOTAL DICTIONARY								 #acumulating the read depth of each nucleotide
	for i in d2_1:
		try: 
			di_total[i] =  1 + di_total[i]
		
		except KeyError:
			di_total[i] = 1
			
	for i in d2_2:
		try: 
			di_total[i] =  1 + di_total[i]
		except KeyError:
			di_total[i] = 1
	
	#LEFF AND RIGHT DICTIONARIES
	for i in d2_1:
		try: 
			di_left[i] = 1 + di_left[i]
		
		except KeyError:
			di_left[i] = 1

	for i in d2_2:
		try: 
			di_right[i] =  1 + di_right[i]
		except KeyError:
			di_right[i] = 1
		
	#Writting in the output file

	for key,value in sorted(di_left.items(), key=lambda i: int(i[0])):
		f2.write('LOCAL\t' + c + '\t' + str(key) + '\t'+ str(value) + '\tLEFT\n')
	for key,value in sorted(di_right.items(), key=lambda i: int(i[0])):
		f2.write('LOCAL\t' + c + '\t' + str(key) + '\t'+ str(value) + '\tRIGHT\n')
	for key,value in sorted(di_total.items(), key=lambda i: int(i[0])):
		f2.write('LOCAL\t' + c + '\t' + str(key) + '\t'+ str(value) + '\tTOTAL\n')

f2.close()