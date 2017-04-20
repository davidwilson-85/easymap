#!/usr/bin/python
import math
import argparse
from string import maketrans
#Script used in order to obtain primers.

parser = argparse.ArgumentParser()
parser.add_argument('-file', action="store", dest = 'File', required = "True")
parser.add_argument('-fasta', action = "store", dest = "genome", required = "True")

args = parser.parse_args()

genome = args.genome

#Tm calculation of an oligo. Based on biophp script of Joseba Bikandi https://www.biophp.org/minitools/melting_temperature/demo.php
def Tm_calculation(oligo):
	primer = float(400) #400 nM are supposed as a standard [primer]	
	mg = float(2) #2 mM are supposed as a standard [Mg2+]
	salt = float(40)  #40 mM are supposed as a standard salt concentration
	s = 0
	h = 0
	#Enthalpy and entrophy values from http://www.ncbi.nlm.nih.gov/pmc/articles/PMC19045/table/T2/ (SantaLucia, 1998)
	dic = {

	"AA": [-7.9,-22.2], 
	"AC": [-8.4,-22.4],
	"AG": [-7.8, -21.0],
	"AT": [-7.2,-20.4],	
	"CA": [-8.5,-22.7],
    "CC": [-8.0, -19.9],
    "CG": [-10.6,-27.2],
    "CT": [-7.8,-21.0],
    "GA": [-8.2,-22.2],
    "GC": [-9.8, -24.4],
    "GG": [-8.0, -19.9],
    "GT": [-8.4, -22.4],
    "TA": [-7.2,-21.3],
    "TC": [-8.2,-22.2],
    "TG": [-8.5,-22.7],
    "TT": [-7.9,-22.2]}

	#Effect on entropy by salt correction; von Ahsen et al 1999
	#Increase of stability due to presence of Mg
	salt_effect = (salt/1000)+((mg/1000)*140)
	#effect on entropy
	s+=0.368 * (len(oligo)-1)* math.log(salt_effect)
	#terminal corrections. Santalucia 1998
	firstnucleotide= oligo[0]
	if firstnucleotide=="G" or firstnucleotide=="C": h+=0.1; s+=-2.8
	if firstnucleotide=="A" or firstnucleotide=="T": h+=2.3; s+=4.1

	lastnucleotide= oligo[-1]
	if lastnucleotide=="G" or lastnucleotide=="C": h+=0.1; s+=-2.8
	if lastnucleotide=="A" or lastnucleotide=="T": h+=2.3; s+=4.1
	#compute new H and s based on sequence. Santalucia 1998
	for i in range(0,len(oligo)-1):
		f = i+ 2
		substring = oligo[i:f]

		h = h + float(dic[substring][0])
		s =s + float(dic[substring][1])

	
   	tm=((1000*h)/(s+(1.987*math.log(primer/2000000000))))-273.15
   	return tm

def obtain_values(File):
	compilation= []
	positions = open(File,"r")
	n = 0
	for line in positions.readlines():
			line = line.split("\t")
			if line[5] != "nh" and n != 0:
				data = [line[1], line[2]]
				mode = line[0]
				compilation.append(data)
			n += 1
	return compilation,mode

def reverse_complementary(oligo):
	revcomp = oligo.translate(maketrans('ACGT', 'TGCA'))[::-1]
	return revcomp	


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



data = obtain_values(args.File)[0]
mode = obtain_values(args.File)[1]

with open(genome) as fp:
	for name_contig, seq_contig in read_fasta(fp):
		if name_contig[1:].lower() == data[0][0]:
			genome = seq_contig


if mode == "snp":
	Tm_min = 60
	Tm_max = 64
	oligos = []
	for snp in data:
		position = snp[1]
		o_position1 = int(position) - 300
		o_position2 = int(position) + 500
		take_look = 0
		actual_position = 0

		while True: 
			actual_position = actual_position + o_position1+take_look+1 #position in the genome 
			window = genome[actual_position: int(position)+1] 
			g_search =window.find("G")
			c_search = window.find("C")
			if g_search== -1 or c_search == -1:
				if g_search == -1: take_look = c_search
				if c_search == -1: take_look = g_search
				else:
					oligos.append([position, "no oligo found"])
					break
			else:  take_look = min([g_search,c_search])
			
			oligo = genome[o_position1+ 1+ take_look - 21 : o_position1+ 1+take_look +1] #Python no coge la ultima posicion
			Tm = Tm_calculation(oligo)
			print take_look
			if Tm >= Tm_min and Tm <= Tm_max:
				result = [position, oligo, o_position1+ 1+ take_look]
				oligos.append(result)
				break
			if  take_look > 99:
				oligos.append([position, "no oligo found"])
				break
			else: continue 

print oligos 

	
