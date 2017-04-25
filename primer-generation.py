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
		try:
			h = h + float(dic[substring][0])
			s =s + float(dic[substring][1])
		except:
			return 0

	
   	tm=((1000*h)/(s+(1.987*math.log(primer/2000000000))))-273.15
   	return tm

def obtain_values(File):
	compilation= set()
	ordered_list = []
	positions = open(File,"r")
	n = 0
	for line in positions.readlines():
		line = line.split("\t")
		if line[5] != "nh" and n != 0:
			data = (line[1], line[2])
			ordered_list.append(data)
			mode = line[0]
			compilation.add(data)
		n += 1
	compilation = sorted(compilation, key=lambda x: ordered_list.index(x)) #sets doesn't keep the order, this function sorts the set as it was made in the list, set is necessary in order to remove the redundant mutations
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


inputt= obtain_values(args.File)
datos = inputt[0]
mode = inputt[1]
data = []
for items in datos:
	data.append(items)


with open(genome) as fp:
	for name_contig, seq_contig in read_fasta(fp):
		if name_contig[1:].lower() == data[0][0]:
			genome = seq_contig

def rule_1(oligo,location):
	last_element = len(oligo)
	if location == "downstream" : oligo == reverse_complementary(oligo)
	while True:
		end_of_primer = 21
		begin_of_primer = 0
		oligo = oligo[:last_element]
		guanine = oligo.rfind("G")
		cytosine = oligo.rfind("C")
		
		#Checking wheter there are both G and C in the oligo
		if guanine != -1 and cytosine != -1:
			last_element =  max([guanine,cytosine])
		else:
			if guanine == -1 and cytosine == -1:
				found = "no"
				return found, found, "-"
				break
			elif guanine == -1 and cytosine != -1:
				last_element = cytosine
			elif guanine != -1 and cytosine == -1:
				last_element = guanine
		begin = last_element - end_of_primer 
		end = last_element 
		while end - begin < 26 and end - begin > 16:
			if begin < 0:
				break
			primer = oligo[begin:end+1]
			Tm = Tm_calculation(primer)
			if Tm > 60 and Tm < 64:
				found = "yes"
				return found, primer, Tm
			elif Tm < 60:
				begin -= 1
			elif Tm > 64:
				begin += 1
def rule_2(oligo,location):
	begin_of_primer = 0
	if location == "dowsntream" : oligo = reverse_complementary(oligo)
	while True:
		end_of_primer =21 + begin_of_primer

		if end_of_primer > len(oligo):
			found = "no"
			return found, found, "-"
		while end_of_primer - begin_of_primer < 26 and end_of_primer - begin_of_primer > 16:
			primer = oligo[begin_of_primer:end_of_primer+1]
			Tm = Tm_calculation(primer)
			if Tm >= 60 and Tm <= 64:
				found = "yes"
				return found, primer, Tm
			elif Tm < 60:
				end_of_primer += 1
			elif Tm > 64:
				end_of_primer -= 1		
		begin_of_primer +=1


#Total size of the amplification 
size_i = 300
size_f = 500
if mode == "snp":
	dic = {}
	dic_Tm= {}
	for snp in data:
		dic[snp[1]] = []
		dic_Tm[snp[1]] = []
		#upstream primer
		up_primer_pos = int(snp[1]) - size_i
		oligo = genome[up_primer_pos-1 : up_primer_pos + 100]
		result = rule_1(oligo, "upstream")
		if result[0] == "no":
			result = rule_2(oligo, "upstream")
			if result[0] == "no":
				dic[snp[1]].append("not found")
				dic[snp[1]].append("-")
				dic_Tm[snp[1]].append("-")
				dic_Tm[snp[1]].append("-")
		if result[0] == "yes":
			dic[snp[1]].append(result[1])
			dic_Tm[snp[1]].append(str(result[2]))
			#downstream primer
			down_primer_pos = int(snp[1]) + size_f
			oligo = genome[down_primer_pos-1 : down_primer_pos + 100]
			result = rule_1(oligo,"downstream")
			if result[0] == "no":
				result = rule_2(oligo, "downstream")
				if result[0] == "no":
					dic[snp[1]].append("not found")
					dic_Tm[snp[1]].append("-")
			if result[0] == "yes":
				dic[snp[1]].append(result[1])
				dic_Tm[snp[1]].append(str(result[2]))

result = open("primers_file", "w")
m = 0
for snp in data:
	if m == 0:
		result.write("Mode\tchromosome\tposition\tforward primer\tTm forward\treverse primer\tTm reverse\n")
	result.write(mode+"\t"+snp[0]+"\t"+snp[1]+"\t"+dic[snp[1]][0]+"\t"+dic_Tm[snp[1]][0]+"\t"+dic[snp[1]][1]+ "\t"+dic_Tm[snp[1]][1]+"\n")
	m +=1
	print snp[1], dic[snp[1]]
print m
		

	
