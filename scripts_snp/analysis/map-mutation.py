#!/usr/bin/python
#Map-snp.py is a script written in python as a part of Easymap software. It is used in mapping by sequencing of SNP. The script can take as input three different situations,
# a backcross, an outcross when the mutated ecotype is not in the reference background and an outcross when the mutated ecotype is in the reference background.

#par f2wt

##Arguments:
	#"fichero" which is the file that comes as an input.
	#"output" refers to the name of the output file is going to be generated.
	#"window_size" it is a number which indicates the size of the windows that will be created in order to generate allele frequencies averages
	#"window_space" number which indicates the space between the number which is in the middle of a window
	#"fasta" corresponds to the fasta file from which the input has been obtained
	#"mode" can be back meaning backcross or out meaning outcross
	#"width" . . . 
	#"parental_modality" takes as values ref and noref, meaning if the line sequenced is from the reference background or not.
	#"filtering_mode" . . . 
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-fichero', action="store", dest = 'input', required = "True")
parser.add_argument('-output', action="store", dest = 'output', required = "True")
parser.add_argument('-window_size', action="store", dest = 'size', required = "True")
parser.add_argument('-window_space', action="store", dest = 'space', required = "True")
parser.add_argument('-fasta', action="store", dest = 'fasta_input', required = "True")
parser.add_argument('-mode', action="store", dest = 'mode', required = "True")
parser.add_argument('-interval_width', action="store", dest = 'interval_width', required = "True") 
parser.add_argument('-parental_modality', action="store", dest = 'modality') #ref = parental in reference background; noref = parental not in reference background
parser.add_argument("-snp_analysis_mode", action="store", dest = "control") #Depending on which control is being used: par, f2wt

fasta_input = args.fasta_input
size = int(args.size)
space = int(args.space)
mode = args.mode
output = args.output
result= open( output ,"w")
modality = args.modality
interval_width = int(args.interval_width)
control = args.control
if mode == "back":
	modality = "n"
result.write("if outcross:\n-@ lines: window, average, boost, chromosome\n-! line: min_max_window, max_max_window, max_boost, chromosome\n-? lines: chromosome, min_big_window, max_big_window\n-* lines: chromosome, min_window, max_window, boost_value\nif backcross:\n-@ line: window, average, chromosome\n-! line: min_window. max_window, max_average, chromosome\n-? line: chromosome, min_big_window, max_big_window\n-* line: chormosome, min_window, max_window, average\nif control is F2 WT\n-@ lines: window, ratio, chromosome\n-! line: min_max_window, max_max_window, max_ratio, chromosome\n-? lines: chromosome, min_big_window, max_big_window\n-* lines: chromosome, min_window, max_window, ratio\n")
result.close()

#Gets all the parameters from a file. Uses arguments chromosome and input file. Creates a dictionary per chromosome. dic[POSITION-SNP]=[list other values stored] 
def getinfo(chro, inpu):
	n = 0 			#counter n will be used in order not to take into account the header
	dicpos = {} 
	for lines in inpu.readlines():
		if n == 0:
			n +=1
			continue
		indiv = lines.split("\t")
		if indiv[0] != chro:
			continue
		dicpos[indiv[1]] = []
		for n in range(2,len(indiv)):
			dicpos[indiv[1]].append(indiv[n].rstrip())
	return dicpos

#Calculates average of a list of AF in a window.	
def calculation_average(li):
	average_list = sum(li)/len(li)
	return average_list
#From the dictionary generated in getinfo, knowing the chromosome and its lenght, the function generates windows according to the parameters size and space between them.
def chromosomal_position(size,space, SNP,ch, chromosomal_lenght, mode, modality, control): 
	if modality == "ref":
		c = 0.7
	elif modality == "noref" or "n": #if we are dealing with a backcross, eventhough it is in the ref background we are looking for high AF SNP, that's why modality is n
		c = 0.3	
	windowsize = float(size)
	windowspace = float(space)
	i = 0 	#is the value in the middle of the windows and the one will be used in order to identify a concrete window																							 
	chromosomal_size = float(chromosomal_lenght) 										
	dictionary_windows = {} 													
	a = 0 																								
	b = 0 
	while i < chromosomal_size:      								
		a = i - windowsize/2
		b = i + windowsize/2
		snps_window = []
		for key in SNP:	
			s = float(key)											
		 	if s >= a and s <b:		
				AF =float(SNP[key][-1])/(float(SNP[key][-2])+float(SNP[key][-1])) 
				if control == "par":	
					if c == 0.7 and AF < c:
						snps_window.append(AF)
					elif c == 0.3 and AF> c:							
			 			snps_window.append(AF)
				elif control == "f2wt":
					AFwt = float(SNP[key][-3])/(float(SNP[key][-3])+ float(SNP[key][-4]))
					snps_window.append(AF)
					snps_window.append(AFwt)
		if len(snps_window) != 0:
			if control == "par":										
				average_FA = calculation_average(snps_window)
				dictionary_windows[i] = []
				dictionary_windows[i].append(average_FA)
			elif control == "f2wt":
				average_FA = calculation_average(snps_window[0])
				average_FA_WT = calculation_average(snps_window[1]) 
				dictionary_windows[i] = []
				dictionary_windows[i].append(average_FA)
				dictionary_windows[i].append(average_FA_WT)
		elif len(snps_window) == 0 and modality == "ref" and control == "par" :  #in the modality outcross of mutant in the reference background, it is possible that a window will not contain any SNP. We will suppose a value near 0.
			average_FA = 0.01	
			dictionary_windows[i] = []
			dictionary_windows[i].append(average_FA)                       
		i += windowspace
	return dictionary_windows		

#This is a threshold step that will remove windows which do not overpass certain values, which will be chosen depending on the mode. This is meant to make the process faster.
def threshold_step(windows, mini_average, maxi_average):
	refined_dic = {}
	for window in windows:					
		FA= windows[window][0]
		if FA > mini_average and FA <= maxi_average:
			refined_dic[window]= windows[window]
	return refined_dic

#This function is just to get a list of window positions which are in order due to the fact dictionaries are not sorted.	
def union_points(windows):
	position = []
	for window in windows:
		position.append(window)
	postion = position.sort()
	return position



#Data function is different for a backcross or an outcross. This function stores the windows and their attributes boost and average in a file.
#It also saves different values that will be used later for further processiong of the data, as the chromosome with the higher attribute, its value or the dictionary of positions of that chromosome
def data_analysis(window, position, chromosome, maximum_position, best_parameter, best_chromosome,best_dictionary, size, mode, control):  
	result= open(output ,"a")
	dictionary ={} 
	for items in position:			
		average = window[items][0]
		items = int(items)
		if mode =="back" and control == "par" :
			parameter = average
			result.write("@"+"\t"+str(items) + "\t"+ str(average) +"\t"+str(chromosome)+ "\n")   
		elif mode =="out" and control == "par":
			boost = 1/abs(1-(1/max(average, 1-average)))
			parameter = boost
			dictionary[items]= boost
			result.write("@"+"\t"+str(items)+"\t"+ str(average)+"\t"+ str(boost)+ "\t"+str(chromosome)+ "\n")
		elif control == "f2wt":
			average_mut = window[items][0]
			average_WT = window[items][1]
			best_ditionary[items][0]
			ratio = int(average_mut) - int(average_WT)
			parameter = ratio
			result.write("@	"+str(items)+ "\t"+ str(ratio) + "\t" + str(chromosme) + "\n")

		if best_parameter == "n/p" or float(parameter) > float(best_parameter):
			maximum_position = []
			items = float(items)
			maximum_position.append(items - size/2)
			maximum_position.append(items + size/2)
			best_chromosome = chromosome
			best_parameter = parameter
	if chromosome == best_chromosome: 		
		if mode == "out" or control == "f2wt":
			best_dictionary= dictionary
		if mode == "back" and control == "par":
			best_dictionary = window
	return maximum_position, best_parameter, best_chromosome, best_dictionary

def final_processing_A(result, interval_width):       #is the one used in the oc
	great_positions = [] #Creation of a list of values with the maximum parameter (boost/ratio), then calculation of the middle value  and create a bigger window
	for windows in result[-1]:
		b_value = result[-1][windows]
		if b_value == result[1]:
			great_positions.append(int(windows)) 
	peak_value = calculation_average(great_positions)
	#Corrector factor in order to take windows which boost is not as big as the biggest but can contain the mutation
	big_window = []
	big_window.append(int(peak_value) + interval_width/2)
	big_window.append(int(peak_value) - interval_width/2)
	possible_windows= []

	for maximum_position in result[3]:
		if float(maximum_position) <= big_window[0] and float(maximum_position)>= big_window[1]:
			maximum_position = float(maximum_position)
			maximum_position = int(maximum_position)
			possible_windows.append(maximum_position)
	#Creation of a list of all the windows in the interval, creation of a big and conservative window. Save of the big window, windows contained on it and best window.		
	possible_windows = sorted(possible_windows)
	min_i = big_window[1]
	maxi_i = big_window[0]
	r= open(output ,"a")
	r.write("!"+"\t"+ str(result[0][0])+"\t" + str(result[0][1]) + "\t"+ str(result[1])+"\t"+ result[2] + "\n") 
	escribir ="?"+"\t"+result[2]+ "\t"+str(min_i)+"\t"+str(maxi_i)+"\n"
	r.write(escribir)
	for windos in possible_windows:
		windos = str(windos)
		escribir = "*"+ "\t"+result[2] +"\t" + str(int(windos)-size/2)+ "\t" + str(int(windos)+size/2)+"\t" + str(result[3][int(windos)]) + "\n"
		r.write(escribir)  


def final_processing_B(result,interval_width):
	#Creation of a conservative bigger window which will contain different windows with an average bigger than a threshold.
	great_positions = []
	for windows in result[-1]:
		average = result[-1][windows][0]
		if average == result[1]:
			great_positions.append(windows)
	peak_value = calculation_average(great_positions)
	min_i = int(peak_value - interval_width/2) 
	maxi_d = int(peak_value + interval_width/2)
	candidate_list = []
	for window in result[-1]:
		if float(result[-1][window][0]) >= 0.9:
			if window >= min_i and window <= maxi_d:
				candidate_list.append(window)			
	candidate_list = sorted(candidate_list)
	r= open(output ,"a")
	r.write("!"+"\t"+str(result[0][0])+"\t"+str(result[0][1]) + "\t"+ str(result[1])+"\t"+str(result[2]) + "\n") 
	r.write("?"+"\t"+str(result[2]) +"\t"+ str(min_i)+"\t"+ str(maxi_d) +"\n")
	for widos in candidate_list:
		r.write("*"+"\t"+str(result[2]) +"\t" + str(int(widos)-size/2) +"\t"+ str(int(widos)+size/2)+"\t" + str(result[-1][widos][0])+ "\n")


#This function enables to obtain data regarding chromosomes and lenght of the,
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

ch = {}
#From the data of read_fasta, I create a dictionary with the name of the contigs and its lenght
with open(fasta_input) as fp:
	for name_contig, seq_contig in read_fasta(fp):
		ch[name_contig[1:len(name_contig)]] = len(seq_contig)

#The calling of the different functions depends on whether we are working in an outcross or backcross

function_used = [final_processing_B, final_processing_A]
if mode == "back":
	final_processing = function_used[0]
if mode == "out":
	final_processing = function_used[1]
if control == "f2wt":
	final_processing = function_used[1]

#Call of the function chromosome by chromosome, the final result comes from the chromosome with the higher parameters
z = 0

for chromosome in ch:  
	inpu = open(args.input, "r")	
	genome = getinfo(chromosome,inpu)		
	windows = chromosomal_position(size, space, genome,chromosome, ch[chromosome], mode, modality, control) 
	if mode == "out":
		if modality == "noref":
			mini_average = 0.4
			maxi_average = 1
		elif modality == "ref":
			mini_average= 0
			maxi_average= 0.6    
		filtered_windows= threshold_step(windows, mini_average, maxi_average) 
	x_value= union_points(filtered_windows)
	if z == 0:
		result_data = data_analysis(filtered_windows, x_value, chromosome, "n/p", "n/p", "n/p", "n/p", size, mode, control)
		best = result_data[1]
		maximum_position = result_data[0]
		best_chromosome = result_data[2]
		best_dictionary = result_data[3]
	else:
		result_data = data_analysis(filtered_windows, x_value, chromosome,maximum_position,best, best_chromosome, best_dictionary, size, mode, control)
	z += 1
final_processing(result_data, interval_width)


		





