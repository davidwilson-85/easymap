#python analysis_Sergio.py -fichero ejemplo -window_size 500000 -window_space 500000 -fasta at_chr1.fas -mode back 
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-fichero', action="store", dest = 'input', required = "True")
parser.add_argument('-output', action="store", dest = 'output', required = "True")
parser.add_argument('-window_size', action="store", dest = 'size', required = "True")
parser.add_argument('-window_space', action="store", dest = 'space', required = "True")
parser.add_argument('-fasta', action="store", dest = 'fasta_input', required = "True")
parser.add_argument('-mode', action="store", dest = 'mode', required = "True")
parser.add_argument('-correction_factor', action="store", dest = 'FC') 
parser.add_argument('-parental_modality', action="store", dest = 'modality') #same = parental in reference background; different = parental not in reference background
args = parser.parse_args()
fasta_input = args.fasta_input
size = int(args.size)
space = int(args.space)
mode = args.mode
output = args.output
result= open( output ,"w")
modality = args.modality
result.write("if outcross:"+"\n"+"@ lines: window, average, boost, chromosome" +"\n" + "! line: min_max_window, max_max_window, max_boost, chromosome" + "\n" + "? lines: chromosome, min_big_window, max_big_window"+ "\n" + "* lines: chromosome, min_window, max_window, boost_value" + "\n" "if backcross:"+"\n"+ "@ line: window, average, chromosome" +"\n" + "! line: min_window. max_window, max_average, chromosome" + "\n"+ "? line: chromosome, min_big_window, max_big_window" + "\n" + "* line: chormosome, min_window, max_window, average"+ "\n")
result.close()

#Gets all the parameters from a file. Uses arguments chromosome and input file. Creates a dictionary per chromosome. dic[POSITION-SNP]=[list other values stored]
def getinfo(chro, inpu):
	n = 0
	dicpos = {} 
	for lines in inpu.readlines():
		if n == 0:
			n +=1
			continue
		indiv = lines.split("\t")
		if ">"+indiv[0] != chro:
			continue
		dicpos[indiv[1]] = []
		for n in range(2,len(indiv)):
			dicpos[indiv[1]].append(indiv[n].rstrip())
	return dicpos

#Calculates average of a list of AF in a window. Two modes are present, if mutant is in the reference background or if it is in a different background.	
def calculation_average(li, modality):
	if modality == "same":
		c = 0.65
	elif modality == "different":
		c = 0.35
	new_list = []
	for items in li:
		items = float(items)
		if c == 0.65 and items < c:
			new_list.append(items)
		elif c == 0.35 and items> c:
			new_list.append(items)
	average_list = sum(new_list)/len(new_list)
	return average_list
#Calculation of the standard deviation of a list of AF in a position. This parameter is not being used so far.
def calcular_sd(lista):
	import math
	m = sum(lista)/len(lista)
	nuevalista = []
	for items in lista:
		nuevalista.append((float(items)- m)**2)
	result= sum(nuevalista)/len(nuevalista)
	return math.sqrt(result)
#From the dictionary generated in getinfo, knowing the chromosome and its lenght, the function generates windows according to the parameters size and space between them.
def chromosomal_position(size,space, SNP,ch, chromosomal_lenght, mode, modality):
	windowsize = float(size)
	windowspace = float(space)
	i = 0 																								 
	chromosomal_size = float(chromosomal_lenght) 										
	final = {} 													
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
		 		snps_window.append(AF)
		if len(snps_window) != 0:										
			average_FA = calculation_average(snps_window, modality)
			final[i] = []
			final[i].append(average_FA)
			if mode == "out":
				sd_FA = calcular_sd(snps_window)		
				final[i].append(sd_FA)
		elif len(snps_window) == 0 and modality == "same":
			average_FA = 0.01	
			final[i] = []
			final[i].append(average_FA)                       
		i += windowspace
	return final		

#This is a threshold step that will remove windows which do not overpass certain values, which will be chosen depending on the mode. This is meant to make the process faster.
def threshold_step(windows, mini_average, maxi_average, mini_Cv):
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
def data_backcross(window, position, chromosome, best_parameter, size):  
	result= open(output ,"a") 
	for items in position:			
		average = window[items][0]
		result.write("@"+str(items) + "\t"+ str(average) +"\t"+str(chromosome)+ "\n")   
		if best_parameter == "T" or float(average) > float(best_parameter):
			maximum_position = []
			items = float(items)
			maximum_position.append(items - size/2)
			maximum_position.append(items + size/2)
			best_parameter = average
			best_chromosome = chromosome
			best_dictionary = window
	return maximum_position, best_parameter, best_chromosome, best_dictionary
def data_outcross(window, position, chromosome, real_best, size): 
	result= open(output ,"a")
	dic_paramet = {}
	for items in position:			
		average = window[items][0]
		boost = 1/abs(1-(1/max(average, 1-average)))
		items = str(int(float(items)))
		average = str(average)
		boost = str(boost)
		dic_paramet[items]= boost
		result.write("@"+items+"\t"+ average+"\t"+ boost+ "\t"+str(chromosome)+ "\n")   
		if float(boost) > float(real_best) or real_best == 0:   
			best_position = []
			items = float(items)
			best_position.append(items - size/2)
			best_position.append(items + size/2)
			real_best = boost
			best_chromosome = chromosome
			best_dic = dic_paramet
	return best_position, real_best, best_chromosome, dic_paramet

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
z = 0
ch = {}
#From the data of read_fasta, I create a dictionary with the name of the contigs and its lenght
with open(fasta_input) as fp:
	for name_contig, seq_contig in read_fasta(fp):
		ch[name_contig] = len(seq_contig)

#The calling of the different functions depends on whether we are working in an outcross or backcross
if mode == "out":
	#Call of the function chromosome by chromosome, the final result comes from the chromosome with the higher parameters
	for chromosome in ch:  
		inpu = open(args.input, "r")	
		genome = getinfo(chromosome,inpu)		
		windows = chromosomal_position(size, space, genome,chromosome, ch[chromosome], mode, modality)
		if modality == "different":
			mini_average = 0.4
			maxi_average = 1
		elif modality == "same":
			mini_average= 0
			maxi_average= 0.6    
		filtered_windows= threshold_step(windows, mini_average, maxi_average, 0) 
		x_value= union_points(filtered_windows)
		if z == 0:
			result = data_outcross(filtered_windows, x_value, chromosome, 0, size)
			best = result[1]
		else:
			result = data_outcross(filtered_windows, x_value, chromosome, best, size)
		z += 1
	#Corrector factor in order to take windows which boost is not as big as the biggest but can contain the mutation
	FC = float(args.FC)
	threshold = float(result[1]) * FC
	possible_windows= []
	for best_position in result[3]:
		if float(result[3][best_position]) >= threshold:
			best_position = float(best_position)
			best_position = int(best_position)
			possible_windows.append(best_position)
	#Creation of a list of all the windows bigger than the threshold, creation of a big and conservative window. Save of the big window, windows contained on it and best window.		
	possible_windows = sorted(possible_windows)
	a = min(possible_windows)
	c = max(possible_windows)
	b = size/2
	min_i = a - b
	maxi_i = c+b
	r= open(output ,"a")
	r.write("!"+ str(result[0][0])+"\t" + str(result[0][1]) + "\t"+ str(result[1])+"\t"+ result[2] + "\n") 
	escribir ="?"+result[2]+ "\t"+str(min_i)+"\t"+str(maxi_i)+"\n"
	r.write(escribir)
	for windos in possible_windows:
		windos = str(windos)
		escribir = "*"+result[2] +"\t" + str(int(windos)-size/2)+ "\t" + str(int(windos)+size/2)+"\t" + result[3][windos] + "\n"
		r.write(escribir)  
		
#The calling of the different functions depends on whether we are working in an outcross or backcross
elif mode == "back":
	#Call of the function chromosome by chromosome, the final result comes from the chromosome with the higher parameters
	for chromosome in ch: 
		inpu = open(args.input, "r")
		genome = getinfo(chromosome, inpu)		
		windows = chromosomal_position(size, space, genome,chromosome, ch[chromosome], mode)    
		x_value= union_points(windows)
		if z ==0 : 								
			result = data_backcross(windows, x_value, chromosome, "T", size)
			best_parameter = result[1]
		else:
			result = data_backcross(windows, x_value, chromosome, best_parameter, size)
		z +=1
	#Creation of a conservative bigger window which will contain different windows with an average bigger than a threshold. 	
	candidate_list = []
	for window in result[-1]:
		if float(result[-1][window][0]) >= 0.9:
			candidate_list.append(float(window))			
	min_i = min(candidate_list)- size/2
	maxi_d = max(candidate_list) + size/2
	r= open(output ,"a")
	r.write("!"+str(result[0][0])+"\t"+str(result[0][1]) + "\t"+ str(result[1])+"\t"+str(result[2]) + "\n") 
	r.write("?"+result[2] +"\t"+ min_i+"\t"+ maxi_i +"\n")
	for widos in candidate_list:
		r.write("*"+result[2] +"\t" + int(windos)-size/2 +"\t"+ int(windos)+size/2+"\t" + result[-1][windos]+ "\n")