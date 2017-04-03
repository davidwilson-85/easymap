#!/usr/bin/python
import argparse
parser = argparse.ArgumentParser()
from subprocess import call
import os

parser.add_argument('--project_name', "-P", action="store", dest = 'project_name', required = "True")
parser.add_argument('--workflow', "-w" ,action="store", dest = 'workflow', choices=set(('snp','ins')), required = "True")
parser.add_argument('--data_source', "-d", action="store", dest = 'data_source', choices=set(('sim','exp')), required = "True")
parser.add_argument('--lib_type', "-l", action="store", dest = 'lib_type', choices=set(('se','pe')), required = "True")
#parser.add_argument('--ref_seqs', "-rs", action="store", dest = 'ref_seqs', required = "True") 
parser.add_argument('--ins_seq', "-i", action="store",default = 'n/p', dest = 'ins_seq') 
parser.add_argument('--read_s', "-A", action="store",default = 'n/p', dest = 'read_s')
parser.add_argument('--read_f', "-B", action="store",default = 'n/p', dest = 'read_f')
parser.add_argument('--read_r', "-C",action="store",default = 'n/p', dest = 'read_r')
parser.add_argument('--gff_file',"-g",action="store", dest = 'gff_file', required = "True")
parser.add_argument('--ann_file', "-a", action="store",default = 'n/p', dest = 'ann_file')
parser.add_argument('--sim_mut', "-sm", action="store",default = 'n/p', dest = 'sim_mut')
parser.add_argument('--sim_recsel',"-sr", action="store", default = 'n/p', dest = 'sim_recsel')
parser.add_argument('--sim_seq', "-ss", action="store", default = 'n/p', dest = 'sim_seq')
parser.add_argument ("--snp_analysis_type","-s", action="store", default = "n/p",dest ="snp_analysis_type", choices=set(("par","f2wt")))
parser.add_argument('--read_s_control',"-D", action="store",default = 'n/p', dest = 'read_s_control')
parser.add_argument('--read_f_control', "-E", action="store",default = 'n/p', dest = 'read_f_control')
parser.add_argument('--read_r_control', "-F", action="store",default = 'n/p', dest = 'read_r_control')
parser.add_argument('--cross_type', "-c", action="store",default = 'n/p', dest = 'cross_type', choices=set(('oc','bc')))
parser.add_argument('--is_ref_strain',"-r", action="store",default = 'n/p', dest = 'is_ref_strain', choices=set(('ref','noref'))) 
parser.add_argument('--parental_used_as_control', "-p", action="store",default = 'n/p', dest = 'parental_used_as_control', choices=set(('mut','nomut')))


args = parser.parse_args()

project_name = args.project_name
workflow = args.workflow
data_source = args.data_source
lib_type = args.lib_type
ins_seq = args.ins_seq
read_s = args.read_s
read_f= args.read_f
read_r= args.read_r
gff_file = args.gff_file
ann_file = args.ann_file
sim_mut = args.sim_mut
sim_recsel = args.sim_recsel
sim_seq = args.sim_seq
snp_analysis_type = args.snp_analysis_type #Which control is being used: parental one or F2wt?
read_s_control = args.read_s_control
read_f_control = args.read_f_control
read_r_control = args.read_r_control
cross_type = args.cross_type
is_ref_strain = args.is_ref_strain
parental_used_as_control = args.parental_used_as_control


problems = []
error = 0


# Check whether the user has provided a project_name and it has no space:  
is_there_space = project_name.split(" ")
if len(is_there_space) != 1:
	error = 1
	problems.append("Please select a project name without any space")
#Check whether user_data directory exists
input_folder = "./user_data"
if os.path.isdir(input_folder) == False:
	error = 1
	problems.append("Easymap could not find "+ input_folder + ". Please create the folder and place your input files inside. Please mind reference genome must be provided inside the path ./user_data/genome_ref")
else:
	#Check whether there are reference file/s in 0_input/gnm_ref
	try:
		gnm_ref_folder = "./user_data/gnm_ref" #Looks in the path were the file/s should be found
		if not os.listdir(gnm_ref_folder): #This function creates a list of the items in the specified path. Thus, if no items are found, the result of the function will be FALSE
		    error = 1
		    problems.append("The reference file/files should be in the directory " + gnm_ref_folder)
	except:
		error = 1
		problems.append("No folder"+ gnm_ref_folder + " has been found, please create the folder and include the reference genome") 

#This script will not check the existance of user_projects because is supposed to be created in its absence by easymap.sh
input_files = os.listdir(input_folder)
#Check the existance of gff file in the appropiate folder
if gff_file not in input_files:
		error = 1
		problems.append("gff file name provided is not in " + input_folder)

#In case the user decides to give an annotation file, it will be checked whether it is in the corresponding place.
if ann_file != "n/p":
	if ann_file not in input_files:
		error =1
		problems.append("You have specify you want to use an annotation file (ann_file), but the name given has not been succesfully found in project/0_input directory")

#In the analysis mode of snp, cross_type, is_ref_strain and parental_used_as_control are required parameters
if workflow=="snp":
	if snp_analysis_type == "f2wt":
		if cross_type == "oc":
			error = 1
			problems.append("This program do not support an experimental design in which an outcross is performed and the wild type F2 is used as control. Please see documentation for more details.")
		else:
			cross_type = "bc"
			parental_used_as_control = "n/p"
			is_ref_strain = "n/p"

	if snp_analysis_type == "par":
		if cross_type == "n/p" or is_ref_strain== "n/p" or parental_used_as_control == "n/p": #Por el bien de los tontos: dividir en tres 
			error = 1
			problems.append("Snp mode requires value oc/bc in cross_type") if cross_type == "n/p"
			problems.append("Snp mode requires value ref/noref in is_ref_strain") if is_ref_strain== "n/p"
			problems.append("Snp mode requires value mut/nomut in parental_used_as_control") if parental_used_as_control == "n/p"
		if cross_type == "oc" and is_ref_strain == "noref" and parental_used_as_control == "nomut":
			error = 1
			problems.append("Unfortunatelly this software is not thought to allow the process of an outcross which is not in the reference background and the sequenced parental is not the pre-mutagenized one")
			
#If the analysis mode is insertions, a insertions sequence must be given
if workflow == "ins":
	if ins_seq not in input_files: #If ins_seq is not in input_folder it can be due to the fact ins_seq are not provided, in that case its value would be n/p, or due to the fact the name given is not located in such directory.
		error = 1
		problems.append("The mode insertion (ins) requires a insertion file in " + input_folder) #CORREGIR ERRORES, ESPECIFICAR COMO DEBE SER DADO EL ARCHIVO


#If the user is not providing its own data, there are a number of extra parameters to fullfil, related to the simulation process: $sim_mut $sim_recsel $sim_seq

    
if data_source == "sim":
	if workflow == "snp": 
		if sim_mut == "n/p" or sim_recsel == "n/p" or sim_seq == "n/p" or snp_analysis_type == "n/p":
			error = 1
			problems.append("Simulated snp data requires sim_mut parameter. Format: 40+e ; were the first number represents the number of mutations and the second character (e/d) if they are due to EMS or spontaniously") if sim_mut == "n/p"
			problems.append("Simulated snp data requires sim_recsel parameter. Format 0,14;1,31/0,24;1,42+1,10000000+r+50; were the first long value should be between "" and represents the probability of suffering recombinations (first number is the number of times that recombines it is separed from the probability by a ;) and chromosomes are separed by / . The secondo parameter is the chormosome and the position where the causal mutation will be held  . . .") if sim_recsel == "n/p"
			problems.append("Simulated snp data requires sim_seq parameter. Format 1+100,0+500,100+1+50 ") if snp_analysis_type == "n/p"
	else:
		if sim_mut == "n/p" or sim_seq == "n/p":
			error = 1
			problems.append("In order to perform the simulation mode (sim) for inserctions (ins) parameters sim_mut and sim_seq must be given, please find the information regarding then in the software documentation")
	


	#Check sim_mut ex: 40+e 
	try:
		values = sim_mut.split("+")
		if len(values) == 2:
			number = int(values[0])
			options = ["e", "d"]
			if values[1] not in options and workflow == "snp":
				error =1
				problems.append("sim_mut requires a second parameter that has to be <  e  > if the mode is EMS mutations or <  d  > if the mode is drift mutations")
		elif len(values)== 1 and workflow == "ins":
			number = int(values[0])
			sim_mut = sim_mut+"+li"

		else:	
			problems.append("sim_mut requires two/one argument/s number_of_mutations(+mode)*     *Only required in snp mode, in insertion mode it is not need ")
	except:
		error = 1
		problems.append("Please mind the format of sim_mut: number_of_mutations+mode  where number_of_mutations is any positive number and mode has to be chosen between e (EMS) and d (drift")
	if workflow == "snp":
		try: 
			values = sim_recsel.split("+")
			if len(values) == 4:
				#check of values[0]? 
				val1 = values[1].split(",")
				int(val1[0])
				int(val1[1])
				options = ["r", "d"] #MIRAR QUE VALORES PUEDEN SER
				if values[2] not in options:
					error = 1
					problems.append("Please, choose a valid selection mode in sim_recsel: < r > for recesive mutations or < d > dominant mutations.  Ex: recombination_frequences+mutant_chromosome,position_mutation+r/d")	
				int(values[3])
			else:
				error =1
				problems.append("Please, mind the format of sim_recsel: recombination_chromosome1/recombination_chromosome2....+chromosome_with_causal_mutation(number),position_causal_mutation(number)+selection_mode(d/r), ???") 
		except:
			error = 1
			problems.append("Please, mind the format of sim_recsel: recombination_chromosome1/recombination_chromosome2....+chromosome_with_causal_mutation(number),position_causal_mutation(number)+selection_mode(d/r), ???")
	#sim_seq 1+100,0+500,100+1+50
	try:
		values = sim_seq.split("+")
		val1 = values[1].split(",")
		val2 = values[2].split(",")
		int(values[0])
		int(val1[0])
		int(val1[1])
		int(val2[0])
		int(val2[1])
		int(values[3])
		int(values[4])
	except:
		error = 1
		problems.append("Please, mind the format of sim_seq: ")
		#append error


#If the user is providing its own data, reads should be provided:
if data_source == "exp":
	#If single-end reads or paired-end reads are not given
	if lib_type == "pe":
		if read_f not in input_files: 
			error = 1
			problems.append("Please introduce the files containing the forward reads in " + input_folder)
		if read_r not in input_files:
			error = 1
			problems.append("Please introduce the files containing reverse reads in " + input_folder)
		else: 
			if read_f == read_r:
				error = 1
				problems.append("Please notice both  forward and reverse reads files have the same name")
	if lib_type == "se":
		if read_s not in input_files:
			error = 1
			problems.append("Please introduce the file containing the reads in easymap/project/0_input directory")	


	#For snp mode, control reads should be provided for the further procesing. 
	if workflow == "snp":
		if lib_type == "pe":
			if read_f_control not in input_files:
				error = 1
				problems.append("Please introduce the files containing the control forward reads in " + input_folder)
			if read_r_control not in input_files:
				error = 1
				problems.append("Please introduce the files containing the control reverse reads in " + input_folder)
			else:
				if read_f_control == read_r_control:
					error = 1
					problems.append("Please notice both control read files have the same name")
		if lib_type == "se":
			if read_s not in input_files:
				error = 1
				problems.append("Please introduce the file containing the control reads in " + input_folder)	




if error == 1:
	print "Warning: one or more errors have been introduced as input" +"\n"
	for items in problems:
		print "-" + items
	quit()	


master_program_input = project_name + " " + workflow + " " + data_source + " "+ lib_type + " " + ins_seq + " " + read_s + " " + read_f + " " + read_r + " " + gff_file + " " + ann_file + " " + sim_mut + " " + sim_recsel + " " + sim_seq + " " + read_s_control + " " + read_f_control + " " + read_r_control + " " + is_ref_strain+ " "+ cross_type+ " "+ parental_used_as_control+ " "+ snp_analysis_type 

print "Succes!"
quit()
call("./easymap.sh " + master_program_input, shell=True)
