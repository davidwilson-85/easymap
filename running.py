#!/usr/bin/env python
#mutated_positions = [5925329,5864788,5807566,5769043,5637307,5622262,5497045,5405459,5085973,5067118]
import os, subprocess,shutil
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-mut', action="store", dest = 'mutations', required = "True")

args = parser.parse_args()
mutations = args.mutations
mutated_positions= mutations.split(",")



	
def purge():
	for folders in os.walk("./user_projects"):
		for z in folders[1]:
			path = "./user_projects/"+z+"/1_intermediate_files"
			try:
				shutil.rmtree(path+"/sim_data/sim_recsel_output_recessive")
				shutil.rmtree(path +"/sim_data/sim_seq_output")
				shutil.rmtree(path +"/sim_data/sim_seq_output")
			except:
				i = "nothing"
			for i in os.walk(path):
	
				for j in i[2]:
					if j[-4:] == ".sam" or j[-4] == ".bai":
						os.remove(path+"/"+j)

#os.remove("./user_projects/*/1_intermediate_files/*.sam")
#os.mdir("./user_projects/*/1_intermediate_files/sim_data")
#subprocess.call(["rm", "-f" ,"./user_projects/*/1_intermediate_files/*.sam"])
#subprocess.call(["rm", "-rf" ,"./user_projects/*/1_intermediate_files/sim_data"


for i in range(len(mutated_positions)):
	a = str(mutated_positions[i])
	#Case 3

	os.system("./easymap -n take_caso_3-"+str(i+11)+" -w snp -sim -r at4 -g chr1+4.gff -rs ref -cr oc -co par_mut -sm 80+e -sr 0,16-1,34-2,31-3,14-4,4-5,1+1,"+a+"+r+100 -ss 25+100,0+500,100+1+50+se")

	subprocess.call(["rm", "-rf" ,"./user_projects/*/1_intermediate_files/*.sam"])
	subprocess.call(["rm", "-rf" ,"./user_projects/*/1_intermediate_files/*.bam"])
	subprocess.call(["rm", "-rf" ,"./user_projects/*/1_intermediate_files/sim_data"])
	
	purge()
	exit()

	#Case 4

	os.system("./easymap -n take_caso_4"+str(i+11)+" -w snp -sim -r at4 -g chr1+4.gff -rs ref -cr oc -co par_nomut -sm 80+e -sr 0,16-1,34-2,31-3,14-4,4-5,1+1,"+a+"+r+100 -ss 25+100,0+500,100+1+50+se")

	subprocess.call(["rm", "-rf" ,"./user_projects/*/1_intermediate_files/*.sam"])
	subprocess.call(["rm", "-rf" ,"./user_projects/*/1_intermediate_files/*.bam"])
	subprocess.call(["rm", "-rf" ,"./user_projects/*/1_intermediate_files/sim_data"])
	purge()
	exit()
	#Case 6, former 7 

	os.system("./easymap -n take_caso_6"+str(i+11)+" -w snp -sim -r at4 -g chr1+4.gff -rs noref -cr oc -co par_mut -sm 80+e -sr 0,16-1,34-2,31-3,14-4,4-5,1+1,"+a+"+r+100 -ss 25+100,0+500,100+1+50+se")


	subprocess.call(["rm", "-rf" ,"./user_projects/*/1_intermediate_files/*.sam"])
	subprocess.call(["rm", "-rf" ,"./user_projects/*/1_intermediate_files/*.bam"])
	subprocess.call(["rm", "-rf" ,"./user_projects/*/1_intermediate_files/sim_data"])
	purge()