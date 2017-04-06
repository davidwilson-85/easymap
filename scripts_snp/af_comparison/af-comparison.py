import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-f2_mut', action="store", dest = 'input_mut')
parser.add_argument('-f2_wt', action="store", dest = 'input_wt')
parser.add_argument('-out', action="store", dest = 'output')

args = parser.parse_args()

import time
start_time = time.time()

#Input 
input1 = args.input_mut
f1 = open(input1, 'r')
mut_lines = f1.readlines()	

input2 = args.input_wt
f2 = open(input2, 'r')
wt_lines = f2.readlines()	

#Output
output = args.output
f3 = open(output, 'w')

f2mut = list()
f2wt = list()


for i, line in enumerate(mut_lines):
	if not line.startswith('#'):
		sp = line.split('\t')
		innerlist = [sp[0], sp[1], sp[2], sp[3], sp[4], sp[5], sp[6].strip('\n')]
		f2mut.append(innerlist)

for i, line in enumerate(wt_lines):
	if not line.startswith('#'):
		sp = line.split('\t')
		innerlist = [sp[0], sp[1], sp[2], sp[3], sp[4], sp[5], sp[6].strip('\n')]
		f2wt.append(innerlist)

counter = 0
pos = 0

for mut in f2mut: 

	pos = pos + 1

	print "line" , pos
	fa_mut = int(float(float(mut[6])/(float(mut[6])+float(mut[5])))*1000)
	c = 0
	for wt in f2wt[counter:] :
			c = c + 1
			fa_wt = int(float(float(wt[6])/(float(wt[6])+float(wt[5])))*1000)
			if (wt[0].strip()+'-'+wt[1].strip()) == (mut[0].strip()+'-'+mut[1].strip()):
				if not fa_mut in range((fa_wt - 50), (fa_wt + 50)): 
					f3.write(str(mut[0]) + '\t' + str(mut[1]) + '\t' + str(mut[2]) + '\t' + str(mut[3]) + '\t' + str(mut[4]) + '\t' + str(mut[5]) + '\t' +  str(mut[6]) + '\t' + str(wt[5]) + '\t' + str(wt[6]) + '\n')
					counter = c + counter + 1
					continue


print("--- %s seconds ---" % (time.time() - start_time))
