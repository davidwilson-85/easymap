import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-a', action="store", dest = 'input')
parser.add_argument('-b', action="store", dest = 'output')
parser.add_argument('-chr', action="store", dest = 'chr', default = '*', nargs='+')
parser.add_argument('-mut_type', action="store", dest = 'mut_type', default = 'all') # / EMS
parser.add_argument('-qual_min', action="store", dest = 'qual_min', default = 0)
parser.add_argument('-dp_min', action="store", dest = 'dp_min', default = 0)
parser.add_argument('-dp_max', action="store", dest = 'dp_max', default = 8000)
parser.add_argument('-af_min', action="store", dest = 'af_min', default = 0)
parser.add_argument('-af_max', action="store", dest = 'af_max', default = 2)
parser.add_argument('-pos_min', action="store", dest = 'pos_min', default = 0)
parser.add_argument('-pos_max', action="store", dest = 'pos_max', default = 1000000000)

parser.add_argument('-step', action="store", dest = 'step')
parser.add_argument('-cand_reg_file', action="store", dest = 'cand_reg_file')

args = parser.parse_args()


#Input 
input = args.input
f1 = open(input, 'r')
lines = f1.readlines()	

#Output
output = args.output
f2 = open(output, 'w')


#_________________________________CANDIDATE REGION FILTER___________________________________________________________________________________
step = args.step
if step == '1':
	pass
elif step == '2':
	f3 = open(args.cand_reg_file, 'r')
	f3lines = f3.readlines()
	for i, line in enumerate(f3lines):
		if line.startswith('?'):
			sp = line.split()
			args.chr = sp[1].strip('>')
			args.pos_min = int(sp[2].strip())
			args.pos_max = int(sp[3].strip())
#__________________________________________________________________________________________________________________________________________

def limits():
	global selector

	selector = 0
	
	if (
			(int(args.pos_min) < int(sp[1].strip()) < int(args.pos_max))
			and float(sp[4].strip()) > float(args.qual_min) 
			and int(args.dp_min) < (int(sp[6].strip()) + int(sp[5].strip())) < int(args.dp_max)
			and float(args.af_min) < ((float(sp[6].strip())/ ((float(sp[6].strip())) + (float(sp[5].strip()))))) < float(args.af_max)
		
	   ):
	   		selector = 1
	   		
	else:
		selector = 0
	
	return [selector]


#__________________________________________________________________________________________________________________________________________



chromosome = args.chr

for i, line in enumerate(lines):
	if not line.startswith('#'):
		sp = line.split('\t')
		if args.mut_type.strip() == 'EMS' and ((str(sp[0].strip())  in chromosome) or (chromosome[0] == '*')): 
			limits()
			ref_b = sp[2]
			alt_b = sp[3]
			if (
					(selector == 1)
					and ((ref_b.strip() == 'G' and alt_b.strip() == 'A')
					or (ref_b.strip() == 'C' and alt_b.strip() == 'T'))
				):
					f2.write(line)
					
		elif args.mut_type.strip() == 'all' and ((str(sp[0].strip())  in chromosome) or (chromosome[0] == '*')):
			limits()
			if selector == 1: 
				f2.write(line)	
		
		
'''
if step == '3':
	for i, line in enumerate(lines)
		if not line.startswith('#'):
			sp = line.split('\t')
			if sp[] ...
'''