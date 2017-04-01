import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-f2_mut', action="store", dest = 'input_mut')
parser.add_argument('-f2_wt', action="store", dest = 'input_wt')
parser.add_argument('-out', action="store", dest = 'output')

args = parser.parse_args()


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

for i, line in enumerate(mut_lines):
	if not line.startswith('#'):
		line.strip('\n')
		sp = line.split('\t')
		fa_mut = float(float(sp[6])/(float(sp[6])+float(sp[5])))
		fa_mut = int(fa_mut * 1000)
		for j, line2 in enumerate(wt_lines):
			line2.strip('\n')
			if not line2.startswith('#'):
				sp2 = line2.split('\t')
				fa_wt = float(float(sp2[6])/(float(sp2[6])+float(sp2[5])))
				fa_wt = int(fa_wt * 1000)
				if (sp2[0].strip()+'-'+sp2[1].strip()) == (sp[0].strip()+'-'+sp[1].strip()):
					if not fa_mut in range((fa_wt - 50), (fa_wt + 50)): 
						f3.write(line.strip('\n') + '\t' + sp2[5] + '\t' + sp2[6])