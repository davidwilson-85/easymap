
import argparse

#We create the input and output objects
parser = argparse.ArgumentParser()
parser.add_argument('-a', action="store", dest = 'input')
parser.add_argument('-b', action="store", dest = 'output')
parser.add_argument('-m', action="store", dest = 'mode')

args = parser.parse_args()

i = 1
f1 = open(str(args.input), 'r')
f2 = open(str(args.output), 'w')

if args.mode == "split":
	for line in f1:
		if i == 1:
			name = str(line).strip('\n')

		if i == 2: 
			seq1= str(line)[0:126].strip('\n')
			seq2= str(line)[126:252].strip('\n')

		if i == 3:
			pass

		if i == 4:
			qual1= str(line)[0:126].strip('\n')
			qual2= str(line)[126:252].strip('\n')

			f2.write(name+'.1'+  '\n' + seq1 + '\n' + '+' + '\n' + qual1 + '\n' + name+'.2'  + '\n' + seq2 + '\n' + '+' + '\n' + qual2 + '\n')
			i = 0

		i = i + 1


if args.mode == "trim":
	for line in f1:
		if i == 1:
			name = str(line).strip('\n')

		if i == 2: 
			seq1= str(line)[6:51].strip('\n')

		if i == 3:
			pass

		if i == 4:
			qual1= str(line)[6:51].strip('\n')

			f2.write(name+'.1'+  '\n' + seq1 + '\n' + '+' + '\n' + qual1 + '\n')
			i = 0

		i = i + 1
