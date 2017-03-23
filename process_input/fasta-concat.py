#!/usr/bin/python

#
# This script scans a folder that only contains fasta files, reads the content of each one
# and writes it to a single fasta-formatted file. I use this script as helper script for
# the mutation analyzer.
#

import argparse, os, shutil

# Process command arguments
parser = argparse.ArgumentParser()
parser.add_argument('-in_dir', action="store", dest='input_dir')
parser.add_argument('-out_dir', action="store", dest='output_dir')
args = parser.parse_args()

input_dir = args.input_dir
output_dir = args.output_dir
# Create a subdirectory to place the reads. If it already exists, remove it first
if os.path.exists(output_dir):
	shutil.rmtree(output_dir)
os.makedirs(output_dir)



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

# Function to divide a long string ('data') into chunks of defined length ('batch_size')
def batch_gen(data, batch_size):
	for i in range(0, len(data), batch_size):
		yield data[i:i+batch_size]
		
	
# Create list of files in input dir
input_files = sorted(os.listdir('./' + input_dir))

# Create and open output file
output = open(output_dir + '/genome.fa', 'w')

# Get the content of each file and append it to output fasta file
for input_file in input_files:
	with open(input_dir + '/' + input_file) as fp:
		for name_contig, seq_contig in read_fasta(fp):
			split_name_contig = name_contig.split(' ')
			output.write(split_name_contig[0] + '\n')
			
			# Write to file a small chunk of the contig sequence in each line
			for chunk in batch_gen(seq_contig, 80):
				output.write(chunk + '\n')

# Close output file
output.close()
