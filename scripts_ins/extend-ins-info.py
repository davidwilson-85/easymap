#!/usr/bin/python

# This scrip retrieves the flanking sequences of the insertions 
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--project-name', action="store", dest='project_name', required=True)
args = parser.parse_args()
project = args.project_name

# Input files: fasta genome and insertions_output.txt
fp = open(project + '/1_intermediate_files/gnm_ref_merged/genome.fa', 'r')
input_file = open(project + '/3_workflow_output/insertions_output.txt', 'r')

# This function parses the information from the fasta files
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

# We create a list of all the contigs with the format [[contig_name, sequence], [contig_name, sequence] ...]
fastalist = list()
for name_contig, seq_contig in read_fasta(fp):
	fastalist.append([name_contig.lower(), seq_contig])

# We retrieve the upstream and downstream sequences of each insertion from the fastalist. We also create a new list with the complete lines
final_lines = list()
for line in input_file:
	if not line.startswith('@'):
		sp = line.split()
		chromosome = str(sp[1])
		position = str(sp[2])
		for chrm in fastalist:
			if chrm[0].strip() == '>'+chromosome.strip().lower():
				upstream =  chrm[1][int(position)-51:int(position)-1]
				downstream =  chrm[1][int(position):int(position)+50]
				line = line.strip('\n')
				line = line + '\t' + upstream + '\t' + downstream
				final_lines.append(line)

# We re-write the file with the extended information to the output file from the final_lines list
output_file = open(project + '/3_workflow_output/insertions_output.txt', 'w')

output_file.write('@type\tcontig\tposition\tref_base\talt_base\thit\tmrna_start\tmrna_end\tstrand\tgene_model\tgene_element\taa_pos\taa_ref\taa_alt\tgene_funct_annot\tf_primer\ttm_f_primer\tinsertion_primer_5\ttm_insertion_primer_5\tinsertion_primer_3\ttm_insertion_primer_3\tr_primer\ttm_r_primer\tupstream\tdownstream\n')    
for line in final_lines:
	output_file.write(line + '\n')



