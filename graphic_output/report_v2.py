#	python ./report_v2.py -variants ./files/variants2.txt -log ./files/log.log -output_html ./report.html -project Report_development -sorted_insertions ./files/sorted_insertions.txt -mut_type lin#from __future__ import division

import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-mut_type', action="store", dest = 'mut_type')
parser.add_argument('-variants', action="store", dest = 'variants')
parser.add_argument('-sorted_insertions', action="store", dest = 'sorted_insertions')
parser.add_argument('-log', action="store", dest = 'log')
parser.add_argument('-output_html', action="store", dest = 'output_html')
parser.add_argument('-project', action="store", dest = 'project')

args = parser.parse_args()

project_name = args.project

#Log input
input_log = args.log

#VAriants input
input_var = args.variants


#HTML output
output_html = args.output_html
output = open(output_html, 'w')

#Others
mut_type = args.mut_type

#______________________________________________List output files__________________________________________________
from os import listdir
from os.path import isfile, join
files = [f for f in listdir('./files/') if isfile(join('./files/', f))]


#______________________________________________Get list of insertions_____________________________________________
insertions_list = list()
with open(args.sorted_insertions) as f:
	for line in f:
		if not line.startswith('@'):
			sp = line.split()
			if sp[2].strip() not in insertions_list:
				insertions_list.append(sp[2].strip())

insertions_pos_list = list()
with open(input_var) as f1:
	for line in f1:
		if not line.startswith('@'):
			sp = line.split()
			ins_localizer = str(sp[1].strip().lower() + '-' + sp[2].strip())
			with open(args.sorted_insertions) as f2:
				for line in f2:
					if not line.startswith('@'):
						sp2 = line.split()
						ins_localizer_2 = str(sp2[1].strip().lower() + '-' + sp2[3].strip())
						if ins_localizer == ins_localizer_2:
							for insertion in insertions_list:
								if insertion == sp2[2]:
									sublist = [insertion, str(sp2[3]), sp[1]]
									if sublist not in insertions_pos_list:
										insertions_pos_list.append(sublist)


#______________________________________________Log info___________________________________________________________
with open(input_log, 'r') as f1:
	for line in f1:
		if line.startswith('Reference sequence:'):
			sp = line.split()
			ref_files = str(sp[-1])

		if line.startswith('Insertion sequence:'):
			sp = line.split()
			ins_file = str(sp[-1])

		if line.startswith('Library type (problem sample):'):
			sp = line.split()
			reads_type = str(sp[-1])

		if line.startswith('Single-end reads (problem sample):'):
			sp = line.split()
			reads_s = str(sp[-1])

		if line.startswith('Forward reads (problem sample):'):
			sp = line.split()
			reads_f = str(sp[-1])

		if line.startswith('Reverse reads (problem sample):'):
			sp = line.split()
			reads_r = str(sp[-1])

		if line.startswith('Library type (control sample):'):
			sp = line.split()
			reads_type_control = str(sp[-1])

		if line.startswith('Single-end reads (control sample):'):
			sp = line.split()
			reads_s_control = str(sp[-1])

		if line.startswith('Forward reads (control sample):'):
			sp = line.split()
			reads_f_control = str(sp[-1])

		if line.startswith('Reverse reads (control sample):'):
			sp = line.split()
			reads_r_control = str(sp[-1])

		if line.startswith('GFF file:'):
			sp = line.split()
			gff_file = str(sp[-1])

		if line.startswith('Annotation file:'):
			sp = line.split()
			if str(sp[-1]) == "n/p": ann_file = "Not provided"

		if line.startswith('Data source:'):
			sp = line.split()
			data_source = str(sp[-1])

		if line.startswith('Type of cross [bc/oc]:'):
			sp = line.split()
			if str(sp[-1]) == 'bc': cross_type = 'backcross'
			if str(sp[-1]) == 'oc': cross_type = 'outcross'

		if line.startswith('Mutant strain [ref/noref]:'):
			sp = line.split()
			if str(sp[-1]) == 'ref': mut_background = 'reference'
			if str(sp[-1]) == 'noref': mut_background = 'non-reference'

		if line.startswith('SNP analysis type [par/f2wt]:'):
			sp = line.split()
			if str(sp[-1]) == 'par': snp_analysis_type = 'parental'
			if str(sp[-1]) == 'f2wt': snp_analysis_type = 'wild type F2'

		if line.startswith('Parental used as control [mut/nomut/np]:'):
			sp = line.split()
			if str(sp[-1]) == 'mut': parental_used_as_control = 'mutant'
			if str(sp[-1]) == 'nomut': parental_used_as_control = 'wild type'

#SNP mappint control samples
if snp_analysis_type == 'parental' and parental_used_as_control == 'mutant' : control = ' parental of the mutant strain.'
if snp_analysis_type == 'parental' and parental_used_as_control == 'wild type' : control = ' wild type parental of the mapping crossing.'
if snp_analysis_type == 'wild type F2' : control = ' wild type F2 of the mapping cross.'

if data_source == 'sim':
	with open(input_log, 'r') as f1:
		for line in f1:
			if line.startswith('Simulator (sim-mut.py) command:'):
				sp = line.split()
				number_mutations = str(sp[-1].split('+')[0])

			if line.startswith('Simulator (sim-seq.py) command:'):
				sp = line.split()
				read_depth = str(sp[-1].split('+')[0])



#______________________________________________Writting header____________________________________________________

#HTML/CSS stuff
output.write(
'<!DOCTYPE html>' + '\n'
'<html>' + '\n'
'<head>' + '\n'

'	<meta charset="utf-8" />' + '\n'
'	<title>Easymap - report</title>' + '\n'
'	<style>' + '\n'
'		#wrapper {' + '\n'
'			width: 100%;' + '\n'
'			max-width: 1050px;' + '\n'
'			margin: 0 auto;' + '\n'
'			font-family: arial, helvetica;' + '\n'
'		}	' + '\n'
'		.easymap {' + '\n'
'			border: 0;' + '\n'
'			color: rgb(139, 167, 214);' + '\n'
'			background-color: rgb(139, 167, 214);' + '\n'
'			height: 5px;' + '\n'
		
'		}	' + '\n'
'		.img { width: 100%; }' + '\n'
'		.result1 { max-width: 864px; }' + '\n'
'		.result2 { max-width: 1000px; }' + '\n'


'		table {border-collapse:collapse; table-layout:fixed; width:310px;}' + '\n'
'		table td {border:solid 0px #fab; width:100px; word-wrap:break-word; vertical-align:top;}' + '\n'
   
'		#t {  border: 0px solid red; word-wrap:break-word; table-layout:fixed; }' + '\n'
'	</style>' + '\n'
'</head>' + '\n'
)

#Header and run summary
output.write(
'<body>' + '\n'
'	<div id="wrapper">' + '\n'
'		<hr class="easymap">' + '\n'
'		<h1>Poject: ' +  project_name + '</h1>' + '\n'
'		<hr class="easymap">' + '\n'

)

#Exp/sim and read files
output.write(
'		<h2>Run summary</h2>' + '\n'
'		<table id="t">' + '\n'
'		<col width="300">' + '\n'
'		<col width="700">' + '\n'
	)

output.write(

'		<tr>' + '\n'
'			<td> <b>Input genome files:</b></td>' + '\n'
'			<td>' + ref_files + '</td>' + '\n'
'		</tr>' + '\n'

	)


#Read files 

if mut_type == 'lin' and data_source == 'exp': 
	if reads_type == 'pe':
		output.write(
'		<tr>' + '\n'
'			<td> <b>Input reads files:</b></td>' + '\n'
'			<td>' + reads_f + ', &nbsp;&nbsp;&nbsp;&nbsp;' + reads_r + '</td>' + '\n'
'		</tr>' + '\n'
		)

	elif reads_type == 'se':
		output.write(
'		<tr>' + '\n'
'			<td> <b>Input reads files:</b></td>' + '\n'
'			<td>' + reads_s + '</td>' + '\n'
'		</tr>' + '\n'
			)


if mut_type == 'snp' and data_source == 'exp': 
	if reads_type == 'pe':
		output.write(
'		<tr>' + '\n'
'			<td> <b>Input problem reads files:</b></td>' + '\n'
'			<td>' + reads_f + ', &nbsp;&nbsp;&nbsp;&nbsp;' + reads_r + '</td>' + '\n'
'		</tr>' + '\n'
		)

	elif reads_type == 'se':
		output.write(
'		<tr>' + '\n'
'			<td> <b>Input problem reads files:</b></td>' + '\n'
'			<td>' + reads_s + '</td>' + '\n'
'		</tr>' + '\n'
			)

	if reads_type_control == 'pe':
		output.write(
'		<tr>' + '\n'
'			<td> <b>Input control reads files:</b></td>' + '\n'
'			<td>' + reads_f_control + ', &nbsp;&nbsp;&nbsp;&nbsp;' + reads_r_control + '</td>' + '\n'
'		</tr>' + '\n'
		)

	elif reads_type_control == 'se':
		output.write(
'		<tr>' + '\n'
'			<td> <b>Input control reads files:</b></td>' + '\n'
'			<td>' + reads_s_control + '</td>' + '\n'
'		</tr>' + '\n'
			)




if data_source == 'sim':
	if mut_type == 'lin':
		output.write(
'		<tr>' + '\n'
'			<td> <b>Number of insertions (simulation)</b></td>' + '\n'
'			<td>' + number_mutations + '</td>' + '\n'
'		</tr>' + '\n'
'		<tr>' + '\n'
'			<td> <b>Read depth (simulation):</b></td>' + '\n'
'			<td>' + read_depth + '</td>' + '\n'
'		</tr>' + '\n'
		)

	elif mut_type == 'snp':
		output.write(
'		<tr>' + '\n'
'			<td> <b>Number of mutations (simulation)</b></td>' + '\n'
'			<td>' + number_mutations + '</td>' + '\n'
'		</tr>' + '\n'
'		<tr>' + '\n'
'			<td> <b>Read depth (simulation):</b></td>' + '\n'
'			<td>' + read_depth + 'x</td>' + '\n'
'		</tr>' + '\n'
		)
		

if mut_type == 'snp': 
	output.write(

'		<tr>' + '\n'
'			<td> <b>Experimental design: </b></td>' + '\n'
'			<td> Mutation in ' + mut_background + ' genetic background. A ' + cross_type + ' was performed to obtain the mapping population. The control sample is the ' + control +'</td>' + '\n'
'		</tr>' + '\n'

		)

#Gff and ann files
output.write(

'		<tr>' + '\n'
'			<td> <b>Structural anotation file:</b></td>' + '\n'
'			<td>' + gff_file + '</td>' + '\n'
'		</tr>' + '\n'

'		<tr>' + '\n'
'			<td> <b>Functional anotation file:</b></td>' + '\n'
'			<td>' + ann_file + '</td>' + '\n'
'		</tr>' + '\n'
'		</table>' + '\n'

)

#Link to log file
output.write(
'		<a href= ' + input_log + ' target="_blank">Cick to see log file</a>' + '\n'
'		<hr class="easymap">' + '\n'
)


#Input data quality assessment
if data_source == 'exp':
	output.write(
	'		<h2>Input data quality assessment</h2>' + '\n'
		)

	if reads_type == 'pe':
		output.write(

'		<b>Paired end reads quality assessment<br></b>' + '\n'
		)

	if reads_type == 'se':
		output.write(
'		<b>Single end reads quality assessment<br></b>' + '\n'
		)

	output.write(
'		<b>Read depth distribution<br></b>' + '\n'
'		<hr class="easymap">' + '\n'
		)


#__________________________________LIN cartographic report________________________________________________________________
if mut_type == 'lin': 
	#Genome overview
	output.write(
	'		<h2>Genomic overview</h2>' + '\n'
	'		<center> <img src="files/insertions_overview.png" align="middle" >  </center>' + '\n'
	#Table
	'		<center><b>Table goes here<br></b></center>' + '\n'
	'		<hr class="easymap">' + '\n'
	)

	#Insertions
	for ins in insertions_list:
		for f in sorted(files): 
			if '_ins_' + ins[0] in str(f):
				output.write(
				'		<h2> Insertion   ' +  str(ins[0]) +'</h2>' + '\n'
				'		<center> <img src="files/'  +  str(f)  + ' " align="middle" >  </center>' + '\n'
				)
		for f in sorted(files):
			if '_lin_' + ins[0] in str(f):
				output.write(
				'		<center> <img src="files/'  +  str(f)  + ' " align="middle" >  </center>' + '\n'
				)


if mut_type == 'snp':

	#Chromosomes FA vs POS
	output.write(
	'		<h2>Genomic overview</h2>' + '\n'
		) 
	for f in sorted(files):
		if 'img_2_mapping' in str(f):
			output.write(
			'		<left> <img src="files/'  +  str(f)  + ' " align="middle" >  </left>' + '\n'
			)

	#Candidates table:
	output.write(
		#Table
		'		<center><b>Table goes here<br></b></center>' + '\n'
		'		<hr class="easymap">' + '\n'
		)

	#Candidate SNPs
	output.write(
	'		<h2>Candidate genes</h2>' + '\n'
		) 

	for f in sorted(files):
		if 'gene_plot_snp' in str(f):
			sp = str(f).split('_')
			pos = str(sp[3])
			gene = str(sp[5].split('.')[0])
			output.write(
			'		<h3>' + gene + '</h3>' + '\n'
			'		<left> <img src="files/'  +  str(f)  + ' " align="middle" >  </left>' + '\n'
			)






















'''
#_____________________________HTML file


output_html = args.output_html
f3 = open(output_html, 'w')

f3.write('<html>' + '\n')

f3.write('<br><h3>Insertions overview</h3>' + '\n')
f3.write('<br>' + '\n')
f3.write('<img src="insertions_overview.png"></img>' + '\n')
f3.write('<br>' + '\n')

for e in insertions: 
	f3.write('<br><h3>Insertion ' + e + '</h3>' + '\n')
	f3.write('<table border="1">' + '\n')
	f3.write('	<tr>' + '\n')
	f3.write('		<td><img src="local_' + e + '.png"></img></td>' + '\n')
	f3.write('	</tr>' + '\n')
	f3.write('</table>' + '\n')


f3.write('<br><br>Here is your processed data:<br>' + '\n')
f3.write('<a href="sorted_insertions.txt" target="_blank">Sorted insertions</a>' + '\n')

f3.write('</html>')
'''