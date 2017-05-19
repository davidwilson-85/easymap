# python ./report.py -variants ./files/variants.txt -log ./files/log.log -output_html ./report.html -project Report_development -sorted_insertions ./files/sorted_insertions.txt
#from __future__ import division
import argparse

parser = argparse.ArgumentParser()
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
		if line.startswith('ref_seqs:'):
			sp = line.split()
			ref_files = str(sp[-1])

		if line.startswith('ins_seq:'):
			sp = line.split()
			ins_file = str(sp[-1])

		if line.startswith('lib_type_sample'):
			sp = line.split()
			reads_type = str(sp[-1])

		if line.startswith('read_s:'):
			sp = line.split()
			reads_s = str(sp[-1])

		if line.startswith('read_f:'):
			sp = line.split()
			reads_f = str(sp[-1])

		if line.startswith('read_r:'):
			sp = line.split()
			reads_r = str(sp[-1])

		if line.startswith('gff_file:'):
			sp = line.split()
			gff_file = str(sp[-1])

		if line.startswith('ann_file:'):
			sp = line.split()
			ann_file = str(sp[-1])

		if line.startswith('data_source:'):
			sp = line.split()
			data_source = str(sp[-1])

		if line.startswith('sim_mut:'):
			sp = line.split()
			d = str(sp[-1])
			sp2 = d.split('+')
			ins_sim = str(sp2[0])

if reads_type == "se":
	reads_type = 'Single end reads'
elif reads_type == "pe":
	reads_type = 'Paired end reads'
if ann_file == "n/p":
	ann_file = "Not provided"


#______________________________________________Writting header____________________________________________________

#Header and run summary
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
'	</style>' + '\n'
'</head>' + '\n'

'<body>' + '\n'
'	<div id="wrapper">' + '\n'
'		<hr class="easymap">' + '\n'
'		<h1>Poject: ' +  project_name + '</h1>' + '\n'
'		<hr class="easymap">' + '\n'
'		<h2>Run summary</h2>' + '\n'

'		<table border="0">' + '\n'
'		<tr>' + '\n'
'		<td> <b>Reference genome file: <b> </td>' + '\n'
'		<td>&nbsp;&nbsp;&nbsp;&nbsp;</td>' + '\n'
'		<td>' +  ref_files +  '</td>' + '\n'
'		</tr>' + '\n'

'		<td> <b>Insertion sequence file:<b> </td>' + '\n'
'		<td>&nbsp;&nbsp;&nbsp;&nbsp;</td>' + '\n'
'		<td>' +  ins_file +  '</td>' + '\n'
'		</tr>' + '\n'
)

#Exp/sim and read files
if data_source == 'exp':

	if reads_type == 'Paired end reads':
		output.write(

'		<td> <b>Reads type: <b> </td>' + '\n'
'		<td>&nbsp;&nbsp;&nbsp;&nbsp;</td>' + '\n'
'		<td>' +  reads_type +  '</td>' + '\n'
'		<td>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td>' + '\n'

'		<td> <b>Forward reads file: <b> </td>' + '\n'
'		<td>&nbsp;&nbsp;&nbsp;&nbsp;</td>' + '\n'
'		<td>' +  reads_f +  '</td>' + '\n'
'		<td>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td>' + '\n'

'		<td> <b>Reverse reads file: <b> </td>' + '\n'
'		<td>&nbsp;&nbsp;&nbsp;&nbsp;</td>' + '\n'
'		<td>' +  reads_r +  '</td>' + '\n'
'		</tr>' + '\n'
		)

	elif reads_type == 'Single end reads':
		output.write(


		'		<td> <b>Reads type: <b> </td>' + '\n'
		'		<td>&nbsp;&nbsp;&nbsp;&nbsp;</td>' + '\n'
		'		<td>' +  reads_type +  '</td>' + '\n'
		'		</tr>' + '\n'

		'		<td> <b>Reads file: <b> </td>' + '\n'
		'		<td>&nbsp;&nbsp;&nbsp;&nbsp;</td>' + '\n'
		'		<td>' +  reads_s +  '</td>' + '\n'
		'		</tr>' + '\n'
		)

if data_source == 'sim':
	output.write(
	'		<td> <b>Simulated insertions: <b> </td>' + '\n'
	'		<td>&nbsp;&nbsp;&nbsp;&nbsp;</td>' + '\n'
	'		<td>' +  ins_sim +  '</td>' + '\n'
	'		</tr>' + '\n'

	)

#Gff and ann files
output.write(
'		<td> <b>Genome features file: <b> </td>' + '\n'
'		<td>&nbsp;&nbsp;&nbsp;&nbsp;</td>' + '\n'
'		<td>' +  gff_file +  '</td>' + '\n'
'		</tr>' + '\n'

'		<td> <b>Genome anotation file: <b> </td>' + '\n'
'		<td>&nbsp;&nbsp;&nbsp;&nbsp;</td>' + '\n'
'		<td>' +  ann_file +  '</td>' + '\n'
'		</tr>' + '\n'
'		</table>' + '\n'
'		<hr class="easymap">' + '\n'
)

#Genome overview
output.write(
'		<div id="wrapper">' + '\n'
'		<h2>Genome overview</h2>' + '\n'
'		<center> <img src="files/insertions_overview.png" align="middle" >  </center>' + '\n'
'		<hr class="easymap">' + '\n'

)

#Insertions overview
output.write(
'		<h2>Insertions overview</h2>' + '\n'
'		<table border="0" align="center">' + '\n'
'		<tr>' + '\n'
'		<td> <b>Insertion<b> </td>' + '\n'
'		<td>&nbsp;&nbsp;&nbsp;&nbsp;</td>' + '\n'
'		<td> <b>Chromosome<b> </td>' + '\n'
'		<td>&nbsp;&nbsp;&nbsp;&nbsp;</td>' + '\n'
'		<td> <b>Position<b> </td>' + '\n'
'		<td>&nbsp;&nbsp;&nbsp;&nbsp;</td>' + '\n'
'		<td> <b>Gene hit<b> </td>' + '\n'
'		<td>&nbsp;&nbsp;&nbsp;&nbsp;</td>' + '\n'
'		<td> <b>Gene element<b> </td>' + '\n'
'		</tr>' + '\n'
)

with open(input_var) as f:
	for line in f:
		#default values
		t_gene_hit = 'No gene hit'
		if not line.startswith('@'):
			sp = line.split()
			for ins in insertions_pos_list:
				if int(ins[1]) == int(sp[2]):
					t_insertion = str(ins[0])
					t_chr = str(sp[1])
					t_pos = str(ins[1])
					if str(sp[9]) != '-':
						t_gene_hit = str(sp[9])
					t_gene_element = str(sp[10])

			output.write(
			'		<tr>' + '\n'
			'		<td> '+t_insertion+' </td>' + '\n'
			'		<td>&nbsp;&nbsp;&nbsp;&nbsp;</td>' + '\n'
			'		<td> '+t_chr+' </td>' + '\n'
			'		<td>&nbsp;&nbsp;&nbsp;&nbsp;</td>' + '\n'
			'		<td> '+t_pos+' </td>' + '\n'
			'		<td>&nbsp;&nbsp;&nbsp;&nbsp;</td>' + '\n'
			'		<td> '+t_gene_hit+' </td>' + '\n'
			'		<td>&nbsp;&nbsp;&nbsp;&nbsp;</td>' + '\n'
			'		<td> '+t_gene_element+' </td>' + '\n'
			'		</tr>' + '\n'
			)

output.write(
'		</table>' + '\n'
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