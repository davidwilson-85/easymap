#from __future__ import division
import argparse
from PIL import Image, ImageDraw, ImageFont, ImageOps


parser = argparse.ArgumentParser()
parser.add_argument('-a', action="store", dest = 'input')		
parser.add_argument('-b', action="store", dest = 'input_f')		#Fasta genome input
parser.add_argument('-gff', action="store", dest = 'gff')		#Genome feature file
parser.add_argument('-iva', action="store", dest = 'input_va') 	#Output de varanalyzer
parser.add_argument('-rrl', action="store", dest = 'rrl') 	#Regulatory region lenght
parser.add_argument('-f', action="store", dest = 'output_html')
parser.add_argument('-m', action="store", dest = 'mode', default = 'P')
args = parser.parse_args()

#Input 1
input = args.input
f1 = open(input, 'r')
lines = f1.readlines()	

#__________________________________________Insertions overview image_____________________________________________________________
#________________________________________________________________________________________________________________________________

#Input 2
finput = args.input_f
f2 = open(finput, 'r')
flines = f2.readlines()	
#define a superlist with innerlists, each of them containing all the info of each contig 
superlist = list()

#Create a list with all the genome contigs
contigs = []
length = 0
n = 0
#dict_contigs = dict()
lengthlist = list()

#read fasta file to determine number of contigs
for i, line in enumerate(flines):
	if line.startswith('>'): #fasta sequences start with '>'
		sp = line.split(' ')  #because some names have whitespaces and extra info that is not written to sam file
		cont = sp[0].strip()  #strip() is to remove the '\r\n' hidden chars
		cont = cont[1:]       #to remove the first char of string (>)
		if cont not in contigs:
			contigs.append(cont)
			innerlist = list()
			innerlist.append(cont)
			superlist.append(innerlist)

#Calculate the width of the image acording to the number of contigs
num_contigs = 0
for c in superlist: 
	num_contigs+=1
contigs_image_width = 200 * num_contigs

im = Image.new("RGB", (contigs_image_width + 80, 600), (255,255,255))

contig_source = args.input_f
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


# Read contig fasta file
with open(contig_source) as fp:
	fastalist = list()
	for name_contig, seq_contig in read_fasta(fp):
		innerlist = list()
		innerlist.append(name_contig.strip('>'))
		innerlist.append(len(seq_contig))
		fastalist.append(innerlist)
try:
	max_list = list()
	for c in fastalist:
		max_list.append(int(c[1]))
	max_length = max(max_list)

except:
	max_length = fastalist[0][1]


for c in superlist:
	for l in fastalist:
		if c[0] == l[0]:
			c.append(l[1])


contigs_scaling_factor = max_length / 400.0


#translate real chromosome lengths into image i and f coordinates
contig_counter = 1
contig_yi_coord = 100

for c in superlist:
	contig_x_coord = ((contigs_image_width / num_contigs) * contig_counter) - 100
	contig_yf_coord = int(contig_yi_coord + c[1]/contigs_scaling_factor)
	contig_counter +=1
	c.append(contig_x_coord)
	c.append(contig_yi_coord)
	c.append(contig_yf_coord)

#add insertions aproximate position to superlist
if args.mode == 'P':
	tag_list = list()
	for i, line in enumerate(lines):
		if not line.startswith('@'):	
			sp = line.split('\t')
			contig = str(sp[1].strip())
			insertion = str(sp[2])
			if '-' not in insertion:
				tag = contig + '-' + insertion
				if tag not in tag_list: 
					tag_list.append(tag)
	for c in superlist: 
		position_list = list()
		for t in tag_list:
			n_reads = 1
			p_reads = 0 		
			for i, line in enumerate(lines):
				if not line.startswith('@'): 
					sp = line.split('\t')
					contig = str(sp[1].strip())
					insertion = str(sp[2])
					if '-' not in insertion and sp[5].strip() == 'TOTAL' and contig == c[0].strip():
						tag2 = contig + '-' + insertion
						if tag2 == t:
							n_reads = n_reads + 1
							p_reads = p_reads + int(sp[3])											
						else:
							aprox_position = p_reads/n_reads
				aprox_position = p_reads/n_reads
			if aprox_position != 0: 
				position_list.append(aprox_position)	
		c.append(position_list)

elif args.mode == 'S':
	tag_list = list()
	for i, line in enumerate(lines):
		if not line.startswith('@'):	
			sp = line.split('\t')
			contig = str(sp[1].strip())
			insertion = str(sp[2])
			if '-' not in insertion:
				tag = contig + '-' + insertion
				if tag not in tag_list: 
					tag_list.append(tag)
	for c in superlist: 
		position_list = list()
		tag2_list = list()
		for t in tag_list:		
			for i, line in enumerate(lines):
				if not line.startswith('@'): 
					sp = line.split('\t')
					contig = str(sp[1].strip())
					insertion = str(sp[2])
					if '-' not in insertion and sp[5].strip() == 'TOTAL' and contig == c[0].strip():
						tag2 = contig + '-' + insertion
						if tag2 not in tag2_list: 
							if tag2 == t:
								position_list.append(int(sp[3]))
								tag2_list.append(tag2)
						else:
							pass
						
		c.append(position_list)


#initialize draw
draw = ImageDraw.Draw(im)
#get fonts from foler 'fonts'
fnt1 = ImageFont.truetype('fonts/arial.ttf', 30)
fnt2 = ImageFont.truetype('fonts/arial.ttf', 16)
fnt3 = ImageFont.truetype('fonts/arial.ttf', 24)
fnt4 = ImageFont.truetype('fonts/arial.ttf', 20)

tab = 50


number = 1

#Drawing the chromosomes:
for c in superlist:
	draw.line((c[2], c[3]) + (c[2], c[4]), fill=256, width=8)
	draw.ellipse((c[2]-2, c[3]-2, c[2]+4, c[3]+3), fill=256)
	draw.ellipse((c[2]-2, c[4]-2, c[2]+4, c[4]+3), fill=256)
	draw.text(((c[2] - tab), (c[3] - 50)), c[0], font=fnt1, fill=(0,0,0,255))
	for i in c[5]:
		draw.polygon([(c[2]+ 5, (i/contigs_scaling_factor+contig_yi_coord)), (c[2]+15, (i/contigs_scaling_factor+contig_yi_coord)+10), (c[2]+15, (i/contigs_scaling_factor + contig_yi_coord)-10)], fill = (200, 0, 0, 200))
		draw.text(((c[2] + 20), (i/contigs_scaling_factor+contig_yi_coord - 8)), ('Insertion ' + str(number)), font=fnt3, fill=(0,0,0,255))
		number = number + 1

		im.save("3_workflow-output/insertions_overview.png")



#_________________________________________________________Local and paired analysis graphs________________________________________________________________	
#_________________________________________________________________________________________________________________________________________________________
if args.mode == 'P':
	#_____________________________________________________________________________________________________________________________________________________
	insertions = list()

	for i, line in enumerate(lines):
		if not line.startswith('@'):	
			sp = line.split('\t')
			insertion = str(sp[2]).strip()
			if insertion not in insertions and insertion != '-':
				insertions.append(insertion)
	
			
	for e in insertions:
		try:
			del region_min
		except:
			pass
		region_max = 0
		rd_max_paired = 0
		rd_max_local = 0
		for i, line in enumerate(lines):
			if not line.startswith('@'):
				sp = line.split('\t')
				#Max and min for genome region in graphic
				if sp[2] == e: 
					if int(sp[3]) > region_max:
						region_max = int(sp[3])
					else:
						try:
							if sp[3] < region_min: 
								region_min = int(sp[3])
						except:
							region_min = int(sp[3])
		
				#Max and min read depth 
				if sp[2] == e and sp[0] == 'PAIRED': 
					if int(sp[4]) > rd_max_paired:
						rd_max_paired = int(sp[4])		

				if sp[2] == e and sp[0] == 'LOCAL_RD': 
					if int(sp[4]) > rd_max_local:
						rd_max_local = int(sp[4])		



		rd_max_paired = rd_max_paired + 10
		rd_max_local = rd_max_local + 5
		region_max = region_max + 100
		if region_min > 200:
			region_min = region_min - 100

		#Images, axes and title 
		im = Image.new("RGB", (1000, 1000), (255,255,255))
		draw = ImageDraw.Draw(im)
		draw.line((120, 450) + (900, 450), fill=256, width=3) 								#x axis paired
		draw.line((120, 150) + (120, 450), fill=256, width=3)   							#y axis paired
		draw.line((120, 755) + (900, 755), fill=256, width=3) 								#x axis local
		draw.line((120, 455) + (120, 755), fill=256, width=3)   							#y axis local
		
	
		draw.text(((450), (790)), ('Nucleotide'), font=fnt1, fill=(0,0,0,255))
		draw.text(((20), (120)), ('RD'), font=fnt1, fill=(0,0,0,255))
		draw.text(((140), (155)), ('Paired-reads analysis'), font=fnt3, fill=(0,0,0,255))
		draw.text(((140), (460)), ('Single-reads analysis'), font=fnt3, fill=(0,0,0,255))

	
		#Scaling factors 
		nucleotides = region_max - region_min
		scaling_factor_x = nucleotides/780.0
		scaling_factor_y_paired = rd_max_paired/300.0
		scaling_factor_y_local = rd_max_local/300.0


		#LOCAL/PAIRED GRAPHICS
		for i, line in enumerate(lines):
			if not line.startswith('@'):
				sp = line.split('\t')
				ins_contig = sp[1]
				if sp[2] == e and sp[0].strip() == 'PAIRED' and sp[5].strip() != 'TOTAL':
					raw_x_position = int(sp[3])
					img_x_position = int(raw_x_position/scaling_factor_x)
					img_relative_x_position = img_x_position - int(region_min/scaling_factor_x) + 121

					try:
						if img_relative_x_position == img_relative_x_position_2:
							raw_y_position = int(sp[4])
							img_y_position_p = 450 - int(raw_y_position/scaling_factor_y_paired)
					
							if img_y_position_p > img_y_position_p_2: 	
								img_y_position_p_2 = img_y_position_p
					
						img_relative_x_position_2 = img_relative_x_position
					except:
						img_relative_x_position_2 = img_relative_x_position
						raw_y_position = int(sp[4])
						img_y_position_p = 450 - int(raw_y_position/scaling_factor_y_paired)

			
					#draw
					if sp[5].strip() == 'R':
						draw.line((img_relative_x_position, 448) + (img_relative_x_position, img_y_position_p), fill=(255, 0, 0, 200), width=1)	
					
					elif sp[5].strip() == 'F':
						draw.line((img_relative_x_position, 448) + (img_relative_x_position, img_y_position_p), fill=(0, 0, 200, 200), width=1)


		for i, line in enumerate(lines):
			if not line.startswith('@'):
				sp = line.split('\t')
				ins_contig = sp[1]
				if sp[2] == e and sp[0].strip() == 'LOCAL_RD' and sp[5].strip() != 'TOTAL_RD':
					raw_x_position = int(sp[3])
					img_x_position = int(raw_x_position/scaling_factor_x)
					img_relative_x_position = img_x_position - int(region_min/scaling_factor_x) + 121

					raw_y_position = int(sp[4])
					img_y_position_l = int(raw_y_position/scaling_factor_y_local)
					img_relative_y_position = 755 - img_y_position_l 


					#draw
					if sp[5].strip() == 'RIGHT_RD':
						draw.line((img_relative_x_position, 753) + (img_relative_x_position, img_relative_y_position), fill=(255, 0, 0, 200), width=1)

					if sp[5].strip() == 'LEFT_RD':
						draw.line((img_relative_x_position, 753) + (img_relative_x_position, img_relative_y_position), fill=(0, 0, 200, 100), width=1)


		#Candidate regions
		for i, line in enumerate(lines):
			if line.startswith('@#'):
				sp = line.split(',')
				if int(sp[2].strip()) == int(e):
						draw.text(((120), (860)), ('Your candidate region is (' + sp[0].strip('@#') + ', ' + sp[1].strip()+')'), font=fnt3, fill=(0,0,0,255))
						draw.line((((120 +int(sp[0].strip('@#'))/scaling_factor_x - int(region_min/scaling_factor_x)) , 448) + ((120 +int(sp[0].strip('@#'))/scaling_factor_x - int(region_min/scaling_factor_x)) , 150)), fill=(150, 0, 150, 0), width=1)
						draw.line((((120 +int(sp[1].strip())/scaling_factor_x - int(region_min/scaling_factor_x)) , 448) + ((120 +int(sp[1].strip())/scaling_factor_x - int(region_min/scaling_factor_x)) , 150)), fill=(150, 0, 150, 0), width=1)


		#Axis anotations
		#x Axis
		x_p = 120 + int(100/scaling_factor_x)
		ruler = region_min + 100
	
		while x_p in range(120, 900):
			draw.line((x_p, 755) + (x_p, 760), fill=256, width=1)
			draw.text((x_p - 20, 762), (str(ruler)), font=fnt2, fill=(0,0,0,255))  #NOPE
			ruler = ruler + 100
			x_p = int(x_p + (100/scaling_factor_x)) #Ruler with 100 nts separations
	
		#y Axis - paired
		y_p = 450 - int(10/scaling_factor_y_paired)
		counter = 10
		while y_p in range(150, 451): 
			draw.line((120, y_p) + (115, y_p), fill=256, width=3)
			draw.text((80, y_p-10), ( str(counter)), font=fnt4, fill=(0,0,0,255))
			counter = counter + 10
			y_p = int(y_p - (10/scaling_factor_y_paired))


		#y Axis - local
		y_p = 755 - int(10/scaling_factor_y_paired)
		counter = 10
		while y_p in range(455, 751): 
			draw.line((120, y_p) + (115, y_p), fill=256, width=3)
			draw.text((80, y_p-10), ( str(counter)), font=fnt4, fill=(0,0,0,255))
			counter = counter + 10
			y_p = int(y_p - (10/scaling_factor_y_paired))


		#save image, specifying the format with the extension
		im.save('3_workflow-output/img_1_ins_' + str(e) + '.png')


#_________________________________________________________________________________________________________________________________________________________
if args.mode == 'S':
	insertions = list()

	for i, line in enumerate(lines):
		if not line.startswith('@'):	
			sp = line.split('\t')
			insertion = str(sp[2]).strip()
			if insertion not in insertions and insertion != '-':
				insertions.append(insertion)
	
			
	for e in insertions:
		try:
			del region_min
		except:
			pass
		region_max = 0
		rd_max_local = 0
		for i, line in enumerate(lines):
			if not line.startswith('@'):
				sp = line.split('\t')
				#Max and min for genome region in graphic
				if sp[2] == e: 
					if int(sp[3]) > region_max:
						region_max = int(sp[3])
					else:
						try:
							if sp[3] < region_min: 
								region_min = int(sp[3])
						except:
							region_min = int(sp[3])
		
				#Max and min read depth 
				if sp[2] == e and sp[0] == 'LOCAL_RD': 
					if int(sp[4]) > rd_max_local:
						rd_max_local = int(sp[4])		

		region_max = region_max + 100
		if region_min > 200:
			region_min = region_min - 100

		#Images, axes and title 
		im = Image.new("RGB", (1000, 600), (255,255,255))
		draw = ImageDraw.Draw(im)
		draw.line((120, 450) + (900, 450), fill=256, width=3) 								#x axis 
		draw.line((120, 150) + (120, 450), fill=256, width=3)   							#y axis 

			
		draw.text(((450), (500)), ('Nucleotide'), font=fnt1, fill=(0,0,0,255))
		draw.text(((20), (120)), ('RD'), font=fnt1, fill=(0,0,0,255))

	
		#Scaling factors 
		nucleotides = region_max - region_min
		scaling_factor_x = nucleotides/780.0
		scaling_factor_y_local = rd_max_local/300.0


		#LOCAL GRAPHICS


		for i, line in enumerate(lines):
			if not line.startswith('@'):
				sp = line.split('\t')
				ins_contig = sp[1]
				if sp[2] == e and sp[0].strip() == 'LOCAL_RD' and sp[5].strip() != 'TOTAL_RD':
					raw_x_position = int(sp[3])
					img_x_position = int(raw_x_position/scaling_factor_x)
					img_relative_x_position = img_x_position - int(region_min/scaling_factor_x) + 121

					raw_y_position = int(sp[4])
					img_y_position_l = int(raw_y_position/scaling_factor_y_local)
					img_relative_y_position = 450 - img_y_position_l 


					#draw
					if sp[5].strip() == 'RIGHT_RD':
						draw.line((img_relative_x_position, 448) + (img_relative_x_position, img_relative_y_position), fill=(255, 0, 0, 200), width=1)

					if sp[5].strip() == 'LEFT_RD':
						draw.line((img_relative_x_position, 448) + (img_relative_x_position, img_relative_y_position), fill=(0, 0, 200, 100), width=1)


		#Axis anotations
		#x Axis
		x_p = 120 + int(100/scaling_factor_x)
		ruler = region_min + 100
	
		while x_p in range(120, 900):
			draw.line((x_p, 450) + (x_p, 455), fill=256, width=1)
			draw.text((x_p - 20, 457), (str(ruler)), font=fnt2, fill=(0,0,0,255))  #NOPE
			ruler = ruler + 100
			x_p = int(x_p + (100/scaling_factor_x)) #Ruler with 100 nts separations
	
		#y Axis 
		y_p = 450 - int(10/scaling_factor_y_local)
		counter = 10
		while y_p in range(150, 451): 
			draw.line((120, y_p) + (115, y_p), fill=256, width=3)
			draw.text((80, y_p-10), ( str(counter)), font=fnt4, fill=(0,0,0,255))
			counter = counter + 10
			y_p = int(y_p - (10/scaling_factor_y_local))




		#save image, specifying the format with the extension
		im.save('3_workflow-output/img_1_ins_' + str(e) + '.png')




#________________________________________________________Gene Plot________________________________________________________________________________________	
#_________________________________________________________________________________________________________________________________________________________


#Output varanalyzer
input = args.input_va
f3 = open(input, 'r')
lines_va = f3.readlines()	

#Input gff
input = args.gff
f4 = open(input, 'r')
lines_gff = f4.readlines()	

#We create an 'intermediate list', which will contain the information necesary for the gene plot gathered from the output file of the varanalyzer module, the 
#sorted_insertions.txt file and the genome feature file. Intermediate list format:
#	['Chr', 'insertion position', 'insertion number', 'gene name', [list of gene features: 'type', 'start', 'end'], [list of positions(required for calculations)]]

intermediate_list = list()
for i, line in enumerate(lines_va):
	if not line.startswith('@'):
		sp = line.split('\t')
		if sp[5].strip() != 'nh':
			temp_list = list()
			temp_list.append(sp[1])
			temp_list.append(sp[2])
			temp_list.append('0')
			temp_list.append(sp[9])
			intermediate_list.append(temp_list)

for p in intermediate_list: 
	for i, line in enumerate(lines): 
		sp = line.split()
		if p[0] == sp[1] and p[1] == sp[3] and sp[5] == 'TOTAL' and sp[2] not in p[2]:
			p[2] = sp[2]


for p in intermediate_list:
	features = list()
	positions = list()
	for g, line in enumerate(lines_gff):
		sp = line.split('\t')
		if p[3] in sp[8]:
			feature = [sp[2], sp[3], sp[4]]
			positions.append(int(sp[3]))
			positions.append(int(sp[4]))
			features.append(feature)
	p.append(features)
	p.append(positions)

for p in intermediate_list:
	p[4].append(['rr', int(args.rrl)])
	p[5].append(min(p[5]) - int(args.rrl))

for p in intermediate_list:
	wide=1000 #<-------------------------------------------------------------------------------- SET IMAGE SIZE
	height=(35/100.0)*wide
	im = Image.new("RGB", (wide, int(height)), (255,255,255))
	draw = ImageDraw.Draw(im)
	gene_max_raw = max(p[5])
	gene_min_raw = min(p[5])
	gene_length = gene_max_raw - gene_min_raw

	gene_px_length = float((0.7)*wide)

	gene_scaling_factor = gene_length/gene_px_length  #pb/pixel

	#Fonts
	fnt1 = ImageFont.truetype('fonts/arial.ttf', int(0.03*wide))
	fnt2 = ImageFont.truetype('fonts/arial.ttf', int(0.016*wide))
	fnt3 = ImageFont.truetype('fonts/arial.ttf', int(0.024*wide))
	fnt4 = ImageFont.truetype('fonts/arial.ttf', int(0.02*wide))


	draw.line((int(0.15*wide), int(180/350.0*height)) + (int(0.15*wide) + gene_px_length, int(180/350.0*height)), fill=(0, 0, 160, 0), width=int(0.004*wide))
	

	for e in p[4]:
		if e[0].strip() == 'rr':
			inicio = int(0.15*wide)
			fin = int(0.15*wide) + int(e[1])/gene_scaling_factor
			draw.line((inicio, int(180/350.0*height)) + (fin, int(180/350.0*height)), fill=(150, 0, 150, 0), width=int(0.004*wide))

		if e[0].strip() == 'exon':
			inicio = int((int(e[1]) - gene_min_raw)/gene_scaling_factor)  + int(0.15*wide)
			fin = int((int(e[2]) - gene_min_raw)/gene_scaling_factor)  + int(0.15*wide) 
			draw.line((inicio, int(180/350.0*height)) + (fin, int(180/350.0*height)), fill=(0, 50, 150, 0), width=int(0.02*wide))
			draw.line((inicio, int(170/350.0*height)) + (fin, int(170/350.0*height)), fill=(0, 0, 0, 0), width=1)
			draw.line((inicio, int(190/350.0*height)) + (fin, int(190/350.0*height)), fill=(0, 0, 0, 0), width=1)
			draw.line((inicio, int(170/350.0*height)) + (inicio, int(190/350.0*height)), fill=(0, 0, 0, 0), width=1)
			draw.line((fin, int(170/350.0*height)) + (fin, int(190/350.0*height)), fill=(0, 0, 0, 0), width=1)

		
		if 'UTR' in e[0].strip():
			inicio = int((int(e[1]) - gene_min_raw)/gene_scaling_factor)  + int(0.15*wide) 
			fin = int((int(e[2]) - gene_min_raw)/gene_scaling_factor)  + int(0.15*wide) 
			draw.line((inicio, int(180/350.0*height)) + (fin, int(180/350.0*height)), fill=(0, 100, 255, 0), width=int(0.02*wide))
			draw.line((inicio, int(170/350.0*height)) + (fin, int(170/350.0*height)), fill=(0, 0, 0, 0), width=1)
			draw.line((inicio, int(190/350.0*height)) + (fin, int(190/350.0*height)), fill=(0, 0, 0, 0), width=1)
		
		'''
		if 'five_prime' in  e[0].strip():
			print 'OKs'
			temp = list()
			temp.append(int(e[1]))
			temp.append(int(e[2]))
			direct = max(temp)
			direct = int((direct - gene_min_raw)/gene_scaling_factor)  + 150
			draw.polygon([(direct + 20, 180), (direct, 160), (direct, 200)], fill = (200, 200, 0, 200))
		'''
		

	ins_pos = int((int(p[1]) - gene_min_raw)/gene_scaling_factor)  + int(0.15*wide) 

	draw.polygon([(ins_pos, int(170/350.0*height)), (ins_pos - int(0.02*wide), int(170/350.0*height) - int(0.025*wide)), (ins_pos + int(0.02*wide), int(170/350.0*height) - int(0.025*wide))], fill = (200, 0, 0, 200))
	draw.line((ins_pos, int(170/350.0*height)) + (ins_pos - int(0.02*wide), int(170/350.0*height) - int(0.025*wide)), fill=(0, 0, 0, 0), width=1)
	draw.line((ins_pos, int(170/350.0*height)) + (ins_pos + int(0.02*wide), int(170/350.0*height) - int(0.025*wide)), fill=(0, 0, 0, 0), width=1)
	draw.line((ins_pos - int(0.02*wide), int(170/350.0*height) - int(0.025*wide)) + (ins_pos + int(0.02*wide), int(170/350.0*height) - int(0.025*wide)), fill=(0, 0, 0, 0), width=1)

	draw.text((ins_pos - int(0.04*wide), int(0.115*wide)), ('Insertion ' + str(p[2])), font=fnt4, fill=(0,0,0,255))

	draw.text((int(0.05*wide), int(0.02*wide)), (str(p[3])), font=fnt3, fill=(0,0,0,255))


	#save image, specifying the format with the extension
	im.save('3_workflow-output/img_2_ins_' + str(p[2]) + '_gene_' + str(p[3])+ '.png')



#________________________________________________________Sequence Plot____________________________________________________________________________________	
#_________________________________________________________________________________________________________________________________________________________

#Fasta input
finput = args.input_f
f2 = open(finput, 'r')
flines = f2.readlines()

#Output varanalyzer
input = args.input_va
f3 = open(input, 'r')
lines_va = f3.readlines()	

#Input gff
input = args.gff
f4 = open(input, 'r')
lines_gff = f4.readlines()	

#Input 1
input = args.input
f1 = open(input, 'r')
lines = f1.readlines()	


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


intermediate_list_2 = list()						#['Chr', 'insertion position', 'insertion number']
for i, line in enumerate(lines_va):
	if not line.startswith('@'):
		sp = line.split('\t')
		temp_list = list()
		temp_list.append(sp[1])
		temp_list.append(sp[2])
		temp_list.append('0')
		intermediate_list_2.append(temp_list)

for p in intermediate_list_2: 
	for i, line in enumerate(lines): 
		sp = line.split()
		if p[0] == sp[1] and p[1] == sp[3] and sp[5] == 'TOTAL' and sp[2] not in p[2]:
			p[2] = sp[2]


adjacent_sequences = list()					 		#['Chr', 'insertion position', 'insertion number', '20 pbs upstream ', '20 pbs downstream']
for p in intermediate_list_2:
	with open(finput) as fp:
		fastalist = list()
		for name_contig, seq_contig in read_fasta(fp):
			if name_contig.lower() == '>'+p[0].strip():
				upstream =  seq_contig[int(p[1])-20:int(p[1])]
				downstream =  seq_contig[int(p[1]):int(p[1])+20]
				print upstream + '[ins]' + downstream
				p.append(upstream)
				p.append(downstream)

