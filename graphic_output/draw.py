import argparse, math
from PIL import Image, ImageDraw, ImageFont, ImageOps


parser = argparse.ArgumentParser()
parser.add_argument('-my_mut', action="store", dest = 'my_mut')		#snp or lin
parser.add_argument('-asnp', action="store", dest = 'input_snp')		
parser.add_argument('-bsnp', action="store", dest = 'input_f_snp')		#Fasta genome input
#parser.add_argument('-gff', action="store", dest = 'gff')			#Genome feature file
#parser.add_argument('-iva', action="store", dest = 'input_va') 	#Output de varanalyzer
#parser.add_argument('-rrl', action="store", dest = 'rrl') 			#Regulatory region lenght
#parser.add_argument('-f', action="store", dest = 'output_html')


parser.add_argument('-a', action="store", dest = 'input')		
parser.add_argument('-b', action="store", dest = 'input_f')				#Fasta genome input
parser.add_argument('-gff', action="store", dest = 'gff')				#Genome feature file
parser.add_argument('-iva', action="store", dest = 'input_va')	 		#Output de varanalyzer
parser.add_argument('-rrl', action="store", dest = 'rrl') 				#Regulatory region lenght
parser.add_argument('-f', action="store", dest = 'output_html')
parser.add_argument('-m', action="store", dest = 'mode', default = 'P')
parser.add_argument('-pname', action="store", dest='project_name')
parser.add_argument('-cross', action="store", dest='my_cross')
parser.add_argument('-snp_analysis_type', action="store", dest='my_snp_analysis_type')

args = parser.parse_args()

project = args.project_name


def red(p):
	if len(str(p)) <= 3:
		r = p

	elif len(str(p)) == 4:
		r = str(p)[:1]  + ' kb'

	elif len(str(p)) == 5:
		r =  str(p)[:2]  + ' kb'

	elif len(str(p)) == 6:
		r =  '0' + '.' + str(p)[:2] + ' Mb'

	elif len(str(p)) == 7:
		r = str(p)[:1] + '.' + str(p)[1:3] + ' Mb'

	elif len(str(p)) == 8:
		r = str(p)[:2] + ' Mb'

	elif len(str(p)) == 9:
		r = str(p)[:3] + ' Mb'

	return r; 




args = parser.parse_args()


#############################################################################################################
#																											#
# 							SNP - Alelic frequence VS Chromosome position									#
#																											#
#############################################################################################################
def fa_vs_pos():
	#Input 1
	input1 = args.input_snp
	f1 = open(input1, 'r')
	lines = f1.readlines()	

	#Input 2
	input2 = args.input_f_snp
	f2 = open(input2, 'r')
	lines_f = f2.readlines()	
	contig_source = args.input_f_snp


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
	contig_lengths = list()

	with open(contig_source) as fp:
		fastalist = list()
		for name_contig, seq_contig in read_fasta(fp):
			innerlist = list()
			innerlist.append(name_contig.strip('>'))
			innerlist.append(len(seq_contig))
			fastalist.append(innerlist)
			contig_lengths.append(len(seq_contig))

	max_contig_len = 0
	for i in contig_lengths:
		if int(i) > max_contig_len:
			max_contig_len = int(i)


	#FA vs POS graphs 
	for i in fastalist:
		
		wide=int(880*float(i[1])/max_contig_len) + 120								 #<-------------------------------------------------------------------------------- SET IMAGE SIZE
		height=500
		im = Image.new("RGB", (wide, int(height)), (255,255,255))
		draw = ImageDraw.Draw(im)
		
		#get fonts from foler 'fonts'
		fnt1 = ImageFont.truetype('fonts/VeraMono.ttf', 28)
		fnt2 = ImageFont.truetype('fonts/VeraMono.ttf', 14)
		fnt3 = ImageFont.truetype('fonts/VeraMono.ttf', 22)
		fnt4 = ImageFont.truetype('fonts/VeraMono.ttf', 18)

		r = red(int(i[1]))

		if 'Mb' in r:
			max_graph_x = int(math.ceil(int(i[1])/1000000.0))*1000000

		elif 'kb' in r: 
			max_graph_x = i[1]

		#Scaling factors
		scaling_factor_x = (max_graph_x)/(wide - 120)							#nts/pixel         <-----------------------------------------------------------#######################
		scaling_factor_y = (1.1/(63/100.0*height))								#fa/pixels

		#snps
		for l, line in enumerate(lines):
			sp = line.split()
			if i[0].lower() == sp[0].lower():
				fa = float(sp[6])/(float(sp[6])+float(sp[5]))
				fa_img = int(80/100.0*height) - int(fa/scaling_factor_y) - 1
				pos_img = int(int(sp[1])/scaling_factor_x) + 70
				draw.ellipse((pos_img-1, fa_img-1, pos_img+1, fa_img+1), fill=(147, 147, 147))


		if args.my_snp_analysis_type == 'f2wt':
			for l, line in enumerate(lines):
				sp = line.split()
				if i[0].lower() == sp[0].lower():
					fa = float(sp[8])/(float(sp[8])+float(sp[7]))
					fa_img = int(80/100.0*height) - int(fa/scaling_factor_y)
					pos_img = int(int(sp[1])/scaling_factor_x) + int(70)
					draw.ellipse((pos_img-2, fa_img-2, pos_img+2, fa_img+2), fill=(171, 219, 208))


		################################################################################################################################################################################################################
		my_cross = str(args.my_cross)
		#Boost / mm 																						
	
		binput = open(project + '/1_intermediate_files/map_info.txt', 'r')
		blines = binput.readlines()

		#Boost line
		if my_cross == 'oc' :
			for b, bline in enumerate(blines):
				sp = bline.split()
				if bline.startswith('!'):
					boost_max = float(sp[3])
			for b, bline in enumerate(blines):
				sp = bline.split()					
				if bline.startswith('@') and sp[4].lower().strip('>') == i[0].lower():
					boost_value = float(sp[3].strip())/boost_max
					boost_value_img = int(80/100.0*height) - int(boost_value/scaling_factor_y )

					window_position = int(sp[1])
					window_position_img = int(window_position/scaling_factor_x) + 70

					try:
						draw.line(((window_position_img, boost_value_img) + (window_position_img_2, boost_value_img_2)), fill=(255, 0, 0, 0), width=1)	
						window_position_img_2 = window_position_img
						boost_value_img_2 = boost_value_img

					except:
						window_position_img_2 = window_position_img
						boost_value_img_2 = boost_value_img

			window_position_img = None 
			boost_value_img = None
			window_position_img_2 = None 
			boost_value_img_2 = None

		#MM line
		if my_cross == 'oc' :
			for b, bline in enumerate(blines):
				sp = bline.split()					
				if bline.startswith('@') and sp[4].lower().strip('>') == i[0].lower():
					mm_value = float(sp[2].strip())
					mm_value_img = int(80/100.0*height) - int(mm_value/scaling_factor_y )
					window_position = int(sp[1])
					window_position_img = int(window_position/scaling_factor_x) + 70
					try:
						draw.line(((window_position_img, mm_value_img) + (window_position_img_2, mm_value_img_2)), fill=(46, 255, 0), width=1)	
						window_position_img_2 = window_position_img
						mm_value_img_2 = mm_value_img
					except:
						window_position_img_2 = window_position_img
						mm_value_img_2 = mm_value_img

		if my_cross == 'bc' :
			for b, bline in enumerate(blines):
				sp = bline.split()					
				if bline.startswith('@') and sp[3].lower().strip('>') == i[0].lower():
					mm_value = float(sp[2].strip())
					mm_value_img = int(80/100.0*height) - int(mm_value/scaling_factor_y )
					window_position = int(sp[1])
					window_position_img = int(window_position/scaling_factor_x) + 70
					try:
						draw.line(((window_position_img, mm_value_img) + (window_position_img_2, mm_value_img_2)), fill=(46, 255, 0), width=1)	
						window_position_img_2 = window_position_img
						mm_value_img_2 = mm_value_img
					except:
						window_position_img_2 = window_position_img
						mm_value_img_2 = mm_value_img

		window_position_img = None 
		mm_value_img = None
		window_position_img_2 = None 
		mm_value_img_2 = None

		#Axes
		draw.line((68, int(15/100.0*height)) + (68, int(80/100.0*height)), fill=(0, 0, 0, 0), width=1)	#Y axis
		draw.line((68, int(80/100.0*height)) + ((wide - 50), int(80/100.0*height)), fill=(0, 0, 0, 0), width=1)	#X axis

		draw.line(((wide - 50), int(15/100.0*height)) + ((wide - 50), int(80/100.0*height)), fill=(0, 0, 0, 0), width=1)	#-Y axis
		draw.line((68, int(15/100.0*height)) + ((wide - 50), int(15/100.0*height)), fill=(0, 0, 0, 0), width=1)	#-X axis

		#draw.text(((int(3/100.0*wide)), (int(7.5/100.0*height))), ('AF'), font=fnt3, fill=(0,0,0,255))
		
		#Axis rulers_____________________
		#X Axis 																		
		r = red(int(i[1]))
		mb_max = (max_graph_x/1000000)
		axis_px = (wide - 120)
		if 'Mb' in r:
			mb = 1
			for mb in range(1, mb_max + 1):	
				pos_x = int((mb*(float(axis_px)/mb_max)) + 70)
				draw.line((pos_x, int(81/100.0*height) ) + (pos_x, int(80/100.0*height)), fill=(0, 0, 0, 0), width=1)	
				
				if len(str(mb)) == 1:
					draw.text(((pos_x - 4), (int(81.5/100.0*height))), (str(mb).strip()), font=fnt2, fill=(0,0,0,255))
				elif len(str(mb)) == 2: 
					draw.text(((pos_x - 8), (int(81.5/100.0*height))), (str(mb).strip()), font=fnt2, fill=(0,0,0,255))


		#Y axis
		draw.line(( 68 , int(22.7/100.0*height) ) + ( 63 , int(22.7/100.0*height) ), fill=(0, 0, 0, 0), width=1)	
		draw.line(( 68 , int((28.65 + 22.7)/100.0*height) ) + ( 63 , int((28.65 + 22.7)/100.0*height) ), fill=(0, 0, 0, 0), width=1)	
		draw.line(( 68 , int((28.65 + 22.7 - 28.65/2)/100.0*height) ) + ( 65 , int((28.65 + 22.7 - 28.65/2)/100.0*height) ), fill=(0, 0, 0, 0), width=1)	
		draw.line(( 68 , int((28.65 + 22.7 + 28.65/2)/100.0*height) ) + ( 65 , int((28.65 + 22.7 + 28.65/2)/100.0*height) ), fill=(0, 0, 0, 0), width=1)	


		draw.text(((32), (int(21/100.0*height))), ( '1.0' ), font=fnt2, fill=(0,0,0,255))
		draw.text(((32), (int((29 + 20.7)/100.0*height))), ( '0.5' ), font=fnt2, fill=(0,0,0,255))
		
		#draw.text(((32), (int((28.65 + 22.7 - 28.65/2)/100.0*height))), ( '0.75' ), font=fnt2, fill=(0,0,0,255))
		#draw.text(((32), (int((28.65 + 22.7 + 28.65/2)/100.0*height))), ( '0.25' ), font=fnt2, fill=(0,0,0,255))

		#______________________________________

		#Y axis label

		label = Image.new("RGB", (140, 20), (255,255,255))
		draw2 = ImageDraw.Draw(label)
		draw2.text((1, 1), "Allele frequency", font=fnt2, fill=(0,0,0))
		label.rotate(90)
		im.paste(label.rotate(90), (2, int(30/100.0*height)))

		#X axis label
		x_title = str(i[0]) + ' (Mb)'
		w, h = draw.textsize(str(x_title))
		draw.text((( (wide-120)/2- w/2 +70), (int(87/100.0*height))), (x_title), font=fnt2, fill=(0,0,0,255))

		#save image, specifying the format with the extension
		im.save(project + '/3_workflow_output/img_2_' + str(i[0]) + '.png')


#############################################################################################################
#																											#
# 									LIN - INSERTIONS OVERVIEW & HISTOGRAMS									#
#																											#
#############################################################################################################

def insertions_overview_and_histograms():

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
	if args.mode == 'pe':
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

	elif args.mode == 'se':
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
	fnt1 = ImageFont.truetype('fonts/VeraMono.ttf', 30)
	fnt2 = ImageFont.truetype('fonts/VeraMono.ttf', 16)
	fnt3 = ImageFont.truetype('fonts/VeraMono.ttf', 14)
	fnt4 = ImageFont.truetype('fonts/VeraMono.ttf', 20)
	fnt5 = ImageFont.truetype('fonts/VeraMono.ttf', 18)

	tab = 50

	number = 1

	#Drawing the chromosomes:
	for c in superlist:
		draw.line((c[2], c[3]) + (c[2], c[4]), fill=(4, 14, 73), width=18)
		draw.ellipse((c[2]-8, c[3]-8, c[2]+8, c[3]+7), fill=(4, 14, 73))
		draw.ellipse((c[2]-8, c[4]-8, c[2]+8, c[4]+7), fill=(4, 14, 73))
		draw.text(((c[2] - tab), (c[3] - 60)), c[0], font=fnt5, fill=(0,0,0,255))
		for i in c[5]:
			draw.polygon([(c[2]+ 11, (i/contigs_scaling_factor+contig_yi_coord)), (c[2]+21, (i/contigs_scaling_factor+contig_yi_coord)+10), (c[2]+21, (i/contigs_scaling_factor + contig_yi_coord)-10)], fill = (200, 0, 0, 200))
			draw.line((c[2]+ 11, (i/contigs_scaling_factor+contig_yi_coord)) + (c[2]+21, (i/contigs_scaling_factor+contig_yi_coord)+10), fill=256, width=1) 
			draw.line((c[2]+ 11, (i/contigs_scaling_factor+contig_yi_coord)) + (c[2]+21, (i/contigs_scaling_factor + contig_yi_coord)-10), fill=256, width=1) 
			draw.line((c[2]+21, (i/contigs_scaling_factor + contig_yi_coord)-10) + (c[2]+21, (i/contigs_scaling_factor+contig_yi_coord)+10), fill=256, width=1) 
			draw.text(((c[2] + 30), (i/contigs_scaling_factor+contig_yi_coord - 8)), ('Insertion ' + str(number)), font=fnt3, fill=(0,0,0,255))
			number = number + 1

			im.save(project + "/3_workflow_output/insertions_overview.png")



	#_________________________________________________________Local and paired analysis graphs________________________________________________________________	
	#_________________________________________________________________________________________________________________________________________________________
	if args.mode == 'pe':
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
					if sp[2] == e and sp[0] == 'PAIRED' : 
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
			draw.line((120, 449) + (900, 449), fill=256, width=1) 								#x axis paired
			draw.line((120, 150) + (120, 449), fill=256, width=1)   							#y axis paired
			draw.line((120, 754) + (900, 754), fill=256, width=1) 								#x axis local
			draw.line((120, 455) + (120, 754), fill=256, width=1)   							#y axis local
			
			draw.line((120, 150) + (900, 150), fill=256, width=1) 								#-x axis paired
			draw.line((900, 150) + (900, 449), fill=256, width=1)   							#-y axis paired
			draw.line((120, 455) + (900, 455), fill=256, width=1) 								#-x axis local
			draw.line((900, 455) + (900, 754), fill=256, width=1)   							#-y axis local


			draw.text(((450), (795)), ('Nucleotide'), font=fnt3, fill=(0,0,0,255))
			draw.text(((140), (155)), ('Paired-reads analysis'), font=fnt3, fill=(0,0,0,255))
			draw.text(((140), (460)), ('Single-reads analysis'), font=fnt3, fill=(0,0,0,255))

			#Y axis label
			label = Image.new("RGB", (150, 30), (255,255,255))
			draw2 = ImageDraw.Draw(label)
			draw2.text((1, 1), "Read depth (x)", font=fnt3, fill=(0,0,0))
			label.rotate(90)
			im.paste(label.rotate(90), (35, 350))


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
							draw.line((img_relative_x_position, 448) + (img_relative_x_position, img_y_position_p), fill=(64, 159, 65, 100), width=1)	
						
						elif sp[5].strip() == 'F':
							draw.line((img_relative_x_position, 448) + (img_relative_x_position, img_y_position_p), fill=(31, 120, 180, 100), width=1)


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
							draw.line((img_relative_x_position, 753) + (img_relative_x_position, img_relative_y_position), fill=(64, 159, 65, 100), width=1)

						if sp[5].strip() == 'LEFT_RD':
							draw.line((img_relative_x_position, 753) + (img_relative_x_position, img_relative_y_position), fill=(31, 120, 180, 100), width=1)


			#Candidate regions
			for i, line in enumerate(lines):
				if line.startswith('@#'):
					sp = line.split(',')
					if int(sp[2].strip()) == int(e):
							cr = [int(sp[0].strip('@#')), int(sp[1].strip())]
							cr_min = min(cr)
							cr_max = max(cr)
							draw.text(((120), (840)), ('Your candidate region is (' + str(cr_min) + ', ' + str(cr_max) + ')'), font=fnt3, fill=(0,0,0,255))
							draw.line((((120 +int(sp[0].strip('@#'))/scaling_factor_x - int(region_min/scaling_factor_x)) , 448) + ((120 +int(sp[0].strip('@#'))/scaling_factor_x - int(region_min/scaling_factor_x)) , 190)), fill=(150, 0, 150, 0), width=1)
							draw.line((((120 +int(sp[1].strip())/scaling_factor_x - int(region_min/scaling_factor_x)) , 448) + ((120 +int(sp[1].strip())/scaling_factor_x - int(region_min/scaling_factor_x)) , 190)), fill=(150, 0, 150, 0), width=1)


			#Axis anotations
			#x Axis
			x_p = 120 + int(100/scaling_factor_x)
			x_p_2 = 120 + int(200/scaling_factor_x)
			ruler = region_min + 100

			while x_p in range(120, 900):
				draw.line((x_p, 755) + (x_p, 762), fill=256, width=1)
				w, h = draw.textsize(str(ruler))
				draw.text((x_p - w/2 - 5, 766), (str(ruler)), font=fnt3, fill=(0,0,0,255))  
				ruler = ruler + 200
				x_p = int(x_p + (200/scaling_factor_x)) #Ruler with 200 nts separations
			
			while x_p_2 in range(120, 900):
				draw.line((x_p_2, 755) + (x_p_2, 758), fill=256, width=1)
				x_p_2 = int(x_p_2 + (200/scaling_factor_x)) #Ruler with 100 nts separations


			#y Axis - paired
			y_p = 450 - int(10/scaling_factor_y_paired)
			counter = 10
			while y_p in range(150, 451): 
				draw.line((120, y_p) + (115, y_p), fill=256, width=1)
				draw.text((90, y_p-8), ( str(counter)), font=fnt3, fill=(0,0,0,255))
				counter = counter + 10
				y_p = int(y_p - (10/scaling_factor_y_paired))


			#y Axis - local
			y_p = 755 - int(10/scaling_factor_y_paired)
			counter = 10
			while y_p in range(455, 751): 
				draw.line((120, y_p) + (115, y_p), fill=256, width=1)
				draw.text((90, y_p-8), ( str(counter)), font=fnt3, fill=(0,0,0,255))
				counter = counter + 10
				y_p = int(y_p - (10/scaling_factor_y_paired))


			#save image, specifying the format with the extension
			w, h = im.size
			im.crop((0, 100, w, h-100)).save(project + '/3_workflow_output/img_1_ins_' + str(e) + '.png')

	#_________________________________________________________________________________________________________________________________________________________
	if args.mode == 'se':
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
			draw.line((120, 450) + (900, 450), fill=256, width=1) 								#x axis 
			draw.line((120, 150) + (120, 450), fill=256, width=1)   							#y axis 
			draw.line((120, 150) + (900, 150), fill=256, width=1) 								#-x axis 
			draw.line((900, 150) + (900, 450), fill=256, width=1)   							#-y axis 

				
			draw.text(((450), (500)), ('Nucleotide'), font=fnt3, fill=(0,0,0,255))
			#draw.text(((20), (120)), ('RD'), font=fnt3, fill=(0,0,0,255))

			#Y axis label
			label = Image.new("RGB", (150, 30), (255,255,255))
			draw2 = ImageDraw.Draw(label)
			draw2.text((1, 1), "Read depth (x)", font=fnt3, fill=(0,0,0))
			label.rotate(90)
			im.paste(label.rotate(90), (35, 195))

		
			#Scaling factors 
			nucleotides = region_max - region_min
			scaling_factor_x = nucleotides/780.0
			scaling_factor_y_local = rd_max_local/280.0


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
							draw.line((img_relative_x_position, 449) + (img_relative_x_position, img_relative_y_position), fill=(64, 159, 65, 100), width=3)

						if sp[5].strip() == 'LEFT_RD':
							draw.line((img_relative_x_position, 449) + (img_relative_x_position, img_relative_y_position), fill=(31, 120, 180, 100), width=3)


			#Axis anotations
			#x Axis
			x_p = 120 + int(25/scaling_factor_x)
			x_p_2 = 120 + int(50/scaling_factor_x)
			ruler = region_min + 25
		
			while x_p in range(120, 900):

				draw.line((x_p, 450) + (x_p, 457), fill=256, width=1)

				w, h = draw.textsize(str(ruler))
				draw.text((x_p - w/2 - 5, 457), (str(ruler)), font=fnt3, fill=(0,0,0,255))  

				ruler = ruler + 50

				x_p = int(x_p + (50/scaling_factor_x)) #Ruler with 200 nts separations

			while x_p_2 in range(120, 900):
				draw.line((x_p_2, 450) + (x_p_2, 455), fill=256, width=1)
				x_p_2 = int(x_p_2 + (50/scaling_factor_x)) #Ruler with 100 nts separations

			#y Axis 
			y_p = 450 - int(5/scaling_factor_y_local)
			counter = 5
			while y_p in range(150, 451): 
				draw.line((120, y_p) + (115, y_p), fill=256, width=1)
				draw.text((90, y_p-8), ( str(counter)), font=fnt3, fill=(0,0,0,255))
				counter = counter + 5
				y_p = int(y_p - (5/scaling_factor_y_local))


			#save image, specifying the format with the extension
			im.save(project + '/3_workflow_output/img_1_ins_' + str(e) + '.png')


#############################################################################################################
#																											#
# 											GENE PLOT														#
#																											#
#############################################################################################################

def gene_plot(): 

	if args.my_mut == 'lin':
		#Input 1
		input = args.input
		f1 = open(input, 'r')
		lines = f1.readlines()	

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
	#	['Chr', 'insertion/snp position', 'insertion number / -', 'gene name', [ref-aa, alt-aa, pos-aa, ref-base, alt-base, strand], [list of gene features: 'type', 'start', 'end'], [list of positions(required for calculations)]]


	intermediate_list = list()
	for i, line in enumerate(lines_va):
		if not line.startswith('@'):
			sp = line.split('\t')
			if sp[5].strip() != 'nh':
				temp_list = list()
				temp_list.append(sp[1])
				temp_list.append(sp[2])
				temp_list.append('-')
				temp_list.append(sp[9])
				refalt = [sp[12].strip(), sp[13].strip(), sp[11].strip(), sp[3].strip(), sp[4].strip(), sp[8].strip()]
				temp_list.append(refalt)
				intermediate_list.append(temp_list)


	if args.my_mut == 'lin':
		for p in intermediate_list: 
			for i, line in enumerate(lines): 
				sp = line.split()
				if p[0].lower().strip() == sp[1].lower().strip() and p[1].strip() == sp[3].strip() and sp[5] == 'TOTAL' and sp[2] not in p[2]:
					p[2] = sp[2]


	for p in intermediate_list:
		features = list()
		positions = list()
		refalt = list()
		for g, line in enumerate(lines_gff):
			if not line.startswith('#'):
				sp = line.split('\t')
				if p[3] in sp[8]:
					feature = [sp[2], sp[3], sp[4]]
					positions.append(int(sp[3]))
					positions.append(int(sp[4]))
					features.append(feature)
					p.append(features)
					p.append(positions)
			

	for p in intermediate_list:
		p[5].append(['rr', int(args.rrl)])
		if p[4][5].strip() == '+':
			p[6].append(min(p[6]) - int(args.rrl))
		if p[4][5].strip() == '-':
			p[6].append(max(p[6]) +  int(args.rrl))

	for p in intermediate_list:
		wide=1000 #<-------------------------------------------------------------------------------- SET IMAGE SIZE
		height=(35/100.0)*wide
		im = Image.new("RGB", (wide, int(height)), (255,255,255))
		draw = ImageDraw.Draw(im)
		gene_max_raw = max(p[6])
		gene_min_raw = min(p[6])
		gene_length = gene_max_raw - gene_min_raw

		gene_px_length = float((0.7)*wide)

		gene_scaling_factor = gene_length/gene_px_length  #pb/pixel

		#Fonts
		fnt1 = ImageFont.truetype('fonts/arial.ttf', int(0.03*wide))
		fnt2 = ImageFont.truetype('fonts/arial.ttf', int(0.016*wide))
		fnt3 = ImageFont.truetype('fonts/arial.ttf', int(0.024*wide))
		fnt4 = ImageFont.truetype('fonts/arial.ttf', int(0.02*wide))


		#Gene name
		draw.text((int(0.05*wide), int(0.03*wide)), (str(p[3])), font=fnt3, fill=(0,0,0,255))

		#Gene line
		if p[4][5] == '+':
			draw.line((int(0.15*wide) + int(int(args.rrl)/gene_scaling_factor), int(180/350.0*height)) + (int(0.15*wide) + gene_px_length, int(180/350.0*height)), fill=(14, 54, 119), width=int(0.004*wide))
		if p[4][5] == '-':
			draw.line((int(0.15*wide), int(180/350.0*height)) + (int(0.15*wide) + gene_px_length - int(int(args.rrl)/gene_scaling_factor), int(180/350.0*height)), fill=(14, 54, 119), width=int(0.004*wide))
		
		#Gene features
		atg_list = list()
		if p[4][5] == '+':
			for e in p[5]:		
				if e[0].strip() == 'rr':
					inicio = int(0.15*wide)
					fin = int(0.15*wide) + int(e[1])/gene_scaling_factor
					step = int(0.005*wide)
					s = inicio
					while s in range(inicio, int(fin)):
						draw.line((s, int(180/350.0*height)) + (s + step, int(180/350.0*height)), fill=(14, 54, 119), width=int(0.004*wide))
						s = s + step*2


		if p[4][5] == '-':
			for e in p[5]:		
				if e[0].strip() == 'rr':
					fin = int(0.85*wide)
					inicio = int(0.85*wide) - int(e[1])/gene_scaling_factor
					step = int(0.005*wide)
					s = int(inicio)
					while s in range(int(inicio), int(fin)):
						draw.line((s, int(180/350.0*height)) + (s + step, int(180/350.0*height)), fill=(14, 54, 119), width=int(0.004*wide))
						s = s + step*2



		for e in p[5]:

			if e[0].strip() == 'exon':
				inicio = int((int(e[1]) - gene_min_raw)/gene_scaling_factor)  + int(0.15*wide)
				fin = int((int(e[2]) - gene_min_raw)/gene_scaling_factor)  + int(0.15*wide) 
				draw.line((inicio, int(180/350.0*height)) + (fin, int(180/350.0*height)), fill=(59, 119, 214), width=int(0.02*wide))
				draw.line((inicio, int(170/350.0*height)) + (fin, int(170/350.0*height)), fill=(0, 4, 71, 0), width=2)	
				draw.line((inicio, int(190/350.0*height)) + (fin+1, int(190/350.0*height)), fill=(0, 4, 71, 0), width=2)
				draw.line((inicio, int(170/350.0*height)) + (inicio, int(190/350.0*height)), fill=(0, 4, 71, 0), width=2)
				draw.line((fin, int(170/350.0*height)) + (fin, int(190/350.0*height)), fill=(0, 4, 71, 0), width=2)
				if p[4][5] == '+':
					draw.line((inicio, int(170/350.0*height)) + (inicio, int(190/350.0*height)), fill=(0, 4, 71, 0), width=2)
				if p[4][5] == '-':
					draw.line((fin, int(170/350.0*height)) + (fin, int(190/350.0*height)), fill=(0, 4, 71, 0), width=2)

			if 'UTR' in e[0].strip():																														# Backup UTR drawing
				inicio = int((int(e[1]) - gene_min_raw)/gene_scaling_factor)  + int(0.15*wide) 
				fin = int((int(e[2]) - gene_min_raw)/gene_scaling_factor)  + int(0.15*wide) 
				draw.line((inicio, int(180/350.0*height)) + (fin, int(180/350.0*height)), fill=(188, 209, 242), width=int(0.02*wide))
				draw.line((inicio, int(170/350.0*height)) + (fin, int(170/350.0*height)), fill=(0, 4, 71, 0), width=2)
				draw.line((inicio, int(190/350.0*height)) + (fin, int(190/350.0*height)), fill=(0, 4, 71, 0), width=2)
				if p[4][5] == '+':
					draw.line((inicio, int(170/350.0*height)) + (inicio, int(190/350.0*height)), fill=(0, 4, 71, 0), width=2)
				if p[4][5] == '-':
					draw.line((fin, int(170/350.0*height)) + (fin, int(190/350.0*height)), fill=(0, 4, 71, 0), width=2)

			if 'utr' in (e[0].strip()).lower() and 'five' in (e[0].strip()).lower():																		# Backup UTR drawing
				inicio = int((int(e[1]) - gene_min_raw)/gene_scaling_factor)  + int(0.15*wide) 
				fin = int((int(e[2]) - gene_min_raw)/gene_scaling_factor)  + int(0.15*wide) 
				draw.line((inicio, int(180/350.0*height)) + (fin, int(180/350.0*height)), fill=(188, 209, 242), width=int(0.02*wide))
				draw.line((inicio, int(170/350.0*height)) + (fin, int(170/350.0*height)), fill=(0, 4, 71, 0), width=2)
				draw.line((inicio, int(190/350.0*height)) + (fin, int(190/350.0*height)), fill=(0, 4, 71, 0), width=2)
				draw.line((inicio, int(170/350.0*height)) + (inicio, int(190/350.0*height)), fill=(0, 4, 71, 0), width=2)
				draw.line((fin, int(170/350.0*height)) + (fin, int(190/350.0*height)), fill=(0, 4, 71, 0), width=2)


			if 'utr' in (e[0].strip()).lower() and 'three' in (e[0].strip()).lower():																		# Backup UTR drawing
				inicio = int((int(e[1]) - gene_min_raw)/gene_scaling_factor)  + int(0.15*wide) 
				fin = int((int(e[2]) - gene_min_raw)/gene_scaling_factor)  + int(0.15*wide) 
				draw.line((inicio, int(180/350.0*height)) + (fin, int(180/350.0*height)), fill=(188, 209, 242), width=int(0.02*wide))
				draw.line((inicio, int(170/350.0*height)) + (fin, int(170/350.0*height)), fill=(0, 4, 71, 0), width=2)
				draw.line((inicio, int(190/350.0*height)) + (fin, int(190/350.0*height)), fill=(0, 4, 71, 0), width=2)
				if p[4][5] == '+':
					draw.line((inicio, int(170/350.0*height)) + (inicio, int(190/350.0*height)), fill=(0, 4, 71, 0), width=2)
				if p[4][5] == '-':
					draw.line((fin, int(170/350.0*height)) + (fin, int(190/350.0*height)), fill=(0, 4, 71, 0), width=2)




		#Gene direction
		if p[4][5] == '+':
			draw.polygon([(int(0.84*wide), int(168/350.0*height)), (int(0.851*wide) , int(168/350.0*height)), (int(0.851*wide), int(181/350.0*height))], fill = (255, 255, 255, 0))
			draw.polygon([(int(0.84*wide), int(192/350.0*height)), (int(0.851*wide) , int(192/350.0*height)), (int(0.851*wide), int(180/350.0*height))], fill = (255, 255, 255, 0))
			draw.line((int(0.841*wide), int(170/350.0*height)) + (int(0.85*wide), int(180/350.0*height)), fill=(0, 4, 71, 0), width=2)
			draw.line((int(0.841*wide), int(190/350.0*height)) + (int(0.851*wide), int(180/350.0*height)), fill=(0, 4, 71, 0), width=2)


		
		if p[4][5] == '-':
			draw.polygon([(int(0.16*wide), int(168/350.0*height)), (int(0.149*wide) , int(168/350.0*height)), (int(0.149*wide), int(180/350.0*height))], fill = (255, 255, 255, 0))
			draw.polygon([(int(0.16*wide), int(192/350.0*height)), (int(0.149*wide) , int(192/350.0*height)), (int(0.149*wide), int(180/350.0*height))], fill = (255, 255, 255, 0))
			draw.line((int(0.158*wide), int(170.5/350.0*height)) + (int(0.148*wide), int(180.5/350.0*height)), fill=(0, 4, 71, 0), width=2)
			draw.line((int(0.159*wide), int(190/350.0*height)) + (int(0.149*wide), int(180/350.0*height)), fill=(0, 4, 71, 0), width=2)

			draw.line((int(0.15*wide), int(169/350.0*height)) + (int(0.16*wide), int(169/350.0*height)), fill=(255, 255, 255, 0), width=1)


		#Scale bar

		scale = 100
		px_scale = float(scale/gene_scaling_factor)
		draw.line((int(0.95*wide) - int(px_scale), int(290/350.0*height)) + (int(0.95*wide), int(290/350.0*height)), fill=(0, 0, 0, 0), width=int(0.002*wide))
		draw.text((int(0.9*wide), int(300.8/350.0*height)), ('100 pb'), font=fnt2, fill=(0,0,0,255))


		#Insertion triangle and info
		if args.my_mut == 'lin':
			ins_pos = int((int(p[1]) - gene_min_raw)/gene_scaling_factor)  + int(0.15*wide) 
			draw.polygon([(ins_pos, int(170/350.0*height)), (ins_pos - int(0.02*wide), int(170/350.0*height) - int(0.025*wide)), (ins_pos + int(0.02*wide), int(170/350.0*height) - int(0.025*wide))], fill = (200, 0, 0, 200))
			draw.line((ins_pos, int(170/350.0*height)) + (ins_pos - int(0.02*wide), int(170/350.0*height) - int(0.025*wide)), fill=(0, 0, 0, 0), width=1)
			draw.line((ins_pos, int(170/350.0*height)) + (ins_pos + int(0.02*wide), int(170/350.0*height) - int(0.025*wide)), fill=(0, 0, 0, 0), width=1)
			draw.line((ins_pos - int(0.02*wide), int(170/350.0*height) - int(0.025*wide)) + (ins_pos + int(0.02*wide), int(170/350.0*height) - int(0.025*wide)), fill=(0, 0, 0, 0), width=1)
			draw.text((ins_pos - int(0.04*wide), int(0.115*wide)), ('Insertion ' + str(p[2])), font=fnt4, fill=(0,0,0,255))


		#SNP arrow and info
		if args.my_mut == 'snp':
			snp_pos = int((int(p[1]) - gene_min_raw)/gene_scaling_factor)  + int(0.15*wide) 
			draw.line((snp_pos, int(194/350.0*height)) + (snp_pos , int(194/350.0*height) + int(0.03*wide)), fill=(180, 0, 0, 0), width=int(0.005*wide))						
			draw.polygon([(snp_pos, int(191/350.0*height)), (snp_pos - int(0.01*wide), int(191/350.0*height) + int(0.01*wide)), (snp_pos + int(0.01*wide), int(191/350.0*height) + int(0.01*wide))], fill = (200, 0, 0, 200))

			#Aa change
			if p[4][0].strip() != '-' : 																#		<---------------------------- DELETE THIS ARROW IF NO ERRORS COME FROM THIS LINE 
				
				draw.text((int(snp_pos - int(0.092*wide)), int(0.75*height)), (
					str(p[4][0])+ ' (' + str(p[4][2]) +')' +  '    >    '  +
					str(p[4][1])), font=fnt4, fill=(0,0,0,255))   

			#Base change
			draw.text((int(snp_pos - int(0.036*wide)), int(0.67*height)), (
				str(p[4][3]) +   '    >    '  +
				str(p[4][4])), font=fnt4, fill=(0,0,0,255))   



		#save image, specifying the format with the extension
		if args.my_mut == 'lin':
			im.save(project + '/3_workflow_output/gene_plot_' + str(args.my_mut) + '_' + str(p[2]) + '_gene_' + str(p[3])+ '.png')

		if args.my_mut == 'snp':
			im.save(project + '/3_workflow_output/gene_plot_' + str(args.my_mut) + '_' + str(p[1]) + '_gene_' + str(p[3])+ '.png')

#############################################################################################################
#																											#
# 											GENE PLOT - PRUEBAS												#
#																											#
#############################################################################################################

def gene_plot2(): 

	W=1000 #<-------------------------------------------------------------------------------- SET IMAGE SIZE
	H=(350)
	im = Image.new("RGB", (W, H), (255,255,255))
	draw = ImageDraw.Draw(im)



	text = 'potato'



	w, h = draw.textsize(text)


	point_x = 600
	point_y = 350/2

	draw.text(((point_x - w/2), (point_y-h/2)), text, fill="black")


	draw.ellipse((point_x-2, point_y-2, point_x+2, point_y+2), fill=(255, 25, 216))



	im.save('project/3_workflow_output/prueba' + '.png')

