#!/usr/bin/python
import argparse
from PIL import Image, ImageDraw, ImageFont, ImageOps

parser = argparse.ArgumentParser()
parser.add_argument('-fasq', action="store", dest = 'File', required = "True")

args = parser.parse_args()
reverse = "no"
files = args.File.split(",")
if len(files) == 2:
	forward = files[0]
	reverse = files[1]
if len(files) == 1:
	forward = files[0]

#phred_not="J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,[,\,],^,_,`,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,{,|,},~"

def character_to_ASCII(string):
	st = []
	for items in string:
		ascii = ord(items)
		if ascii > 73: ascii = ascii-64
		else: 	ascii = ascii-33
		st.append(str(ascii))
	return st


def fasq_process(fil):
	fastaq_dic = {}
	i = 0
	with open(fil) as fastaq:
		n = 0
		for lines in fastaq:
			n += 1
			lines = lines.rstrip()
			if lines.startswith("@"):
				pos = lines
				fastaq_dic[pos]= ""
				m = 0
			elif lines.startswith("+"):
				m = 1
				continue
			elif m == 1:
				fastaq_dic[pos]+=lines
				i += 1
			if i >= 50000:
				break

		#number_reads = float(n)/4
		#if number_reads < 1000000:
		#	representative_amount = number_reads
		#else:
		#	representative_amount = 1000000
	qual_dic= {}
	n = 0	
	for position in fastaq_dic:
		qual_dic[position] = character_to_ASCII(fastaq_dic[position])
		n += 1
		#if n >= representative_amount:
		#	break
	return qual_dic
def average(lista):
	n = 0 
	l = []
	for items in lista:
		l.append(float(items))
		n += 1
	average = sum(l)/n
	return average 
def boxplot_stats(lista):
	lista = sorted(lista)
	lenght = len(lista)
	
	position1 = float(lenght/2)
	if position1.is_integer():
		per_50 = lista[int(position1)]
	else:
		position1 = position1.split(".")[0]
		position2 = position1 + 1
		position1 =  lista[int(position1)]
		position2 = lista[int(position2)]
		per_50 = (postion1 + position2)/2 

	position2 = float((lenght+1)/4)
	if position2.is_integer():
		per_25 = lista[int(position2)]
	else:
		position2 = position2.split(".")[0]
		position3 = position2 + 1
		position2 =  lista[int(position2)]
		position3 = lista[int(position3)]
		position4 = position2 + position3*position2.split(".")[1]*((lenght+1)/4)
		per_25 = int(position4) 

	position3 = float((lenght+1)*3/4)
	if position3.is_integer():
		per_75 = lista[int(position3)]
	else:
		position3 = position3.split(".")[0]
		position4 = position4 + 1
		position3 =  lista[int(position3)]
		position4 = lista[int(position4)]
		position5 = position3 + position4*position3.split(".")[1]*((lenght+1)*3/4)
		per_75 = int(position5) 
	
	IQR = int(per_75) - int(per_25)
	Extreme_positive = 1.5*IQR + float(per_75)
	Extreme_negative =  float(per_25)-1.5*IQR 
	if Extreme_negative < 0:
		Extreme_negative = 0
	if float(min(lista))> Extreme_negative: Extreme_negative = min(lista)
	if float(max(lista))< Extreme_positive: Extreme_positive = max(lista)
	
	return Extreme_negative, per_25, per_50, per_75, Extreme_positive

def calculations(dic_pos):
	value_list = []
	biggest = 0
	for reads in dic_pos:	
		if len(dic_pos[reads]) > biggest:
			biggest = len(dic_pos[reads]) 

	dic_qual =  {}		
	for position in range(biggest):
		dic_qual[position] = []
		for reads in dic_pos:
			dic_qual[position].append(dic_pos[reads][position])
			
	dic_pos.clear()
	dic_final = {}
	m = 0 
	for values in dic_qual:
		if m == 0:
			val = []	
		val.extend(dic_qual[values])
		m += 1
		if m == 5:
			intervalo = int(values) - m
			intervalo = str(intervalo)+"-"+str(values)
			if intervalo not in dic_final: dic_final[intervalo] = []
			dic_final[intervalo].append(average(val))
			stats = boxplot_stats(val)
			
			dic_final[intervalo].extend([stats[0],stats[1],stats[2],stats[3],stats[4]])
			m = 0

	return dic_final 


def Draw_box_plot(table):
	fnt1 = ImageFont.truetype('easymap-server/fonts/arial.ttf', 7)
	#Size of the window
	a = 20
	b = 20
	c = 20
	x_window = len(table)*10+ 10+ a + b 
	y_window = 400
	#Generation of the file
	im = Image.new("RGB", (x_window, y_window), (255,255,255))	
	draw = ImageDraw.Draw(im)
	#Creation of the axes: exes will start with an indentation of 20 above, below and beside. Each Phred quality score will be in a X% proportion of the y_window  px
	size_y_exe = y_window-40 #Total size minus above and below indentations
	position_y_exe= size_y_exe+20
	draw.line(((a, c) + (a, position_y_exe)), fill=(0, 0, 0, 0), width=1)
	size_x_exe = len(table)*10 +10 #number of positions*10 pxls which is one will take + 10 for the position 1. 
	draw.line(((a, position_y_exe) + (a+size_x_exe, position_y_exe)), fill=(0, 0, 0, 0), width=1) 
	
	#Vertical values
	step = float(size_y_exe)/42
	for values in range(42,-1,-1):
		#pos = str(values+1)
		draw.line(((a,20+abs(values-42)*step) + (a-4,20+abs(values-42)*step)), fill=(0, 0, 0, 0), width=1)
		#draw.text((a-8,20+values*step) + (a-8,20+values*step), pos, font=fnt1, fill=(0,0,0,0))
		if values%5 == 0:
			draw.line(((a,20+abs(values-42)*step) + (a-5,20+abs(values-42)*step)), fill=(0, 0, 0, 0), width=1)

	i = 10 + a #indentation + space for the first box (same space as in size_x_exe)
	for position in table:
		#write the position in the x axe
		draw.line(((i, position_y_exe) + (i, position_y_exe+5)), fill=(0, 0, 0, 0), width=1)

		#Create a line from the begining to the end of the parameters
		beg = float(table[position][1]) * step
		end = float(table[position][-1]) * step
		draw.line(((i, position_y_exe-beg) + (i, position_y_exe-end)), fill=(0, 0, 0, 0), width=1)

		#Create the boxplot using CORREGIR!!!!!!!!!!
		beg = float(table[position][2]) * step
		end = float(table[position][-2]) * step
		draw.rectangle([(i-3, position_y_exe-beg), (i+3, position_y_exe-end)], fill=(24, 56, 214), outline= None)

		#Draw the average and the MEDIANA?
		av = float(table[position][0]) * step
		med = float(table[position][3]) * step
		draw.line(((i-3, position_y_exe-med) + (i+3, position_y_exe-med)), fill=(191, 17, 54), width=1)
		draw.line(((i, position_y_exe-av) + (i, position_y_exe-av)), fill=(50, 214, 25), width=1)
		i +=10


	#save image, specifying the format with the extension
	im.save("result.png")






forward_reads = fasq_process(forward)
forward_table = calculations(forward_reads)
reverse_reads = "no"
if reverse != "no": reverse_reads = fasq_process(reverse); reverse_table = calculations(reverse_reads)

Draw_box_plot(forward_table)


