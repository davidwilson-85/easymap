#!/usr/bin/python
import argparse

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

		number_reads = float(n)/4
		if number_reads < 1000000:
			representative_amount = number_reads
		else:
			representative_amount = 1000000
	qual_dic= {}
	n = 0	
	for position in fastaq_dic:
		qual_dic[position] = character_to_ASCII(fastaq_dic[position])
		n += 1
		if n >= representative_amount:
			break
	return qual_dic
def average(lista):
	n = 0 
	l = []
	for items in lista:
		l.append(float(items))
		n += 1
	average = sum(l)/n
	return average 


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
	dic_final = {}
	for values in dic_qual:
		if values not in dic_final: dic_final[values] = []
		dic_final[values].append(average(dic_qual[values]))
		#append more statistics values

	return dic_final 

forward_reads = fasq_process(forward)
forward_table = calculations(forward_reads)
reverse_reads = "no"
if reverse != "no": reverse_reads = fasq_process(reverse); reverse_table = calculations(reverse_reads)

print forward_table