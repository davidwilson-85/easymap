#!/usr/bin/python3





file = open("./dev/complete.gff", 'r')



lista = list()

for line in file:
	sp = line.split()
	if sp[0].strip() not in lista:
		lista.append(sp[0].strip())
print lista