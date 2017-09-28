file = open("./user_data/complete.gff", 'r')


target = "AT1G24706"


for line in file:
	if target in line:
		print line