#!/bin/bash

mutation=2137705,6306495,12385219,17669884
#6751663,2336221,5710945,14171828,12310030,4246170,4916889,16553990,5163558,7100534,10003051,6231302,17546556,2185860,2384375,1788833,3575452,11864850,14295161,16875834,9133823,13774804,13501110,672934,1987346,2979912
#2137704
python running.py -mut $mutation


######################
for i in ./user_projects/*take_*/
do
	#mv $i ../easy_map-boost_essay/user_projects/
	
	INPUT=$i
	SUBSTRING=$(echo $INPUT| cut -d'_' -f 5)
	SUBSTRING2=$(echo $SUBSTRING| cut -d'/' -f 1)

	#mv ../easy_map-boost_essay/user_projects/$i ../easy_map-boost_essay/user_projects/$SUBSTRING
	m=$(echo $SUBSTRING2| cut -d'-' -f 1)
	k=$(echo $SUBSTRING2| cut -d'-' -f 2)
	
	mv $i ../easy_map-boost_essay/user_projects/$SUBSTRING2


	bash ../easy_map-boost_essay/run.sh $m $k
	


done 

#find -name variants.txt > lista-rutas
#Ordenar la lista de rutas de alguna forma
#Usar Root.py
######################

