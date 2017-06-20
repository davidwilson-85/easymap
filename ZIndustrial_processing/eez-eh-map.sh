#!/bin/bash

mutation=2137704
 #,6306480,12385209,17669884,6751662,2336217,5710941,14171827,12310030,4246170,4916888,16553985,5163557,7100533,10003051,6231301,17546550,2185860,2384373,1788833,3575452,11864848,14295161,16875832,9133820,13774801,13501109,672933,1987346,2979912
echo $mutation
echo hola $mutation
#python running.py -mut $mutation
exit




######################
for i in ../user_projects/*take*:
do
	mv $i ../easy_map_boost_essay/user_projects
	
	INPUT=../easy_map_boost_essay/user_projects/$i
	SUBSTRING=$(echo $INPUT| cut -d'_' -f 4)
	mv ../easy_map_boost_essay/user_projects/$i ../easy_map_boost_essay/user_projects/$SUBSTRING


done 
######################