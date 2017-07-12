#IPT9
#echo IPT9 Mapping/us
#time ./easymap -n IPT9 -w snp -sim -r at -g complete.gff -rs ref -cr bc -co par_mut -sm 1000+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,26-1,42-2,25-3,5-4,1-5,1/0,20-1,39-3,28-4,2-5,1/0,24-1,43-2,25-3,6-4,1-5,1/0,16-1,34-2,31-3,14-4,4-5,1+5,6769500+r+100 -ss 25+100,0+500,100+1+50+se
#rm -rf .er_projects/*/1_intermediate_files/*.sam
#rm -rf ./user_projects/*/1_intermediate_files/sim_data/

#CPX
#echo CPX Mapping
#time ./easymap -n CPX -w snp -sim -r at -g complete.gff -rs ref -cr bc -co par_mut -sm 1000+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,26-1,42-2,25-3,5-4,1-5,1/0,20-1,39-3,28-4,2-5,1/0,24-1,43-2,25-3,6-4,1-5,1/0,16-1,34-2,31-3,14-4,4-5,1+1,18012623+r+100 -ss 25+100,0+500,100+1+50+se
#rm -rf ./user_projects/*/1_intermediate_files/*.sam
#rm -rf ./user_projects/*/1_intermediate_files/sim_data/

#ABC
#echo ABC Mapping
#time ./easymap -n ABC -w snp -sim -r at -g complete.gff -rs ref -cr bc -co par_mut -sm 1000+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,26-1,42-2,25-3,5-4,1-5,1/0,20-1,39-3,28-4,2-5,1/0,24-1,43-2,25-3,6-4,1-5,1/0,16-1,34-2,31-3,14-4,4-5,1+3,4460890+r+100 -ss 25+100,0+500,100+1+50+se
#rm -rf ./user_projects/*/1_intermediate_files/*.sam
#rm -rf ./user_projects/*/1_intermediate_files/sim_data/


#Caso 2, caso 5, seminario
time ./easymap -n caso5 -w snp -sim -r chr -g complete.gff -rs noref -cr bc -co f2wt -sm 300+e -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,4000000+r+100 -ss 30+100,0+500,100+1+50+se
rm -rf ./user_projects/*/1_intermediate_files/sim_data/
time ./easymap -n caso2 -w snp -sim -r chr -g complete.gff -rs ref -cr bc -co f2wt -sm 300+e -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,4000000+r+100 -ss 30+100,0+500,100+1+50+se
rm -rf ./user_projects/*/1_intermediate_files/sim_data/

exit 

#Insertions
./easymap -P TFM -w ins -rs at -d sim -g complete.gff -sm 6+li -i pbinprok2.fa -ss 40+100,0+500,100+1+50+pe
rm -rf ./user_projects/*/1_intermediate_files

./easymap -P TFM -w ins -rs at -d sim -g complete.gff -sm 7+li -i pbinprok2.fa -ss 40+100,0+500,100+1+50+pe
rm -rf ./user_projects/*/1_intermediate_files

./easymap -P TFM -w ins -rs at -d sim -g complete.gff -sm 5+li -i pbinprok2.fa -ss 40+100,0+500,100+1+50+pe
rm -rf ./user_projects/*/1_intermediate_files

./easymap -P TFM -w ins -rs at -d sim -g complete.gff -sm 5+li -i pbinprok2.fa -ss 40+100,0+500,100+1+50+pe
rm -rf ./user_projects/*/1_intermediate_files

./easymap -P TFM -w ins -rs at -d sim -g complete.gff -sm 6+li -i pbinprok2.fa -ss 40+100,0+500,100+1+50+pe
rm -rf ./user_projects/*/1_intermediate_files

./easymap -P TFM -w ins -rs at -d sim -g complete.gff -sm 7+li -i pbinprok2.fa -ss 40+100,0+500,100+1+50+pe
rm -rf ./user_projects/*/1_intermediate_files

./easymap -P TFM -w ins -rs at -d sim -g complete.gff -sm 5+li -i pbinprok2.fa -ss 40+100,0+500,100+1+50+pe
rm -rf ./user_projects/*/1_intermediate_files

./easymap -P TFM -w ins -rs at -d sim -g complete.gff -sm 5+li -i pbinprok2.fa -ss 40+100,0+500,100+1+50+pe
rm -rf ./user_projects/*/1_intermediate_files

time ./easymap -n caso2 -w snp -sim -r mini -g complete.gff -rs ref -cr bc -co f2wt -sm 300+e -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,4000000+r+100 -ss 30+100,0+500,100+1+50+se
