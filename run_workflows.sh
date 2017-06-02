#!/bin/bash



./easymap -P caso_1_2_mut -w snp -rs at -d sim -g chr1+4.gff -sm 400+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,5840000-2,1000000+dr+200 -ss 25+100,0+500,100+1+50+se -c bc -r ref -t par_mut
rm -rf ./user_projects/*/1_intermediate_files

./easymap -P caso_1_2_mut -w snp -rs at -d sim -g chr1+4.gff -sm 400+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,5840000-1,9236758+dr+200 -ss 25+100,0+500,100+1+50+se -c bc -r ref -t par_mut
rm -rf ./user_projects/*/1_intermediate_files

./easymap -P caso_3_2_mut -w snp -rs at -d sim -g chr1+4.gff -sm 400+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,5840000-2,1000000+dr+200 -ss 25+100,0+500,100+1+50+se -c oc -r ref -t par_mut
rm -rf ./user_projects/*/1_intermediate_files

./easymap -P caso_3_2_mut -w snp -rs at -d sim -g chr1+4.gff -sm 400+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,5840000-1,9236758+dr+200 -ss 25+100,0+500,100+1+50+se -c oc -r ref -t par_mut
rm -rf ./user_projects/*/1_intermediate_files


#./easymap -P caso_1_2_mut -w snp -rs at -d sim -g chr1+4.gff -sm 400+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,5840000-2,1000000+dr+200 -ss 25+100,0+500,100+1+50+se -c bc -r ref -t par_mut
#rm -rf ./user_projects/*/1_intermediate_files

#./easymap -P caso_2_2_mut -w snp -rs at -d sim -g chr1+4.gff -sm 400+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,5840000-2,1000000+dr+200 -ss 25+100,0+500,100+1+50+se -c bc -r ref -t f2wt
#rm -rf ./user_projects/*/1_intermediate_files

#./easymap -P caso_3_2_mut -w snp -rs at -d sim -g chr1+4.gff -sm 400+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,5840000-2,1000000+dr+200 -ss 25+100,0+500,100+1+50+se -c oc -r ref -t par_mut
#rm -rf ./user_projects/*/1_intermediate_files

#./easymap -P caso_4_2_mut -w snp -rs at -d sim -g chr1+4.gff -sm 400+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,5840000-2,1000000+dr+200 -ss 25+100,0+500,100+1+50+se -c oc -r ref -t par_nomut
#rm -rf ./user_projects/*/1_intermediate_files


#./easymap -P caso_7_2_mut -w snp -rs at -d sim -g chr1+4.gff -sm 400+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,584000000000-2,1000000+dr+200 -ss 25+100,0+500,100+1+50+se -c bc -r noref -t f2wt
#rm -rf ./user_projects/*/1_intermediate_files

#./easymap -P caso_8_2_mut -w snp -rs at -d sim -g chr1+4.gff -sm 400+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,5840000-2,1000000+dr+200 -ss 25+100,0+500,100+1+50+se -c oc -r noref -t par_mut
#rm -rf ./user_projects/*/1_intermediate_files



#****************************************************************************************************

#./easymap -P caso_1 -w snp -rs at -d sim -g chr1+4.gff -sm 400+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,5840000+r+200 -ss 25+100,0+500,100+1+50+se -c bc -r ref -t par_mut
#rm -rf ./user_projects/*/1_intermediate_files

#./easymap -P caso_2 -w snp -rs at -d sim -g chr1+4.gff -sm 400+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,5840000+r+200 -ss 25+100,0+500,100+1+50+se -c bc -r ref -t f2wt
#rm -rf ./user_projects/*/1_intermediate_files

#./easymap -P caso_3 -w snp -rs at -d sim -g chr1+4.gff -sm 400+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,5840000+r+200 -ss 25+100,0+500,100+1+50+se -c oc -r ref -t par_mut
#rm -rf ./user_projects/*/1_intermediate_files

#./easymap -P caso_4 -w snp -rs at -d sim -g chr1+4.gff -sm 400+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,5840000+r+200 -ss 25+100,0+500,100+1+50+se -c oc -r ref -t par_nomut
#rm -rf ./user_projects/*/1_intermediate_files


#./easymap -P caso_7 -w snp -rs at -d sim -g chr1+4.gff -sm 400+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,5840000+r+200 -ss 25+100,0+500,100+1+50+se -c bc -r noref -t f2wt
#rm -rf ./user_projects/*/1_intermediate_files

#./easymap -P caso_8 -w snp -rs at -d sim -g chr1+4.gff -sm 400+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,5840000+r+200 -ss 25+100,0+500,100+1+50+se -c oc -r noref -t par_mut
#rm -rf ./user_projects/*/1_intermediate_files



