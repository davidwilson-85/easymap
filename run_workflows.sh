#!/bin/bash



./easymap -P caso_1 -w snp -d sim -g chr1+4.gff -sm 40+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,1000000+r+50 -ss 25+100,0+500,100+1+50+se -c bc -r ref --control_type par_mut

rm -rf ./user_projects/*/1_intermediate_files

./easymap -P caso_2 -w snp -d sim -g chr1+4.gff -sm 40+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,1000000+r+50 -ss 25+100,0+500,100+1+50+se -c bc -r ref --control_type f2wt

rm -rf ./user_projects/*/1_intermediate_files

./easymap -P caso_3 -w snp -d sim -g chr1+4.gff -sm 40+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,1000000+r+50 -ss 25+100,0+500,100+1+50+se -c oc -r ref --control_type par_mut

rm -rf ./user_projects/*/1_intermediate_files

./easymap -P caso_4 -w snp -d sim -g chr1+4.gff -sm 40+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,1000000+r+50 -ss 25+100,0+500,100+1+50+se -c bc -r ref --control_type par_nomut

rm -rf ./user_projects/*/1_intermediate_files


./easymap -P caso_7 -w snp -d sim -g chr1+4.gff -sm 40+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,1000000+r+50 -ss 25+100,0+500,100+1+50+se -c bc -r noref --control_type f2wt

rm -rf ./user_projects/*/1_intermediate_files

./easymap -P caso_8 -w snp -d sim -g chr1+4.gff -sm 40+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,1000000+r+50 -ss 25+100,0+500,100+1+50+se -c oc -r noref --control_type par_mut


rm -rf ./user_projects/*/1_intermediate_files



#./easymap -P caso_4 -w snp -d sim -g chr1+4.gff -sm 40+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,1000000+r+50 -ss 25+100,0+500,100+1+50+se -c oc -r noref --control_type par_nomut




