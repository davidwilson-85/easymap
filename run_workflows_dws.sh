#!/bin/bash


# INSERCIONES
#./easymap -n ins-test -w ins -sim -r at2 -i at_chloroplast.fa -g TAIR10_GFF3_genes.gff -sm 10 -ss 15+100,0+500,100+1+50+pe


#./easymap -n caso1 -w snp -sim -r at -g TAIR10_GFF3_genes.gff -rs ref -cr bc -co par_mut -sm 150+e -sr 0,24-1,43-2,25-3,6-4,1-5,1+1,5860210+r+100 -ss 25+100,0+500,100+1+50+se


#./easymap -n caso1 -w snp -sim -r at -g TAIR10_GFF3_genes.gff -rs ref -cr bc -co par_mut -sm 75+e -sr 0,14-1,31-2,33-3,15-4,5-5,2+1,5860210+r+100 -ss 25+100,0+500,100+1+50+se
#rm -rf ./user_projects/*/1_intermediate_files

#./easymap -n caso1 -w snp -sim -r at -g TAIR10_GFF3_genes.gff -rs ref -cr bc -co par_mut -sm 150+e -sr 0,14-1,31-2,33-3,15-4,5-5,2+1,5860210+r+100 -ss 25+100,0+500,100+1+50+se
#rm -rf ./user_projects/*/1_intermediate_files

#./easymap -n caso1 -w snp -sim -r at -g TAIR10_GFF3_genes.gff -rs ref -cr bc -co par_mut -sm 225+e -sr 0,14-1,31-2,33-3,15-4,5-5,2+1,5860210+r+100 -ss 25+100,0+500,100+1+50+se
#rm -rf ./user_projects/*/1_intermediate_files

#./easymap -n caso1 -w snp -sim -r at -g TAIR10_GFF3_genes.gff -rs ref -cr bc -co par_mut -sm 300+e -sr 0,14-1,31-2,33-3,15-4,5-5,2+1,5860210+r+100 -ss 25+100,0+500,100+1+50+se
#rm -rf ./user_projects/*/1_intermediate_files












./easymap -n caso1 -w snp -sim -r at -g TAIR10_GFF3_genes.gff -rs ref -cr bc -co par_mut -sm 400+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,5487500+r+200 -ss 25+100,0+500,100+1+50+se
rm -rf ./user_projects/*/1_intermediate_files

./easymap -n caso2 -w snp -sim -r at -g TAIR10_GFF3_genes.gff -rs ref -cr bc -co f2wt -sm 400+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,5487500+r+200 -ss 25+100,0+500,100+1+50+se
rm -rf ./user_projects/*/1_intermediate_files

./easymap -n caso3 -w snp -sim -r at -g TAIR10_GFF3_genes.gff -rs ref -cr oc -co par_mut -sm 400+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,5860210+r+100 -ss 25+100,0+500,100+1+50+se
rm -rf ./user_projects/*/1_intermediate_files

./easymap -n caso4 -w snp -sim -r at -g TAIR10_GFF3_genes.gff -rs ref -cr oc -co par_nomut -sm 400+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,5860210+r+100 -ss 25+100,0+500,100+1+50+se
rm -rf ./user_projects/*/1_intermediate_files

./easymap -n caso6 -w snp -sim -r at -g TAIR10_GFF3_genes.gff -rs noref -cr oc -co par_mut -sm 400+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,5860210+r+100 -ss 25+100,0+500,100+1+50+se
rm -rf ./user_projects/*/1_intermediate_files










# ref bc par_mut. SIMULADOR FUNCIONA
#./easymap -P caso_1 -w snp -sim -rs at -g TAIR10_GFF3_genes.gff -sm 400+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,5487500+r+200 -ss 25+100,0+500,100+1+50+se -r ref -c bc -t par_mut
#rm -rf ./user_projects/*/1_intermediate_files

# ref bc f2wt. SIMULADOR FUNCIONA
#./easymap -P caso_2 -w snp -sim -rs at -g TAIR10_GFF3_genes.gff -sm 400+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,5487500+r+200 -ss 25+100,0+500,100+1+50+se -r ref -c bc -t f2wt
#rm -rf ./user_projects/*/1_intermediate_files

# ref oc par_mut. SIMULADOR FUNCIONA
#./easymap -P caso_3 -w snp -sim -rs at -g TAIR10_GFF3_genes.gff -sm 400+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,5487500+r+200 -ss 25+100,0+500,100+1+50+se -r ref -c oc -t par_mut
#rm -rf ./user_projects/*/1_intermediate_files

# ref oc par_nomut. SIMULADOR FUNCIONA
#./easymap -P caso_4 -w snp -sim -rs at -g TAIR10_GFF3_genes.gff -sm 400+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,5487500+r+200 -ss 25+100,0+500,100+1+50+se -r ref -c oc -t par_nomut
#rm -rf ./user_projects/*/1_intermediate_files

# ref oc f2wt. NO ME HA DEJADO easymap.py
#./easymap -P caso_5 -w snp -sim -rs at -g TAIR10_GFF3_genes.gff -sm 400+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,5487500+r+200 -ss 25+100,0+500,100+1+50+se -r ref -c oc -t f2wt
#rm -rf ./user_projects/*/1_intermediate_files

# noref bc par_mut. NO ME HA DEJADO easymap.py
#./easymap -P caso_6 -w snp -sim -rs at -g TAIR10_GFF3_genes.gff -sm 600+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,5487500+r+200 -ss 25+100,0+500,100+1+50+se -r noref -c bc -t par_mut
#rm -rf ./user_projects/*/1_intermediate_files

# noref bc f2wt
#./easymap -P caso_7 -w snp -sim -rs at -g TAIR10_GFF3_genes.gff -sm 600+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,5487500+r+200 -ss 25+100,0+500,100+1+50+se -r noref -c bc -t f2wt
#rm -rf ./user_projects/*/1_intermediate_files

# noref oc par_mut. SIMULADOR FUNCIONA
#./easymap -P caso_8 -w snp -sim -rs at -g TAIR10_GFF3_genes.gff -sm 600+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,5487500+r+200 -ss 25+100,0+500,100+1+50+se -r noref -c oc -t par_mut
#rm -rf ./user_projects/*/1_intermediate_files

# noref oc par_nomut. NO ME HA DEJADO easymap.py
#./easymap -P caso_9 -w snp -sim -rs at -g TAIR10_GFF3_genes.gff -sm 600+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,5487500+r+200 -ss 25+100,0+500,100+1+50+se -r noref -c oc -t par_nomut
#rm -rf ./user_projects/*/1_intermediate_files

# noref oc par_f2wt. NO ME HA DEJADO easymap.py
#./easymap -P caso_10 -w snp -sim -rs at -g TAIR10_GFF3_genes.gff -sm 600+e -sr 0,14-1,31-2,33-3,15-4,5-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,5487500+r+200 -ss 25+100,0+500,100+1+50+se -r noref -c oc -t f2wt
#rm -rf ./user_projects/*/1_intermediate_files


sleep 60s

#shutdown -h +5
#shutdown -h now
