


echo INS___________________________________________________________________________________________
./easymap -n BBI -w ins -r mid -sim -g complete.gff -sm 8 -i pbinprok2.fa -ss 5+100,0+500,100+1+50+pe -a TAIR10_gene_info.txt
rm -rf ./user_projects/*/1_intermediate_files/sim_data/

./easymap -n BBI -w ins -r mid -sim -g complete.gff -sm 8 -i pbinprok2.fa -ss 15+100,0+500,100+1+50+pe -a TAIR10_gene_info.txt
rm -rf ./user_projects/*/1_intermediate_files/sim_data/

./easymap -n BBI -w ins -r mid -sim -g complete.gff -sm 8 -i pbinprok2.fa -ss 5+100,0+500,100+1+50+se -a TAIR10_gene_info.txt
rm -rf ./user_projects/*/1_intermediate_files/sim_data/

./easymap -n BBI -w ins -r mid -sim -g complete.gff -sm 8 -i pbinprok2.fa -ss 15+100,0+500,100+1+50+se -a TAIR10_gene_info.txt
rm -rf ./user_projects/*/1_intermediate_files/sim_data/




echo CASO 5___________________________________________________________________________________________
./easymap -n caso5 -w snp -sim -r mid -g complete.gff -mb noref -cr bc -co par_mut -sm 90 -sr 0,24-1,42-2,25-3,6-4,1-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,3000000+100 -ss 25+100,0+500,100+1+50+se -a TAIR10_gene_info.txt
rm -rf ./user_projects/*/1_intermediate_files/sim_data/


echo CASO 7___________________________________________________________________________________________
./easymap -n caso7 -w snp -sim -r mid -g complete.gff -mb noref -cr oc -co par_mut -sm 90 -sr 0,24-1,42-2,25-3,6-4,1-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,3000000+100 -ss 25+100,0+500,100+1+50+se -a TAIR10_gene_info.txt
rm -rf ./user_projects/*/1_intermediate_files/sim_data/

echo CASO 6___________________________________________________________________________________________
./easymap -n caso6 -w snp -sim -r mid -g complete.gff -mb noref -cr bc -co f2wt -sm 90 -sr 0,24-1,42-2,25-3,6-4,1-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,3000000+100 -ss 25+100,0+500,100+1+50+se -a TAIR10_gene_info.txt

echo CASO 2___________________________________________________________________________________________
./easymap -n caso2 -w snp -sim -r mid -g complete.gff -mb ref -cr bc -co f2wt -sm 90 -sr 0,24-1,42-2,25-3,6-4,1-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,3000000+100 -ss 25+100,0+500,100+1+50+se -a TAIR10_gene_info.txt
rm -rf ./user_projects/*/1_intermediate_files/sim_data/

echo CASO 4___________________________________________________________________________________________
./easymap -n caso4 -w snp -sim -r mid -g complete.gff -mb ref -cr oc -co par_nomut -sm 90 -sr 0,24-1,42-2,25-3,6-4,1-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,3000000+100 -ss 25+100,0+500,100+1+50+se -a TAIR10_gene_info.txt
rm -rf ./user_projects/*/1_intermediate_files/sim_data/

echo CASO 1___________________________________________________________________________________________
./easymap -n caso1 -w snp -sim -r mid -g complete.gff -mb ref -cr bc -co par_mut -sm 90 -sr 0,24-1,42-2,25-3,6-4,1-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,3000000+100 -ss 25+100,0+500,100+1+50+se -a TAIR10_gene_info.txt

echo CASO 3___________________________________________________________________________________________
./easymap -n caso3 -w snp -sim -r mid -g complete.gff -mb ref -cr oc -co par_mut -sm 90 -sr 0,24-1,42-2,25-3,6-4,1-5,2/0,24-1,42-2,25-3,6-4,1-5,2+1,3000000+100 -ss 25+100,0+500,100+1+50+se -a TAIR10_gene_info.txt
rm -rf ./user_projects/*/1_intermediate_files/sim_data/