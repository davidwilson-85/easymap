echo CASO 1___________________________________________________________________________________________
./easymap -n minicaso1 -w snp -sim -r insim -g complete.gff -mb ref -cr bc -co par_mut -sm 90 -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,30000+100 -ss 30+100,0+500,100+1+50+se -a TAIR10_gene_info.txt
rm -rf ./user_projects/*/1_intermediate_files/sim_data/

echo CASO 2___________________________________________________________________________________________
./easymap -n minicaso2 -w snp -sim -r insim -g complete.gff -mb ref -cr bc -co f2wt -sm 90 -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,30000+100 -ss 30+100,0+500,100+1+50+se -a TAIR10_gene_info.txt
rm -rf ./user_projects/*/1_intermediate_files/sim_data/

echo CASO 3___________________________________________________________________________________________
./easymap -n minicaso3 -w snp -sim -r insim -g complete.gff -mb ref -cr oc -co par_mut -sm 90 -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,30000+100 -ss 30+100,0+500,100+1+50+se -a TAIR10_gene_info.txt
rm -rf ./user_projects/*/1_intermediate_files/sim_data/

echo CASO 4___________________________________________________________________________________________
./easymap -n minicaso4 -w snp -sim -r insim -g complete.gff -mb ref -cr oc -co par_nomut -sm 90 -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,30000+100 -ss 30+100,0+500,100+1+50+se -a TAIR10_gene_info.txt
rm -rf ./user_projects/*/1_intermediate_files/sim_data/

echo CASO 5___________________________________________________________________________________________
./easymap -n minicaso5 -w snp -sim -r insim -g complete.gff -mb noref -cr bc -co par_mut -sm 90 -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,30000+100 -ss 30+100,0+500,100+1+50+se -a TAIR10_gene_info.txt
rm -rf ./user_projects/*/1_intermediate_files/sim_data/

echo CASO 6___________________________________________________________________________________________
./easymap -n minicaso6 -w snp -sim -r insim -g complete.gff -mb noref -cr bc -co f2wt -sm 90 -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,30000+100 -ss 30+100,0+500,100+1+50+se -a TAIR10_gene_info.txt
rm -rf ./user_projects/*/1_intermediate_files/sim_data/

echo CASO 7___________________________________________________________________________________________
./easymap -n minicaso7 -w snp -sim -r insim -g complete.gff -mb noref -cr oc -co par_mut -sm 90 -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,30000+100 -ss 30+100,0+500,100+1+50+se -a TAIR10_gene_info.txt
rm -rf ./user_projects/*/1_intermediate_files/sim_data/













