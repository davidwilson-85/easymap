
echo CASO 2___________________________________________________________________________________________
./easymap -n minicaso2_lou -w snp -sim -r mid -g complete.gff -ed ref_bc_f2wt  -a TAIR10_gene_info.txt --low-stringency
#rm -rf ./user_projects/*/1_intermediate_files/sim_data/


echo CASO 2___________________________________________________________________________________________
./easymap -n minicaso2_hai -w snp -sim -r mid -g complete.gff -ed ref_bc_f2wt  -a TAIR10_gene_info.txt
#rm -rf ./user_projects/*/1_intermediate_files/sim_data/








exit

echo CASO 4___________________________________________________________________________________________
./easymap -n minicaso4 -w snp -sim -r insim -g complete.gff -ed ref_oc_parpol -a TAIR10_gene_info.txt
#rm -rf ./user_projects/*/1_intermediate_files/sim_data/

echo CASO 6___________________________________________________________________________________________
./easymap -n minicaso6 -w snp -sim -r insim -g complete.gff -ed noref_bc_f2wt -a TAIR10_gene_info.txt
#rm -rf ./user_projects/*/1_intermediate_files/sim_data/

echo CASO 7___________________________________________________________________________________________
./easymap -n minicaso7 -w snp -sim -r insim -g complete.gff -ed noref_oc_parmut -a TAIR10_gene_info.txt
#rm -rf ./user_projects/*/1_intermediate_files/sim_data/

echo CASO 1___________________________________________________________________________________________
./easymap -n minicaso1 -w snp -sim -r insim -g complete.gff -ed ref_bc_parmut  -a TAIR10_gene_info.txt
#rm -rf ./user_projects/*/1_intermediate_files/sim_data/

echo CASO 2___________________________________________________________________________________________
./easymap -n minicaso2 -w snp -sim -r insim -g complete.gff -ed ref_bc_f2wt  -a TAIR10_gene_info.txt
#rm -rf ./user_projects/*/1_intermediate_files/sim_data/

echo CASO 3___________________________________________________________________________________________
./easymap -n minicaso3 -w snp -sim -r insim -g complete.gff -ed ref_oc_parmut -a TAIR10_gene_info.txt
#rm -rf ./user_projects/*/1_intermediate_files/sim_data/


exit

echo CASO 4___________________________________________________________________________________________
./easymap -n minicaso4 -w snp -sim -r insim -g complete.gff -ed ref_oc_parpol -sm 50 -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,400000+100 -ss 30+100,0+500,100+1+50+se -a TAIR10_gene_info.txt
#rm -rf ./user_projects/*/1_intermediate_files/sim_data/

echo CASO 6___________________________________________________________________________________________
./easymap -n minicaso6 -w snp -sim -r insim -g complete.gff -ed noref_bc_f2wt -sm 50 -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,400000+100 -ss 30+100,0+500,100+1+50+se -a TAIR10_gene_info.txt
#rm -rf ./user_projects/*/1_intermediate_files/sim_data/

echo CASO 7___________________________________________________________________________________________
./easymap -n minicaso7 -w snp -sim -r insim -g complete.gff -ed noref_oc_parmut -sm 50 -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,400000+100 -ss 30+100,0+500,100+1+50+se -a TAIR10_gene_info.txt
#rm -rf ./user_projects/*/1_intermediate_files/sim_data/

echo CASO 1___________________________________________________________________________________________
./easymap -n minicaso1 -w snp -sim -r insim -g complete.gff -ed ref_bc_parmut  -sm 50 -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,400000+100 -ss 30+100,0+500,100+1+50+se -a TAIR10_gene_info.txt
#rm -rf ./user_projects/*/1_intermediate_files/sim_data/

echo CASO 2___________________________________________________________________________________________
./easymap -n minicaso2 -w snp -sim -r insim -g complete.gff -ed ref_bc_f2wt  -sm 50 -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,400000+100 -ss 30+100,0+500,100+1+50+se -a TAIR10_gene_info.txt
#rm -rf ./user_projects/*/1_intermediate_files/sim_data/

echo CASO 3___________________________________________________________________________________________
./easymap -n minicaso3 -w snp -sim -r insim -g complete.gff -ed ref_oc_parmut -sm 50 -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,400000+100 -ss 30+100,0+500,100+1+50+se -a TAIR10_gene_info.txt
#rm -rf ./user_projects/*/1_intermediate_files/sim_data/
