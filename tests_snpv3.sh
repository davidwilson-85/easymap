

echo 633___________________________________________________________________________________________
./easymap -n 633 -w snp -r at -g complete.gff -ed ref_bc_f2wt -a TAIR10_gene_info.txt -C Wachsman_633wt_SRR5029633_1.fq,Wachsman_633wt_SRR5029633_2.fq  -P Wachsman_633mt_SRR5029631_1.fq,Wachsman_633mt_SRR5029631_2.fq
rm -rf ./user_projects/*/1_intermediate_files/*.sam







exit










echo Schn-back___________________________________________________________________________________________
./easymap -n bc_schn_new_all_filters -w snp -r at -g complete.gff -mb ref -cr bc -co par_mut -a TAIR10_gene_info.txt -P BC.fg.reads1.fq,BC.fg.reads2.fq -C BC.bg.reads2.fq,BC.bg.reads1.fq
rm -rf ./user_projects/*/1_intermediate_files/*.sam


rm -rf ./user_projects/*/1_intermediate_files/BC.bg*
rm -rf ./user_projects/*/1_intermediate_files/BC.fg*



















echo Schn-back___________________________________________________________________________________________
./easymap -n BCschn -w snp -P BC.fg.reads1.fq,BC.fg.reads2.fq -C BC.bg.reads1.fq,BC.bg.reads2.fq -r at -g complete.gff -mb ref -cr bc -co par_mut  -a TAIR10_gene_info.txt
rm -rf ./user_projects/*/1_intermediate_files/*.sam

echo Schn-out____________________________________________________________________________________________
./easymap -n OCschn -w snp -P OC.fg.reads1.fq,OC.fg.reads2.fq -C OC.bg.reads1.fq,OC.bg.reads2.fq -r at -g complete.gff -mb ref -cr oc -co par_nomut  -a TAIR10_gene_info.txt




echo CASO 1___________________________________________________________________________________________
./easymap -n caso1 -w snp -sim -r insim -g complete.gff -mb ref -cr bc -co par_mut -sm 90 -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,30000+100 -ss 25+100,0+500,100+1+50+se -a TAIR10_gene_info.txt

echo CASO 2___________________________________________________________________________________________
./easymap -n caso2 -w snp -sim -r insim -g complete.gff -mb ref -cr bc -co f2wt -sm 90 -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,30000+100 -ss 25+100,0+500,100+1+50+se -a TAIR10_gene_info.txt

echo CASO 3___________________________________________________________________________________________
./easymap -n caso3 -w snp -sim -r insim -g complete.gff -mb ref -cr oc -co par_mut -sm 90 -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,30000+100 -ss 25+100,0+500,100+1+50+se -a TAIR10_gene_info.txt

echo CASO 4___________________________________________________________________________________________
./easymap -n caso4 -w snp -sim -r insim -g complete.gff -mb ref -cr oc -co par_nomut -sm 90 -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,30000+100 -ss 25+100,0+500,100+1+50+se -a TAIR10_gene_info.txt

echo CASO 5___________________________________________________________________________________________
./easymap -n caso5 -w snp -sim -r insim -g complete.gff -mb noref -cr bc -co par_mut -sm 90 -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,30000+100 -ss 25+100,0+500,100+1+50+se -a TAIR10_gene_info.txt

echo CASO 6___________________________________________________________________________________________
./easymap -n caso6 -w snp -sim -r insim -g complete.gff -mb noref -cr bc -co f2wt -sm 90 -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,30000+100 -ss 25+100,0+500,100+1+50+se -a TAIR10_gene_info.txt

echo CASO 7___________________________________________________________________________________________
./easymap -n caso7 -w snp -sim -r insim -g complete.gff -mb noref -cr oc -co par_mut -sm 90 -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,30000+100 -ss 25+100,0+500,100+1+50+se -a TAIR10_gene_info.txt




echo INS___________________________________________________________________________________________
./easymap -n BBI -w ins -r insim -sim -g complete.gff -sm 8 -i pbinprok2.fa -ss 5+100,0+500,100+1+50+pe -a TAIR10_gene_info.txt
rm -rf ./user_projects/*/1_intermediate_files/sim_data/

./easymap -n BBI -w ins -r mid -sim -g complete.gff -sm 7 -i pbinprok2.fa -ss 7+100,0+500,100+1+50+pe -a TAIR10_gene_info.txt
rm -rf ./user_projects/*/1_intermediate_files/sim_data/

./easymap -n BBI -w ins -r insim -sim -g complete.gff -sm 8 -i pbinprok2.fa -ss 5+100,0+500,100+1+50+se -a TAIR10_gene_info.txt
rm -rf ./user_projects/*/1_intermediate_files/sim_data/

./easymap -n BBI -w ins -r insim -sim -g complete.gff -sm 8 -i pbinprok2.fa -ss 15+100,0+500,100+1+50+se -a TAIR10_gene_info.txt
rm -rf ./user_projects/*/1_intermediate_files/sim_data/






exit



echo CASO 1___________________________________________________________________________________________
./easymap -n caso1 -w snp -sim -r insim -g complete.gff -mb ref -cr bc -co par_mut -sm 90 -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,30000+100 -ss 25+100,0+500,100+1+50+se -a TAIR10_gene_info.txt

echo CASO 2___________________________________________________________________________________________
./easymap -n caso2 -w snp -sim -r insim -g complete.gff -mb ref -cr bc -co f2wt -sm 90 -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,30000+100 -ss 25+100,0+500,100+1+50+se -a TAIR10_gene_info.txt

echo CASO 3___________________________________________________________________________________________
./easymap -n caso3 -w snp -sim -r insim -g complete.gff -mb ref -cr oc -co par_mut -sm 90 -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,30000+100 -ss 25+100,0+500,100+1+50+se -a TAIR10_gene_info.txt

echo CASO 4___________________________________________________________________________________________
./easymap -n caso4 -w snp -sim -r insim -g complete.gff -mb ref -cr oc -co par_nomut -sm 90 -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,30000+100 -ss 25+100,0+500,100+1+50+se -a TAIR10_gene_info.txt

echo CASO 5___________________________________________________________________________________________
./easymap -n caso5 -w snp -sim -r insim -g complete.gff -mb noref -cr bc -co par_mut -sm 90 -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,30000+100 -ss 25+100,0+500,100+1+50+se -a TAIR10_gene_info.txt

echo CASO 6___________________________________________________________________________________________
./easymap -n caso6 -w snp -sim -r insim -g complete.gff -mb noref -cr bc -co f2wt -sm 90 -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,30000+100 -ss 25+100,0+500,100+1+50+se -a TAIR10_gene_info.txt

echo CASO 7___________________________________________________________________________________________
./easymap -n caso7 -w snp -sim -r insim -g complete.gff -mb noref -cr oc -co par_mut -sm 90 -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,30000+100 -ss 25+100,0+500,100+1+50+se -a TAIR10_gene_info.txt















echo CASO 1___________________________________________________________________________________________
./easymap -n caso1 -w snp -P 1_se_reads_sample.fq -C 1_se_reads_control.fq -r insim -g complete.gff -mb ref -cr bc -co par_mut  -a TAIR10_gene_info.txt

echo CASO 2___________________________________________________________________________________________
./easymap -n caso2 -w snp -P 2_se_reads_sample.fq -C 2_se_reads_control.fq -r insim -g complete.gff -mb ref -cr bc -co f2wt  -a TAIR10_gene_info.txt

echo CASO 3___________________________________________________________________________________________
./easymap -n caso3 -w snp -P 3_se_reads_sample.fq -C 3_se_reads_control.fq -r insim -g complete.gff -mb ref -cr oc -co par_mut  -a TAIR10_gene_info.txt

echo CASO 4___________________________________________________________________________________________
./easymap -n caso4 -w snp -P 4_se_reads_sample.fq -C 4_se_reads_control.fq -r insim -g complete.gff -mb ref -cr oc -co par_nomut  -a TAIR10_gene_info.txt

echo CASO 5___________________________________________________________________________________________
./easymap -n caso5 -w snp -P 5_se_reads_sample.fq -C 5_se_reads_control.fq -r insim -g complete.gff -mb noref -cr bc -co par_mut  -a TAIR10_gene_info.txt

echo CASO 6___________________________________________________________________________________________
./easymap -n caso6 -w snp -P 6_se_reads_sample.fq -C 6_se_reads_control.fq -r insim -g complete.gff -mb noref -cr bc -co f2wt  -a TAIR10_gene_info.txt

echo CASO 7___________________________________________________________________________________________
./easymap -n caso7 -w snp -P 7_se_reads_sample.fq -C 7_se_reads_control.fq -r insim -g complete.gff -mb noref -cr oc -co par_mut  -a TAIR10_gene_info.txt














