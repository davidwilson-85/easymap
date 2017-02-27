start_time=`date +%s`

#Create input variables
my_ps= #[p, s], paired/single
my_rd=IonProton_reads_30X.fq #reads (single)
my_rf=ref_10kb_f_30X.fq #forward reads
my_rr=ref_10kb_r_30X.fq #reverse reads
my_gs=ref_10kb.fasta #genome sequence
my_ix=${my_gs%.*}  #name of index (the construct retains the part before the dot)


#Save path to bowtie2-build and bowtie2 in variable BT2
export location="$PWD" #to do: use pwd to make this automatic


#Create folder variables
f0=0_user-input
f1=1_intermediate-files
f2=2_logs
f3=3_workflow-output


#Execute vcf filter
python $location/filter.py -a $f1/clean_vcf.vcf -b $f1/filtered_variants.vcf 
