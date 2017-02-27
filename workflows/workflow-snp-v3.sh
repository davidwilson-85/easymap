#!/bin/bash
# This is the command sent by 'master.sh':
# ./workflow-x.sh $my_log_file $project_name $workflow $data_source $lib_type $ins_seq $read_s $read_f $read_r $gff_file $ann_file
#
# This command is always the same, regardless of the workflow
#
# Fields that can be equal to 'n/p' (= not provided) are the following: $ins_seq, $read_s, $read_f, $read_r, $ann_file
# If $ins_seq = n/p, that is because $workflow = snp. Therefore, no insertion seq is needed and $ins_seq is ignored by hte program.
# If $read_s = n/p, that is because $lib_type = pe, so it is ignored by the program.
# If $read_f and $read_r = n/p, that is because $lib_type = se, so it is ignored by the program.
# If $ann_file = n/p, this is because user did no have it. This program has the deal with this: if data not provided, simply do not
# include gene annotated info to the report.
#
# $my_log_file		>	$1
# $project_name		>	$2
# $workflow			>	$3
# $data_source		>	$4
# $lib_type			>	$5
# $ins_seq			>	$6
# $read_s			>	$7
# $read_f			>	$8
# $read_r			>	$9
# $gff_file			>	${10}
# $ann_file			>	${11}
#


# Set 'exit_code' (flag variable) to 0
exit_code=0

# For tests:
my_log_file=$1
echo $0 $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} >> $my_log_file


start_time=`date +%s`


#Create input variables
my_log_file=$1
project_name=$2
my_mode=$5 														#[P, S], paired/single  <------------------------------------------------------------------------
my_rd=$7											 			#reads (single)
my_rf=$8 														#forward reads
my_rr=$9												 		#reverse reads
my_is=$6		 												#insertion sequence
my_gs=gnm_ref_merged/genome.fa 									#genome sequence
my_ix=genome_index 							
my_gff=${10}													#Genome feature file
my_ann=${11}													
my_rrl=250 														#Regulatory region length
my_log_file=$1

my_vcf2=myvcf2.vcf 												#empty file 
my_mut=snp  													#my_mut takes the values 'snp' in this workflow and 'lin' in the large insertions workflow, for the execution of the graphic output module


#Define the folders in the easymap directory 
f0=$project_name/0_input
f1=$project_name/1_intermediate_files
f2=$project_name/2_logs
f3=$project_name/3_workflow_output


#Save path to bowtie2-build and bowtie2 in variable BT2
export location="$PWD" 


#Execute bowtie2-build on genome sequence 
{
	$location/bowtie2/bowtie2-build $f0/$my_gs $f1/$my_ix 1> $f2/bowtie2-build_std1.txt 2> $f2/bowtie2-build_std2.txt

} || {
	echo 'bowtie2-build on genome sequence returned an error. See log files.' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo Bowtie2-build finished >> $my_log_file

if [ $my_mode == se ] 
then
	#Execute bowtie2 unpaired to align raw reads to insertion sequence
	{
		$location/bowtie2/bowtie2 --very-sensitive --mp 3,2 -x $f1/$my_ix -U $f0/$my_rd -S $f1/alignment1.sam 2> $f2/bowtie2_std2.txt

	} || {
		echo 'bowtie2 returned an error. See log files.' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo Bowtie2 finished >> $my_log_file
fi


if [ $my_mode == pe ] 
then
	#Execute bowtie2 paired to align raw reads to insertion sequence
	{
		$location/bowtie2/bowtie2 --very-sensitive -X 1000 --mp 3,2 -x $f1/$my_ix -1 $f0/$my_rf -2 $f0/$my_rr -S $f1/alignment1.sam 2> $f2/bowtie2_std2.txt

	} || {
		echo 'bowtie2 returned an error. See log files.' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo Bowtie2 finished >> $my_log_file
fi

#SAM to BAM
{
	$location/samtools1/samtools sort $f1/alignment1.sam > $f1/alignment1.bam

} || {
	echo 'Error: SAM to BAM' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo SAM to BAM finished >> $my_log_file


#Variant calling
{
	$location/samtools1/samtools mpileup  -B -t DP,ADF,ADR -vuo $f1/raw_variants_temp.vcf  -f $f0/$my_gs $f1/alignment1.bam  2> $f2/mpileup_std.txt
	$location/bcftools-1.3.1/bcftools call -vmO v -vo $f1/raw_variants.vcf $f1/raw_variants_temp.vcf 2> $f2/bcf_std.txt
	
} || {
	echo 'Error during variant-calling' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo Variant calling finished >> $my_log_file


#Groom vcf
{
	python $location/scripts_snp/groomer/vcf_groomer.py -a $f1/raw_variants.vcf -b $f1/clean_vcf.vcf 

} || {
	echo 'error: vcf_groomer.py' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo VCF grooming finished >> $my_log_file


#Execute vcf filter
{
	python $location/scripts_snp/filter/filter.py -a $f1/clean_vcf.vcf -b $f1/filtered_variants.vcf 

} || {
	echo 'error: filter.py' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo VCF filter finished >> $my_log_file


#Execute vcf operations
{
	python $location/scripts_snp/operations/vcf_operations.py -a $f1/filtered_variants.vcf -b $f0/$my_vcf2 -c $f3/final_variants.vcf -mode N -primary 1  

} || {
	echo 'error: operations.py' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo VCF operations finished >> $my_log_file


#Execute vcf analysis 


#varanalyzer
# {
#	python $location/varanalyzer/varanalyzer_v1.py -itp lim -con $f0/$my_gs -gff $f0/$my_gff -var $f1/INPUT.txt -rrl $my_rrl    #	<------------------- Input para varanalyzer
	
#} || {
#	echo 'error: varanalyzer_v1.py' >> $my_log_file
#	exit_code=1
#	echo $exit_code
#	exit
#}
#echo finished. >> $my_log_file


#Graphic output
{
	python $location/graphic_output/graphic-output-v3.py -my_mut $my_mut -asnp $f3/final_variants.vcf -bsnp $f0/$my_gs -rrl $my_rrl -iva $2/3_workflow_output/varanalyzer_output.txt -gff $f0/$my_gff -pname $2          #<--------------------------------- temporal
	
} || {
	echo 'error: graphic-output-v3.py' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo Graphic output created. >> $my_log_file

echo $exit_code