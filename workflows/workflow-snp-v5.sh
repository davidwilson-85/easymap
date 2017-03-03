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
my_mode=$5 														#[pe, se], paired/single  
my_p_mode=se														#						<------------------------------------------------------------------------
my_rd=$7											 			#reads (single)
my_rf=$8 														#forward reads
my_rr=$9												 		#reverse reads 			
my_p_rd=col-lab.fq											 			#reads (single) parent	<------------------------------------------------------------------------
my_p_rf=none 														#forward reads parent	<------------------------------------------------------------------------
my_p_rr=none												 		#reverse reads parent	<------------------------------------------------------------------------
my_gs=gnm_ref_merged/genome.fa 									#genome sequence
my_ix=genome_index 							
my_gff=${10}													#Genome feature file
my_ann=${11}													
my_rrl=250 														#Regulatory region length
my_log_file=$1
my_mut=snp  													#my_mut takes the values 'snp' in this workflow and 'lin' in the large insertions workflow, for the execution of the graphic output module



my_mutbackgroud=ref 											#ref / noref : genetic background of the mutation									<------------------------------------------------------------------------
my_cross=oc														#oc / bc : f2 obtained by outcross or backcross 									<------------------------------------------------------------------------
my_pseq=mut  													#mut / nomut : sequenced parental provided is the mutagenized one or the other		<------------------------------------------------------------------------





#Define the folders in the easymap directory 
f0=$project_name/0_input
f1=$project_name/1_intermediate_files
f2=$project_name/2_logs
f3=$project_name/3_workflow_output


#Save path to bowtie2-build and bowtie2 in variable BT2
export location="$PWD" 


'''

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


##################################################################################################################################################################################
#																																												 #
#																																												 #
#																	F2 FQ PROCESSING 																							 #
#																																												 #
#																																												 #
##################################################################################################################################################################################

if [ $my_mode == se ] 
then
	#Execute bowtie2 unpaired to align raw F2 reads to genome 
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
	#Execute bowtie2 paired to align raw F2 reads to genome 
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
	python $location/scripts_snp/groomer/vcf-groomer.py -a $f1/raw_variants.vcf -b $f1/F2_raw.va 

} || {
	echo 'error: vcf-groomer.py' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo VCF grooming finished >> $my_log_file
'''
'''
#Execute vcf filter
{
	python $location/scripts_snp/filter/variants-filter.py -a $f1/F2_raw.va -b $f1/F2_filtered.va -step 1

} || {
	echo 'error: variants-filter.py' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo VCF filter finished >> $my_log_file
'''

'''
##################################################################################################################################################################################
#																																												 #
#																																												 #
#																	PARENTAL FQ PROCESSING																						 #
#																																												 #
#																																												 #
##################################################################################################################################################################################

if [ $my_mode == se ] 
then
	#Execute bowtie2 unpaired to align raw F2 reads to genome 
	{
		$location/bowtie2/bowtie2 --very-sensitive --mp 3,2 -x $f1/$my_ix -U $f0/$my_p_rd -S $f1/alignment1P.sam 2> $f2/bowtie2_std2.txt

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
	#Execute bowtie2 paired to align raw F2 reads to genome 
	{
		$location/bowtie2/bowtie2 --very-sensitive -X 1000 --mp 3,2 -x $f1/$my_ix -1 $f0/$my_p_rf -2 $f0/$my_p_rr -S $f1/alignment1P.sam 2> $f2/bowtie2_std2.txt

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
	$location/samtools1/samtools sort $f1/alignment1P.sam > $f1/alignment1P.bam

} || {
	echo 'Error: SAM to BAM' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo SAM to BAM finished >> $my_log_file


#Variant calling
{
	$location/samtools1/samtools mpileup  -B -t DP,ADF,ADR -vuo $f1/raw_p_variants_temp.vcf  -f $f0/$my_gs $f1/alignment1P.bam  2> $f2/mpileup_std.txt
	$location/bcftools-1.3.1/bcftools call -vmO v -vo $f1/raw_p_variants.vcf $f1/raw_p_variants_temp.vcf 2> $f2/bcf_std.txt
	
} || {
	echo 'Error during variant-calling' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo Variant calling finished >> $my_log_file


#Groom vcf
{
	python $location/scripts_snp/groomer/vcf-groomer.py -a $f1/raw_p_variants.vcf -b $f1/parental_raw.va 

} || {
	echo 'error: vcf-groomer.py' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo VCF grooming finished >> $my_log_file


#Execute vcf filter
{
	python $location/scripts_snp/filter/variants-filter.py -a $f1/parental_raw.va -b $f1/parental_filtered.va -step 1

} || {
	echo 'error: variants-filter.py' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo VCF filter finished >> $my_log_file
'''


'''
##################################################################################################################################################################################
#																																												 #
#																																												 #
#																			DATA ANALYSIS 																						 #
#																																												 #
#																																												 #
##################################################################################################################################################################################

#_____________________________________________________________________________OPERATIONS__________________________________________________________________________________________
#Setting up operation mode: 

#Mutation in refference backgroud, outcross with non-refference background, sequencing refference parental
if [ $my_mutbackgroud == ref ] && [ $my_pseq == mut ] && [ $my_cross == oc ]
then
	my_operation_mode=A
fi

#Mutation in refference backgroud, outcross with non-refference background, sequencing non-refference parental
if [ $my_mutbackgroud == ref ] && [ $my_pseq == nomut ] && [ $my_cross == oc ]
then
	my_operation_mode=I
fi

#Mutation in refference backgroud, backcross with non-refference background, sequencing refference parental
if [ $my_mutbackgroud == ref ] && [ $my_pseq == mut ] && [ $my_cross == bc ]
then
	my_operation_mode=A
fi

#Mutant in non-reference background, outcross with reference, sequencing non-refference parental
if [ $my_mutbackgroud == noref ] && [ $my_pseq == mut ] && [ $my_cross == oc ]
then
	my_operation_mode=I
fi

#Execute vcf operations
{
	python $location/scripts_snp/operations/variants-operations.py -a $f1/F2_filtered.va -b $f1/parental_filtered.va -c $f1/F2_parental_comparison.va -mode $my_operation_mode -primary 1  

} || {
	echo 'error: operations.py' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo VCF operations finished >> $my_log_file


'''

#_________________________________________________________________________________ANALYSIS_____________________________________________________________________________________________

#Setting up analysis mode:

#Mutation in refference backgroud, outcross with non-refference background, sequencing refference parental
if [ $my_cross == oc ] 
then
	my_analysis_mode=out
fi

#Mutation in refference backgroud, outcross with non-refference background, sequencing non-refference parental
if [ $my_cross == bc ] 
then
	my_analysis_mode=back
fi

#Execute vcf analysis 
{
	python $location/scripts_snp/analysis/map-mutation.py -fichero $f1/F2_parental_comparison.va -fasta $f0/$my_gs -mode $my_analysis_mode -window_size 500000 -window_space 500000 -output $f1/map_info.txt -parental_modality $my_mutbackgroud -correction_factor 0.6


} || {
	echo 'error: map-mutation.py' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo VCF analysis finished >> $my_log_file

#__________________________________________________________________________________FILTER____________________________________________________________________________________________


#Execute vcf filter, selecting snps in the candidate region defined by map-mutation.py, with an alelic frequency > 0.8 and corresponding to EMS mutations
{
	python $location/scripts_snp/filter/variants-filter.py -a $f1/F2_filtered.va -b $f1/F2_filtered2.va -step 2 -cand_reg_file $f1/map_info.txt -af_min 0.8 -mut_type EMS

} || {
	echo 'error: variants-filter.py' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo VCF filter finished >> $my_log_file



#__________________________________________________________________________________OPERATIONS________________________________________________________________________________________


#Setting up operation mode: 

#Mutation in refference backgroud, outcross with non-refference background, sequencing refference parental
if [ $my_mutbackgroud == ref ] && [ $my_pseq == mut ] && [ $my_cross == oc ]
then
	my_operation_mode=A
fi

#Mutation in refference backgroud, outcross with non-refference background, sequencing non-refference parental
if [ $my_mutbackgroud == ref ] && [ $my_pseq == nomut ] && [ $my_cross == oc ]
then
	my_operation_mode=A
fi

#Mutation in refference backgroud, backcross with non-refference background, sequencing refference parental
if [ $my_mutbackgroud == ref ] && [ $my_pseq == mut ] && [ $my_cross == bc ]
then
	my_operation_mode=A
fi

#Mutant in non-reference background, outcross with reference, sequencing non-refference parental
if [ $my_mutbackgroud == noref ] && [ $my_pseq == mut ] && [ $my_cross == oc ]
then
	my_operation_mode=A
fi


#Execute vcf operations
{
	python $location/scripts_snp/operations/variants-operations.py -a $f1/F2_filtered2.va -b $f1/parental_filtered.va -c $f1/final_variants.txt -mode $my_operation_mode -primary 1  

} || {
	echo 'error: operations.py' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo vcf operations finished >> $my_log_file

#__________________________________________________________________________________VARANALYZER INPUT_________________________________________________________________________________


#snp-to-varanalyzer.py
{
	python $location/scripts_snp/snp-to-varanalyzer.py -a $f1/final_variants.txt -b $f1/final_variants2.txt	
	
} || {
	echo 'error: snp-to-varanalyzer.py' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo snp-to-varanalyzer.py finished. >> $my_log_file



#__________________________________________________________________________________VARANALYZER_______________________________________________________________________________________


#varanalyzer
{
	python $location/varanalyzer/varanalyzer_v1.py -itp lim -con $f0/$my_gs -gff $f0/$my_gff -var $f1/final_variants2.txt -rrl $my_rrl   -pname $project_name

} || {
	echo 'error: varanalyzer_v1.py' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo Varanalyzer finished. >> $my_log_file


##################################################################################################################################################################################
#																																												 #
#																																												 #
#																				REPORT 	 																						 #
#																																												 #
#																																												 #
##################################################################################################################################################################################


#__________________________________________________________________________________GRAPHIC OUTPUT_____________________________________________________________________________________


#Graphic output
{
	python $location/graphic_output/graphic-output-v3.py -my_mut $my_mut -asnp $f1/F2_parental_comparison.va -bsnp $f0/$my_gs -rrl $my_rrl -iva $2/3_workflow_output/variants.txt -gff $f0/$my_gff -pname $2  -cross $my_cross
	
} || {
	echo 'error: graphic-output-v3.py' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo Graphic output created. >> $my_log_file

echo $exit_code


#HTML file