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

# Set location of log file
my_log_file=$1



start_time=`date +%s`


my_mut=lin 			#my_mut takes the values 'lin' in this workflow and 'snp' in the snp workflow, for the execution of the graphic output module

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
my_ix2=insertion_index 						
my_gff=${10}													#Genome feature file
my_ann=${11}													
my_rrl=250 														#Regulatory region length


#Define the folders in the easymap directory 
f0=user_data
f1=$project_name/1_intermediate_files
f2=$project_name/2_logs
f3=$project_name/3_workflow_output

# Write PID to status file
my_status_file=$project_name/$f2/status
echo 'pid workflow '$BASHPID >> $my_status_file

#Save path to bowtie2-build and bowtie2 in variable BT2
export location="$PWD" 

#Execute bowtie2-build on insertion and genome sequence 
{
	$location/bowtie2/bowtie2-build $f0/$my_is $f1/$my_ix2 1> $f2/bowtie2-build_std1.txt 2> $f2/bowtie2-build_std2.txt
	
} || {
	echo 'bowtie2-build on insertion sequence returned an error. See log files.' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo bowtie2-build insertion index finished. >> $my_log_file

{
	$location/bowtie2/bowtie2-build $f1/$my_gs $f1/$my_ix 1> $f2/bowtie2-build2_std3.txt 2> $f2/bowtie2-build2_std4.txt
	
} || {
	echo 'bowtie2-build on genome sequence returned an error. See log files.' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo bowtie2-build genome index finished. >> $my_log_file


#Execute bowtie2 paired to align raw reads to insertion sequence
if [ $my_mode == 'pe' ]
then  
	{
		$location/bowtie2/bowtie2 -x $f1/$my_ix2 -1 $my_rf -2 $my_rr -S $f1/alignment1.sam 2> $f2/bowtie2_std2.txt
	
	} || {
		echo 'bowtie2 on the insertion sequence returned an error. See log files.' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}

	echo bowtie2 paired finished. >> $my_log_file
fi

if [ $my_mode == 'se' ] 
then
	{
		$location/bowtie2/bowtie2 --very-sensitive --mp 3,2 -x $f1/$my_ix -U $my_rd -S $f1/alignment1.sam 2> $f2/bowtie2_std2.txt
	
	} || {
		echo 'bowtie2 on the insertion sequence returned an error. See log files.' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}

	echo bowtie2 unpaired finished. >> $my_log_file
fi

#Check output .sam file format
{
	python $location/scripts_ins/sam_file_check/sam-file-check.py -a $f1/alignment1.sam 
	
} || {
	echo 'intermediate .sam file incorrect. See log files.' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo Correct .sam file. >> $my_log_file


if [ $my_mode == 'pe' ]
then  
	#Execute filter1
	{
		python $location/scripts_ins/filter_1/filter1.py -a $f1/alignment1.sam -b $f1/output_F1.fq
	
	} || {
		echo 'error: filter1.py' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo First filter finished. >> $my_log_file


	#Execute bowtie2 to align filtered reads to genome sequence
	{
		$location/bowtie2/bowtie2 -x $f1/$my_ix -U $f1/output_F1.fq -S $f1/alignment2.sam 2> $f2/bowtie2_std4.txt
	
	} || {
		echo 'bowtie2 on the genome sequence returned an error. See log files.' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo bowtie2 finished. >> $my_log_file


	#Check output .sam file format
	{
		python $location/scripts_ins/sam_file_check/sam-file-check.py -a $f1/alignment2.sam 
	
	} || {
		echo 'intermediate .sam file incorrect. See log files.' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo Correct .sam file. >> $my_log_file
fi


if [ $my_mode == 'pe' ]
then  
	#Execute bowtie2 to make a local aligment of the reads with the insertion
	{
		$location/bowtie2/bowtie2 --local -x $f1/$my_ix2 -1 $my_rf -2 $my_rr -S $f1/alignment3.sam 2> $f2/bowtie2_local_std1.txt
	
	} || {
		echo 'bowtie2 local alignment to the insertion sequence returned an error. See log files.' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo bowtie2 finished. >> $my_log_file
fi


if [ $my_mode == 'se' ]
then  	
	#Execute bowtie2 to make a local aligment of the reads with the insertion
	{
		$location/bowtie2/bowtie2 --local -x $f1/$my_ix2 -U $my_rd -S $f1/alignment3.sam 2> $f2/bowtie2_local_std1.txt
	
	} || {
		echo 'bowtie2 local alignment to the insertion sequence returned an error. See log files.' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo bowtie2 finished. >> $my_log_file
fi


#Execute filter2
{
	python $location/scripts_ins/filter_2/filter2.py -a $f1/alignment3.sam -b $f1/output_F2.fq

} || {
	echo 'error: filter2.py' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo Second filter finished. >> $my_log_file


#Execute bowtie2 to align filtered reads to genome sequence
{
	$location/bowtie2/bowtie2 --local -x $f1/$my_ix -U $f1/output_F2.fq -S $f1/alignment4.sam 2> $f2/bowtie2_local_std2.txt

} || {
	echo 'bowtie2 local alignment to the genome sequence returned an error. See log files.' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo bowtie2 finished. >> $my_log_file




#Count read depth and find candidate region
if [ $my_mode == 'pe' ]
then  
	{
		python $location/scripts_ins/analysis/paired-analysis.py -a $f1/alignment2.sam -b $f1/output_analysis.txt -c $f1/$my_gs 

	} || {
		echo 'error: paired-analysis.py' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	{
		python $location/scripts_ins/analysis/local-analysis.py -a $f1/alignment4.sam -b $f1/output_analysis.txt -c $f1/$my_gs -m $my_mode

	} || {
		echo 'error: local-analysis.py' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}

fi	

if [ $my_mode == 'se' ]
then  
	{
		python $location/scripts_ins/analysis/local-analysis.py -a $f1/alignment4.sam -b $f1/output_analysis.txt -c $f1/$my_gs -m $my_mode

	} || {
		echo 'error: local-analysis.py' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
fi
echo Analysis finished. >> $my_log_file



#Sort insertions
{
	python $location/scripts_ins/sort_insertions/sort.py -a $f1/output_analysis.txt -b $f1/$my_gs -c $f1/output_ordered.csv -d $f3/sorted_insertions.txt -m $my_mode
	
} || {
	echo 'error: sort.py' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo Insertions sorted. >> $my_log_file


#ma-input.py
{
	python $location/scripts_ins/ins_to_varanalyzer/ins-to-varanalyzer.py -a $f3/sorted_insertions.txt -b $f1/ins-to-varanalyzer.txt
	
} || {
	echo 'error: ins-to-varanalyzer.py' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo Done. >> $my_log_file


#varanalyzer
{
	python $location/varanalyzer/varanalyzer_v1.py -itp lim -con $f1/$my_gs -gff $f0/$my_gff -var $f1/ins-to-varanalyzer.txt -rrl $my_rrl -pname $project_name
	
} || {
	echo 'error: varanalyzer_v1.py' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo Done. >> $my_log_file


#Graphic output
{
	python $location/graphic_output/graphic-output.py -my_mut $my_mut -a $f3/sorted_insertions.txt -b $f1/$my_gs -f $f3/report.html -m $my_mode	-gff $f0/$my_gff  -iva $f3/variants.txt -rrl $my_rrl -pname $project_name
	
} || {
	echo 'error:graphic-output.py' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo Graphic output created. >> $my_log_file

echo run time is $(expr `date +%s` - $start_time) s >> $my_log_file


#Report
#firefox ./project/3_workflow_output/report.html


echo $exit_code
