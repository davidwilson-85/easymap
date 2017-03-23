#!/bin/bash


# Check all the files and report to log if they are correct or not, but always continue
# ultil all input files checked. This way, the user can know in a single run all the files
# that are not correct.
# At the same time, if any file is incorrect, set exit_code=1. master.sh retrieves this value and
# stops the workflow. Then, user can see in log stream which file was incorrect.
#
# TO DO: Maybe add error checks after the execution of each python program, as in the rest of the
# workflows.
#
#


#######################################################################################################
# This block of code does some steps previous to the analysis
# 
# Description:
# 
# 
# 
#######################################################################################################


# Set 'exit_code' (flag variable) to it's initial value (0)
exit_code=0

# This is the command sent by 'master.sh':
# ./process-input.sh $my_log_file $project_name $analysis_type $data_source $lib_type $ins_seq $read_s $read_f $read_r $gff_file $ann_file`
# example: ./process-input.sh project/log.log project ins sim se ins.fa fq_se.fq fq_1.fq fq_2.fq gff.gff ann.ann

# Store the location of each folder in a variable
f0=0_input
f1=1_intermediate_files

# Get command arguments and assign them to variables
# The reference to the template genome is not passed from 'master.sh' because it is static
my_log_file=$1
project_name=$2
analysis_type=$3
data_source=$4
lib_type=$5
ins_seq=$project_name/$f0/$6
read_s=$project_name/$f0/$7
read_f=$project_name/$f0/$8
read_r=$project_name/$f0/$9
gff_file=$project_name/$f0/${10}
ann_file=$project_name/$f0/${11}
ann_option=${11}


# Establish locations of reference genome
ref_seqs_dir=$project_name/$f0/gnm_ref
ref_seqs_merged_dir=$project_name/$f1/gnm_ref_merged
ref_seqs_merged_file=$project_name/$f1/gnm_ref_merged/genome.fa


#######################################################################################################
# This block of code does the initial processing of the input data provided by the user
#                                                                                        
# Description:
# This workflow calls several times 'process-input/verify-input.py', each time with a single option.
# When called this way, 'process-input/verify-input.py' prints a single value. In most cases: 
# 0=input is ok, 1=input is not ok. In other cases: 0=input ok, 1=input has problem,
# 2=input has another problem. This bash program interprets the value and decides whether
# to continue or to exit prematurely.
# Options are: -gnm, -ins, -fq, -gff, -ann, -match
# This workflow also calls 'fasta-merger.py'. The script gets all the fasta files in the input
# folder and concatenates their content, which is written to a new file. The purpose of creating
# an all-contigs-in-one-place file is multiple: First, it eases comparing contig names between
# fasta and gff inputs; second, this file s required for varanalyzer, which is called downstream;
# third, it is easier to run bowtie2-build.
#######################################################################################################


# Check fasta input(s)
fa=`python process_input/verify-input.py -gnm $ref_seqs_dir`

if [ $fa == 0 ]
then
	{
		echo $(date)": Genome fasta input check passed." >> $my_log_file
	}
else
	{
		echo $(date)": Genome fasta input check failed. One or more genome fasta inputs are empty or have an incorrect format. Please provide new file(s)." >> $my_log_file
		exit_code=1
		echo $exit_code
		exit # Exit because fasta files are required for fasta-concat.py and for fasta-gff comparison
		     # In any other check I do not exit prematurely of process-input.sh so the program checks
		     # all the files and therefore several incorrect files can be revealed in a single run.
	}
fi

# Concatenate fasta input and store in a single file

{
	python process_input/fasta-concat.py -in_dir $ref_seqs_dir -out_dir $ref_seqs_merged_dir
} || {
	echo echo $(date)": Processing of genome fasta input failed: fasta-concat.py could not concatenate fasta files into one file." >> $my_log_file
	exit_code=1
}
echo $(date)": Processing of genome fasta input completed." >> $my_log_file


# Check fasta input with insertion sequence
# Do this only if analyzing insertions

if [ $analysis_type == 'ins' ]
then
	{
		fa=`python process_input/verify-input.py -ins $ins_seq`
		
		if [ $fa == 0 ]
		then
			{
				echo $(date)": Insertion fasta input check passed." >> $my_log_file
			}
		else
			{
				echo $(date)": Insertion fasta input check failed. File is empty or has an incorrect format. Please provide a new file." >> $my_log_file
				exit_code=1
			}
		fi
	}
fi

# Check fastq input file(s)
# Do this only if data source is exp (reads provided by the user, not simulated)

if [ $data_source == 'exp' ]
then
	{
		if [ $lib_type == 'se' ]
		then  
			{
				fq=`python process_input/verify-input.py -fq $read_s`
				
				if [ $fq == 0 ]
				then
					{
						echo $(date)": Single-end fastq input passed." >> $my_log_file
					}
				else
					{
						echo $(date)": Single-end fastq input failed. File is empty or has an incorrect format. Please provide a new file." >> $my_log_file
						exit_code=1
					}
				fi
			}
		fi

		if [ $lib_type == 'pe' ]
		then  
			{
				fq=`python process_input/verify-input.py -fq $read_f`
				
				if [ $fq == 0 ]
				then
					{
						echo $(date)": Paired-end forward fastq input passed." >> $my_log_file
					}
				else
					{
						echo $(date)": Paired-end forward fastq input failed. File is empty or has an incorrect format. Please provide a new file." >> $my_log_file
						exit_code=1
					}
				fi

				fq=`python process_input/verify-input.py -fq $read_r`
				
				if [ $fq == 0 ]
				then
					{
						echo $(date)": Paired-end reverse fastq input passed." >> $my_log_file
					}
				else
					{
						echo $(date)": Paired-end forward fastq input failed. File is empty or has an incorrect format. Please provide a new file." >> $my_log_file
						exit_code=1
					}
				fi
			}
		fi
	}
fi


# Check gff input
gff=`python process_input/verify-input.py -gff $gff_file`

if [ $gff == 0 ]
then
	{
		echo $(date)": GFF3 input check passed." >> $my_log_file
	}
else
	{
		echo $(date)": GFF3 input check failed. File is empty or has an incorrect format. Please provide a new file." >> $my_log_file
		exit_code=1
	}
fi


# Check gene funtional annotation input
# Do this only if user has provided file (it's optional)

if [ $ann_option != 'n/p' ]
then
	{

		ann=`python process_input/verify-input.py -ann $ann_file`
		echo 'Checking gene functional annotation input...' >> $my_log_file
		if [ $ann == 0 ]
		then
			{
				echo $(date)": Gene annotation file check passed." >> $my_log_file
			}
		else
			{
				echo $(date)": Gene annotation file check failed. File is empty or has an incorrect format. Please replace input gene functional annotation file or turn off the gene annotation option." >> $my_log_file
				exit_code=1
			}
		fi
	}
fi


# Check contigs match between fasta and gff3 files
match=`python process_input/verify-input.py -fa_match $ref_seqs_merged_file -gff_match $gff_file`

if [ $match == 0 ]
then
	{
		echo $(date)": Contigs match check between FASTA and GFF3 inputs passed." >> $my_log_file
	}
elif [ $match == 1 ]
then
	{
		echo $(date)": Contigs match check between FASTA and GFF3 inputs failed. FASTA input, GFF3 input, or both are empty. Please provide new files." >> $my_log_file
		exit_code=1
	}
elif [ $match == 2 ]
then
	{
		echo $(date)": Contigs match check between FASTA and GFF3 inputs failed. The contig names in the FASTA are not in GFF3. Please provide new files." >> $my_log_file
		exit_code=1
	}
elif [ $match == 3 ]
then
	{
		echo $(date)": Contigs match check between FASTA and GFF3 inputs failed. Some contig names in the GFF3 are not in FASTA. The process will proceed, please mind the existance of such differences." >> $my_log_file
	}	
fi

echo $exit_code
