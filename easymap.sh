#!/bin/bash

# Command structure:
#                                verify-input.py
#  [0] ./easymap.sh               .
#  [1] $project_name             .
#  [2] $workflow[ins/snp]   .                      Maybe add a 3rd workflow: Analysis of SNPs
#  [3] $data_source[exp/sim]     .
#  [4] $lib_type[se/pe]          .
#  [5] $ref_seqs                 *
#  [6] $ins_seq                  *
#  [7] $read_s                   *
#  [8] $read_f                   *
#  [9] $read_r                   *
# [10] $gff_file                 *
# [11] $ann_file                 *
# [12] $sim_mut                  .                      nbr+mod
# [13] $sim_recsel               .                      rfd+pos+mod+nre
# [14] $sim_seq                  .                      rd+rl+fl+ber+gbs
# [15] $read_s_par               TO DO
# [16] $read_f_par               TO DO
# [17] $read_r_par               TO DO
# [18] $is_ref_strain            .                      Only for linkage analysis mapping				ref/noref
# [19] $cross_type               .                      Only for linkage analysis mapping				oc/bc
# [20] $parental_reads           .                      Only for linkage analysis mapping				mut/nomut
# [21] 
# [22]
#
# sim-mut.py
# nbr:		${12}[0]
# mod:		${12}[1]
# con:		[5]
# out:		constant
#
# sim-recsel.py
# rfd:		${13}[0]
# pos:		${13}[1]
# mod:		${13}[2]
# nre:		${13}[3]
# mut:		constant
# pol:		constant
# out:		constant
#
# sim-seq.py
# if:		constant
# out:		$f0
# mod:		$4
# rd:		${14}[0]
# rl:		${14}[1]
# fl:		${14}[2]
# ber:		${14}[3]
# gbs:		${14}[4]
#
# The command has three levels of arguments:
# 1st: Main arguments are separated by whitespace.
# 2nd: Each simulator script receives an argument composed of a string of arguments reparated
#      by the '+' character.
# 3rd: Some of the 2nd level arguments are in turn composed of a string of values separated by
#      the ',' character.
#
# ./easymap.sh $project_name $workflow $data_source $ref_seq $ins_seq $read_s
# $reads_f $reads_r $gff_file $ann_file $sim-mut $sim-recsel $sim-seq
#
# example: ./easymap.sh project ins sim pe genome.fa pbinprok2.fa n/p n/p n/p chr1.gff n/p 10+li n/p 10+100,0+500,100+1+100 n/p n/p n/p n/p n/p n/p
#
# Command to do mapping by sequencing with simulated data
# ./easymap.sh project snp sim pe genome.fa n/p n/p n/p n/p TAIR10_GFF3_genes_transposons-2c.gff n/p 5000+e "0,14;1,31;2,33;3,15;4,5;5,2/0,24;1,42;2,25;3,6;4,1;5,1"+1,10000000+r+50 25+100,0+500,100+1+50 n/p n/p n/p oc ref mut
#
# Simulated MbS with just chromosome 1
# ./easymap.sh project snp sim se genome.fa n/p n/p n/p n/p chr1.gff n/p 150+e "0,14;1,31;2,33;3,15;4,5;5,2"+1,10000000+r+50 25+200,40+0,0+1+100 n/p n/p n/p oc ref mut

############################################################
# Get command arguments and assign them to variables

project_name=$1
workflow=$2
data_source=$3
lib_type=$4
ref_seqs=$5
ins_seq=$6
read_s=$7
read_f=$8
read_r=$9
gff_file=${10}
ann_file=${11}
sim_mut=${12}
sim_recsel=${13}
sim_seq=${14}
read_s_par=${15}
read_f_par=${16}
read_r_par=${17}
cross_type=${18}
is_ref_strain=${19}
parental_used_as_control=${20}


############################################################
# Preparation steps

# Declare a flag variable that will be used as exit code, and set it to 0 (no error)
exit_code=0

# Store the location of each folder in variables
f0=0_input
f1=1_intermediate_files				
f2=2_logs
f3=3_workflow_output

# Check if some folders exist. If no, create them
[ -d $project_name/$f1 ] || mkdir $project_name/$f1
[ -d $project_name/$f2 ] || mkdir $project_name/$f2
[ -d $project_name/$f3 ] || mkdir $project_name/$f3

# Delete intermediate and final files from previous executions, For that, check whether dirs have any
# content (folders or files) and, if so, remove it
[ "$(ls -A $project_name/$f1)" ] && rm --recursive $project_name/$f1/*
[ "$(ls -A $project_name/$f2)" ] && rm --recursive $project_name/$f2/*
[ "$(ls -A $project_name/$f3)" ] && rm --recursive $project_name/$f3/*

# Define path of log file and create it
my_log_file=$project_name/$f3/log.log
touch $my_log_file

# Check that the folders /project/0_input and /project/0_input/gnm_ref exist. If they do not, 
if ! [ -d $project_name/$f0 ]
then
	{
		echo $(date)": Execution could not start because folder /project/0_input could not be found. Please, create the folder and use it to place the files to analyze." > $my_log_file
		exit_code=1
		echo $exit_code
		exit	
	}
fi
if ! [ -d $project_name/$f0/gnm_ref ]
then
	{
		echo $(date)": Execution could not start because folder /project/0_input/gnm_ref could not be found. Please, create the folder and use it to place the your reference genome." > $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
fi

############################################################
# Start easymap
echo $(date)": Execution of project {" $project_name "} started." > $my_log_file
echo $(date)": Project data directories created." >> $my_log_file


############################################################
# Overwrite read_s, read_f and read_r if use chose to simulate data
if [ $data_source == 'sim' ]
then
	{
		if [ $lib_type == 'se' ]
		then
			{
				read_s=$project_name/$f1/sim_data/sim_seq_output/sample/se_reads.fq
				read_s_par=$project_name/$f1/sim_data/sim_seq_output/control/se_reads.fq
			}
		else
			{
				read_f=$project_name/$f1/sim_data/sim_seq_output/sample/pe-for_reads.fq
				read_r=$project_name/$f1/sim_data/sim_seq_output/sample/pe-rev_reads.fq
				read_f_par=$project_name/$f1/sim_data/sim_seq_output/control/pe-for_reads.fq
				read_r_par=$project_name/$f1/sim_data/sim_seq_output/control/pe-rev_reads.fq
			}
		fi
	}
fi


############################################################
# Run 'process-input.sh'

echo $(date)": STARTING INPUT PROCESSING..." >> $my_log_file

process_input=`./process_input/process-input.sh $my_log_file $project_name $workflow $data_source $lib_type $ins_seq $read_s $read_f $read_r $gff_file $ann_file`


if [ $process_input == 0 ]
then
	{
		echo $(date)": All inputs correct." >> $my_log_file
	}
else 
	{
		echo $(date)": One or more inputs incorrect (see details above in this log). Quit." >> $my_log_file
		exit
	}
fi


############################################################
# Run 'simulator.sh'

if [ $data_source == 'sim' ]
then
	{
		echo $(date)": STARTING DATA SIMULATION..." >> $my_log_file
		
		simulator=`./simulator/simulator.sh $my_log_file $project_name $workflow $lib_type $ins_seq $sim_mut $sim_recsel $sim_seq $cross_type $is_ref_strain $parental_used_as_control`

		if [ $simulator == 0 ]
		then
			{
				echo $(date)": Simulation completed." >> $my_log_file
			}
		else 
			{
				echo $(date)": Simulation failed (see details above in this log). Quit." >> $my_log_file
				exit
			}
		fi
	}
fi


############################################################
# Run the chosen analysis workflow

if [  $workflow == 'ins' ]
then
	{		
		workflow_result=`./workflows/workflow-ins-v3.sh $my_log_file $project_name $workflow $data_source $lib_type $ins_seq $read_s $read_f $read_r $gff_file $ann_file`


		if [ $workflow_result == 0 ]
		then
			{
				echo $(date)": Analysis workflow finished correctly." >> $my_log_file
			}
		else 
			{
				echo $(date)": Analysis workflow failed (see details above in this log)." >> $my_log_file
				exit
			}
		fi
	}
fi

if [  $workflow == 'snp' ]
then
	{		
		workflow_result=`./workflows/workflow-snp-v5.sh $my_log_file $project_name $workflow $data_source $lib_type $ins_seq $read_s $read_f $read_r $gff_file $ann_file $read_s_par $read_f_par $read_r_par $cross_type $is_ref_strain $parental_used_as_control`

		if [ $workflow_result == 0 ]
		then
			{
				echo $(date)": Analysis workflow finished correctly." >> $my_log_file
			}
		else 
			{
				echo $(date)": Analysis workflow failed (see details above in this log)." >> $my_log_file
				exit
			}
		fi
	}
fi

echo $(date)": Execution of project {" $project_name "} finished." >> $my_log_file

echo $exit_code
