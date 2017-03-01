#!/bin/bash

# Command structure:
#                                verify-input.py
#  [0] ./master.sh               .
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
# [13] $sim_recsel               .                      {rfd}+pos+mod+nre
# [14] $sim_seq                  .                      rd+rl+fl+ber+gbs
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
# if:			constant
# out:		$f0
# mod:		$4
# rd:			${14}[0]
# rl:			${14}[1]
# fl:			${14}[2]
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
# ./master.sh $project_name $workflow $data_source $ref_seq $ins_seq $read_s
# $reads_f $reads_r $gff_file $ann_file $sim-mut $sim-recsel $sim-seq
#
# example: ./master.sh project ins sim se genome.fa ins.fa n/p n/p n/p gff.gff n/p 1+ins n/p 10+30,0+0,0+0.1+50
#
#
# ./master.sh project ins sim pe genome.fa pbinprok2.fa n/p n/p n/p TAIR10_GFF3_genes_transposons-2c.gff n/p 10+ins n/p 30+100,0+0,0+1+50


# ./master.sh project ins exp pe genome.fa pbinprok2.fa n/p pe-for_reads_20170105123045.fq pe-rev_reads_20170105123045.fq TAIR10_GFF3_genes_transposons-2c.gff n/p n/p n/p n/p

#
#
# ./master.sh project snp exp se 34k_genome.fa n/p se_reads.fq n/p n/p TAIR10_GFF3_genes_transposons-2c.gff n/p n/p n/p n/p
#
#
#
#

start_time=`date +%s`
	
# Store the location of each folder in a variables				<--------------------------Nombres correctos ? para que son ? 
f0=0_input
f1=1_intermediate_files				
f2=2_logs
f3=3_workflow_output



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


# Overwrite read_s, read_f and read_r if use chose to simulate data
if [ $data_source == 'sim' ]
then
	{
		if [ $lib_type == 'se' ]
		then
			{
				read_s=$project_name/$f0/sim_data/sim_reads/se_reads.fq
			}
		else
			{
				read_f=$project_name/$f0/sim_data/sim_reads/pe-for_reads.fq
				read_r=$project_name/$f0/sim_data/sim_reads/pe-rev_reads.fq
			}
		fi
	}
fi

# Define project location, folders, and log file
#now=$(date +"%T")
my_log_file=$project_name/log.log

# Remove all projects from server: '(sudo) rm -r project*'

# Create project folder and data folders
#mkdir $project_name
#mkdir $project_name/$f0
#mkdir $project_name/$f0/fa_input
#mkdir $project_name/$f1
#mkdir $project_name/$f2
#mkdir $project_name/$f3


# Browser has already uploaded files to generic folder


# Create log file
echo $(date)": Execution of project {" $project_name "} started." > $my_log_file
echo $(date)": Project data directories created." >> $my_log_file


############################################################
# Run 'process-input.sh'

echo $(date)": STARTING INPUT PROCESSING..." >> $my_log_file

process_input=`./process-input/process-input.sh $my_log_file $project_name $workflow $data_source $lib_type $ins_seq $read_s $read_f $read_r $gff_file $ann_file`


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
		
		simulator=`./simulator/simulator.sh $my_log_file $project_name $workflow $lib_type $ins_seq $sim_mut $sim_recsel $sim_seq`

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
		workflow_result=`./workflows/workflow-snp-v5.sh $my_log_file $project_name $workflow $data_source $lib_type $ins_seq $read_s $read_f $read_r $gff_file $ann_file`

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

echo run time is $(expr `date +%s` - $start_time) s
