#!/bin/bash

# When user has selected to simualte data, this bash script is run, which in turn runs
# python programs: sim-mut.py, sim-recsel.py, sim-seq.py.
# If the user chooses to simulate data, hse is not asked to provide reads files.
# If the user chooses workflow insertion analysis, only sim-mut.py and sim-seq.py are
# required because this analysis does not involve a mapping population.
# If the user chooses snp mapping through a recombinant mapping population, in this case
# sim-recsel.py is also executed to create the mapping population from where the reads are
# generated.


#######################################################################################################
# This block of code does some steps previous to the analysis
# 
# Description:
# 
# 
# 
#######################################################################################################


# This is the command sent by 'master.sh':
# ./simulator.sh $my_log_file $project_name $analysis_type $lib_type $ins_seq $sim-mut{nbr+mod} $sim-recsel{} $sim-seq{rd+rl+fl+ber+gbs}
# example: ./simulator.sh project/log.log project ins se ins.fa 1+ins {} 10+30,0+0,0+1+50

# Set 'exit_code' (flag variable) to it's initial value (0)
exit_code=0

# Store the location of each folder in a variable
f0=0_input

# Get command arguments and assign them to variables
my_log_file=$1
project_name=$2
analysis_type=$3
lib_type=$4
ins_seq=$project_name/$f0/$5

# Get the string taht contains the parameters for sim-mut.py and extract them by splitting the string by the '+' character
sim_mut_statement=$6
IFS='+' read -ra sim_mut_array <<< "$sim_mut_statement"
nbr_muts=${sim_mut_array[0]}
mut_mode=${sim_mut_array[1]}

# Get the string taht contains the parameters for sim-recsel.py and extract them by splitting the string by the '+' character
sim_recsel_statement=$7
IFS='+' read -ra sim_recsel_array <<< "$sim_recsel_statement"
rec_freq_distr=${sim_recsel_array[0]} # Recombination frequency distribution. Pass it ty program as a string and analyze it with python
mut_pos=${sim_recsel_array[1]}
sel_mode=${sim_recsel_array[2]}
nbr_rec_chrs=${sim_recsel_array[3]}
mutant_parental=$project_name/$f0/sim_data/sim_mut_output/mutated_genome/mutated_genome.fa
polymorphic_parental=$project_name/$f0/gnm_ref_merged/genome.fa


# Get the string taht contains the parameters for sim-seq.py and extract them by splitting the string by the '+' character
sim_seq_statement=$8
IFS='+' read -ra sim_seq_array <<< "$sim_seq_statement"
read_depth=${sim_seq_array[0]}
basecalling_error_rate=${sim_seq_array[3]}
gc_bias_strength=${sim_seq_array[4]}

# Split '$rl' (= read length) by the ',' character because it contains a pair of values (mean, sd)
rl=${sim_seq_array[1]}
IFS=',' read -ra rl_array <<< "$rl"
read_length_mean=${rl_array[0]}
read_length_sd=${rl_array[1]}

# Split '$fl' (= fragment length) by the ',' character because it contains a pair of values (mean, sd)
fl=${sim_seq_array[2]}
IFS=',' read -ra fl_array <<< "$fl"
fragment_length_mean=${fl_array[0]}
fragment_length_sd=${fl_array[1]}

# For development:
#echo nbr_muts: $nbr_muts >> $my_log_file
#echo mut_mode: $mut_mode >> $my_log_file
#echo rec_freq_distr: $rec_freq_distr >> $my_log_file
#echo mut_pos: $mut_pos >> $my_log_file
#echo sel_mode: $sel_mode >> $my_log_file
#echo nbr_rec_chrs: $nbr_rec_chrs >> $my_log_file
#echo mutant_parental: $mutant_parental >> $my_log_file
#echo polymorphic_parental: $polymorphic_parental >> $my_log_file
#echo read_depth: $read_depth >> $my_log_file
#echo read_length_mean: $read_length_mean >> $my_log_file
#echo read_length_sd: $read_length_sd >> $my_log_file
#echo fragment_length_mean: $fragment_length_mean >> $my_log_file
#echo fragment_length_sd: $fragment_length_sd >> $my_log_file
#echo basecalling_error_rate: $basecalling_error_rate >> $my_log_file
#echo gc_bias_strength: $gc_bias_strength >> $my_log_file

# Establish the location of the reference sequence
ref_seqs_merged_file=$project_name/$f0/gnm_ref_merged/genome.fa

# Establish location of some folders required for the simulation
sim_mut_output_folder=$project_name/$f0/sim_data/sim_mut_output
# mutated sequence is in: $project_name/$f0/sim_data/sim_mut_output/mutated_genome/mutated_genome.fa
# info is in: $project_name/$f0/sim_data/sim_mut_output/info/info_all_mutations.txt
sim_recsel_output_folder=$project_name/$f0/sim_data/sim_recsel_output
sim_seq_output_folder=$project_name/$f0/sim_data/sim_reads


if [ $analysis_type == 'ins' ]
then
	{
		# Run sim-mut.py
		{
			python simulator/sim-mut.py -nbr $nbr_muts -mod li -con $ref_seqs_merged_file -ins $ins_seq -out $sim_mut_output_folder
	
		} || {
			echo $(date)": Simulation of mutagenesis failed. Quit." >> $my_log_file
			exit_code=1
			echo exit_code
			exit
		}
		echo $(date)": Simulation of mutagenesis completed." >> $my_log_file
		
		# Run sim-seq.py. The input is a folder becasuse the program works with all the fasta files that finds in a folder. Thos is necessary to simulate the sequencing of bulked DNA.
		{
			python simulator/sim-seq.py -if $sim_mut_output_folder/mutated_genome -out $sim_seq_output_folder -mod $lib_type -rd $read_depth -rlm $read_length_mean -rls $read_length_sd -flm $fragment_length_mean -fls $fragment_length_sd -ber $basecalling_error_rate -gbs $gc_bias_strength
	
		} || {
			echo $(date)": Simulation of high-throughput sequencing failed. Quit." >> $my_log_file
			exit_code=1
			echo exit_code
			exit
		}
		echo $(date)": Simulation of high-throughput sequencing completed." >> $my_log_file
	}
fi


if [ $analysis_type == 'snp' ]
then
	{	
		# Run sim-mut.py
		{
			python simulator/sim-mut.py -nbr $nbr_muts -mod $mut_mode -con $ref_seqs_merged_file -out $sim_mut_output_folder

		} || {
			echo $(date)": Simulation of mutagenesis failed. Quit." >> $my_log_f
			exit_code=1
			echo exit_code
			exit
		}
		echo $(date)": Simulation of mutagenesis completed." >> $my_log_file

		# Run sim-recsel.py
		{
			python simulator/sim-recsel.py -outdir $sim_recsel_output_folder -recombination_frequency $rec_freq_distr -parmut $mutant_parental -parpol $polymorphic_parental -mutapos $mut_pos -smod $sel_mode -nrec $nbr_rec_chrs 
		} || {
			echo $(date)": Simulation of recombination and phenotype selection failed. Quit." >> $my_log_file
			exit_code=1
			echo exit_code
			exit
		}
		echo $(date)": Simulation of recombination and phenotype selection completed." >> $my_log_file

		# Run sim-seq.py. The input is a folder becasuse the program works with all the fasta files that finds in a folder. Thos is necessary to simulate the sequencing of bulked DNA.
		{
			python simulator/sim-seq.py -if $sim_recsel_output_folder -out $sim_seq_output_folder -mod $lib_type -rd $read_depth -rlm $read_length_mean -rls $read_length_sd -flm $fragment_length_mean -fls $fragment_length_sd -ber $basecalling_error_rate -gbs $gc_bias_strength
	
		} || {
			echo $(date)": Simulation of high-throughput sequencing failed. Quit." >> $my_log_file
			exit_code=1
			echo exit_code
			exit
		}
		echo $(date)": Simulation of high-throughput sequencing completed." >> $my_log_file
	}
fi

echo $exit_code
