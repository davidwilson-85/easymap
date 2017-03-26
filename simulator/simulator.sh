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
f0=user_data
f1=1_intermediate_files

# Get command arguments and assign them to variables
my_log_file=$1
project_name=$2
analysis_type=$3
lib_type=$4
ins_seq=$f0/$5

# Get the string that contains the parameters for sim-mut.py and extract them by splitting the string by the '+' character
sim_mut_statement=$6
IFS='+' read -ra sim_mut_array <<< "$sim_mut_statement"
nbr_muts=${sim_mut_array[0]}
mut_mode=${sim_mut_array[1]}

# Get the string that contains the parameters for sim-recsel.py and extract them by splitting the string by the '+' character
sim_recsel_statement=$7
IFS='+' read -ra sim_recsel_array <<< "$sim_recsel_statement"
rec_freq_distr=${sim_recsel_array[0]} # Recombination frequency distribution. Pass it ty program as a string and analyze it with python
mut_pos=${sim_recsel_array[1]} #This parameter will be used in sim_mut as well, due to the fact the mutations will be generated previously.
sel_mode=${sim_recsel_array[2]}
nbr_rec_chrs=${sim_recsel_array[3]}
mutant_parental=$project_name/$f1/sim_data/sim_mut_output/mutated_genome/mutated_genome.fa
polymorphic_parental=$project_name/$f1/gnm_ref_merged/genome.fa


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
#echo mut_pos: $mut_pos >> $my_log_file This value will be used both for sim-mut in order to create the appropiate mutation and in sim-rec because the recombination should not happen in this position
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
ref_seqs_merged_file=$project_name/$f1/gnm_ref_merged/genome.fa

# Establish location of some folders required for the simulation
sim_mut_output_folder_polymorphicstrain=$project_name/$f1/sim_data/sim_mut_output/polymorphic_strain
sim_mut_output_folder_mutantstrain=$project_name/$f1/sim_data/sim_mut_output/mutant_strain
# mutated sequence is in: $project_name/$f1/sim_data/sim_mut_output/mutated_genome/mutated_genome.fa
# info is in: $project_name/$f1/sim_data/sim_mut_output/info/info_all_mutations.txt
sim_recsel_output_folder=$project_name/$f1/sim_data/sim_recsel_output
sim_seq_output_folder_sample=$project_name/$f1/sim_data/sim_seq_output/sample
sim_seq_output_folder_control=$project_name/$f1/sim_data/sim_seq_output/control

# Determine if user wants backcross or outcross, the background of the mutant (ref/noref), and the parental to sequence and use as control in 'snp' mode
cross_type=${9}
is_ref_strain=${10}
parental_genome_to_sequence=${11}

# Define folders with input files for outcross simulations
if [ $cross_type == 'oc' ]
then
	{

		if [ $is_ref_strain == 'ref' ]
		then
			{
				seq_to_mutate=$project_name/$f1/gnm_ref_merged/genome.fa
				outcross_polymorphic_parental=$sim_mut_output_folder_polymorphicstrain/mutated_genome/mutated_genome.fa

				if [ $parental_genome_to_sequence == 'mut' ]
				then
					{
						parental_genome_location=$project_name/$f1/gnm_ref_merged
					}
				else
					{
						parental_genome_location=$sim_mut_output_folder_polymorphicstrain/mutated_genome
					}
				fi
			}
		else
			{
				seq_to_mutate=$sim_mut_output_folder_polymorphicstrain/mutated_genome/mutated_genome.fa
				outcross_polymorphic_parental=$ref_seqs_merged_file

				if [ $parental_genome_to_sequence == 'mut' ]
				then
					{
						parental_genome_location=$sim_mut_output_folder_polymorphicstrain/mutated_genome
					}
				else
					{
						parental_genome_location=$project_name/$f1/gnm_ref_merged
					}
				fi
			}
		fi

	}
fi






if [ $analysis_type == 'ins' ]
then
	{
		# Run sim-mut.py
		{
			python simulator/sim-mut.py -nbr $nbr_muts -mod $mut_mode -con $ref_seqs_merged_file -ins $ins_seq -out $sim_mut_output_folder_mutantstrain
	
		} || {
			echo $(date)": Simulation of mutagenesis failed. Quit." >> $my_log_file
			exit_code=1
			echo exit_code
			exit
		}
		echo $(date)": Simulation of mutagenesis completed." >> $my_log_file
		
		# Run sim-seq.py. The input is a folder becasuse the program works with all the fasta files that finds in a folder. This is necessary to simulate the sequencing of bulked DNA.
		{
			python simulator/sim-seq.py -if $sim_mut_output_folder_mutantstrain/mutated_genome -out $sim_seq_output_folder_sample -mod $lib_type -rd $read_depth -rlm $read_length_mean -rls $read_length_sd -flm $fragment_length_mean -fls $fragment_length_sd -ber $basecalling_error_rate -gbs $gc_bias_strength
	
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
		# Two different types of workflows depending on the type of cross 
		if [ $cross_type == 'bc' ]
		then
			{
				# Run sim-mut.py to create mutant strain
				{
					python simulator/sim-mut.py -nbr $nbr_muts -mod $mut_mode -con $ref_seqs_merged_file -out $sim_mut_output_folder_mutantstrain -causal_mut $mut_pos

				} || {
					echo $(date)": Simulation of mutagenesis failed. Quit." >> $my_log_file
					exit_code=1
					echo exit_code
					exit
				}
				echo $(date)": Simulation of mutagenesis completed." >> $my_log_file

				# Run sim-recsel.py to create recombinant chromosomes 
				{
					python simulator/sim-recsel.py -outdir $sim_recsel_output_folder -rec_freq_distr $rec_freq_distr  -parmut $sim_mut_output_folder_mutantstrain/mutated_genome/mutated_genome.fa -parpol $ref_seqs_merged_file -mutapos $mut_pos -smod $sel_mode -nrec $nbr_rec_chrs 
				} || {
					echo $(date)": Simulation of recombination and phenotype selection failed. Quit." >> $my_log_file
					exit_code=1
					echo exit_code
					exit
				}
				echo $(date)": Simulation of recombination and phenotype selection completed." >> $my_log_file

				# Run sim-seq.py on parental genome. The input is a folder becasuse the program works with all the fasta files that finds in a folder. This is necessary to simulate the sequencing of bulked DNA.
				{
					python simulator/sim-seq.py -if $project_name/$f1/gnm_ref_merged -out $sim_seq_output_folder_control -mod $lib_type -rd $read_depth -rlm $read_length_mean -rls $read_length_sd -flm $fragment_length_mean -fls $fragment_length_sd -ber $basecalling_error_rate -gbs $gc_bias_strength
	
				} || {
					echo $(date)": Simulation of high-throughput sequencing failed. Quit." >> $my_log_file
					exit_code=1
					echo exit_code
					exit 
				} 
				echo $(date)": Simulation of high-throughput sequencing reads on parental genome completed." >> $my_log_file

				# Run sim-seq.py on F2 recombinant population. The input is a folder becasuse the program works with all the fasta files that finds in a folder. This is necessary to simulate the sequencing of bulked DNA.
				{
					python simulator/sim-seq.py -if $sim_recsel_output_folder -out $sim_seq_output_folder_sample -mod $lib_type -rd $read_depth -rlm $read_length_mean -rls $read_length_sd -flm $fragment_length_mean -fls $fragment_length_sd -ber $basecalling_error_rate -gbs $gc_bias_strength
	
				} || {
					echo $(date)": Simulation of high-throughput sequencing failed. Quit." >> $my_log_file
					exit_code=1
					echo exit_code
					exit
				}
				echo $(date)": Simulation of high-throughput sequencing reads on F2 recombinant population completed." >> $my_log_file
			}	
		else
			{
				# Run sim-mut.py to create polymorphic strain
				{
					python simulator/sim-mut.py -nbr 100000 -mod d -con $ref_seqs_merged_file -out $sim_mut_output_folder_polymorphicstrain

				} || {
					echo $(date)": Simulation of mutagenesis to create polymorphic strain failed. Quit." >> $my_log_file
					exit_code=1
					echo exit_code
					exit
				}
				echo $(date)": Simulation of mutagenesis to create polymorphic strain completed." >> $my_log_file

				# Run sim-mut.py to create mutant strain
				{
					python simulator/sim-mut.py -nbr $nbr_muts -mod $mut_mode -con $seq_to_mutate -out $sim_mut_output_folder_mutantstrain -causal_mut $mut_pos

				} || {
					echo $(date)": Simulation of mutagenesis to create mutant strain failed. Quit." >> $my_log_file
					exit_code=1
					echo exit_code
					exit
				}
				echo $(date)": Simulation of mutagenesis to create mutant strain completed." >> $my_log_file		
	
				# Run sim-recsel.py
				{
					python simulator/sim-recsel.py -outdir $sim_recsel_output_folder -rec_freq_distr $rec_freq_distr -parmut $sim_mut_output_folder_mutantstrain/mutated_genome/mutated_genome.fa -parpol $outcross_polymorphic_parental -mutapos $mut_pos -smod $sel_mode -nrec $nbr_rec_chrs 
				} || {
					echo $(date)": Simulation of recombination and phenotype selection failed. Quit." >> $my_log_file
					exit_code=1
					echo exit_code
					exit
				}
				echo $(date)": Simulation of recombination and phenotype selection completed." >> $my_log_file

				# Run sim-seq.py on parental genome. The input is a folder becasuse the program works with all the fasta files that finds in a folder. Thos is necessary to simulate the sequencing of bulked DNA.
				{
					python simulator/sim-seq.py -if $parental_genome_location -out $sim_seq_output_folder_control -mod $lib_type -rd $read_depth -rlm $read_length_mean -rls $read_length_sd -flm $fragment_length_mean -fls $fragment_length_sd -ber $basecalling_error_rate -gbs $gc_bias_strength
	
				} || {
					echo $(date)": Simulation of high-throughput sequencing reads on parental genome failed. Quit." >> $my_log_file
					exit_code=1
					echo exit_code
					exit
				}
				echo $(date)": Simulation of high-throughput sequencing reads on parental genome completed." >> $my_log_file

				# Run sim-seq.py on F2 recombinant population. The input is a folder becasuse the program works with all the fasta files that finds in a folder. This is necessary to simulate the sequencing of bulked DNA.
				{
					python simulator/sim-seq.py -if $sim_recsel_output_folder -out $sim_seq_output_folder_sample -mod $lib_type -rd $read_depth -rlm $read_length_mean -rls $read_length_sd -flm $fragment_length_mean -fls $fragment_length_sd -ber $basecalling_error_rate -gbs $gc_bias_strength
	
				} || {
					echo $(date)": Simulation of high-throughput sequencing reads on F2 recombinant population failed. Quit." >> $my_log_file
					exit_code=1
					echo exit_code
					exit
				}
				echo $(date)": Simulation of high-throughput sequencing reads on F2 recombinant population completed." >> $my_log_file

			}
		fi
	}
fi

echo $exit_code
