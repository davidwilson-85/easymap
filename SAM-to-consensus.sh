#!/bin/bash

#Generation of a variable with the path were the SAM files of each insertion will be held
directory=./PRIMERS/
#Loop through all files in the directory and DO the SAM to BAM conversion and the genereation of consensus sequence from each insertion consensus file. 
for i in $directory* 
do
    if test -f "$i" 
    then
    #SAM to BAM
    substring=${i%.*}
    #Check whether the number of lines that are not starting with @ to be > 0, if it is, do the rest: we might have a program to do this
       {
			./samtools1/samtools sort $i  > $substring.bam 

		}

		{
			./samtools1/samtools mpileup -uf pbinprok2.fa $substring.bam | ./bcftools-1.3.1/bcftools call -c | ./bcftools-1.3.1/vcfutils.pl vcf2fq > cns.fq
		} || {
			echo $(date) ': Error ' 

		}
		echo $(date) ': finished.'
		

		{
			#sed -i "s/pbinprok2/$substring/g" ./cns.fq
			tail -n +2 cns.fq > cns.fq.temp && mv cns.fq.temp cns.fq
			echo @"$substring" | cat - cns.fq > temp && mv temp cns.fq
		} || { 
			echo Error
		}	



    fi
   #Concatenate all the fastaq files into one big fq file, which will be given as an input for the primer generation script
	cat cns.fq >> all_insertions_cns.fq
done
sed -i "s/n//g" all_insertions_cns.fq
./primer-generation.py ./variants.txt ./genome.fa ./all_insertions_cns.fq 
