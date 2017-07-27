#!/usr/bin/env python
#echo -e "umh2" | sudo -S python changing.py
import subprocess
from tempfile import mkstemp
from shutil import move
import os
from os import remove, close


def replace(file_path, pattern, subst):
    #Create temp file
    fh, abs_path = mkstemp()
    with open(abs_path,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                new_file.write(line.replace(pattern, subst))
    close(fh)
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)



list_x= [250000,200000,150000]
list_y = [25000,20000,15000]


n = 0

for x in list_x:
    z = 0
    for y in list_y:
        if z != 0 and n != 0:
            replace("./workflows/workflow-snp.sh","python $location/scripts_snp/analysis/map-mutation.py -fichero $f1/F2_control_comparison.va -fasta $f1/$my_gs -mode $my_analysis_mode -window_size"+ str(list_x[n-1]) +" -window_space"+ str(list_x[z-1])+ " -output $f1/map_info.txt -control_modality $my_mutbackgroud -interval_width 4000000 -snp_analysis_type $snp_analysis_type","python $location/scripts_snp/analysis/map-mutation.py -fichero $f1/F2_control_comparison.va -fasta $f1/$my_gs -mode $my_analysis_mode -window_size " +str(x)+ " -window_space "+str(y)+" -output $f1/map_info.txt -control_modality $my_mutbackgroud -interval_width 4000000  -snp_analysis_type $snp_analysis_type")
        elif z == 0 and n != 0:
            replace("./workflows/workflow-snp.sh","python $location/scripts_snp/analysis/map-mutation.py -fichero $f1/F2_control_comparison.va -fasta $f1/$my_gs -mode $my_analysis_mode -window_size"+ str(list_x[n]) +" -window_space"+ str(list_x[z-1])+ " -output $f1/map_info.txt -control_modality $my_mutbackgroud -interval_width 4000000 -snp_analysis_type $snp_analysis_type","python $location/scripts_snp/analysis/map-mutation.py -fichero $f1/F2_control_comparison.va -fasta $f1/$my_gs -mode $my_analysis_mode -window_size " +str(x)+ " -window_space "+str(y)+" -output $f1/map_info.txt -control_modality $my_mutbackgroud -interval_width 4000000  -snp_analysis_type $snp_analysis_type")            
        elif z != 0 and n == 0:
            replace("./workflows/workflow-snp.sh","python $location/scripts_snp/analysis/map-mutation.py -fichero $f1/F2_control_comparison.va -fasta $f1/$my_gs -mode $my_analysis_mode -window_size"+ str(list_x[n-1]) +" -window_space"+ str(list_x[z])+ " -output $f1/map_info.txt -control_modality $my_mutbackgroud -interval_width 4000000 -snp_analysis_type $snp_analysis_type","python $location/scripts_snp/analysis/map-mutation.py -fichero $f1/F2_control_comparison.va -fasta $f1/$my_gs -mode $my_analysis_mode -window_size " +str(x)+ " -window_space "+str(y)+" -output $f1/map_info.txt -control_modality $my_mutbackgroud -interval_width 4000000  -snp_analysis_type $snp_analysis_type")
        elif z == 0 and n == 0:
            replace("./workflows/workflow-snp.sh","python $location/scripts_snp/analysis/map-mutation.py -fichero $f1/F2_control_comparison.va -fasta $f1/$my_gs -mode $my_analysis_mode -window_size"+ str(list_x[n]) +" -window_space"+ str(list_x[z])+ " -output $f1/map_info.txt -control_modality $my_mutbackgroud -interval_width 4000000 -snp_analysis_type $snp_analysis_type","python $location/scripts_snp/analysis/map-mutation.py -fichero $f1/F2_control_comparison.va -fasta $f1/$my_gs -mode $my_analysis_mode -window_size " +str(x)+ " -window_space "+str(y)+" -output $f1/map_info.txt -control_modality $my_mutbackgroud -interval_width 4000000  -snp_analysis_type $snp_analysis_type")

        #subprocess.call(['./easymap', "-n caso_7_z_n -w snp -d sim -r at4 -g chr1+4.gff -rs noref -cr oc -co par_mut -sm 80+e -sr 0,16-1-34-2,31-3,14-4,4-5,1+1,7000000+r+100 -ss 25+100,0+500,100+1+50+se")
        #mut = rand_input
        command = subprocess.call(["sudo","chmod", "777", "workflows/workflow-snp.sh"]) 

        os.system("./easymap -n caso_7_z_n -w snp -d sim -r at4 -g chr1+4.gff -rs noref -cr oc -co par_mut -sm 80+e -sr 0,16-1,34-2,31-3,14-4,4-5,1+1,7000000+r+100 -ss 25+100,0+500,100+1+50+se")


        subprocess.call(["rm", "-rf" ,"./user_projects/*/1_intermediate_files"])

        z +=1 
    n+= 1

replace("./workflows/workflow-snp.sh","python $location/scripts_snp/analysis/map-mutation.py -fichero $f1/F2_control_comparison.va -fasta $f1/$my_gs -mode $my_analysis_mode -window_size"+ str(list_x[-1]) +" -window_space"+ str(list_x[-1])+ " -output $f1/map_info.txt -control_modality $my_mutbackgroud -interval_width 4000000 -snp_analysis_type $snp_analysis_type","python $location/scripts_snp/analysis/map-mutation.py -fichero $f1/F2_control_comparison.va -fasta $f1/$my_gs -mode $my_analysis_mode -window_size " +str(list_x[0])+ " -window_space "+str(list_y[0])+" -output $f1/map_info.txt -control_modality $my_mutbackgroud -interval_width 4000000  -snp_analysis_type $snp_analysis_type")
