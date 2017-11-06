#!/bin/bash

################################################################################
#
# This script automates some the steps required after cloning or downloading
# in order to make easymap ready for execution.
#
################################################################################

################################################################################
#
# REQUIREMENTS:
#	
#	- To use Easymap through the command line
#		- ...
#
#	- To use Easymap through the web interface
#		- Web server that runs PHP
#
################################################################################

# Create some folders not present in GitHub repo (e.g. 'user_data' and 'user_projects')

[ -d user_data ] || mkdir user_data
[ -d user_projects ] || mkdir user_projects

################################################################################

# Compile bcftools, bowtie and samtools

cd bcftools-1.3.1 
make clean
make

cd ../bowtie2 
make clean
make

cd ../samtools1 
make clean
make

cd ..

################################################################################

# Get Python-2.7.12
[ -d src ] || mkdir src
cd src
wget https://www.python.org/ftp/python/2.7.12/Python-2.7.12.tgz
tar -zxvf Python-2.7.12.tgz

# Install Python-2.7.12
cd Python-2.7.12
[ -d .localpython ] || mkdir .localpython
./configure --prefix=$PWD/.localpython
make
make install
cd ..

# Get virtualenv-15.1.0
wget https://pypi.python.org/packages/d4/0c/9840c08189e030873387a73b90ada981885010dd9aea134d6de30cd24cb8/virtualenv-15.1.0.tar.gz#md5=44e19f4134906fe2d75124427dc9b716
tar -zxvf virtualenv-15.1.0.tar.gz

# Install virtualenv-15.1.0
cd virtualenv-15.1.0/
../Python-2.7.12/.localpython/bin/python setup.py install

# Create virtual environment "easymap-env"
../Python-2.7.12/.localpython/bin/python virtualenv.py easymap-env -p ../Python-2.7.12/.localpython/bin/python

# Install Pillow with pip
[ -d cache ] || mkdir cache
easymap-env/bin/pip install Pillow --cache-dir cache

cd ../..

################################################################################

# Change permissions to the easymap folder and subfolders so Easymap can be used both from the
# web interface (server user -- e.g. www-data) and the command line of any user
sudo chmod 777 . --recursive


# In file easymap, set absolute path to the Python binaries of the virtual environment
# Rest of Pyhton scripts don't need this because are executed after easymap.sh activates the virtual environment
sed -i -e "s~ABS_PATH_ENV_PYTHON~${PWD}/src/Python-2.7.12/.localpython/bin/python2~g" easymap



# Finally we check if Easymap functions properly by running a small project: 
cp fonts/1.check.fa user_data/
cp fonts/check.gff user_data/

run_result=`./easymap -n setup -w snp -sim -r check -g check.gff -ed ref_bc_parmut  -sm 50 -sr 0,24-1,42-2,25-3,6-4,1-5,2+1,400000+100 -ss 30+100,0+500,100+1+50+se`

# Cleanup
rm -rf user_data/check.gff
rm -rf user_data/1.check.gff 
rm -rf user_projects/*

if [ "$run_result" == "Easymap analysis properly completed." ]; then
	echo "Easymap installation susccesfuly completed" 
else 
	echo "Easymap installation failed" 
	exit
fi