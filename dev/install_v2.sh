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
#		- Python2
#		- pip >= 9.0.1
#
#	- To use Easymap through the web interface
#		- Python2
#		- pip >= 9.0.1
#		- Web server
#		- PHP
#
################################################################################


# install.sh requires one argument, wich has two options [ cli | server ]
# server: does all the regular installation steps plus others required to 
# deal with server particularities

# Deal with argument provided by user
if ! [ $1 ]; then
	echo 'Please provide an argument specifying the type of installation: "cli" or "server". Example: "./install.sh server"'
	exit
fi

if ! [ $1 == server ] && ! [ $1 == cli ]
then
	echo 'Please choose between "cli" and "server". Example: "./install.sh server"'
	exit
fi


# Create some folders not present in GitHub repo (e.g. 'user_data' and 'user_projects')

[ -d user_data ] || mkdir user_data
[ -d user_projects ] || mkdir user_projects

################################################################################

# Compile bcftools, bowtie and samtools

# Install libraries libncurses5-dev and libncursesw5-dev if not already installed
sudo apt-get install libncurses5-dev libncursesw5-dev zlib1g-dev
# zlib1g-dev htslib/bgzf.h:34:18: fatal error: zlib.h: No such file or directory
#sudo yum (for other distros)

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

# Install Pillow (PIL fork easier to maintain and actively developed)
# http://pillow.readthedocs.io/en/3.1.x/installation.html

#sudo apt-get install python-pip    # this can be outdated and cause many problems downstream

# Get and install pip (https://pip.pypa.io/en/stable/installing/)
sudo wget https://bootstrap.pypa.io/get-pip.py
sudo python get-pip.py

# Intall Pillow
sudo pip install Pillow

# If necessary to uninstall Pillow, use this command
#sudo pip uninstall Pillow

################################################################################

# Grant full permissions to logged in user and its group
sudo chmod 770 . --recursive




################################################################################
################################################################################
################################################################################
################################################################################

# The following commands install Apache and PHP.
# Not needed if Apache and PHP already installed in machine
# Some AWS EC2 AMIs already come with Apache and PHP by default
# or can be set when creating the AMI

# Update apt
#sudo apt-get update

# Install Apache
#sudo apt-get install apache2

# Install PHP including component for Apache
#sudo apt-get install php libapache2-mod-php

# Check PHP installation
#php --version

# Restart apache
#sudo service apache2 restart

################################################################################
#
# The following commands make possible to use Easymap from both the terminal of the logged in user
# and by a server (all examples with Apache2 - www-data user)
#
################################################################################

if [ $1 == server ]
then

	# Create tmp directory to host Python unzipped eggs
	#[ -d tmp ] || mkdir tmp
	# Pillow's will be unziped in ./tmp at runtime if Easymap is being run from web interface (www-data user)
	# This directory must be owned by user and have 700 permissions
	#sudo chmod -R 700 tmp

	# Change ownership of easymap folder (recursively) to server user (www-data)
	sudo chown -R www-data:www-data .

	# Add logged in user to server group (www-data)
	user="$(id -u -n)"
	sudo usermod -a -G www-data $user

	# Make .python-eggs directory executable by any user (includes www-data or any other server)
	# [Explanation: When Easymap is executed from the web interface, Pillow's .egg files are extracted in 
	# folder .python-eggs inside the www-data root directory (/var/www). For this, you need to give 
	# www-data user write permission to that folder]
	#[ -d /var/www/.python-eggs ] || mkdir /var/www/.python-eggs
	#sudo chown -R www-data:www-data /var/www/.python-eggs
	#sudo chmod -R 740 /var/www/.python-eggs

fi