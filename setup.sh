#!/bin/bash

# This script automates some of the steps required to perform after cloning from GitHub
# repository in order to make easymap ready for execution.

################################################################################

# Create some folders not present in GitHub repo (e.g. 'user_data' and 'user_projects')

[ -d user_data ] || mkdir user_data
[ -d user_projects ] || mkdir user_projects

# Allow server users like www-data (Apache user) to execute programs, to work with 
# data in "user_data" and "user_projects", and to read configuration in "config/config" 
sudo chmod 777 . --recursive

################################################################################

# Compile bcftools, bowtie and samtools

# Install libraries libncurses5-dev and libncursesw5-dev if not already installed
sudo apt-get install libncurses5-dev libncursesw5-dev
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

# Install several libraries that are Pillow dependencies
sudo apt-get install libtiff5-dev libjpeg8-dev zlib1g-dev libfreetype6-dev liblcms2-dev libwebp-dev tcl8.6-dev tk8.6-dev python-tk
#sudo yum (for other distros)

# Install Pillow
sudo python ./graphic_output/Pillow-4.2.1/setup.py install

# In Ubuntu AMI (python-minimal installed), setup.py needs the package setuptools, which
# is not found and a receive an error.

################################################################################



# Make files executable and give logged user execution permissions (.sh, .py, binaries)



################################################################################

# The following commands install Apache and PHP.
# Not needed if Apache and PHP already installed in machine
# AWS EC2 AMIs already come with Apache nad PHP by default
# or can be set when creating the AMI

# Update apt
#sudo apt-get update

# Install Apache
#sudo apt-get install apache2

# Install PHP including componenet for Apache
#sudo apt-get install php libapache2-mod-php

# Check PHP installation
#php --version

# Restart apache
#sudo service apache2 restart

################################################################################
