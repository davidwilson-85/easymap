#!/bin/bash

# This script automates some of the steps required to perform after cloning from GitHub
# repository in order to make easymap ready for execution.
#


# Command-line version

# Create folders not present in GitHub repo (e.g. 'user_data' and 'user_projects')
mkdir user_data
mkdir user_projects

# Pillow is PIL fork easier to maintain and actively developed
# http://pillow.readthedocs.io/en/3.1.x/installation.html

# Install several libraries that are Pillow dependencies
sudo apt-get install libtiff5-dev libjpeg8-dev zlib1g-dev libfreetype6-dev liblcms2-dev libwebp-dev tcl8.6-dev tk8.6-dev python-tk

# Install pip because with pip I can easily install Pillow afterwards
# pip gets packages from PyPI (Python Package Index), which is official Python library repository
sudo apt install python-pip

# Install Pillow
pip install Pillow

# Uninstall Pillow using pip
#pip uninstall Pillow

# Compile Bowtie2, Samtools, Bcftools
#############

# Make files executable and give logged user execution permissions (.sh, .py, binaries)

# 
