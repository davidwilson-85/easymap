#!/bin/bash

################################################################################
#
# This script automates some the steps required after cloning or downloading
# in order to make easymap ready for execution.
#
################################################################################


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

################################################################################

# Installing some general tools 

yum groupinstall "Development Tools"
yum groupinstall "Development Libraries"

################################################################################

# Create some folders not present in GitHub repo (e.g. 'user_data' and 'user_projects')

[ -d user_data ] || mkdir user_data
[ -d user_projects ] || mkdir user_projects

################################################################################


# Compile bcftools, bowtie and samtools

# Install libraries libncurses5-dev and libncursesw5-dev if not already installed
# sudo apt-get install libncurses5-dev libncursesw5-dev
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

sudo pip install Pillow

################################################################################

# Grant full permissions to logged in user and its group
sudo chmod 770 . --recursive


################################################################################
#
# The following commands make possible to use Easymap from both the terminal of the logged in user
# and by a server (all examples with Apache2 - www-data user)
#
################################################################################

if [ $1 == server ]
then

	# Change ownership of easymap folder (recursively) to server user (www-data)
	sudo chown -R www-data:www-data .

	# Add logged in user to server group (www-data)
	user="$(id -u -n)"
	sudo usermod -a -G www-data $user


fi