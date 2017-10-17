#!/bin/bash


# Deal with argument provided by user
if ! [ $1 ]; then
	echo 'Please provide an argument specifying Linux repo. Example: "./dependencies.sh apt"'
	exit
fi

if ! [ $1 == apt ] && ! [ $1 == yum ]
then
	echo 'Please choose between "apt", "yum". Example: "./dependencies.sh apt"'
	exit
fi


cd ..

if [ $1 == apt ]
then

	# Update repositories list
	sudo apt-get update

	# Install libraries required by Linux, Python or Easymap
	sudo apt-get install build-essential libncurses5-dev libncursesw5-dev zlib1g-dev libssl-dev

	#Install several programs needed downstream
	sudo apt-get install wget, tar, zip, git


	#If Apache and PHP are not installed, install them
	sudo apt-get install apache2
	sudo apt-get install php libapache2-mod-php
	
	#Check PHP installation
	php --version
	
	#Restart apache
	sudo service apache2 restart
	
fi

if [ $1 == yum ]
then

	# Update repositories list

	# Install libraries required by Linux, Python or Easymap
	yum groupinstall "Development Tools"
	yum groupinstall "Development Libraries"
	sudo yum install ncurses-devel

	#Install several programs needed downstream
	sudo yum install git
	sudo yum install wget

	# Install LAMP
	sudo yum install -y httpd24 php70 mysql56-server php70-mysqlnd

fi