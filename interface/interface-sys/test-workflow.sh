#!/bin/bash


arg1=$1
arg2=$2
arg3=$3

mkdir files

echo 'First message' > log.log
echo $arg1 $arg2 $arg3 >> log.log
sleep 3s

#mkdir files/test_dir1
echo $arg1 $arg2 $arg3 >> log.log
sleep 3s

#mkdir files/test_dir2
echo $arg1 $arg2 $arg3 >> log.log
sleep 3s

echo $arg1 $arg2 $arg3 >> log.log
sleep 3s

echo $arg1 $arg2 $arg3 >> log.log
sleep 3s

echo $arg1 $arg2 $arg3 >> log.log
sleep 3s

#rmdir files/test_dir1
#rmdir files/test_dir2
#rmdir files
