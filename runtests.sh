#!/bin/bash

if [ $# != 1 ]
then
	echo "Provide a directory with tests"
	exit
fi


for dir in $1/*
do
    echo $dir
	./bin/heat $dir/settings $dir/u0.bin
done
