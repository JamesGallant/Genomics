#!/bin/bash
##this is a build file##
master_dir="$(pwd "$0")"
#mkdir programs

#URLs for wget

URL=(	"http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip"
	)

for i in ${URL[@]};
do
	echo "downloading packages"
	wget -P $master_dir/programs/  $i
done
