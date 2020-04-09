#!/bin/bash
##this is a build file##
#to do: SRA toolkits 
#novocraft has to be downloaded by user

#grab passwords
if [ $? -eq 0 ]; then
    echo "Online"
    read -s -p "Please enter your super user password: " sudoPW
    
	#set dirs
	master_dir="$(pwd "$0")"
	mkdir programs
	program_dir="$master_dir/programs"

	source "${master_dir}/scripts/install_functions.sh"
	#wget downloads
	echo "downloading packages"

	URL=(	"http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip"
		"https://github.com/broadinstitute/picard/releases/download/2.12.2/picard.jar"
		)

	for i in ${URL[@]};
	do
		wget -P $master_dir/programs/  $i
	done

	mv "${master_dir}/databases/GenomeAnalsysTK.jar" "${program_dir}/GenomeAnalsysTK.jar"


	#pip installations
	sudoInstaller "update"
	sudoInstaller "libboost-all-dev"
	sudoInstaller "gcc"
	sudoInstaller "make"
	sudoInstaller "libbz2-dev"
	sudoInstaller "zlib1g-dev"
	sudoInstaller "libncurses5-dev" 
	sudoInstaller "libncursesw5-dev"
	sudoInstaller "liblzma-dev"
	sudoInstaller "cython"
	sudoInstaller "cython3"
	sudoInstaller "python3-pip"
	sudoInstaller "python-pip"
	sudoInstaller "autoconf"

	pip install pysam
	pip install numpy
	pip3 install pandas
	pip3 install pysam

	#multithread processes
	sudoInstaller "samtools"
	sudoInstaller "bedtools"
	sudoInstaller "libcurl4-openssl-dev"
	sudoInstaller "libssl-dev"
	sudoInstaller "vcftools" 
	sudoInstaller "bwa" 
	sudoInstaller "bcftools" 
	sudoInstaller "abacas" 
	sudoInstaller "fastqc" 
	sudoInstaller "tabix" 
	sudoInstaller "sambamba" 
	sudoInstaller "yad" 
	sudoInstaller "g++" 
	sudoInstaller "cmake" 
	sudoInstaller "gawk"
	sudoInstaller "gnuplot-nox"
	sudoInstaller "soapdenovo2"

	samblaster
	export PATH=$PATH:"${master_dir}/programs/samblaster"
	svtyper
	lumpy
	delly 
	svprops
	zips


	echo "Build finished"
	exit 0

else
    echo "internet connection required"
fi
