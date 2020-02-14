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

source "${master_dir}/install_functions.sh"
#wget downloads
echo "downloading packages"

URL=(	"http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip"
	"https://github.com/broadinstitute/picard/releases/download/2.12.2/picard.jar"
	)

for i in ${URL[@]};
do
	wget -P $master_dir/programs/  $i
done

mv "$master_dir/GenomeAnalysisTK.jar" "${program_dir}/GenomeAnalsysTK.jar"


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
sudoInstaller "python3-pip"
sudoInstaller "python-pip"

pip install svtyper
pip install pysam
pip install numpy

#multithread processes
sudoInstaller "samtools"
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
lumpy
delly 
svprops
zips

echo "changing write permissions"
declare -a arr2=("lumpy-sv" "picard.jar" "trimmomatic-0.36/trimmomatic-0.36.jar" "delly"  "svprops" "samblaster")

#change permissions
for a in "${arr2[@]}";
do
	echo "changing permissions"
	chmod 777 "$program_dir/${a}"
done

echo "Build finished"

startup=$(yad --item-separator="," --separator="\t" \
	--title="Pegasus v1.0" \
	--form \
	--text="Thank you for downloading Pegasus" \
	--field="Launch GUI(Nightly)":CHK \
	--field="Launch Terminal":CHK)

echo "${startup}" > temp.txt

GUI=$(cat temp.txt | awk '{print $1}')
echo "${GUI}"

term=$(cat temp.txt | awk '{print $2}')
echo "${term}"

chmod 755 Pegasus.sh
chmod 755 PegasusGUI.sh

if [[ "${GUI}" == TRUE ]];
then
	echo "Launching gui"
	bash PegasusGUI.sh
fi

if [[ "${term}" == TRUE ]];
then
	echo "launching terminal"
	bash Pegasus.sh
fi

rm temp.txt

exit 0

else
    echo "internet connection required"
fi
