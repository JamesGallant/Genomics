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

#global package installer
sudoInstaller() {
	if ! [ -x "$(command -v $1)" ];
	then
	    "installing $1"
	    echo $sudoPW | sudo -S apt -y install $1
	else
	    echo "$1 found"
	fi
}

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
#git installation functions
#samblaster works
samblaster() {
	echo "Installing samblaster"
	mkdir $program_dir/samblaster
	git clone "https://github.com/GregoryFaust/samblaster.git" $program_dir/samblaster
	cd $program_dir/samblaster
	make
	cd $master_dir
	
	echo "Sanity checks"
	export PATH=$PATH:"${master_dir}/programs/samblaster"
	if ! [ -x "$(command -v samblaster)" ];
	then 
	    echo "samblaster installation failed"
	    echo "Follow these docs: https://github.com/GregoryFaust/samblaster.git"
	    echo "The files need to be in the programs directory"
	    
	    read -s "Would you like to continue the installation? [y/n]: " contInstall
	    
	   if [[ $contInstall == "n" || $contInstall == "N" ]];
	    then
	    	echo "aborting installation"
	    	exit 1
    	    else 
    	    	echo "moving on to next installation" 
    	    fi
	else
	    echo "samblaster found"
	fi

}


#lumpy-SV
lumpy() {
	echo "installing lumpy"
	mkdir $program_dir/lumpy-sv
	git clone --recursive "https://github.com/arq5x/lumpy-sv.git" $program_dir/lumpy-sv/
	cd $program_dir/lumpy-sv
	export ZLIB_PATH="/usr/lib/x86_64-linux-gnu/"
	make
	cd $master_dir
	echo "checking sanity of Lumpy"
	export PATH=$PATH:"${master_dir}/programs/lumpy-sv/bin/"
	if ! [ -x "$(command -v lumpyexpress)" ];
	then 
	    echo "lumpy installation failed"
	    echo "Follow these docs: https://github.com/arq5x/lumpy-sv.git"
	    echo "The files need to be in the programs directory"
	    
	    read -s  "Would you like to continue the installation? [y/n]: " contInstall
	    
	   if [[ $contInstall == "n" || $contInstall == "N" ]];
	    then
	    	echo "aborting installation"
	    	exit 1
    	    else 
    	    	echo "moving on to next installation" 
    	    fi
    	    
	else
	    echo "lumpy found"
	fi
}

#SOAPdenovo2				


#delly works
delly() {
	mkdir $program_dir/delly
	git clone --recursive "https://github.com/dellytools/delly.git" $program_dir/delly
	cd $program_dir/delly
	make all
	cd $master_dir
	echo "checking sanity of delly"
	export PATH=$PATH:"${master_dir}/programs/delly/src/"
	if ! [ -x "$(command -v delly)" ];
	then 
	    echo "Delly installation failed"
	    echo "Follow these docs: https://github.com/dellytools/delly.git"
	    echo "The files need to be in the programs directory"
	     
	    read -s  "Would you like to continue the installation? [y/n]: " contInstall
	    
	    if [[ $contInstall == "n" || $contInstall == "N" ]];
	    then
	    	echo "aborting installation"
	    	exit 1
    	    else 
    	    	echo "moving on to next installation" 
    	    fi
    	    
	else
	    echo "delly found"
	fi
}

#svprops works
svprops() {
	echo "installing svprops"
	mkdir $program_dir/svprops
	git clone --recursive "https://github.com/dellytools/svprops.git" $program_dir/svprops
	cd $program_dir/svprops
	make 
	cd $master_dir
	echo "Sanity checks"
	export PATH=$PATH:"${master_dir}/programs/svprops/src/"
	if ! [ -x "$(command -v svprops)" ];
	then 
	    echo "svprops installation failed"
	    echo "Follow these docs: https://github.com/dellytools/svprops.git"
	    echo "The files need to be in the programs directory"
	     
	    read -s  "Would you like to continue the installation? [y/n]: " contInstall
	    
	    if [[ $contInstall == "n" || $contInstall == "N" ]];
	    then
	    	echo "aborting installation"
	    	exit 1
    	    else 
    	    	echo "moving on to next installation" 
    	    fi
    	    
	else
	    echo "svprops found"
	fi
}

HTSlib() {
	cd /usr/bin
	wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
	tar -vxjf htslib-1.9.tar.bz2
	cd htslib-1.9
	make
	cd $master_dir
}

zips() {
	echo "installing zipped files"
	cd $program_dir
	zipFiles=$(ls | grep "zip")
	for i in ${zipFiles[@]};
	do
		echo $i
		unzip $i
		rm -r $i
	done
}


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
if ! [ -x "$(command -v svtyper)" ];
then 
    "installing svtyper"  
     
    pip3 install git+https://github.com/hall-lab/svtyper.git
else
    echo "svtyper found"
fi
#multithread processes
sudoInstaller "samtools" 
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
	--field="Launch GUI":CHK \
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
	PegasusGUI.sh
fi

if [[ "${term}" == TRUE ]];
then
	echo "launching terminal"
	Pegasus.sh
fi

rm temp.txt

exit 0

else
    echo "internet connection required"
fi
