#!/bin/bash
##this is a build file##
#to do: SRA toolkits, soapdenovo2 is a issue, samblaster and lumpy not checked, 
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
	"https://github.com/broadinstitute/picard/releases/download/2.21.6/picard.jar"
	"https://github.com/broadinstitute/gatk/releases/download/4.1.4.1/gatk-4.1.4.1.zip"
	)

for i in ${URL[@]};
do
	wget -P $master_dir/programs/  $i
done

#git installation functions
#samblaster
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
	    
	    read -s -p "Would you like to continue the installation? [y/n]: " contInstall
	    
	    if [$contInstall == "n" || $contInstall == "N"];
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
	make
	cd $master_dir
	echo "checking sanity of Lumpy"
	export PATH=$PATH:"${master_dir}/programs/lumpy-sv/"
	if ! [ -x "$(command -v lumpyexpress)" ];
	then 
	    echo "lumpy installation failed"
	    echo "Follow these docs: https://github.com/arq5x/lumpy-sv.git"
	    echo "The files need to be in the programs directory"
	    
	    read -s -p "Would you like to continue the installation? [y/n]: " contInstall
	    
	    if [$contInstall == "n" || $contInstall == "N"];
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
soap_install () {
	cd "${programs}"
	echo "Installing SOAPdenovo2"
	git clone --recursive "https://github.com/aquaskyline/SOAPdenovo2.git"
	cd SOAPdenovo2
	make 
	cd $master_dir
	echo "Sanity checks"
	export PATH=$PATH:"${master_dir}/programs/SOAPdenovo2"
	if ! [ -x "$(command -v soapdenovo-63mer)" ];
	then 
	    echo "soapdenovo installation failed"
	    echo "Follow these docs: https://github.com/aquaskyline/SOAPdenovo2.git"
	    echo "The files need to be in the programs directory"
	     
	    read -s -p "Would you like to continue the installation? [y/n]: " contInstall
	    
	    if [$contInstall == "n" || $contInstall == "N"];
	    then
	    	echo "aborting installation"
	    	exit 1
    	    else 
    	    	echo "moving on to next installation" 
    	    fi
    	    
	else
	    echo "SOAPdenovo2 found"
	fi					
}

#delly works
delly() {
	mkdir $program_dir/delly
	git clone --recursive "https://github.com/dellytools/delly.git" $program_dir/delly
	cd $program_dir/delly
	make 
	cd $master_dir
	echo "checking sanity of delly"
	export PATH=$PATH:"${master_dir}/programs/delly/src/"
	if ! [ -x "$(command -v delly)" ];
	then 
	    echo "Delly installation failed"
	    echo "Follow these docs: https://github.com/dellytools/delly.git"
	    echo "The files need to be in the programs directory"
	     
	    read -s -p "Would you like to continue the installation? [y/n]: " contInstall
	    
	    if [$contInstall == "n" || $contInstall == "N"];
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
	     
	    read -s -p "Would you like to continue the installation? [y/n]: " contInstall
	    
	    if [$contInstall == "n" || $contInstall == "N"];
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
if ! [ -x "$(command -v svtyper)" ];
then 
    "installing svtyper"   
    pip install git+https://github.com/hall-lab/svtyper.git
else
    echo "svtyper found"
fi
#multithread processes
sudoInstaller "samtools" &
sudoInstaller "bwa" &
sudoInstaller "bcftools" && fg
wait
sudoInstaller "abacas" &
sudoInstaller "svtyper" &
sudoInstaller "fastqc" && fg
wait
sudoInstaller "tabix" &
sudoInstaller "sambamba" &
sudoInstaller "yad" && fg
wait
sudoInstaller "g++" &
sudoInstaller "cnake" &
sudoInstaller "gawk" && fg
wait
samblaster &
lumpy && fg
wait 
delly &
svprops && fg
wait
soap_install &
zips && fg

echo "changing write permissions"
declare -a arr2=("lumpy-sv" "picard" "trimmomatic-0.36" "delly" "bwa" "SOADPdenovo2" "svprops" "samblaster")

#change permissions
for a in "${arr2[@]}";
do
	echo "changing permissions"
	chmod 777 "${a}"
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

chmod 755 pegasus.sh
chmod 755 pegasusGUI.sh

if [[ "${GUI}" == TRUE ]];
then
	echo "Launching gui"
	./PegasusGUI.sh
fi

if [[ "${term}" == TRUE ]];
then
	echo "launching terminal"
	./main.sh
fi

rm ./temp.txt

exit 0

else
    echo "internet connection required"
fi
