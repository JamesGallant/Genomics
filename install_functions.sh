######InstallFunctions#####
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

######ENDofInstalls###########

######ReferenceFunctions######

H37Rv() {
	cd ${H37Rv_dir}
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz
	ls ${H37Rv_dir}
	gunzip "${H37Rv_dir}/GCF_000195955.2_ASM19595v2_genomic.fna.gz"
	mv "${H37Rv_dir}/GCF_000195955.2_ASM19595v2_genomic.fna" "${H37Rv_dir}/H37Rv.fasta"
	
	samtools faidx "${H37Rv_dir}/H37Rv.fasta"

	bwa index "${H37Rv_dir}/H37Rv.fasta"
	
	java -jar ${picard} CreateSequenceDictionary R="${H37Rv_dir}/H37Rv.fasta" O="${H37Rv_dir}/H37Rv.dict" 
	   
	novoindex "${H37Rv_dir}/H37Rv.ndx" "${H37Rv_dir}/H37Rv.fasta"
	cd ${master_dir}
} 

CDC1551() {
	cd ${CDC1551_dir}
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/008/585/GCA_000008585.1_ASM858v1/GCA_000008585.1_ASM858v1_genomic.fna.gz
	
	gunzip "${CDC1551_dir}/GCA_000008585.1_ASM858v1_genomic.fna.gz"
	mv "${CDC1551_dir}/GCA_000008585.1_ASM858v1_genomic.fna" "${CDC1551_dir}/CDC1551.fasta"
	samtools faidx "${CDC1551_dir}/CDC1551.fasta"
	bwa index "${CDC1551_dir}/CDC1551.fasta"
	java -jar ${picard} CreateSequenceDictionary R="${CDC1551_dir}/CDC1551.fasta" O="${CDC1551_dir}/CDC1551.dict" 

	novoindex "${CDC1551_dir}/CDC1551.ndx" "${CDC1551_dir}/CDC1551.fasta"
	cd ${master_dir}
}

Marinum() {
	cd ${Marinum_dir}
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/723/425/GCF_000723425.2_E11/GCF_000723425.2_E11_genomic.fna.gz
	
	gunzip "${Marinum_dir}/GCF_000723425.2_E11_genomic.fna.gz"
	mv "${Marinum_dir}/GCF_000723425.2_E11_genomic.fna" "${Marinum_dir}/Marinum.fasta"
	samtools faidx "${Marinum_dir}/Marinum.fasta"
	bwa index "${Marinum_dir}/Marinum.fasta"
	java -jar ${picard} CreateSequenceDictionary R="${Marinum_dir}/Marinum.fasta" O="${Marinum_dir}/Marinum.dict" 

	novoindex "${Marinum_dir}/Marinum.ndx" "${Marinum_dir}/Marinum.fasta"
	cd ${master_dir}
}
######End#######
