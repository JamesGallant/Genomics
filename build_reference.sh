#!/bin/bash
master_dir="$(pwd "$0")"
export PATH=$PATH:"${master_dir}/programs/novocraft"
programs="${master_dir}/programs"
picard="${programs}/picard.jar"
ref_dir="${master_dir}/references"
mkdir $ref_dir
source "${master_dir}/scripts/install_functions.sh"

getRefData=$(yad --width=400 --title="Pegasus reference settings" --item-separator="," --separator="\t" \
	    --text="Choose installation settings"\
	    --form \
	    --field="Install M. tuberculosis references":CHK \
  	    --field="Install your custom reference":FL\
	    --field="Add custom annotation file:":FL \
	    --field="Name of your reference")
	   
defaultRef=$(echo $getRefData | awk '{print $1}')  
addUserRef=$(echo $getRefData | awk '{print $2}')
addUserAnno=$(echo $getRefData | awk '{print $3}')
userRefName=$(echo $getRefData | awk '{print $4}')


if [[ $defaultRef == "TRUE" ]]; then
	H37Rv_dir="${ref_dir}/H37Rv"
	CDC1551_dir="${ref_dir}/CDC1551"
	
	mkdir ${H37Rv_dir} & mkdir ${CDC1551_dir}
	H37Rv
	CDC1551
fi

if [[ $addUserRef == "" || $addUserAnno == "" ]]; then
	echo "You did not choose a custom reference or files are missing"
else
	if [[ $userRefName == "" ]]; then
		userRefName="customReference"
	fi

	echo "adding your reference $userRefName"

	if [[ $addUserRef =~ \.gz$ ]]; then
		echo "Your file is zipped, please unzip and reupload"
		exit 1
	fi

	refSanity=$(grep ">" $addUserRef)

	if [[ $refSanity == "" ]]; then
		echo "$addUserRef is not a fasta file, aborting"
		exit 1
	fi

	userRef_dir="${ref_dir}/$userRefName"
	mkdir $userRef_dir
	cp $addUserRef $userRef_dir
	customRefTemp=$( ls $userRef_dir | grep ".fasta")
	userRefOld=$(basename "${userRef_dir}/${customRefTemp}" .fasta)
	mv "${userRef_dir}/${userRefOld}.fasta" "${userRef_dir}/${userRefName}.fasta"

	sed -i "s/^.*>.*$/>$userRefName/" "${userRef_dir}/${userRefName}.fasta"

	#Indexing
	samtools faidx "${userRef_dir}/${userRefName}.fasta"
	bwa index "${userRef_dir}/${userRefName}.fasta"
	java -jar ${picard} CreateSequenceDictionary R="${userRef_dir}/${userRefName}.fasta" O="${userRef_dir}/${userRefName}.dict" 
	novoindex "${userRef_dir}/${userRefName}.ndx" "${userRef_dir}/${userRefName}.fasta"

	#build annotation DB
	python3 /scripts/annoHandler.py --file "${addUserAnno}"
	
	
	
fi


startup=$(yad --item-separator="," --separator="\t" --width=400 \
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




