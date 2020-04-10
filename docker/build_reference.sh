#!/bin/bash
master_dir="$(pwd "$0")"
export PATH=$PATH:"${master_dir}/programs/novocraft"
database_dir="${master_dir}/databases"
if [[ "${1}" == -h || "${1}" == --help || "${1}" == "" ]];
	then
	cat "${database_dir}/help_refs.txt"
exit 0
fi

addUserRef="$1"
addUserAnno="$2"
userRefName="$3"

programs="${master_dir}/programs"
picard="${programs}/picard.jar"
ref_dir="${master_dir}/references"

source "${master_dir}/scripts/install_functions.sh"

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
	
	if [[ $addUserRef =~ \.fa$ || $addUserRef =~ \.fna$ ]]; then
		echo "your file, $addUserRef, has a .fa/.fna extension. Please rename to .fasta"
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
	database_dir="${master_dir}/databases"
	python3 "${master_dir}/scripts/annotationHandler.py" --file "${addUserAnno}" --out "${database_dir}/${userRefName}_DB.txt"
	echo "${database_dir}/${userRefName}_DB.txt"
	
	
	
fi

exit 0




