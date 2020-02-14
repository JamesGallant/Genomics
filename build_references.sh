#!/bin/bash
master_dir="$(pwd "$0")"
export PATH=$PATH:"${master_dir}/programs/novocraft"
programs="${master_dir}/programs"
picard="${programs}/picard.jar"

ref_dir="${master_dir}/references"
mkdir $ref_dir

H37Rv_dir="${ref_dir}/H37Rv"
CDC1551_dir="${ref_dir}/CDC1551"
#Marinum_dir="${ref_dir}/Marinum"

mkdir ${H37Rv_dir} & mkdir ${CDC1551_dir} 
#mkdir ${Marinum_dir}


source "${master_dir}/install_functions.sh"

H37Rv
CDC1551

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


