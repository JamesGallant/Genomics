#!/bin/bash
master_dir="$(pwd "$0")"
export PATH=$PATH:"${master_dir}/programs/novocraft"
programs="${master_dir}/programs"
picard="${programs}/picard.jar"

ref_dir="${master_dir}/references"
mkdir $ref_dir

H37Rv_dir="${ref_dir}/H37Rv"
CDC1551_dir="${ref_dir}/CDC1551"
Marinum_dir="${ref_dir}/Marinum"

mkdir ${H37Rv_dir} & mkdir ${CDC1551_dir} 
mkdir ${Marinum_dir}


source "${master_dir}/Pegasus_functions.sh"

H37Rv
CDC1551
Marinum

exit 0


