#!/bin/bash

#frontend GUI
#Im gonna keep the original and modify this for Fusion calling specifically

master_dir=$(pwd "$0")

input=$(yad --width=400 --title="Pegasus" --item-separator="," --separator="\t" --image="${master_dir}/misc/logo.png" --image-on-top  \
	--form \
	--field="Reference":CB  \
	--field="Illumina raw files":DIR \
	--field="sample list":FL \
	--field="output directory":DIR \
	--field="Threads" \
	--field="Ram" \
	--field="Target regions":FL \
	--field="Call structural variants":CHK  \
	--field="Search gene fusions":CHK  \
	--field="Debug mode":CHK \
	  'H37Rv,CDC1551')

#grap variables
echo ${input} > "${master_dir}/usr.txt"
ref_user=$(cat "${master_dir}/usr.txt" | awk '{print $1}')
echo "${ref_user}"

raw_files=$(cat "${master_dir}/usr.txt" | awk '{print $2}')
echo  "${raw_files}"

files_in=$(cat "${master_dir}/usr.txt" | awk '{print $3}')
echo "${files_in}"

out_dir=$(cat "${master_dir}/usr.txt" | awk '{print $4}')
echo "${out_dir}"

threads=$(cat "${master_dir}/usr.txt" | awk '{print $5}')
echo "${threads}"

ram=$(cat "${master_dir}/usr.txt" | awk '{print $6}')
echo "${ram}"

regions_of_interest=$(cat "${master_dir}/usr.txt" | awk '{print $7}')
echo "${regions_of_interest}"

sv_calling=$(cat "${master_dir}/usr.txt" | awk '{print $8}')
echo "${sv_calling}"

gene_fusions=$(cat "${master_dir}/usr.txt" | awk '{print $9}')
echo "${gene_fusions}"

debug=$(cat "${master_dir}/usr.txt" | awk '{print $10}')
echo "${debug}"

rm "${master_dir}/usr.txt"

#start pipeline here
echo "Executing pipeline" 

if [[ "${debug}" == "TRUE" ]];
	then 
	set -euxo pipefail
	
fi

#requires python 2.7
#requires java
trim="${master_dir}/programs/Trimmomatic-0.36/trimmomatic-0.36.jar"											#works
trim_PE="${master_dir}/programs/Trimmomatic-0.36/adapters/TruSeq2-PE.fa"									#works									
picard="${master_dir}/programs/build/libs/picard.jar"														#this file path will change NB
gatk="${master_dir}/programs/GenomeAnalysisTK.jar"															#still a problem
lumpy_extract_splitters="${master_dir}/programs/lumpy-sv/scripts/extractSplitReads_BwaMem"					#double check path			
lumpy_vcftobed="${master_dir}/programs/lumpy-sv/scripts/vcfToBedpe"											#works
export PATH=$PATH:"${master_dir}/programs/bwa"																#works																						
export PATH=$PATH:"${master_dir}/programs/novocraft"														#drequired manual download
export PATH=$PATH:"${master_dir}/programs/lumpy-sv"															#issue with make
export PATH=$PATH:"${master_dir}/programs/svtyper-master"													#required manual download
export PATH=$PATH:"${master_dir}/programs/FastQC"															#required manual download
export PATH=$PATH:"${master_dir}/programs/SOAPdenovo2"														#works
export PATH=$PATH:"${master_dir}/programs/svprops/src/"														#works
export PATH=$PATH:"${master_dir}/programs/samblaster"														#works

#create directory structures and pre analysi checks
#mkdir "${out_dir}"
log="${out_dir}/log.txt"

#checking input files 
# testing input file integrity

if [[ -s "$files_in" ]];
   then
   echo "data found in sample file" >> "$log"
   else
       echo "your sample file is empty, now exiting" >> "$log"
       exit 1
fi 

if [[ -s "$regions_of_interest" ]];
   then
   echo "data found in regions file" >> "$log"
   else
       echo "your regions file is empty, assuming discovery option" >> "$log"
      
       
fi

#set the reference
if [[ "${ref_user}" == H37Rv ]];  #I want to add a user based reference but will do that later
   then
   echo "H37Rv was chosen as a reference" >> "$log"
   ref="${master_dir}/references/H37Rv/H37Rv.fasta"
   ref_novo="${master_dir}/references/H37Rv/H37Rv.ndx"	
fi

if [[ "${ref_user}" == CDC1551 ]];  
   then
   echo "H37Rv was chosen as a reference" >> "$log"
   ref="${master_dir}/references/CDC1551/CDC1551.fasta"
   ref_novo="${master_dir}/references/CDC1551/CDC1551.ndx"	
fi
  
if [[ "${sv_calling}" == "TRUE" ]];
  then
  SV="${out_dir}/SV"
  mkdir "${SV}"
fi

#get functions for pipe
source "${master_dir}/scripts/main_functions.sh"      

#start the mainloop
echo "Entering Pegasus mainloop"
while IFS='' read -r sample || [[ -n "$sample" ]];  
do
	echo "Starting analysis of ${sample}" >> "$log" 
    	mkdir "${out_dir}/${sample}"
    	out_dir_2="${out_dir}/${sample}"
   	data="${out_dir_2}/data"
   	mkdir "${data}"
   	temp="${out_dir_2}/temp"
   	mkdir "${temp}"  
   	targets="${out_dir_2}/targets"  
   	mkdir "${targets}"
	fasta_header=$(cat ${ref} | awk 'sub(/^>/, "")')
	
	if [[ $sv_calling == "TRUE" ]];
  	then 
		sv_dir="${out_dir_2}/structural_variants"
		mkdir "${sv_dir}"
	fi
	
	if [[ $gene_fusions == "TRUE" ]];
	then
		gfdir="${out_dir_2}/fusion_genes"
	fi
	
	#set raw files vals
	raw_1="${raw_files}/${sample}_R1_001.fastq.gz"
   	raw_2="${raw_files}/${sample}_R2_001.fastq.gz"
	
	#stats
	fastqc_1
	#trim
	trimmomatic
	#stats
	fastqc_2
	
	#done with initial analysis, move on to alingments
	#use multiprocessing
	echo "Initialising best practises workflow" >> "${log}"
	BWA
  
	novo 
	wait
	
	echo "Best practises is done moving on to variant calling" >> "${log}"

#	#Structural variant calling: we will intersect the whole thing but lets first check the output
	if [[ $sv_calling == "TRUE" ]];
	then
		lumpy_bwa &
		lumpy_novo &
		delly_bwa &
		wait 
		lumpy_delly_isec 
	fi
	
	if [[ "${gene_fusions}" == "TRUE" ]]; # made gene fusions optional here. Will still do a first pass over the SV calling and call GF if asked for
	then 
		lumpy_bwa &
		lumpy_novo &
		delly_bwa &
		wait 
		lumpy_delly_isec 
		wait
		gene_fusion_calling
	fi
	
	#targets
	#Have to make this optional but still tie in with GF calling
	if [[ $regions_of_interest == "" ]];
	then
		echo "Target regions file not detected" >> "${log}"
		echo "Assuming only discovery based detection" >> "${log}"
		
		else
			while read name start end
			do 
				echo "starting analysis on region ${name} from ${sample}"
				regions="${targets}/${name}"
				mkdir "$regions"
				target_soap="${regions}/denovo"
				mkdir "${target_soap}"
				target_ordering="${regions}/ordered_contig"
				mkdir "${target_ordering}"
	
				target_region_process

			done< <(tr -d '\r' < "$regions_of_interest")
	fi
	rm -r "${temp}"	
done<$files_in

yad --title="Pegasus" \
	--form \
	--text "Done with Analysis"
exit 0
