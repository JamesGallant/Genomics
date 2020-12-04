#!/bin/bash
#to do: Make structural variant calling optional, i e add new var to the menu also change this in the help menu. Will make life easier. 
#		Could also consider adding a postpipe comeback script. Which will only analyse targets etc. 
# master pipeline
time_start="$(date +%H.%m.%s.%N)"
master_dir="$(pwd "$0")"
echo "${master_dir}"
database_dir="${master_dir}/databases"

if  [ "${1}" == manual ];            
    then
    less "${master_dir}/README.md"
exit 0
fi

if [[ "${1}" == -h || "${1}" == --help || "${1}" == "" ]];
	then
	cat "${database_dir}/help.txt"
exit 0
fi
#set environment variables

raw_files="$1"    				#Path to raw file directory
files_in="$2" 					#path to sample names
out_dir="$3" 					#path to output directory
threads="$4"					#amount of cores
ram="$5"						#amount of ram
gene_fusions="${6}"				#call gene fusions
ref_user="${7}"				#reference file
debug="${8}"					#run in verbose mode

if [[ "${debug}" == "TRUE" ]];
	then 
	set -euxo pipefail
	
fi

#set environment
#NB docker imange 

trim="${master_dir}/programs/Trimmomatic-0.36/trimmomatic-0.36.jar"
trim_PE="${master_dir}/programs/Trimmomatic-0.36/adapters/TruSeq2-PE.fa"
picard="${master_dir}/programs/picard.jar"
#this works on java 1.8 we need gatk-3.8
gatk="${master_dir}/programs/GenomeAnalsysTK.jar" 
lumpy_extract_splitters="${master_dir}/programs/lumpy-sv/scripts/extractSplitReads_BwaMem"
lumpy_vcftobed="${master_dir}/programs/lumpy-sv/scripts/vcfToBedpe"
export PATH=$PATH:"${master_dir}/programs/novocraft"
export PATH=$PATH:"${master_dir}/programs/delly/src/"
export PATH=$PATH:"${master_dir}/programs/lumpy-sv/bin/"
export PATH=$PATH:"${master_dir}/programs/svtyper"
export PATH=$PATH:"${master_dir}/programs/svprops/src/"
export PATH=$PATH:"${master_dir}/programs/samblaster/"


#create directory structures and pre analysi checks
#mkdir "${out_dir}"
log="${out_dir}/log.txt"

#Environment check
echo "Checking user required tools...." >> "$log"
#samtools
if ! [ -x "$(command -v samtools)" ];
then
    echo "samtools is not found in your system, Please install samtools in your ./bin directory" >> "$log"
    echo "Aborting..."
    exit 1
else
    echo "samtools found" >> "$log"
fi
#bcftools
if ! [ -x "$(command -v bcftools)" ];
then
    echo "bcftools is not found in your system, Please install bcftools in your ./bin directory" >> "$log"
    echo "Aborting..."
    exit 1
else
    echo "bcftools found" >> "$log"
fi

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


#set the reference and annotation
echo "$ref_user was chosen as the reference" >> "${log}"
ref="${master_dir}/references/$ref_user/${ref_user}.fasta"
ref_novo="${master_dir}/references/$ref_user/${ref_user}.ndx"
anno_DB="${database_dir}/${ref_user}_DB.txt"

    
#set up directories of common fils NB ln-s (link) [[[I think i want to remove this]]]
if [[ "${sv_calling}" == "TRUE" ]];
  then
  SV="${out_dir}/SV"
  mkdir "${SV}"
fi

#final results
results_dir="${out_dir}/results"
mkdir "${results_dir}"

#source the functions

source "${master_dir}/scripts/Pegasus_functions.sh"

#start the mainloop

while IFS='' read -r sample || [[ -n "$sample" ]];  
do
	echo "Starting analysis of ${sample}" >> "$log" 
    	mkdir "${out_dir}/${sample}"
    	out_dir_2="${out_dir}/${sample}"
   	data="${out_dir_2}/data"
   	mkdir "${data}"
   	temp="${out_dir_2}/temp"
   	mkdir "${temp}"  
	results_sample="${results_dir}/${sample}"
	mkdir "${results_sample}"
	
	
	fasta_header=$(cat ${ref} | awk 'sub(/^>/, "")')

	if [[ ${gene_fusions} == "TRUE" ]];
  	then 
		sv_dir="${out_dir_2}/structural_variants"
		mkdir "${sv_dir}"
		chimera_dir="${out_dir_2}/chimeric_genes"
		mkdir "${chimera_dir}"
	fi

    	raw_1="${raw_files}/${sample}_R1_001.fastq.gz"
   	raw_2="${raw_files}/${sample}_R2_001.fastq.gz"


	#fastqc_1
#trim
	trimmomatic
#stats
	#fastqc_2

#done with initial analysis, move on to alingments
#use multiprocessing
	echo "Initialising best practises workflow" >> "${log}"
	BWA
  
	novo 
	wait
	
	echo "Best practises is done moving on to variant calling" >> "${log}"


	
#variant calling

#Structural variant calling: we will intersect the whole thing but lets first check the output

	if [[ "${gene_fusions}" == "TRUE" ]];
		then 
			lumpy_bwa &
			lumpy_novo &
			delly_bwa &
			wait 
			lumpy_delly_isec 
			gene_fusion_calling
	fi

#moving on to targeted so long, we can get back later
	
	
rm -r "${temp}"	
done < $files_in
exit 0
