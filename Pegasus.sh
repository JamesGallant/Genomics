#!/bin/bash
#to do: Make structural variant calling optional, i e add new var to the menu also change this in the help menu. Will make life easier. 
#		Could also consider adding a postpipe comeback script. Which will only analyse targets etc. 
# master pipeline
time_start="$(date +%H.%m.%s.%N)"
master_dir="$(pwd "$0")"
echo "${master_dir}"

if  [ "${1}" == manual ];            
    then
    less "${master_dir}/databases/manual.txt"
exit 0
fi

if [[ "${1}" == -h || "${1}" == --help || "${1}" == "" ]];
	then
	cat "${master_dir}/databases/help.txt"
exit 0
fi
#set environment variables

raw_files="$1"    				#Path to raw file directory
files_in="$2" 					#path to sample names
out_dir="$3" 					#path to output directory
threads="$4"					#amount of cores
ram="$5"						#amount of ram
snp_calling="$6" 				#snp calling yes or no
sv_calling="$7"					#structural variant calling
lineage_calling="$8"			#only compatible with H37Rv
gene_fusions="${9}"				#call gene fusions
tree="${$10}"					#get files for tree
ref_user="${11}"				#reference file
#ref_novo=${10}					#reference index file for novo.ndx
regions_of_interest="${12}"		#file to target list
debug="${13}"					#run in verbose mode

if [[ "${debug}" == "TRUE" ]];
	then 
	set -euxo pipefail
	
fi

#set environment
#NB docker imange 

trim="${master_dir}/programs/Trimmomatic-0.36/trimmomatic-0.36.jar"
trim_PE="${master_dir}/programs/Trimmomatic-0.36/adapters/TruSeq2-PE.fa"
TBprofiler_path="${master_dir}/programs/TBProfiler/tb-profiler"
picard="${master_dir}/programs/picard.jar"
gatk="${master_dir}/programs/GenomeAnalysisTK.jar" #this works on java 1.8 we need gatk-3.8
lumpy_extract_splitters="${master_dir}/programs/lumpy-sv/scripts/extractSplitReads_BwaMem"
lumpy_vcftobed="${master_dir}/programs/lumpy-sv/scripts/vcfToBedpe"
#export PATH=$PATH:"${master_dir}/programs/bwa-0.7.12/"
export PATH=$PATH:"${master_dir}/programs/novocraft"
export PATH=$PATH:"${master_dir}/programs/delly/src/"
export PATH=$PATH:"${master_dir}/programs/lumpy-sv/bin/"
#export PATH=$PATH:"${master_dir}/programs/svtyper-master"
#export PATH=$PATH:"${master_dir}/programs/FastQC"
#export PATH=$PATH:"${master_dir}/programs/SOAPdenovo2"
export PATH=$PATH:"${master_dir}/programs/svprops/src/"
export PATH=$PATH:"${master_dir}/programs/samblaster/"
#samtools,bcftools,abacas

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
  
if [[ "${ref_user}" == Marinum ]];  
   then
   echo "Marinum M was chosen as a reference" >> "$log"
   ref="${master_dir}/reference/Marinum/M_marinumM.fasta"
   ref_novo="${master_dir}/reference/Marinum/M_marinumM.ndx"	
fi      
       


#set up directories of common fils NB ln-s (link) [[[I think i want to remove this]]]
if [[ "${snp_calling}" == "TRUE" ]];
  then
  vcf="${out_dir}/vcf"
  mkdir "${vcf}"
fi

if [[ "${sv_calling}" == "TRUE" ]];
  then
  SV="${out_dir}/SV"
  mkdir "${SV}"
fi

#source the functions

source "${master_dir}/Pegasus_functions.sh"

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
   	targets="${out_dir_2}/targets"  
   	mkdir "${targets}"
	fasta_header=$(cat ${ref} | awk 'sub(/^>/, "")')

	if [[ $sv_calling == "TRUE" ]];
  	then 
		sv_dir="${out_dir_2}/structural_variants"
		mkdir "${sv_dir}"
		chimera_dir="${out_dir_2}/chimeric_genes"
		mkdir "${chimera_dir}"
	fi
	
	if [[ $snp_calling == "TRUE" ]];
  	then 
		snp_dir="${out_dir_2}/snps"
		mkdir "${snp_dir}"
	fi

	if [[ $lineage_calling == "TRUE" ]];
  	then 
		lineage_dir="${out_dir_2}/TBprofiler"
		mkdir "${lineage_dir}"
	fi

	if [[ "${tree}" == "TRUE" ]];
	then
		tree_dir="${out_dir_2}/Tree"
		mkdir "${tree_dir}"
	fi

    	raw_1="${raw_files}/${sample}_R1_001.fastq.gz"
   	raw_2="${raw_files}/${sample}_R2_001.fastq.gz"

#start with analysis
#TBprofiler if lineage was requested
	if [[ $lineage_calling == "TRUE" ]];
	then 
		TBprofiler
	else
		echo "TBprofiler not initiated, moving on" "${log}"

	fi
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

#GATK best practises
	
	
#variant calling
	if [[ "${snp_calling}" == "TRUE" ]];
	then 
		snp_gatk_novo &
		snp_gatk_bwa &
		snps_pileup_novo &
		snps_pileup_bwa &
		wait
		intersect_snps_indel
	
		echo "snp calling done" >> "${log}"
		#tree functions?
		if [[ "${tree}" == "TRUE" ]];
		then 
			echo "Drawing files for trees, send through snpTREE" >> "${log}"
			tree_draw_files
		fi
	fi

#Structural variant calling: we will intersect the whole thing but lets first check the output
	if [[ $sv_calling == "TRUE" ]];
	then
		lumpy_bwa &
		lumpy_novo &
		delly_bwa &
		wait 
		lumpy_delly_isec 
			if [[ "${gene_fusions}" == "TRUE" ]];
			then 
				gene_fusion_calling
			fi
	fi

#moving on to targeted so long, we can get back later
	if [[ -s "$regions_of_interest" ]];
   	then
   	echo "starting on targeted analysis" >> "$log"
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
   		
		else
       	 		echo "skipping region specific analysis" >> "$log"	
	fi

rm -r "${temp}"	
done < $files_in
time_end=$(date +%H.%m.%s.%N)
time_diff=${echo $time_start - $time_end | bc}
echo "it takes ${time_diff} to execute the pipeline" >> "${log}"
exit 0
