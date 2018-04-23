
#trim reads
trimmomatic () {
		echo "Trimming reads..." >> "$log"
		java -Xmx"${ram}"g -jar $trim PE \
    		-phred33 \
    		-threads "$threads" \
    		"$raw_1" "$raw_2" \
    		"${temp}/${sample}_forward_paired.fq.gz" "${temp}/${sample}_forward_unpaired.fq.gz" \
    		"${temp}/${sample}_reverse_paired.fq.gz" "${temp}/${sample}_reverse_unpaired.fq.gz" \
    		ILLUMINACLIP:"${trim_PE}":2:30:15 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:30
		echo "done trimming" >> "${log}"
}

#TBprofiler
TBprofiler () {
		echo "Gathering drug resistance and lineage data" >>"${log}"
		cd "${lineage_dir}"
		python "${TBprofiler_path}" full -1 "${raw_1}" -2 "${raw_2}" -p "${sample}"
		echo "Lineage determination done" >> "${log}"
		cd "${master_dir}"
}
		
#fastqc
fastqc_1 () {
	   	fastqc "${raw_1}" -o "${data}"
		fastqc "${raw_2}" -o "${data}" 
}
fastqc_2 () {	
		
    		fastqc "${temp}/${sample}_forward_paired.fq.gz" -o "${data}"
    		fastqc "${temp}/${sample}_reverse_paired.fq.gz" -o "${data}"
		echo "fastqc done" >> "${log}"
}

#align with bwa
BWA () {
	echo "Initialising bwa alignment" >> "${log}"
	#BWA alignment
	 bwa mem -t "${threads}" \
		 -c 100 -R "@RG\\tID:${sample}\\tSM:${sample}\\tPL:Illumina" \
		 -M -T 50 \
		 "${ref}" \
		 "${temp}/${sample}_forward_paired.fq.gz" "${temp}/${sample}_reverse_paired.fq.gz" > "${temp}/${sample}_bwa_paired.sam"
	echo "bwa alingment done" >> "${log}"
	#validation of SAM file
        	java -Xmx"${ram}"g -jar "${picard}" ValidateSamFile \
       	        INPUT= "${temp}/${sample}_bwa_paired.sam" \
        	OUTPUT="${data}/${sample}_bwa_sam_validate.txt"

		if grep -q "No" "${data}/${sample}_bwa_sam_validate.txt" == "No";
        		then
           		echo "Sam file is error free" >> "$log"
        		else
           		 echo "SAM file has errors" >> "$log"
        fi

		samblaster \
			--addMateTags \
			-i "${temp}/${sample}_bwa_paired.sam" \
			-o "${temp}/${sample}_bwa_paired_sb.sam"
 
		#NBexport ram to java 		
		#convert SAM into BAM and sort
       	 	java -Xmx"${ram}"g -jar "${picard}" SortSam \
         	INPUT="${temp}/${sample}_bwa_paired_sb.sam" \
         	OUTPUT="${temp}/${sample}_bwa_sorted.bam" \
         	SORT_ORDER=coordinate \
         	VALIDATION_STRINGENCY=LENIENT
      
        	samtools index "${temp}/${sample}_bwa_sorted.bam" 

		#Mark PCR Duplicates 
        	java -Xmx"${ram}"g -jar "${picard}" MarkDuplicates \
                INPUT="${temp}/${sample}_bwa_sorted.bam" \
				OUTPUT="${temp}/${sample}_bwa_sorted_dedup.bam" \
				VALIDATION_STRINGENCY=LENIENT \
				REMOVE_DUPLICATES=TRUE \
				ASSUME_SORTED=TRUE \
				M="${temp}/${sample}_bwa_sorted_dedup.bam.txt"
         
          	#generate index  with picard
        	java -Xmx"${ram}"g -jar "${picard}" BuildBamIndex \
				I="${temp}/${sample}_bwa_sorted_dedup.bam" \
				VALIDATION_STRINGENCY=LENIENT
		
			#Create a target list of intervals to be realigned with GATK
        	java -jar "$gatk" -T  RealignerTargetCreator -R "${ref}" \
				  -I "${temp}/${sample}_bwa_sorted_dedup.bam" \
				  -o "${temp}/${sample}_bwa_sorted_dedup.bam.list"
        
        	#Realignment of target intervals        
        	java -Xmx"${ram}"g -jar "$gatk" \
				-T IndelRealigner \
				-R "${ref}" \
				-I "${temp}/${sample}_bwa_sorted_dedup.bam" \
				-targetIntervals "${temp}/${sample}_bwa_sorted_dedup.bam.list" \
				-o "${temp}/${sample}_bwa_sorted_dedup_realigned.bam"
         
         	samtools sort -@ $threads -m "${ram}"g -O "bam" -T "working" -o "${temp}/${sample}_bwa_sorted_dedup_realigned_sorted.bam" "${temp}/${sample}_bwa_sorted_dedup_realigned.bam"
		
		if [[ -s "${temp}/${sample}_bwa_sorted_dedup_realigned_sorted.bam" ]];
   		then
   		echo "data found in your BAM file" >> "$log"
   		else
       			echo "your BAM file is empty, now exiting" >> "$log"
      			exit 1
		fi	      
		#New index, sort with samtools 
         	samtools index "${temp}/${sample}_bwa_sorted_dedup_realigned_sorted.bam" 
		#Clean up
         	mv "${temp}/${sample}_bwa_sorted_dedup_realigned_sorted.bam" "${data}/${sample}_bwa_sorted_dedup_realigned_sorted.bam"  
           	mv "${temp}/${sample}_bwa_sorted_dedup_realigned_sorted.bam.bai" "${data}/${sample}_bwa_sorted_dedup_realigned_sorted.bam.bai"
			
		echo "Done with GATK best practises for bwa aligned files" >> "${log}"
} 

#align with novoalign
novo () {
	echo "Initialising NOVOalign" >> "${log}"
	gunzip "${temp}/${sample}_forward_paired.fq.gz"
    gunzip "${temp}/${sample}_reverse_paired.fq.gz"
	novoalign -d "${ref_novo}" \
		  -f "${temp}/${sample}_forward_paired.fq" "${temp}/${sample}_reverse_paired.fq" \
		  -o SAM "@RG\\tID:${sample}\\tSM:${sample}\\tPL:Illumina" \
		  2> "${data}/${sample}_novo_stats.novodist"  > "${temp}/${sample}_novo_paired.sam"

	if grep -q "error" "${data}/${sample}_novo_stats.novodist" == "error";
        then
           echo "Novoalign failed, now aborting..." >> "$log"
           exit 1
        else
           echo "Novoalign completed successfully" >> "$log"
	fi

	#validation of SAM file
        	java -Xmx"${ram}"g -jar "${picard}" ValidateSamFile \
       	        INPUT= "${temp}/${sample}_novo_paired.sam" \
				OUTPUT="${data}/${sample}_novo_sam_validate.txt"

		if grep -q "No" "${data}/${sample}_novo_sam_validate.txt" == "No";
        		then
           		echo "Sam file is error free" >> "$log"
        		else
           		 echo "SAM file has errors" >> "$log"
        	fi
		
			samblaster \
				--addMateTags \
				-i "${temp}/${sample}_novo_paired.sam" \
				-o "${temp}/${sample}_novo_paired_sb.sam" 


			#convert SAM into BAM and sort
       	 	java -Xmx"${ram}"g -jar "${picard}" SortSam \
				INPUT="${temp}/${sample}_novo_paired_sb.sam" \
				OUTPUT="${temp}/${sample}_novo_sorted.bam" \
				SORT_ORDER=coordinate \
				VALIDATION_STRINGENCY=LENIENT
      
        	samtools index "${temp}/${sample}_novo_sorted.bam" 

		
			#Mark PCR Duplicates 
        	java -Xmx"${ram}"g -jar "${picard}" MarkDuplicates \
                INPUT="${temp}/${sample}_novo_sorted.bam" \
				OUTPUT="${temp}/${sample}_novo_sorted_dedup.bam" \
				VALIDATION_STRINGENCY=LENIENT \
				REMOVE_DUPLICATES=TRUE \
				ASSUME_SORTED=TRUE \
				M="${temp}/${sample}_novo_sorted_dedup.bam.txt"
         
          	#generate index  with picard
        	java -Xmx"${ram}"g -jar "${picard}" BuildBamIndex \
				I="${temp}/${sample}_novo_sorted_dedup.bam" \
				VALIDATION_STRINGENCY=LENIENT
		
			#Create a target list of intervals to be realigned with GATK
        	java -jar "$gatk" \
				-T  RealignerTargetCreator \
				-R "${ref}" \
				-I "${temp}/${sample}_novo_sorted_dedup.bam" \
				-o "${temp}/${sample}_novo_sorted_dedup.bam.list"
        
        	#Realignment of target intervals        
        	java -Xmx"${ram}"g -jar "$gatk" \
				-T IndelRealigner \
				-R "${ref}" \
				-I "${temp}/${sample}_novo_sorted_dedup.bam" \
				-targetIntervals "${temp}/${sample}_novo_sorted_dedup.bam.list" \
				-o "${temp}/${sample}_novo_sorted_dedup_realigned.bam"
         
         	samtools sort -@ $threads -m "${ram}"g -O "bam" -T "working" -o "${temp}/${sample}_novo_sorted_dedup_realigned_sorted.bam" "${temp}/${sample}_novo_sorted_dedup_realigned.bam"

		
		if [[ -s "${temp}/${sample}_novo_sorted_dedup_realigned_sorted.bam" ]];
   		then
   		echo "data found in your BAM file" >> "$log"
   		else
       			echo "your BAM file is empty, now exiting" >> "$log"
      			exit 1
		fi
		#New index, sort with samtools 
         	samtools index "${temp}/${sample}_novo_sorted_dedup_realigned_sorted.bam" 

		#Clean up
         	mv "${temp}/${sample}_novo_sorted_dedup_realigned_sorted.bam" "${data}/${sample}_novo_sorted_dedup_realigned_sorted.bam"  
           	mv "${temp}/${sample}_novo_sorted_dedup_realigned_sorted.bam.bai" "${data}/${sample}_novo_sorted_dedup_realigned_sorted.bam.bai"
		echo "Done with GATK best practises for Novo aligned files" >> "${log}"

		#require insert length from PICARD

		java -Xmx"${ram}"g -jar "${picard}" CollectInsertSizeMetrics \
     			I="${data}/${sample}_novo_sorted_dedup_realigned_sorted.bam" \
      		 	O="${temp}/${sample}.picard.insert.metrics.tab" \
      			HISTOGRAM_FILE="${data}/${sample}_picard_insert_metrics.pdf"

		#get the insert length
		mean=$(cat "${temp}/${sample}.picard.insert.metrics.tab" | awk '{print $5}' | sed -n '8p' | cut -f1 -d".")
}

snp_gatk_novo () {
	     	 echo "snp calling has been chosen, initallising haplotype caller" >> "$log"  
		 #HaplotypeCaller - better at calling indels B
            	 java -Xmx"${ram}"g -jar "$gatk" \
                 -T HaplotypeCaller \
                 -ploidy 1 \
                 -R "${ref}" \
                 -I "${data}/${sample}_novo_sorted_dedup_realigned_sorted.bam" \
                 -stand_emit_conf 10 \
                 -stand_call_conf 30 \
                 --genotyping_mode DISCOVERY \
                 -o "${temp}/${sample}_novo_raw_variants_GATK_haplotype_caller.vcf"
	         
		 #extracting snps
           	 java -Xmx"${ram}"g -jar "$gatk" \
                -T SelectVariants \
                -R "${ref}" \
                -V "${temp}/${sample}_novo_raw_variants_GATK_haplotype_caller.vcf" \
                -selectType SNP \
                -o "${temp}/${sample}_novo_raw_snps.vcf"
             
             	#hard filtering of snps
           	 java -Xmx"${ram}"g -jar "$gatk" \
                -T VariantFiltration \
                -R "${ref}" \
                -V "${temp}/${sample}_novo_raw_snps.vcf" \
                --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
                --filterName "${sample}" \
                -o "${temp}/${sample}_gatk_novo_filt_snps.vcf"
             
             	#extract indels
             	java -Xmx"${ram}"g -jar "$gatk" \
                -T SelectVariants \
                -R "${ref}" \
                -V "${temp}/${sample}_novo_raw_variants_GATK_haplotype_caller.vcf" \
                -selectType INDEL \
                -o "${temp}/${sample}_novo_raw_INDEL.vcf"
            
           	 #filter
           	 java -Xmx"${ram}"g -jar "$gatk" \
                -T VariantFiltration \
                -R "${ref}" \
                -V "${temp}/${sample}_novo_raw_INDEL.vcf" \
                --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
                --filterName "${sample}" \
                -o "${temp}/${sample}_gatk_novo_filt_INDELs.vcf"

		  cp "${temp}/${sample}_gatk_novo_filt_INDELs.vcf" "${snp_dir}/${sample}_gatk_novo_filt_INDELs.vcf" \
		| cp "${temp}/${sample}_gatk_novo_filt_INDELs.vcf" "${data}/${sample}_gatk_novo_filt_INDELs.vcf" \
		| cp "${temp}/${sample}_gatk_novo_filt_INDELs.vcf" "${vcf}/${sample}_gatk_novo_filt_INDELs.vcf" \
		| cp "${temp}/${sample}_gatk_novo_filt_snps.vcf" "${snp_dir}/${sample}_gatk_novo_filt_snps.vcf" \
		| cp "${temp}/${sample}_gatk_novo_filt_snps.vcf" "${data}/${sample}_gatk_bwa_novo_snps.vcf" \
		| cp "${temp}/${sample}_gatk_novo_filt_snps.vcf" "${vcf}/${sample}_gatk_bwa_novo_snps.vcf"
}

snp_gatk_bwa () {
		#HaplotypeCaller - better at calling indels B
            	java -Xmx"${ram}"g -jar "$gatk" \
                -T HaplotypeCaller \
                -ploidy 1 \
                -R "${ref}" \
                -I "${data}/${sample}_bwa_sorted_dedup_realigned_sorted.bam" \
                -stand_emit_conf 10 \
                -stand_call_conf 30 \
                --genotyping_mode DISCOVERY \
                -o "${temp}/${sample}_bwa_raw_variants_GATK_haplotype_caller.vcf"
           	 # -L targets.interval_list for subset
            
           	 #extracting snps
            	java -Xmx"${ram}"g -jar "$gatk" \
                -T SelectVariants \
                -R "${ref}" \
                -V "${temp}/${sample}_bwa_raw_variants_GATK_haplotype_caller.vcf" \
                -selectType SNP \
                -o "${temp}/${sample}_bwa_raw_snps.vcf"
             
             	#hard filtering of snps
            	java -Xmx"${ram}"g -jar "$gatk" \
                -T VariantFiltration \
                -R "${ref}" \
                -V "${temp}/${sample}_bwa_raw_snps.vcf" \
                --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
                --filterName "${sample}" \
                -o "${temp}/${sample}_gatk_bwa_filt_snps.vcf"
             
             	#extract indels
             	java -Xmx"${ram}"g -jar "$gatk" \
                -T SelectVariants \
                -R "${ref}" \
                -V "${temp}/${sample}_bwa_raw_variants_GATK_haplotype_caller.vcf" \
                -selectType INDEL \
                -o "${temp}/${sample}_bwa_raw_INDEL.vcf"
            
            	#filter
           	 java -Xmx"${ram}"g -jar "$gatk" \
                -T VariantFiltration \
                -R "${ref}" \
                -V "${temp}/${sample}_bwa_raw_INDEL.vcf" \
                --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
                --filterName "${sample}" \
                -o "${temp}/${sample}_gatk_bwa_filt_INDELs.vcf"

		cp "${temp}/${sample}_gatk_bwa_filt_INDELs.vcf" "${snp_dir}/${sample}_gatk_bwa_filt_INDELs.vcf" \
		| cp "${temp}/${sample}_gatk_bwa_filt_INDELs.vcf" "${data}/${sample}_gatk_bwa_filt_INDELs.vcf" \
		| cp "${temp}/${sample}_gatk_bwa_filt_INDELs.vcf" "${vcf}/${sample}_gatk_bwa_filt_INDELs.vcf" \
		| cp "${temp}/${sample}_gatk_bwa_filt_snps.vcf" "${snp_dir}/${sample}_gatk_bwa_filt_snps.vcf" \
		| cp "${temp}/${sample}_gatk_bwa_filt_snps.vcf" "${data}/${sample}_gatk_bwa_filt_snps.vcf" \
		| cp "${temp}/${sample}_gatk_bwa_filt_snps.vcf" "${vcf}/${sample}_gatk_bwa_filt_snps.vcf"
		 
}

snps_pileup_novo () {
			echo "Initialising mpileup of ${sample} and novo" >> "${log}"

			samtools depth "${data}/${sample}_novo_sorted_dedup_realigned_sorted.bam" > "${temp}/${sample}_novo_sorted_dedup_realigned_sorted.coverage"
	             	depth_cov=$(cat "${temp}/${sample}_novo_sorted_dedup_realigned_sorted.coverage" |  awk '{sum+=$3} END { print sum/NR}')
        		depth_cov=${depth_cov%.*}
        		samtools mpileup -B -Q 20 -d $((depth_cov*2)) \
		        -C 50 \
		        -ugf "$ref" \
		        "${data}/${sample}_novo_sorted_dedup_realigned_sorted.bam" \
			2> "${data}/${sample}_novo_pileup.log" | bcftools call -vc \
		        -Ob > "${temp}/${sample}_novo.raw.bcf"

        		bcftools view "${temp}/${sample}_novo.raw.bcf" | vcfutils.pl varFilter \
		       -d 10 \
		       -D $depth_cov > "${temp}/${sample}_novo_raw_variants_Samtools_cons_caller.vcf"


			bgzip "${temp}/${sample}_novo_raw_variants_Samtools_cons_caller.vcf"
      			tabix "${temp}/${sample}_novo_raw_variants_Samtools_cons_caller.vcf.gz"

			bcftools filter -e "QUAL<30 && MIN(DP)<10" \
			-O z \
 			-o "${temp}/${sample}_novo_filt_pilepup.vcf.gz" "${temp}/${sample}_novo_raw_variants_Samtools_cons_caller.vcf.gz"

			cp "${temp}/${sample}_novo_filt_pilepup.vcf.gz" "${data}/${sample}_novo_filt_pilepup.vcf.gz" \
			| cp "${temp}/${sample}_novo_filt_pilepup.vcf.gz" "${snp_dir}/${sample}_novo_filt_pilepup.vcf.gz"
			cp "${temp}/${sample}_novo_filt_pilepup.vcf.gz" "${vcf}/${sample}_novo_filt_pilepup.vcf.gz"
}

snps_pileup_bwa () {	#NB I seem to have two versions installed, this needs to be checked when porting
			echo "Initialising mpileup of ${sample} and bwa" >> "${log}"

			samtools depth "${data}/${sample}_bwa_sorted_dedup_realigned_sorted.bam" > "${temp}/${sample}_bwa_sorted_dedup_realigned_sorted.coverage"

	             	depth_cov=$(cat "${temp}/${sample}_bwa_sorted_dedup_realigned_sorted.coverage" |  awk '{sum+=$3} END { print sum/NR}')
        		depth_cov=${depth_cov%.*}
        		samtools mpileup -B -Q 20 -d $((depth_cov*2)) \
		        -C 50 \
		        -ugf "$ref" \
		        "${data}/${sample}_bwa_sorted_dedup_realigned_sorted.bam" \
			2> "${data}/${sample}_bwa_pileup.log" | bcftools call -vc \
		        -Ob > "${temp}/${sample}_bwa.raw.bcf"

        		bcftools view "${temp}/${sample}_bwa.raw.bcf" | vcfutils.pl varFilter \
		       -d 10 \
		       -D $depth_cov > "${temp}/${sample}_bwa_raw_variants_Samtools_cons_caller.vcf"


			bgzip "${temp}/${sample}_bwa_raw_variants_Samtools_cons_caller.vcf"
      			tabix "${temp}/${sample}_bwa_raw_variants_Samtools_cons_caller.vcf.gz"

			bcftools filter -e "QUAL<30 && MIN(DP)<10" \
			-O z \
 			-o "${temp}/${sample}_bwa_filt_pilepup.vcf.gz" "${temp}/${sample}_bwa_raw_variants_Samtools_cons_caller.vcf.gz"

			cp "${temp}/${sample}_bwa_filt_pilepup.vcf.gz" "${data}/${sample}_bwa_filt_pilepup.vcf.gz" \
			| cp "${temp}/${sample}_bwa_filt_pilepup.vcf.gz" "${snp_dir}/${sample}_bwa_filt_pilepup.vcf.gz"
			cp "${temp}/${sample}_bwa_filt_pilepup.vcf.gz" "${vcf}/${sample}_bwa_filt_pilepup.vcf.gz"
}

intersect_snps_indel () {
			echo "Intersecting snp calls" >> "${log}"
			bgzip "${temp}/${sample}_gatk_novo_filt_snps.vcf"
			tabix "${temp}/${sample}_gatk_novo_filt_snps.vcf.gz"
			bgzip "${temp}/${sample}_gatk_bwa_filt_snps.vcf"
      		tabix "${temp}/${sample}_gatk_bwa_filt_snps.vcf.gz"
			tabix "${temp}/${sample}_novo_filt_pilepup.vcf.gz"
			tabix "${temp}/${sample}_bwa_filt_pilepup.vcf.gz"	

			bcftools isec -n 4 -O z \
		 	"${temp}/${sample}_gatk_novo_filt_snps.vcf.gz" \
		    	"${temp}/${sample}_gatk_bwa_filt_snps.vcf.gz" \
		 	"${temp}/${sample}_novo_filt_pilepup.vcf.gz" \
			"${temp}/${sample}_bwa_filt_pilepup.vcf.gz" \
			-o "${data}/${sample}_high_conf_snps.txt"

			echo "Intersecitng indels" >> "${log}"
		
			bgzip "${temp}/${sample}_gatk_bwa_filt_INDELs.vcf"
			tabix "${temp}/${sample}_gatk_bwa_filt_INDELs.vcf.gz"
			bgzip "${temp}/${sample}_gatk_novo_filt_INDELs.vcf"
			tabix "${temp}/${sample}_gatk_novo_filt_INDELs.vcf.gz"
			
			bcftools isec -n 2 -O z \
			"${temp}/${sample}_gatk_bwa_filt_INDELs.vcf.gz" \
			"${temp}/${sample}_gatk_novo_filt_INDELs.vcf.gz" \
			-o "${data}/${sample}_high_conf_INDEL.txt"

			#annotate indels with gene names
			cut -f2 "${data}/${sample}_high_conf_INDEL.txt" > "${temp}/${sample}_high_conf_INDEL_temp1.txt" 
			awk 'NR == FNR { x[$1] = $1+0; next; } { for (i in x) { if (x[i] > $1+0 && x[i] < $2+0) { print $3 "\t" $4 "\t" $5 } } }' "${temp}/${sample}_high_conf_INDEL_temp1.txt" "${master_dir}/databases/fusiondb_anno_sorted.txt" > "${temp}/${sample}_high_conf_INDEL_annotatedpos_temp2.txt"  

			paste "${data}/${sample}_high_conf_INDEL.txt" "${temp}/${sample}_high_conf_INDEL_annotatedpos_temp2.txt" > "${data}/${sample}_high_conf_INDEL_annotated.txt" 

			  cp "${data}/${sample}_high_conf_snps.txt" "${vcf}" \
			| cp "${data}/${sample}_high_conf_snps.txt" "${snp_dir}"
			  cp "${data}/${sample}_high_conf_INDEL.txt" "${vcf}" \
			| cp  "${data}/${sample}_high_conf_INDEL.txt" "${snp_dir}" \
			| cp  "${data}/${sample}_high_conf_INDEL_annotated.txt" "${vcf}" \
			| cp  "${data}/${sample}_high_conf_INDEL_annotated.txt" "${snp_dir}"
}   

#if we decide to draw tree we have to incorporate tree drawing software   

lumpy_bwa () {
		# Extract the split-read alignments
		echo "Initialising lumpy structural variant caller for bwa" >> "${log}"

  		samtools view -h "${data}/${sample}_bwa_sorted_dedup_realigned_sorted.bam" \
		| "${lumpy_extract_splitters}" -i stdin \
    		| samtools view -Sb - \
    		> "${temp}/${sample}_bwa_sample.splitters.unsorted.bam"  

		# Extract the discordant paired-end alignments.
  		samtools view -b -F 1294 "${data}/${sample}_bwa_sorted_dedup_realigned_sorted.bam" > "${temp}/${sample}_bwa_sample.discordants.unsorted.bam"

		lumpyexpress \
    		-B "${data}/${sample}_bwa_sorted_dedup_realigned_sorted.bam" \
    		-S "${temp}/${sample}_bwa_sample.splitters.unsorted.bam" \
    		-D "${temp}/${sample}_bwa_sample.discordants.unsorted.bam" \
    		-o "${temp}/${sample}_bwa_lumpy.vcf" 

		vcf-sort "${temp}/${sample}_bwa_lumpy.vcf" > "${temp}/${sample}_bwa_lumpy_sorted.vcf"

		#svtyper
  		bgzip "${temp}/${sample}_bwa_lumpy_sorted.vcf"
  		samtools index "${temp}/${sample}_bwa_sample.splitters.unsorted.bam"
  
  		zcat "${temp}/${sample}_bwa_lumpy_sorted.vcf" \
		| svtyper -B "${data}/${sample}_bwa_sorted_dedup_realigned_sorted.bam" \
		-S "${temp}/${sample}_bwa_sample.splitters.unsorted.bam" >  "${temp}/${sample}_bwa_lumpy_typed.vcf"
  		
  		vcf-sort "${temp}/${sample}_bwa_lumpy_typed.vcf" > "${temp}/${sample}_bwa_lumpy_typed_sorted.vcf"
		cp "${temp}/${sample}_bwa_lumpy_typed_sorted.vcf" "${data}/${sample}_bwa_lumpy_typed_sorted.vcf"
		"${lumpy_vcftobed}" \
		-i "${temp}/${sample}_bwa_lumpy_typed_sorted.vcf" \
		-o "${temp}/${sample}_bwa_lumpybed.txt"

		awk -F"\t" '{print $2 "\t" $6 "\t" $8 "\t" $11}' "${temp}/${sample}_bwa_lumpybed.txt" > "${temp}/${sample}_lumpy_bwa_temp1.txt" 
		
		#remove header lines
		sed -i '1d' "${temp}/${sample}_lumpy_bwa_temp1.txt"
		
		#filter for qual   
		awk '{if ($3 > 0) {print $1 "\t" $2 "\t" $3 "\t" $4}}' "${temp}/${sample}_lumpy_bwa_temp1.txt"  > "${temp}/${sample}_lumpy_bwa_temp2.txt"

		#add the sample names to not get confused also add S for the strain name will make life easy
		awk 'BEGIN{OFS="\t"}{print $0, "'"${sample}"'"}' "${temp}/${sample}_lumpy_bwa_temp2.txt" > "${temp}/${sample}_lumpy_bwa_temp3.txt"

		awk '{if ($1 > 0 && $2 > 0) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' "${temp}/${sample}_lumpy_bwa_temp3.txt"  > "${temp}/${sample}_lumpy_bwa_temp4.txt"

		#now we get the lengths
		awk 'BEGIN { OFS = "\t" } { $6 = $2 - $1 } 1' "${temp}/${sample}_lumpy_bwa_temp4.txt" > "${temp}/${sample}_lumpy_bwa_temp5.txt"

		#add the caller
		awk 'BEGIN{OFS="\t"}{print $0, "Lumpy_bwa"}' "${temp}/${sample}_lumpy_bwa_temp5.txt" > "${temp}/${sample}_lumpy_bwa_temp6.txt"
		
		#filter for length
		awk '{if ($6 < 50000 ) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' "${temp}/${sample}_lumpy_bwa_temp6.txt" > "${temp}/${sample}_lumpy_bwa_temp7.txt" 

		#deduplicate
		awk '!a[$1]++ && !b[$2]++' "${temp}/${sample}_lumpy_bwa_temp7.txt" > "${temp}/${sample}_lumpy_bwa_temp8.txt"  
		#annotate
		cut -f1 "${temp}/${sample}_lumpy_bwa_temp8.txt" > "${temp}/${sample}_lumpy_bwa_tempstart.txt"
		cut -f2 "${temp}/${sample}_lumpy_bwa_temp8.txt" > "${temp}/${sample}_lumpy_bwa_tempend.txt"
		
		awk 'NR == FNR { x[$1] = $1+0; next; } { for (i in x) { if (x[i] > $1+0 && x[i] < $2+0) { print x[i] "\t" $3 "\t" $4 "\t" $5 } } }' "${temp}/${sample}_lumpy_bwa_tempstart.txt" "${master_dir}/databases/fusiondb_anno_sorted.txt" > "${temp}/${sample}_lumpy_bwa_anno_tempstart.txt"

		awk 'NR == FNR { x[$1] = $1+0; next; } { for (i in x) { if (x[i] > $1+0 && x[i] < $2+0) { print x[i] "\t" $3 "\t" $4 "\t" $5 } } }' "${temp}/${sample}_lumpy_bwa_tempend.txt" "${master_dir}/databases/fusiondb_anno_sorted.txt" > "${temp}/${sample}_lumpy_bwa_anno_tempend.txt"

		awk '{print $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' "${temp}/${sample}_lumpy_bwa_temp8.txt" > "${temp}/${sample}_lumpy_bwa_temp9.txt" 

		paste "${temp}/${sample}_lumpy_bwa_anno_tempstart.txt" "${temp}/${sample}_lumpy_bwa_anno_tempend.txt" "${temp}/${sample}_lumpy_bwa_temp9.txt"  > "${data}/${sample}_bwa_novo.txt"
		#1start 2cds 3anno 4end 5cds 6anno 7qual 8sv 9strain 10length 11caller  
		 
                if [[ -s "${data}/${sample}_lumpy_bwa.txt" ]];
                    then
                    echo "Lumpy found structural variants from bwa algned files, now searching novoaligned files" >> "${log}"
		    cp "${data}/${sample}_lumpy_bwa.txt" "${sv_dir}" | cp "${data}/${sample}_lumpy_bwa.txt" "${SV}"
                        else 
                        echo "No variants could be found, moving on to novoaligned files" >> "${log}"
                fi
}	

lumpy_novo () {
		# Extract the split-read alignments
		echo "Initialising lumpy structural variant caller for novo" >> "${log}"

  		samtools view -h "${data}/${sample}_novo_sorted_dedup_realigned_sorted.bam" \
		| "${lumpy_extract_splitters}" -i stdin \
    		| samtools view -Sb - \
    		> "${temp}/${sample}_novo_sample.splitters.unsorted.bam"  

		# Extract the discordant paired-end alignments.
  		samtools view -b -F 1294 "${data}/${sample}_novo_sorted_dedup_realigned_sorted.bam" > "${temp}/${sample}_novo_sample.discordants.unsorted.bam"

		lumpyexpress \
    		-B "${data}/${sample}_novo_sorted_dedup_realigned_sorted.bam" \
    		-S "${temp}/${sample}_novo_sample.splitters.unsorted.bam" \
    		-D "${temp}/${sample}_novo_sample.discordants.unsorted.bam" \
    		-o "${temp}/${sample}_novo_lumpy.vcf" 

		vcf-sort "${temp}/${sample}_novo_lumpy.vcf" > "${temp}/${sample}_novo_lumpy_sorted.vcf"

		#svtyper
  		bgzip "${temp}/${sample}_novo_lumpy_sorted.vcf"
  		samtools index "${temp}/${sample}_novo_sample.splitters.unsorted.bam"
  
  		zcat "${temp}/${sample}_novo_lumpy_sorted.vcf" \
		| svtyper -B "${data}/${sample}_novo_sorted_dedup_realigned_sorted.bam" \
		-S "${temp}/${sample}_novo_sample.splitters.unsorted.bam" >  "${temp}/${sample}_novo_lumpy_typed.vcf"
  
  		vcf-sort "${temp}/${sample}_novo_lumpy_typed.vcf" > "${temp}/${sample}_novo_lumpy_typed_sorted.vcf"

		cp "${temp}/${sample}_novo_lumpy_typed_sorted.vcf" "${data}/${sample}_novo_lumpy_typed_sorted.vcf"
		
		"${lumpy_vcftobed}" \
		-i "${temp}/${sample}_novo_lumpy_typed_sorted.vcf" \
		-o "${temp}/${sample}_novo_lumpybed.txt"
		

		awk -F"\t" '{print $2 "\t" $6 "\t" $8 "\t" $11}' "${temp}/${sample}_novo_lumpybed.txt" > "${temp}/${sample}_lumpy_novo_temp1.txt" 

		# remove header
		#remove header lines
		sed -i '1d' "${temp}/${sample}_lumpy_novo_temp1.txt"  
		#filter for qual   
		awk '{if ($3 > 0) {print $1 "\t" $2 "\t" $3 "\t" $4}}' "${temp}/${sample}_lumpy_novo_temp1.txt"  > "${temp}/${sample}_lumpy_novo_temp2.txt"

		#add the sample names to not get confused also add S for the strain name will make life easy
		awk 'BEGIN{OFS="\t"}{print $0, "'"${sample}"'"}' "${temp}/${sample}_lumpy_novo_temp2.txt" > "${temp}/${sample}_lumpy_novo_temp3.txt"

		awk '{if ($1 > 0 && $2 > 0) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' "${temp}/${sample}_lumpy_novo_temp3.txt"  > "${temp}/${sample}_lumpy_novo_temp4.txt"

		#now we get the lengths
		awk 'BEGIN { OFS = "\t" } { $6 = $2 - $1 } 1' "${temp}/${sample}_lumpy_novo_temp4.txt" > "${temp}/${sample}_lumpy_novo_temp5.txt"
		
		#add caller name
		awk 'BEGIN{OFS="\t"}{print $0, "Lumpy_novo"}' "${temp}/${sample}_lumpy_novo_temp5.txt" > "${temp}/${sample}_lumpy_novo_temp6.txt"
		
		#filter for length
		awk '{if ($6 < 50000 ) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' "${temp}/${sample}_lumpy_novo_temp6.txt" > "${temp}/${sample}_lumpy_novo_temp7.txt" 

		#deduplicate
		awk '!a[$1]++ && !b[$2]++' "${temp}/${sample}_lumpy_novo_temp7.txt" > "${temp}/${sample}_lumpy_novo_temp8.txt"  

		#annotate
		cut -f1 "${temp}/${sample}_lumpy_novo_temp8.txt" > "${temp}/${sample}_lumpy_novo_tempstart.txt"
		cut -f2 "${temp}/${sample}_lumpy_novo_temp8.txt" > "${temp}/${sample}_lumpy_novo_tempend.txt"
		
		awk 'NR == FNR { x[$1] = $1+0; next; } { for (i in x) { if (x[i] > $1+0 && x[i] < $2+0) { print x[i] "\t" $3 "\t" $4 "\t" $5 } } }' "${temp}/${sample}_lumpy_novo_tempstart.txt" "${master_dir}/databases/fusiondb_anno_sorted.txt" > "${temp}/${sample}_lumpy_novo_anno_tempstart.txt"

		awk 'NR == FNR { x[$1] = $1+0; next; } { for (i in x) { if (x[i] > $1+0 && x[i] < $2+0) { print x[i] "\t" $3 "\t" $4 "\t" $5 } } }' "${temp}/${sample}_lumpy_novo_tempend.txt" "${master_dir}/databases/fusiondb_anno_sorted.txt" > "${temp}/${sample}_lumpy_novo_anno_tempend.txt"

		awk '{print $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' "${temp}/${sample}_lumpy_novo_temp8.txt" > "${temp}/${sample}_lumpy_novo_temp9.txt" 

		paste "${temp}/${sample}_lumpy_novo_anno_tempstart.txt" "${temp}/${sample}_lumpy_novo_anno_tempend.txt" "${temp}/${sample}_lumpy_novo_temp9.txt"  > "${data}/${sample}_lumpy_novo.txt"
		#1start 2cds 3anno 4end 5cds 6anno 7qual 8sv 9strain 10length 11caller   
                
                if [[ -s "${data}/${sample}_lumpy_novo.txt" ]]
                   then 
                   echo "Lumpy successfully called variants from novoaligned files, now intersecting" >> "${log}"
		   cp "${data}/${sample}_lumpy_novo.txt" "${sv_dir}" |  cp "${data}/${sample}_lumpy_novo.txt" "${SV}"
                    else
                    echo "Lumpy could not find vairants for novoalingned files" >> "${log}"
               fi 
                
}

delly_bwa () {
                #this is still untested, we need to check the format of the output
                #it may be better to annotate at the end of the final intersection
                echo "Scanning for split reads using delly V2" >> "${log}"
                export OMP_NUM_THREADS=$threads
                samtools view -f 2 -b "${data}/${sample}_bwa_sorted_dedup_realigned_sorted.bam" -o "${temp}/${sample}_bwa_sorted_dedup_realigned_sorted_delly.bam"
                samtools index "${temp}/${sample}_bwa_sorted_dedup_realigned_sorted_delly.bam"
                #NB export and put delly on path
                delly call -t DEL -q 10 -o "${data}/${sample}_delly.bcf" -g "${ref}" "${temp}/${sample}_bwa_sorted_dedup_realigned_sorted_delly.bam"
                #convert to vcf
                bcftools view "${data}/${sample}_delly.bcf" > "${data}/${sample}_delly.vcf"
		cp "${data}/${sample}_delly.vcf" "${sv_dir}/${sample}_delly.vcf"
		#processing
                #convert to bed
		svprops "${data}/${sample}_delly.vcf" > "${temp}/${sample}_delly_temp1.tab"
		#remove header
		sed -i '1d' "${temp}/${sample}_delly_temp1.tab"
		#get colums of interest
		awk '{ print $2 "\t" $4 "\t" $26 "\t" $11 }' "${temp}/${sample}_delly_temp1.tab" > "${temp}/${sample}_delly_temp2.tab" #start stop #qualvtype  #sype

		#get sample names
		awk 'BEGIN{OFS="\t"}{print $0, "'"${sample}"'"}' "${temp}/${sample}_delly_temp2.tab" > "${temp}/${sample}_delly_temp3.tab" #start stop #qualvtype  #sype #sample
  
		# filter nonsensical positions
		awk '{if ($1 > 0 && $2 > 0) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' "${temp}/${sample}_delly_temp3.tab"     > "${temp}/${sample}_delly_temp4.tab" #start stop #qualvtype  #sype #sample
 
		# get the lengths
		awk 'BEGIN { OFS = "\t" } { $6 = $2 - $1 } 1' "${temp}/${sample}_delly_temp4.tab" > "${temp}/${sample}_delly_temp5.tab" #start stop qualvtype  sype sample lengths

		#add the caller name
		awk 'BEGIN{OFS="\t"}{print $0, "Delly_bwa"}' "${temp}/${sample}_delly_temp5.tab" > "${temp}/${sample}_delly_temp6.tab" #1start 2stop 3qualvtype  4sype 5sample 6lengths 7caller

		awk '{if ($6 < 50000 ) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' "${temp}/${sample}_delly_temp6.tab" > "${temp}/${sample}_delly_temp7.tab" 

		#filter for quality NB the problem is duplicate lines!!!!!awk will deduplicate during annotation
		awk '{if ($3 > 0) {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7}}' "${temp}/${sample}_delly_temp7.tab" > "${temp}/${sample}_delly_temp8.tab" #1start 2stop 3qualvtype  4sype 5sample 6lengths 7caller

		#attemptdeduplicate otherwise gives problems
		awk '!a[$1]++ && !b[$2]++' "${temp}/${sample}_delly_temp8.tab" > "${temp}/${sample}_delly_temp9.tab" 
		
		#annotate #testing WARN changes made new print plus sort also awk column pull
########################################################################################################################
		cut -f1 "${temp}/${sample}_delly_temp9.tab" > "${temp}/${sample}_delly_tempstart.tab" 
		cut -f2 "${temp}/${sample}_delly_temp9.tab" > "${temp}/${sample}_delly_tempend.tab"

		awk 'NR == FNR { x[$1] = $1+0; next; } { for (i in x) { if (x[i] > $1+0 && x[i] < $2+0) { print x[i] "\t" $3 "\t" $4 "\t" $5 } } }' "${temp}/${sample}_delly_tempstart.tab" "${master_dir}/databases/fusiondb_anno_sorted.txt" > "${temp}/${sample}_delly_annotempstart.tab" 

		awk 'NR == FNR { x[$1] = $1+0; next; } { for (i in x) { if (x[i] > $1+0 && x[i] < $2+0) { print x[i] "\t" $3 "\t" $4 "\t" $5 } } }' "${temp}/${sample}_delly_tempend.tab" "${master_dir}/databases/fusiondb_anno_sorted.txt" > "${temp}/${sample}_delly_annotempend.tab" 
		
		awk '{print $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' "${temp}/${sample}_delly_temp9.tab" > "${temp}/${sample}_delly_temp10.tab" 
		paste "${temp}/${sample}_delly_annotempstart.tab" "${temp}/${sample}_delly_annotempend.tab" "${temp}/${sample}_delly_temp10.tab"  > "${data}/${sample}_delly_bwa_annotated.txt" 

################################################################################################################################3		

                if [[ -s "${data}/${sample}_delly_bwa_annotated.txt" ]];
                    then
                    echo "Delly has completed" >> "${log}"
		    cp "${data}/${sample}_delly_bwa_annotated.txt" "${sv_dir}/${sample}_delly_bwa_annotated.txt" 
                        else
                        echo "Delly has failed or could not find variants" >> "${log}"
                fi
		
		 
}    

lumpy_delly_isec () {	#intersect delly as well figure out format
		cat "${data}/${sample}_lumpy_bwa.txt" "${data}/${sample}_lumpy_novo.txt" "${data}/${sample}_delly_bwa_annotated.txt" > "${data}/${sample}_SV_isec_annotated.txt"
                
                if [[ -s "${data}/${sample}_SV_isec_annotated.txt" ]];
                    then
                    echo "Intersection of lumpy calls has completed successfully" >> "${log}"
		    cp "${data}/${sample}_SV_isec_annotated.txt" "${SV}" 
		    cp "${data}/${sample}_SV_isec_annotated.txt" "${sv_dir}"
		      else 
  		      echo "Intersection of lumpy files has failed" >> "${log}"
		fi
}



# this will be the genefusions section, i need to see delly output first though
#NB parse along with targeted to save space??
gene_fusion_calling () { #genefusions: lumpynovo,dellybwa,lumpybwa,indelsgatk
                        echo "Now searching for gene fusions" >> "${log}"
			#get files from SVisec and filter for CDS
			#now lets filter for CDS
			awk '{if ($3 == "CDS" && $7 == "CDS") print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 }' "${data}/${sample}_SV_isec_annotated.txt" > "${temp}/${sample}_SV_isec_annotated_fusions_temp1.txt"
                        
			#now get regions which are multigene deletions add headers at the end
			awk '{if ($2 != $6) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 }'  "${temp}/${sample}_SV_isec_annotated_fusions_temp1.txt" >  "${temp}/${sample}_SV_isec_annotated_fusions_temp2.txt"

			#remove transposons
			awk '{if ($4 != "insertion_seqs_and_phages" && $8 != "insertion_seqs_and_phages" ) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 }' "${temp}/${sample}_SV_isec_annotated_fusions_temp2.txt" > "${temp}/${sample}_SV_isec_annotated_fusions_temp3.txt" 

			#filter for deletions
			awk '{if ($10 == "DEL") print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 }' "${temp}/${sample}_SV_isec_annotated_fusions_temp3.txt"  > "${temp}/${sample}_SV_isec_annotated_fusions_temp4.txt"

			#deduplicate
			awk -F"\t" '{rows[$2]=$0} END {for (i in rows) print rows[i]}' "${temp}/${sample}_SV_isec_annotated_fusions_temp4.txt" > "${temp}/${sample}_SV_isec_annotated_fusions_temp5.txt"

			#sort 
			sort -k1 -n "${temp}/${sample}_SV_isec_annotated_fusions_temp5.txt" > "${chimera_dir}/${sample}_putative_chimeras.txt"
			

			#extract coordinates from of genes and get actual gene start/stop
			cut -f 2 "${chimera_dir}/${sample}_putative_chimeras.txt" > "${temp}/${sample}_fusions_anno_start.txt"
			cut -f 6 "${chimera_dir}/${sample}_putative_chimeras.txt" > "${temp}/${sample}_fusions_anno_stop.txt" 

			#extract the gene positions
			grep -f "${temp}/${sample}_fusions_anno_start.txt" "${master_dir}/databases/fusiondb_anno_sorted.txt" | awk '{print $1 "\t" $3}' > "${temp}/${sample}_fusions_gene_subset_start.txt" 

			grep -f "${temp}/${sample}_fusions_anno_stop.txt" "${master_dir}/databases/fusiondb_anno_sorted.txt" | awk '{print $1 "\t" $3}' > "${temp}/${sample}_fusions_gene_subset_stop.txt" 

  			paste "${temp}/${sample}_fusions_gene_subset_start.txt" "${temp}/${sample}_fusions_gene_subset_stop.txt" > "${temp}/${sample}_fusionloop.txt" 


			SOAP="${chimera_dir}/denovo"
			abacas="${chimera_dir}/contig_ordering"
			fusion_data="${chimera_dir}/data"
			mkdir "${SOAP}"
			mkdir "${abacas}"
			mkdir "${fusion_data}"

			while read reg_fusionstart reg_fusionannostart reg_fusionstop reg_fusionannostop 
				do
				echo "${sample}"
				echo "${reg_fusionstart}"
				echo "${reg_fusionannostart}"
				echo "${reg_fusionstop}"
				echo "${reg_fusionannostop}"

				echo "starting denovo assembly of ${reg_fusionannostart}-${reg_fusionannostop}" >> "${log}"
				
				#we need to make some dirs
				
				SOAP_targets="${SOAP}/soap_${reg_fusionannostart}-${reg_fusionannostop}"
				abacas_targets="${abacas}/abacas_${reg_fusionannostart}-${reg_fusionannostop}"
				mkdir "${SOAP_targets}"
				mkdir "${abacas_targets}"
				
				
				#extract region from bam files

				samtools view  "${data}/${sample}_novo_sorted_dedup_realigned_sorted.bam" ${fasta_header}:$(($reg_fusionstart-1000))-$(($reg_fusionstop+1000)) > "${temp}/${sample}_fusion_${reg_fusionannostart}-${reg_fusionannostop}.bam" -b -h

				samtools index "${temp}/${sample}_fusion_${reg_fusionannostart}-${reg_fusionannostop}.bam" 
				
				samtools faidx "${ref}" ${fasta_header}:$(($reg_fusionstart-1000))-$(($reg_fusionstop+1000)) > "${abacas_targets}/${sample}_ref_${reg_fusionannostart}-${reg_fusionannostop}.fasta"
				
				samtools faidx "${abacas_targets}/${sample}_ref_${reg_fusionannostart}-${reg_fusionannostop}.fasta" 

				#now lets de novo assemble
				echo "region files are ready, initiating SOAP" >> "${log}"
				echo "making config file" >> "${log}"

				#making config file
				{
           			 echo "max_rd_len=100"
            			 echo "[LIB]"
            			 echo "avg_ins=$mean"
            			 echo "reverse_seq=0"
            			 echo "asm_flags=3"
            			 echo "rank=1"
            			 echo "pair_num_cutoff=3"
            			 echo "map_len=32"
            			 echo "b=${temp}/${sample}_fusion_${reg_fusionannostart}-${reg_fusionannostop}.bam"
				} > "${SOAP_targets}/${sample}_${reg_fusionannostart}-${reg_fusionannostop}_SOAP_subset_config.txt"

				SOAPdenovo-63mer all \
	       	 		-s "${SOAP_targets}/${sample}_${reg_fusionannostart}-${reg_fusionannostop}_SOAP_subset_config.txt" \
		 		-o "${SOAP_targets}/${sample}_${reg_fusionannostart}-${reg_fusionannostop}_SOAP" \
		 		-K 31 \
		 		-F \
		 		-R \
		 		-N \
		 	 	1>"${SOAP_targets}/${sample}_ass.log" 2>"${SOAP_targets}/${sample}_ass.err"

				#order contigs
				abacas -r "${abacas_targets}/${sample}_ref_${reg_fusionannostart}-${reg_fusionannostop}.fasta" \
				       -q "${SOAP_targets}/${sample}_${reg_fusionannostart}-${reg_fusionannostop}_SOAP.scafSeq" \
		               -p nucmer \
		               -b  \
					   -d \
					   -a \
					   -m \
					   -N \
					   -o "${abacas_targets}/${sample}_fusion_${reg_fusionannostart}-${reg_fusionannostop}"

				cp "${temp}/${sample}_fusion_${reg_fusionannostart}-${reg_fusionannostop}.bam" "${fusion_data}"
		
				cp  "${temp}/${sample}_fusion_${reg_fusionannostart}-${reg_fusionannostop}.bam.bai" "${fusion_data}"

				samtools depth "${fusion_data}/${sample}_fusion_${reg_fusionannostart}-${reg_fusionannostop}.bam" > "${fusion_data}/${sample}_fusion_${reg_fusionannostart}-${reg_fusionannostop}_DoC.txt"

				samtools depth  "${fusion_data}/${sample}_fusion_${reg_fusionannostart}-${reg_fusionannostop}.bam" | cut -f 2,3 | gnuplot -e "set terminal png; set output '${fusion_data}/${sample}_fusion_${reg_fusionannostart}-${reg_fusionannostop}.png'; set title 'Coverage'; set xlabel 'position' ; set xtics rotate by 90 right ; set ylabel 'depth' ; plot '-'  using 2:xticlabels(1)  with lines notitle"

			done< <(tr -d '\r' < "${temp}/${sample}_fusionloop.txt")
}




target_region_process () {  #we could add a choice variable to combine fusion finding with targeted and give user a choice

			if [[ "${regions_of_interest}" == "FALSE" ]];
			then
				echo "Skipping targeted analysis" >> "${log}"
					else
                           	 		echo "analysing your target regions" >> "${log}"
                            	 		echo "Analysis is starting with the following params, ${sample}, ${name}, ${start}, ${end}" >> "${log}"
                            	 		#NB we have to extract the fasta header hopefully samtools will parse it, check this though
                           
                            	 		# NB check denovo dir test with ncbi referene 
                           	 		samtools view  "${data}/${sample}_novo_sorted_dedup_realigned_sorted.bam" ${fasta_header}:$(($start-2000))-$(($end+2000)) > "${regions}/${sample}_sorted_dedup_realigned_sorted_subset_${name}.bam" -b -h

                           			samtools index "${regions}/${sample}_sorted_dedup_realigned_sorted_subset_${name}.bam"
                           	 		samtools flagstat "${regions}/${sample}_sorted_dedup_realigned_sorted_subset_${name}.bam" > "${regions}/${sample}_${name}_subseted_bam_stats.txt"  
                           	 		#get subsetted reference file 
                           	 		samtools faidx "${ref}" ${fasta_header}:$(($start-2000))-$(($end+2000)) > "${regions}/${name}_ref_area_intrest.fasta"
                          	 		samtools faidx "${regions}/${name}_ref_area_intrest.fasta"
                           	    
                            			#draw a graph of coverage once i figure out how to call Rscript i'll port
                            			samtools depth  "${regions}/${sample}_sorted_dedup_realigned_sorted_subset_${name}.bam" | cut -f 2,3 | gnuplot -e "set terminal png; set output '${regions}/${sample}_${name}.png'; set title 'Coverage'; set xlabel 'position' ; set xtics rotate by 90 right ; set ylabel 'depth' ; plot '-'  using 2:xticlabels(1)  with lines notitle"
	                    
	                    			#get full per base coverage
	                    	 		samtools view -b "${regions}/${sample}_sorted_dedup_realigned_sorted_subset_${name}.bam" | \
	                    			genomeCoverageBed -bga -ibam stdin -g "${regions}/${name}_ref_area_intrest.fasta" > "${regions}/${sample}_subset_all_coverage.txt"
	                     
	                    			#get zero coverage coords
	                    			samtools view -b "${regions}/${sample}_sorted_dedup_realigned_sorted_subset_${name}.bam" | \
	                    			genomeCoverageBed -bga -ibam stdin -g "${regions}/${name}_ref_area_intrest.fasta" | \
	                    			awk '{if ($4 < 1) {printf ("%d\t%d\n",$3,$4);}}' > "${temp}/${sample}_${name}_zero_coverage_temp.txt" 
	                   			awk 'BEGIN { OFS = "\t" } { $3 = $2 - $1 } 1' "${temp}/${sample}_${name}_zero_coverage_temp.txt" > "${temp}/${sample}_${name}_zero_coverage_temp1.txt" 
	                    			awk 'BEGIN{OFS="\t"}{print $0, "'"${sample}"'"}' "${temp}/${sample}_${name}_zero_coverage_temp1.txt" > "${temp}/${sample}_${name}_zero_coverage_temp2.txt"                                                                                                               
	                    			awk 'BEGIN{OFS="\t"}{print $0, "'"${name}"'"}' "${temp}/${sample}_${name}_zero_coverage_temp2.txt" > "${regions}/${sample}_${name}_zero_coverage.txt" 
	                    
	                   			#get depth as well 
	                    			samtools depth  "${regions}/${sample}_sorted_dedup_realigned_sorted_subset_${name}.bam" | awk '{print $2 "\t" $3}' > "${temp}/${sample}_${name}_Doc_temp1.txt"
			    			awk 'BEGIN{OFS="\t"}{print $0, "'"${sample}"'"}' "${temp}/${sample}_${name}_Doc_temp1.txt" > "${temp}/${sample}_${name}_DoC_temp2.txt" 
			    			awk 'BEGIN{OFS="\t"}{print $0, "'"${name}"'"}' "${temp}/${sample}_${name}_DoC_temp2.txt" > "${regions}/${sample}_${name}_DoC.txt"

			    
			    			#add denovo assemblies
			    			#making config file
						{
           			 			echo "max_rd_len=100"
            			 			echo "[LIB]"
            			 			echo "avg_ins=$mean"
            			 			echo "reverse_seq=0"
            			 			echo "asm_flags=3"
            			 			echo "rank=1"
            			 			echo "pair_num_cutoff=3"
            			 			echo "map_len=32"
            			 			echo "b=${regions}/${sample}_sorted_dedup_realigned_sorted_subset_${name}.bam"
						} > "${target_soap}/${sample}_config_${name}.txt" 


			    			SOAPdenovo-63mer all \
	       	 				-s "${target_soap}/${sample}_config_${name}.txt" \
		 				-o "${target_soap}/${sample}_${name}_SOAP" \
		 				-K 31 \
		 				-F \
		 				-R \
		 				-N \
		 	 			1>"${target_soap}/${sample}_assembly.log" 2>"${target_soap}/${sample}_assembly.err"

		

			    			#order contigs
						abacas -r "${regions}/${name}_ref_area_intrest.fasta" \
				       		-q "${target_soap}/${sample}_${name}_SOAP.scafSeq" \
		                       		-p nucmer \
		                       		-b \
						-d \
						-a \
						-m \
						-N \
						-o "${target_ordering}/${sample}_ordered_contigs_${name}"

			fi
			    
	                    
}	                    
Tree_draw_files () {
					if [[ "${tree}" == "TRUE" ]];
					then
						echo "Initialising alternative reference for snps" >> "${log}"
						gunzip "${temp}/${sample}_gatk_novo_filt_snps.vcf.gz"
					
						java -Xmx"${ram}"g -jar "$gatk" \
							-T FastaAlternateReferenceMaker \
							-R "${ref}" \
							-o "${tree_dir}/${sample}_altref_tree.fasta" \
							-V "${temp}/${sample}_gatk_novo_filt_snps.vcf"
					fi

}
