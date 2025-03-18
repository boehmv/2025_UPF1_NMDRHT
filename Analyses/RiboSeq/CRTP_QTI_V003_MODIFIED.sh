#!/bin/bash

# Complete Ribo-TISH pipeline [CRTP] (Ribo-TISH PMID: 29170441) (Version 0.0.3)

# Usage:
# CRTP_XXXX_V003.sh -h

# Version 0.0.3 (10.11.22):
# - Use Gencode.v42.SIRVomeERCCome as reference in STAR step

# Version 0.0.2 (12.08.20):
# - Use bowtie2 for mapping to rRNA, tRNA, miRNA and snoRNA (as done in 'ORFquant'; https://github.com/lcalviell/ORFquant/issues/6)

####################################
#
# Perform multiple analyses steps: please set the parameters accordingly
#
####################################

# Change these global variables for each run

####################################
# Give absolute paths to the experiment file and to the project folder

# Give this a meaningful name (e.g. PSAP, NMD, ...)
myname="UPF1_RiboSeq_3rd"

# Give this a meaningful name (e.g. PSAP, NMD, ...) or use existing folders
srvdir="/srv/2024_08_UPF1_RiboSeq_3rd"

# Give this a meaningful name (e.g. PSAP, NMD, ...)
mydir="/home/volker/Riboseq/2024_08_UPF1_RiboSeq_3rd"

# This experiment.txt file has three columns: sample name, condition and method. IMPORTANT: methods can be: QTI, riboseq or RNAseq
experiment_file="${mydir}/experiment.txt"

# Please give this information to the best of your knowledge, if unsure check FastQC output files (MultiQC can help for an overview)
# https://dnatech.genomecenter.ucdavis.edu/wp-content/uploads/2019/03/illumina-adapter-sequences-2019-1000000002694-10.pdf
# TruSeq Ribo Profile - Adapter Trimming
adapter="AGATCGGAAGAGCACACGTCT"


####################################

# Define and check for options
if [ $# -eq 0 ]
then
	echo "Missing options!"
	echo "(run $0 -h for help)"
	echo ""
	exit 0
fi

ECHO="false"

while getopts "aTUFRMdcQPpDh" OPTION; do
	case $OPTION in

		a)
		ECHO="true"
		trim_enable="true"
		UMI_enable="true"
		Filter_Length_enable="true"
		rRNA_enable="true"
		STAR_enable="true"
		UMI_dedup_enable="true"
		count_enable="true"
		quality_enable="true"
		predict_riboseqOnly_enable="true"
		diff_enable="true"
		;;

		T)
		ECHO="true"
		trim_enable="true"
		;;
		
		U)
		ECHO="true"
		UMI_enable="true"
		;;
		
		F)
		ECHO="true"
		Filter_Length_enable="true"
		;;

		R)
		ECHO="true"
		rRNA_enable="true"
		;;

		M)
		ECHO="true"
		STAR_enable="true"
		;;
		
		d)
		ECHO="true"
		UMI_dedup_enable="true"
		;;
		
		c)
		ECHO="true"
		count_enable="true"
		;;

		Q)
		ECHO="true"
		quality_enable="true"
		;;

		P)
		ECHO="true"
		predict_enable="true"
		;;

		p)
		ECHO="true"
		predict_riboseqOnly_enable="true"
		;;

		D)
		ECHO="true"
		diff_enable="true"
		;;
		h)
		echo "Usage:"
		echo "CRTP_XXX_V003.sh [-OPTION] "
		echo ""
		echo "   -a     to execute the full analysis: trimming (-T), UMI extraction (-U), Filter length (-F), rRNA, tRNA, miRNA, snoRNA  removal (-R), STAR (-M), deduplication (-d), featureCounts (-c), Quality check (-Q), ORF prediction with riboseq only (-p), Differential analysis (-D)"
		echo "   -h     help (this output)"
		echo "   -T     to execute trimming with cutadapt"
		echo "   -U     to execute UMI extraction"
		echo "   -F     Filter length"
		echo "   -R     to execute rRNA, tRNA, miRNA, snoRNA  removal"
		echo "   -M     to execute STAR"
		echo "   -d     to execute deduplication"
		echo "   -c     to execute featureCounts"
		echo "   -Q     to execute Quality check"
		echo "   -P     to execute ORF prediction"
		echo "   -p     to execute ORF prediction with Riboseq data only"
		echo "   -D     to execute Differential analysis"
		exit 0
		;;

	esac
done

if [ $ECHO = "true" ]
then

	# Get date
	date=$(date "+%Y-%m-%d")

	#Define log file and redirect stdout and stderr to this file
	if [ ! -d "${mydir}/Logs/" ]; then
		mytime=$(date "+%Y-%m-%d %H:%M:%S")
		echo "$mytime Make directory ${mydir}/Logs/"
		mkdir ${mydir}/Logs/
	fi
	log_file="${mydir}/Logs/log_$date"
	exec &> >(tee -a "$log_file")

	# Give absolute paths to local folders
	bowtie_rRNA_dir="/home/volker/reference/Bowtie2/rRNA_tRNA_miRNA_snoRNA"
	STAR_dir="/home/volker/reference"
	signalp_dir="/home/volker/Tools/signalp-5.0b/bin"
	targetp_dir="/home/volker/Tools/targetp-2.0/bin"

	# State Script Version
	echo "###########################################"
	echo "## Complete Ribo-TISH pipeline [CRTP] Version 0.0.3"
	echo "###########################################"
	echo ""

	# Show paths, etc. in log
	echo "###################"
	echo "## General settings"
	echo "###################"
	echo ""
	echo "Name of the study: $myname"
	echo "Location of srv2 files: $srvdir"
	echo "Location of study files: $mydir"
	echo "Experiment file used for this analysis: $experiment_file"
	echo ""

	# Load experiment file, parse for samples and save them into the array
	let p=0
	while read -r f1 f2 f3; do
		samples[p]="${f1}"
		((++p))
	done < $experiment_file

	# Declare condition arrays
	declare -a cond

	# Load experiment file, parse for conditions and save unique conditions into array. Note: Condition1 is always "control"
	let i=1
	cond[0]="control"
	while read -r f1 f2 f3; do
		if [[ " ${cond[*]} " == *"$f2"* ]];then
			continue
		else
			cond[i]="${f2}"
			((++i))
		fi
	done < $experiment_file

	# Declare individual condition QTI and riboseq arrays
	arr_length="$((${#cond[@]}-1))"
	for i in $( eval echo {0..${arr_length}} );do
		declare -a cond${i}_QTI
		declare -a cond${i}_riboseq
		declare -a cond${i}_RNAseq
	done

	# Load experiment file again, parse for conditions and save filenames into condition-specific arrays.
	while read -r f1 f2 f3 f4; do
		for i in $( eval echo {0..${arr_length}} );do
			if [[ "$f2" == "${cond[i]}" ]];then
				if [[ "$f3" == "QTI" ]];then
					eval cond${i}_QTI[cond${i}_QTI_count]="${f1}"
					((++cond${i}_QTI_count))
				elif [[ "$f3" == "riboseq" ]];then
					eval cond${i}_riboseq[cond${i}_riboseq_count]="${f1}"
					((++cond${i}_riboseq_count))
				elif [[ "$f3" == "RNAseq" ]];then
					eval cond${i}_RNAseq[cond${i}_RNAseq_count]="${f1}"
					((++cond${i}_RNAseq_count))
				fi
			fi
		done
	done < $experiment_file

	# State the conditions and samples for this analysis
	echo "###################"
	echo "## Conditions and samples"
	echo "###################"
	echo ""

	arr_length="$((${#cond[@]}-1))"
	for i in $( eval echo {0..${arr_length}} );do
		echo -e "cond${i} \t ${cond[i]} \t QTI $(eval echo \${cond${i}_QTI[*]}) \t riboseq $(eval echo \${cond${i}_riboseq[*]}) \t RNAseq $(eval echo \${cond${i}_RNAseq[*]})"
	done


	# Generate experiment folder in srvdir
	if [ ! -d "${srvdir}/${myname}/" ]; then
		mytime=$(date "+%Y-%m-%d %H:%M:%S")
		echo "$mytime Make directory ${srvdir}/${myname}/"
		mkdir ${srvdir}/${myname}/
	fi

	# Generate srv BAM folder
	if [ ! -d "${srvdir}/${myname}/BAM/" ]; then
		mytime=$(date "+%Y-%m-%d %H:%M:%S")
		echo "$mytime Make directory ${srvdir}/${myname}/BAM/"
		mkdir ${srvdir}/${myname}/BAM/
	fi

	# Generate home BAM folder
	if [ ! -d "${mydir}/BAM/" ]; then
		mytime=$(date "+%Y-%m-%d %H:%M:%S")
		echo "$mytime Make directory ${mydir}/BAM/"
		mkdir ${mydir}/BAM/
	fi

	# Generate Tracks folder for easy visualization
	if [ ! -d "${mydir}/Tracks/" ]; then
		mytime=$(date "+%Y-%m-%d %H:%M:%S")
		echo "$mytime Make directory ${mydir}/Tracks/"
		mkdir ${mydir}/Tracks/
	fi

	# Return state of analysis options
	echo ""
	echo "###################"
	echo "## Analysis options"
	echo "###################"
	echo ""
	echo "trim_enable = ${trim_enable}"
	echo "UMI_enable = ${UMI_enable}"
	echo "Filter_Length_enable= ${Filter_Length_enable}"
	echo "rRNA_enable = ${rRNA_enable}"
	echo "STAR_enable = ${STAR_enable}"
	echo "UMI_dedup_enable = ${UMI_dedup_enable}"
	echo "quality_enable = ${quality_enable}"
	echo "predict_enable = ${predict_enable}"
	echo "predict_riboseqOnly_enable = ${predict_riboseqOnly_enable}"
	echo "diff_enable = ${diff_enable}"
	echo ""

	if  [ "$trim_enable" = "true" ]; then
		######
		######################
		#
		# Use cutadapt to trim
		#
		######################
		######

		echo ""
		echo "###################"
		echo "## Trim reads with cutadapt"
		echo "###################"
		echo ""

		source /home/volker/miniconda3/bin/activate cutadaptenv

		for i in $( eval echo {0..${arr_length}} );do
			for p in $( eval echo \${cond${i}_QTI[@]}); do
				if [ -e "${srvdir}/FASTQ/${p}_trimmed.fq.gz" ]; then
					mytime=$(date "+%Y-%m-%d %H:%M:%S")
					echo "$mytime ${p} FASTQ file already trimmed"
				else
					echo "Process ${p} with cutadapt"
					cutadapt -a $adapter --overlap=10 --minimum-length=0 --discard-untrimmed --compression-level=8 --cores 0 -o ${srvdir}/FASTQ/${p}_trimmed.fq.gz ${srvdir}/FASTQ/${p}.fq.gz > $mydir/FASTQ/${p}.cutadapt_stat.txt
					ln -s ${srvdir}/FASTQ/${p}_trimmed.fq.gz $mydir/FASTQ/${p}_trimmed.fq.gz
					echo "Process ${p} with cutadapt finished"
				fi
			done
			for p in $( eval echo \${cond${i}_riboseq[@]}); do
				if [ -e "${srvdir}/FASTQ/${p}_trimmed.fq.gz" ]; then
					mytime=$(date "+%Y-%m-%d %H:%M:%S")
					echo "$mytime ${p} FASTQ file already trimmed"
								
				else
					echo "Process ${p} with cutadapt"
					cutadapt -a $adapter --overlap=10 --minimum-length=0 --discard-untrimmed --compression-level=8 --cores 0 -o ${srvdir}/FASTQ/${p}_trimmed.fq.gz ${srvdir}/FASTQ/${p}.fq.gz > $mydir/FASTQ/${p}.cutadapt_stat.txt
					ln -s ${srvdir}/FASTQ/${p}_trimmed.fq.gz $mydir/FASTQ/${p}_trimmed.fq.gz
					echo "Process ${p} with cutadapt finished"
				fi
			done
		done
		conda deactivate
	fi	
	
	if  [ "$UMI_enable" = "true" ]; then
		######
		######################
		#
		# Use UMI_tools
		#
		######################
		######

		echo ""
		echo "###################"
		echo "## Extract UMIs with UMI-tools"
		echo "###################"
		echo ""

		source /home/volker/miniconda3/bin/activate umi_tools

		for i in $( eval echo {0..${arr_length}} );do
			for p in $( eval echo \${cond${i}_QTI[@]}); do
				if [ -e "${srvdir}/FASTQ/${p}_trimmed_UMI.fq.gz" ]; then
					mytime=$(date "+%Y-%m-%d %H:%M:%S")
					echo "$mytime ${p} FASTQ file already UMI extracted"
				else
					echo "Process ${p} with UMI-tools"
					umi_tools extract --extract-method=regex -p "^(?P<umi_1>.{4}).*(?P<umi_2>.{4})$" -I ${srvdir}/FASTQ/${p}_trimmed.fq.gz -S ${srvdir}/FASTQ/${p}_trimmed_UMI.fq.gz --log=${srvdir}/FASTQ/${p}_trimmed_UMI_processed.log 
					ln -s ${srvdir}/FASTQ/${p}_trimmed_UMI.fq.gz $mydir/FASTQ/${p}_trimmed_UMI.fq.gz
					echo "Process ${p} with UMI-tools finished"
				fi
			done
			for p in $( eval echo \${cond${i}_riboseq[@]}); do
				if [ -e "${srvdir}/FASTQ/${p}_trimmed_UMI.fq.gz" ]; then
					mytime=$(date "+%Y-%m-%d %H:%M:%S")
					echo "$mytime ${p} FASTQ file already UMI extracted"
								
				else
					echo "Process ${p} with UMI-tools"
					umi_tools extract --extract-method=regex -p "^(?P<umi_1>.{4}).+(?P<umi_2>.{4})$" -I ${srvdir}/FASTQ/${p}_trimmed.fq.gz -S ${srvdir}/FASTQ/${p}_trimmed_UMI.fq.gz --log=${srvdir}/FASTQ/${p}_trimmed_UMI_processed.log 
					ln -s ${srvdir}/FASTQ/${p}_trimmed_UMI.fq.gz $mydir/FASTQ/${p}_trimmed_UMI.fq.gz
					echo "Process ${p} with UMI-tools finished"
				fi
			done
		done
		conda deactivate
	fi
	
	if  [ "$Filter_Length_enable" = "true" ]; then
		######
		######################
		#
		# Filter length
		#
		######################
		######

		echo ""
		echo "###################"
		echo "## Filter length"
		echo "###################"
		echo ""

		for i in $( eval echo {0..${arr_length}} );do
			for p in $( eval echo \${cond${i}_QTI[@]}); do
				if [ -e "${srvdir}/FASTQ/${p}_trimmed_UMI_filtered.fq.gz" ]; then
					mytime=$(date "+%Y-%m-%d %H:%M:%S")
					echo "$mytime ${p} FASTQ file already filtered"
				else
					echo "Process ${p}"
					
					zcat ${srvdir}/FASTQ/${p}_trimmed_UMI.fq.gz | awk 'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 16 && length(seq) <= 40) {print header, seq, qheader, qseq}}' | gzip -8 > ${srvdir}/FASTQ/${p}_trimmed_UMI_filtered.fq.gz
										
					ln -s ${srvdir}/FASTQ/${p}_trimmed_UMI_filtered.fq.gz $mydir/FASTQ/${p}_trimmed_UMI_filtered.fq.gz
					echo "Process ${p} finished"
				fi
			done
			for p in $( eval echo \${cond${i}_riboseq[@]}); do
				if [ -e "${srvdir}/FASTQ/${p}_trimmed_UMI_filtered.fq.gz" ]; then
					mytime=$(date "+%Y-%m-%d %H:%M:%S")
					echo "$mytime ${p} FASTQ file already filtered"
								
				else
					echo "Process ${p}"
					zcat ${srvdir}/FASTQ/${p}_trimmed_UMI.fq.gz | awk 'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 16 && length(seq) <= 40) {print header, seq, qheader, qseq}}' | gzip -8 > ${srvdir}/FASTQ/${p}_trimmed_UMI_filtered.fq.gz
					
					ln -s ${srvdir}/FASTQ/${p}_trimmed_UMI_filtered.fq.gz $mydir/FASTQ/${p}_trimmed_UMI_filtered.fq.gz
					echo "Process ${p} finished"
				fi
			done
		done
	fi

	if  [ "$rRNA_enable" = "true" ]; then
		######
		######################
		#
		# Use bowtie2 to filter out the rRNA, tRNA, miRNA, snoRNA reads
		#
		######################
		######

		echo ""
		echo "###################"
		echo "## Perform mapping to rRNA, tRNA, miRNA, snoRNA with Bowtie2"
		echo "###################"
		echo ""

		source /home/volker/miniconda3/bin/activate bowtie2

		for i in $( eval echo {0..${arr_length}} );do
			for p in $( eval echo \${cond${i}_QTI[@]}); do
				if [ -e "${srvdir}/FASTQ/${p}_unaligned.fq.gz" ]; then
					mytime=$(date "+%Y-%m-%d %H:%M:%S")
					echo "$mytime ${p} FASTQ file already aligned vs rRNA, tRNA, miRNA, snoRNA"
				else
					echo "Process rRNA, tRNA, miRNA, snoRNA mapping ${p} with bowtie2"
					bowtie2 -p 15 --no-unal --local --un ${srvdir}/FASTQ/${p}_unaligned.fq -x $bowtie_rRNA_dir -U ${srvdir}/FASTQ/${p}_trimmed_UMI_filtered.fq.gz -S ${srvdir}/FASTQ/${p}_rRNA.align.sam
					rm ${srvdir}/FASTQ/${p}_rRNA.align.sam
					pigz -v ${srvdir}/FASTQ/${p}_unaligned.fq
					ln -s ${srvdir}/FASTQ/${p}_unaligned.fq.gz $mydir/FASTQ/${p}_unaligned.fq.gz
				fi
			done
			for p in $( eval echo \${cond${i}_riboseq[@]}); do
				if [ -e "${srvdir}/FASTQ/${p}_unaligned.fq.gz" ]; then
					mytime=$(date "+%Y-%m-%d %H:%M:%S")
					echo "$mytime ${p} FASTQ file already aligned vs rRNA, tRNA, miRNA, snoRNA"
				else
					echo "Process rRNA, tRNA, miRNA, snoRNA mapping ${p} with bowtie2"
					bowtie2 -p 15 --no-unal --local --un ${srvdir}/FASTQ/${p}_unaligned.fq -x $bowtie_rRNA_dir -U ${srvdir}/FASTQ/${p}_trimmed_UMI_filtered.fq.gz -S ${srvdir}/FASTQ/${p}_rRNA.align.sam
					rm ${srvdir}/FASTQ/${p}_rRNA.align.sam
					pigz -v ${srvdir}/FASTQ/${p}_unaligned.fq
					ln -s ${srvdir}/FASTQ/${p}_unaligned.fq.gz $mydir/FASTQ/${p}_unaligned.fq.gz
				fi
			done
		done
		
		conda deactivate
		
		# Execute multiQC
		
	# Generate FastQC folder
	if [ ! -d "${mydir}/FASTQ/FastQC" ]; then
		mytime=$(date "+%Y-%m-%d %H:%M:%S")
		echo "$mytime Make directory ${mydir}/FASTQ/FastQC"
		mkdir ${mydir}/FASTQ/FastQC
	fi
	
	# Generate QC folder
	if [ ! -d "${mydir}/QC" ]; then
		mytime=$(date "+%Y-%m-%d %H:%M:%S")
		echo "$mytime Make directory ${mydir}/QC"
		mkdir ${mydir}/QC
	fi
		
	
	/home/volker/Tools/FastQC/fastqc $mydir/FASTQ/*.fq.gz -o $mydir/FASTQ/FastQC/ -t 15
	
	source /home/volker/miniconda3/bin/activate activate multiqc
	
	multiqc --version

	multiqc ${mydir}/FASTQ -v -f -o ${mydir}/QC/ 
	
	conda deactivate
	
	fi

	if  [ "$STAR_enable" = "true" ]; then
		######
		######################
		#
		# Perform STAR mapping, indexing, on gencode.v42+SIRVomeERCCome annotation/genome
		#
		######################
		######

		echo ""
		echo "###################"
		echo "## Perform STAR mapping, indexing on gencode.v42+SIRVomeERCCome"
		echo "###################"
		echo ""
		echo "STAR directory: $STAR_dir"
		echo -n "STAR version: "
		STAR --version

		for i in "${samples[@]}"; do
			mytime=$(date "+%Y-%m-%d %H:%M:%S")
			echo "$mytime Start STAR mapping of sample ${i}"

			# Generate individual folder
			if [ ! -d "${srvdir}/${myname}/BAM/${i}/" ]; then
				mytime=$(date "+%Y-%m-%d %H:%M:%S")
				echo "$mytime Make directory ${srvdir}/${myname}/BAM/${i}/"
				mkdir ${srvdir}/${myname}/BAM/${i}/
			fi

			# Check if BAM file was already generated
			if [ -e "${srvdir}/${myname}/BAM/${i}/Aligned.sortedByCoord.out.bam" ]; then
				mytime=$(date "+%Y-%m-%d %H:%M:%S")
				echo "$mytime ${i} BAM file already generated"
			else
				STAR \
				--runThreadN 15 \
				--genomeDir $STAR_dir/gencode.v42.SIRVomeERCCome \
				--genomeLoad NoSharedMemory \
				--readFilesIn ${srvdir}/FASTQ/${i}_unaligned.fq.gz \
				--readFilesCommand zcat \
				--outFilterIntronMotifs RemoveNoncanonicalUnannotated \
				--outFilterMismatchNmax 2 \
				--alignEndsType EndToEnd \
				--alignIntronMax 200000 \
				--outMultimapperOrder Random \
				--outSAMmultNmax 1 \
				--outFileNamePrefix ${srvdir}/${myname}/BAM/${i}/ \
				--outSAMattributes All \
				--outSAMtype BAM SortedByCoordinate
			fi	
			
			# Index BAM file if necessary
			if [ -e "${srvdir}/${myname}/BAM/${i}/Aligned.sortedByCoord.out.bam.bai" ]; then
				mytime=$(date "+%Y-%m-%d %H:%M:%S")
				echo "$mytime ${i} BAM file already indexed"
			else
				mytime=$(date "+%Y-%m-%d %H:%M:%S")
				echo "$mytime Index ${srvdir}/${myname}/BAM/${i}"

				samtools index \
				${srvdir}/${myname}/BAM/${i}/Aligned.sortedByCoord.out.bam \
				${srvdir}/${myname}/BAM/${i}/Aligned.sortedByCoord.out.bam.bai
			fi
		done
	fi
	
	if  [ "$UMI_dedup_enable" = "true" ]; then
		######
		######################
		#
		# Use UMI_tools for deduplication
		#
		######################
		######

		echo ""
		echo "###################"
		echo "## Deduplication with UMI-tools"
		echo "###################"
		echo ""

		source /home/volker/miniconda3/bin/activate umi_tools

		for i in $( eval echo {0..${arr_length}} );do
			for p in $( eval echo \${cond${i}_QTI[@]}); do
				if [ -e "${srvdir}/${myname}/BAM/${p}/Aligned.sortedByCoord.out_dedup.bam" ]; then
					mytime=$(date "+%Y-%m-%d %H:%M:%S")
					echo "$mytime ${p} already deduplicated"
				else
					echo "Process ${p} with UMI-tools for deduplication"
					umi_tools dedup -I ${srvdir}/${myname}/BAM/${p}/Aligned.sortedByCoord.out.bam -S ${srvdir}/${myname}/BAM/${p}/Aligned.sortedByCoord.out_dedup.bam --log=${srvdir}/${myname}/BAM/${p}/${p}_dedup.log --extract-umi-method read_id --method unique
					ln -s ${srvdir}/${myname}/BAM/${p}/Aligned.sortedByCoord.out_dedup.bam $mydir/BAM/${p}.bam
					echo "Process ${p} with UMI-tools finished"
				fi
			done
			for p in $( eval echo \${cond${i}_riboseq[@]}); do
				if [ -e "${srvdir}/${myname}/BAM/${p}/Aligned.sortedByCoord.out_dedup.bam" ]; then
					mytime=$(date "+%Y-%m-%d %H:%M:%S")
					echo "$mytime ${p} already deduplicated"
				else
					echo "Process ${p} with UMI-tools for deduplication"
					umi_tools dedup -I ${srvdir}/${myname}/BAM/${p}/Aligned.sortedByCoord.out.bam -S ${srvdir}/${myname}/BAM/${p}/Aligned.sortedByCoord.out_dedup.bam --log=${srvdir}/${myname}/BAM/${p}/${p}_dedup.log --extract-umi-method read_id --method unique
					ln -s ${srvdir}/${myname}/BAM/${p}/Aligned.sortedByCoord.out_dedup.bam $mydir/BAM/${p}.bam
					echo "Process ${p} with UMI-tools finished"
				fi
			done
		done
		conda deactivate
		
		for i in "${samples[@]}"; do
		
			# Index BAM file if necessary
			if [ -e "${srvdir}/${myname}/BAM/${i}/Aligned.sortedByCoord.out_dedup.bam.bai" ]; then
				mytime=$(date "+%Y-%m-%d %H:%M:%S")
				echo "$mytime ${i} BAM file already indexed"
			else
				mytime=$(date "+%Y-%m-%d %H:%M:%S")
				echo "$mytime Index ${srvdir}/${myname}/BAM/${i}"

				samtools index \
				${srvdir}/${myname}/BAM/${i}/Aligned.sortedByCoord.out_dedup.bam \
				${srvdir}/${myname}/BAM/${i}/Aligned.sortedByCoord.out_dedup.bam.bai
			fi
		
			# Generate symbolic index link for easy access
			if [ -e "$mydir/BAM/${i}_dedup.bam.bai" ]; then
				mytime=$(date "+%Y-%m-%d %H:%M:%S")
				echo "$mytime ${i} BAM file index link already created"
			else
				ln -s ${srvdir}/${myname}/BAM/${i}/Aligned.sortedByCoord.out_dedup.bam.bai $mydir/BAM/${i}.bam.bai
			fi

			# Generate bedGraph using Alfred and convert to tdf using IGVtools
			if [ -e "$mydir/Tracks/${i}.tdf" ]; then
				mytime=$(date "+%Y-%m-%d %H:%M:%S")
				echo "$mytime ${i} BAM file already converted to track"
			else
				if [ ! -d "${mydir}/Tracks/BedGraph/" ]; then
					mytime=$(date "+%Y-%m-%d %H:%M:%S")
					echo "$mytime Make directory ${mydir}/Tracks/BedGraph/"
					mkdir ${mydir}/Tracks/BedGraph/
				fi

				# Use Alfred to generate bedGraphs
				alfred tracks -r 1 -o $mydir/Tracks/BedGraph/${i}.bedGraph.gz $mydir/BAM/${i}.bam
				cd /home/volker/Tools/IGV_2.14.1/

				# Convert to tdf
				/home/volker/Tools/IGV_2.14.1/igvtools totdf $mydir/Tracks/BedGraph/${i}.bedGraph.gz $mydir/Tracks/${i}.tdf gencode.v42.SIRVomeERCCome
				# Remove bedGraph files to save disk space
				rm $mydir/Tracks/BedGraph/${i}.bedGraph.gz
			fi
		done
	fi
	
	if  [ "$count_enable" = "true" ]; then
		######
		######################
		#
		# Perform featureCounts analysis
		#
		######################
		######

		echo ""
		echo "###################"
		echo "## Perform featureCounts on gencode.v42+SIRVomeERCCome"
		echo "###################"
		echo ""
		featureCounts -v
		
		# Generate counts folder
	if [ ! -d "${mydir}/counts/" ]; then
		mytime=$(date "+%Y-%m-%d %H:%M:%S")
		echo "$mytime Make directory ${mydir}/counts/"
		mkdir ${mydir}/counts/
	fi
		# Prepare BAM files
		printf -v BAM_samples "${mydir}/BAM/%s.bam \n" "${samples[@]}"
		echo $BAM_samples
		
			mytime=$(date "+%Y-%m-%d %H:%M:%S")
			echo "$mytime Start featureCounts for ${i}"


			# Check if featureCounts was already generated
			if [ -e "${mydir}/counts/${myname}_featureCounts.txt" ]; then
				mytime=$(date "+%Y-%m-%d %H:%M:%S")
				echo "$mytime ${i} featureCounts already generated"
			else
				featureCounts -T 10 -s 1 -M \
  				-a /home/volker/reference/Gencode/gencode.v42.SIRVomeERCCome.annotation.gtf \
  				-o ${mydir}/counts/${myname}_featureCounts.txt \
  				$BAM_samples
			fi
			
			
	fi

	if  [ "$quality_enable" = "true" ]; then
		######
		######################
		#
		# Perform ribo-TISH quality check
		#
		######################
		######

		echo ""
		echo "###################"
		echo "## Perform Ribo-TISH quality check"
		echo "###################"
		echo ""
		echo -n "Ribo-TISH version: "
		ribotish --version

		# GENCODE
		# Loop over conditions and then loop over files in QTI and riboseq arrays
		for i in $( eval echo {0..${arr_length}} );do
			for p in $( eval echo \${cond${i}_QTI[@]}); do
				if [ -e "${mydir}/BAM/${p}_qual.pdf" ]; then
					mytime=$(date "+%Y-%m-%d %H:%M:%S")
					echo "$mytime ${p} Quality check for GENCODE already performed"
				else
					echo "Process ${p} for GENCODE with ribotish-quality"
					ribotish quality -p 10 -v -b ${mydir}/BAM/${p}.bam -g /home/volker/reference/Gencode/gencode.v42.SIRVomeERCCome.annotation.gtf -t
				fi
			done
			for p in $( eval echo \${cond${i}_riboseq[@]}); do
				if [ -e "${mydir}/BAM/${p}_qual.pdf" ]; then
					mytime=$(date "+%Y-%m-%d %H:%M:%S")
					echo "$mytime ${p} Quality check for GENCODE already performed"
				else
					echo "Process ${p} for GENCODE with ribotish-quality"
					ribotish quality -p 10 -v -b ${mydir}/BAM/${p}.bam -g /home/volker/reference/Gencode/gencode.v42.SIRVomeERCCome.annotation.gtf
				fi
			done
		done
		
		# NMDRHT
		# Loop over conditions and then loop over files in QTI and riboseq arrays
		for i in $( eval echo {0..${arr_length}} );do
			for p in $( eval echo \${cond${i}_QTI[@]}); do
				if [ -e "${mydir}/BAM_NMDRegHumanTxome/${p}_qual.pdf" ]; then
					mytime=$(date "+%Y-%m-%d %H:%M:%S")
					echo "$mytime ${p} Quality check for NMDRHT already performed"
				else
					echo "Process ${p} for NMDRHT with ribotish-quality"
					ribotish quality -p 10 -v -b ${mydir}/BAM_NMDRegHumanTxome/${p}.bam -g /home/volker/2023_UPF1project/NMDRegHumanTxome/NMDRegHumanTxome.v1.1.sort.gtf -t
				fi
			done
			for p in $( eval echo \${cond${i}_riboseq[@]}); do
				if [ -e "${mydir}/BAM_NMDRegHumanTxome/${p}_qual.pdf" ]; then
					mytime=$(date "+%Y-%m-%d %H:%M:%S")
					echo "$mytime ${p} Quality check for NMDRHT already performed"
				else
					echo "Process ${p} for NMDRHT with ribotish-quality"
					ribotish quality -p 10 -v -b ${mydir}/BAM_NMDRegHumanTxome/${p}.bam -g /home/volker/2023_UPF1project/NMDRegHumanTxome/NMDRegHumanTxome.v1.1.sort.gtf
				fi
			done
		done
	fi

  if  [ "$predict_enable" = "true" ]; then
    ######
    ######################
    #
    # Perform ribo-TISH predict TIS
    #
    ######################
    ######

    echo ""
    echo "###################"
    echo "## Perform ribo-TISH predict TIS"
    echo "###################"
    echo ""
    echo -n "Ribo-TISH version: "
    ribotish --version

    # Generate predict folder
    if [ ! -d "${mydir}/Predict/" ]; then
      mytime=$(date "+%Y-%m-%d %H:%M:%S")
      echo "$mytime Make directory ${mydir}/Predict/"
      mkdir ${mydir}/Predict/
    fi
   
    cd ${mydir}/BAM/

    # Loop over conditions and then loop over files in QTI and riboseq arrays
    for i in $( eval echo {0..${arr_length}} );do
    
      # Generate Prediction file
      if [ -e "${mydir}/Predict/${cond[i]}.txt" ]; then
        mytime=$(date "+%Y-%m-%d %H:%M:%S")
        echo "$mytime ${cond[i]} Prediction already performed"
      else
        echo "Start ribotish predict on ${cond[i]}"
        QTI=$( eval echo \${cond${i}_QTI[@]})
        riboseq=$( eval echo \${cond${i}_riboseq[@]})
        ribotish predict -v --seq --alt -p 16 --aaseq -t ${mydir}/BAM/${QTI}.bam -b ${mydir}/BAM/${riboseq}.bam -g /home/volker/reference/Gencode/gencode.v42.SIRVomeERCCome.annotation.gtf -f /home/volker/reference/Gencode/GRCh38.primary.SIRVomeERCCome.fa -o ${mydir}/Predict/${cond[i]}.txt
      fi

    done
  fi

if  [ "$predict_riboseqOnly_enable" = "true" ]; then
    ######
    ######################
    #
    # Perform ribo-TISH predict TIS using only riboseq data
    #
    ######################
    ######

    echo ""
    echo "###################"
    echo "## Perform ribo-TISH predict TIS using only riboseq data"
    echo "###################"
    echo ""
    echo -n "Ribo-TISH version: "
    ribotish --version

    # Generate predict folder
    if [ ! -d "${mydir}/Predict/" ]; then
      mytime=$(date "+%Y-%m-%d %H:%M:%S")
      echo "$mytime Make directory ${mydir}/Predict/"
      mkdir ${mydir}/Predict/
    fi
    
        # Generate predict_NMDRegHumanTxome folder
    if [ ! -d "${mydir}/Predict_NMDRegHumanTxome/" ]; then
      mytime=$(date "+%Y-%m-%d %H:%M:%S")
      echo "$mytime Make directory ${mydir}/Predict_NMDRegHumanTxome/"
      mkdir ${mydir}/Predict_NMDRegHumanTxome/
    fi

    # Loop over conditions and then loop over files in riboseq arrays

    # Define filenames comma-separated
    for i in $( eval echo {0..${arr_length}} );do
	delim=""
	riboseq_samples=""
	for p in $( eval echo \${cond${i}_riboseq[@]}); do
		riboseq_samples="${riboseq_samples}${delim}${mydir}/BAM/${p}.bam"
 		delim=","
	done
    
	echo $riboseq_samples

      # Generate Prediction file for GENCODE
      if [ -e "${mydir}/Predict/${cond[i]}_riboseq_only.txt" ]; then
        mytime=$(date "+%Y-%m-%d %H:%M:%S")
        echo "$mytime ${cond[i]} Prediction for GENCODE already performed"
      else
        echo "Start ribotish predict for GENCODE on ${cond[i]}"
        ribotish predict -v --seq -p 16 --aaseq --inframecount -b ${riboseq_samples} -g /home/volker/reference/Gencode/gencode.v42.SIRVomeERCCome.annotation.gtf -f /home/volker/reference/Gencode/GRCh38.primary.SIRVomeERCCome.fa --framebest -o ${mydir}/Predict/${cond[i]}_riboseq_only.txt
      fi
      
            # Generate Prediction file for NMDRHT
      if [ -e "${mydir}/Predict_NMDRegHumanTxome/${cond[i]}_riboseq_only.txt" ]; then
        mytime=$(date "+%Y-%m-%d %H:%M:%S")
        echo "$mytime ${cond[i]} Prediction for NMDRHT already performed"
      else
        echo "Start ribotish predict for NMDRHT on ${cond[i]}"
        ribotish predict -v --seq -p 16 --aaseq --inframecount -b ${riboseq_samples} -g /home/volker/2023_UPF1project/NMDRegHumanTxome/NMDRegHumanTxome.v1.1.sort.gtf -f /home/volker/reference/Gencode/GRCh38.primary_assembly.genome.fa --framebest -o ${mydir}/Predict_NMDRegHumanTxome/${cond[i]}_riboseq_only.txt
      fi

    done
  fi

  if  [ "$diff_enable" = "true" ]; then
    ######
    ######################
    #
    # Perform ribo-TISH differential TIS analysis
    #
    ######################
    ######

    echo ""
    echo "###################"
    echo "## Perform ribo-TISH differential TIS analysis"
    echo "###################"
    echo ""
    echo -n "Ribo-TISH version: "
    ribotish --version

    # Generate DiffTIS folder
    if [ ! -d "${mydir}/DiffTIS/" ]; then
      mytime=$(date "+%Y-%m-%d %H:%M:%S")
      echo "$mytime Make directory ${mydir}/DiffTIS/"
      mkdir ${mydir}/DiffTIS/
    fi

    # Loop over conditions and then loop over files in QTI arrays
    for i in $( eval echo {1..${arr_length}} );do
      p=0
      if [ -e "${mydir}/DiffTIS/diffTIS_${cond[0]}_${cond[i]}.txt" ]; then
        mytime=$(date "+%Y-%m-%d %H:%M:%S")
        echo "$mytime ${cond[i]} differential TIS analysis was already performed"
      else
        echo "Start ribotish diff on ${cond[p]} vs ${cond[i]}"
        QTI_control=$( eval echo \${cond${p}_QTI[@]})
        QTI_sample=$( eval echo \${cond${i}_QTI[@]})
        ribotish tisdiff -1 ${mydir}/Predict/${cond[p]}.txt -2 ${mydir}/Predict/${cond[i]}.txt -a ${mydir}/BAM/${QTI_control}.bam -b ${mydir}/BAM/${QTI_sample}.bam -g /home/volker/reference/Gencode/gencode.v42.SIRVomeERCCome.annotation.gtf -o ${mydir}/DiffTIS/diffTIS_${cond[p]}_${cond[i]}.txt -p 12 -v --plotout ${mydir}/DiffTIS/diff_${cond[p]}_${cond[i]}.pdf
      fi
    done
  fi
fi
