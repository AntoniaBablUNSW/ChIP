#!/bin/sh

# Initialize conda (adjust the path to conda.sh as needed)
source /home/julie/anaconda3/etc/profile.d/conda.sh
conda activate envCM
# $1 = input directory (fastq)
# $2 = DB (fasta)

mkdir -p "$1/../Log"

log_file="$1/../Log/Log.txt"
stat_file="$1/../Log/Stats.txt"

log() {
	echo "$1"
    echo "Align: $(date '+%Y-%m-%d %H:%M:%S') - $1" >> "$log_file"
}
log2() {
	
    echo "Align: $(date '+%Y-%m-%d %H:%M:%S') - $1" >> "$log_file"
}

stat_log() {
	
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" >> "$stat_file"
}

# Log basic infos
log "Starting script Align"
log2 "User: $(whoami)"
log2 "Working Directory: $(pwd)"
log2 "System Info: $(uname -a)"
log2 "Active Conda Environment: $(conda info --envs | grep '*' | awk '{print $1}')"
log2 "Conda Environment Details: $(conda list)"


if test $# -lt 2
then	
	log "Need 2 arguments, input dir and DB location"
	log "Usage: $0 <input_directory> <reference_database>"
	exit 0
fi

#Preparation

if [ ! -d "$1/../aligned" ]
then
	mkdir "$1/../aligned"
	log "Made directory aligned"
fi

if [ ! -d "$1/../BigWig" ]
then
	mkdir "$1/../BigWig"
	log "Made directory BigWig"
fi


#File processing


log "Start Aligning reads"

for f in $(ls "$1")
do

	source /home/julie/anaconda3/etc/profile.d/conda.sh
	conda activate envCM


	if [[ "$f" != "cut_"* ]]; then 

		continue  #go to the next iteration if it's not a cutadapt clipped file


	fi
	#Set up the variable from the file names (first) check if it's PE or SE then align appropriately

	
	


	if [[ "$f" == *"_1"* ]]; then

		
		basename="${f##*/}"  # Remove path
    	basename="${basename%%_1*}"  # Remove '_1*' suffix
		log "START Aligning paired end for $basename"

		#BWA alignment
	
		#bwa mem -t 40 $2 "$1/fastq/${basename}_1.fq.gz" "$1/fastq/${basename}_2.fq.gz" > "$1/aligned/$basename.sam"  #PE alignment original
		bwa mem -t 40 -R "@RG\tID:${basename}\tSM:${basename}\tPL:ILLUMINA" $2 $1/${basename}_1.* $1/${basename}_2.* > "$1/../aligned/$basename.sam" #PE alignment with chatGP help on read groups
		log " END BWA MEM for $basename"

		if [ ! -f "$1/../aligned/$basename.sam" ]; then
		    log  " "
    		log "Error: SAM File  not found!"
    		log " "
		else
    	    log " Succsefull BWA MEM for $basename"
		fi


	elif [[ "$f" == *"_2"* ]]; then

		log "PE data, ignoring read_2 of $basename as variable input"

		continue #break the loop and continue to the next iteration if the file is read 2 from a PE read

	else

		basename="${f##*/}"  # Remove path
    	basename="${basename%%_1*}"  # Remove '_1*' suffix
		log "START Aligning single end read $basename"
	

		#BWA alignment
	
		bwa mem -t 40 -R "@RG\tID:${basename}\tSM:${basename}\tPL:ILLUMINA"  $2 "$1/${basename}.fq.gz" > "$1/../aligned/basename.sam"  #SE alignment, check this works OK
		
		
		if [ ! -f "$1/../aligned/$basename.sam" ]; then
		    log  " "
    		log "Error: SAM File not found!"
    		log " "
		else
    	    log " Succsefull BWA MEM for $basename"
		fi

	fi



		
	#SAM to BAM
	
	log "SAM to BAM for $basename"
	samtools view -S -b "$1/../aligned/$basename.sam" > "$1/../aligned/$basename.bam"

	if [ -f "$1/../aligned/$basename.bam" ]; then
    # If the BAM file exists, delete the SAM file
    rm "$1/../aligned/$basename.sam"
	log " SAM file deleted for $basename"
	log " Successfull SAM to BAM for $basename"
	else
    # If the BAM file does not exist, print an error message

    log " \n Error: BAM file not created. SAM file not deleted.  \n"

	fi
			


	#bamtools stats
	
	stat_log "Bamtools stats for  $basename"
	stat_log "" 
	bamtools stats -in "$1/../aligned/$basename.bam" >> "$1/../Log/Stats.txt"
	
	# sort BAM
	
	log  "Sorting BAM  $basename"
	samtools sort "$1/../aligned/$basename.bam" -o "$1/../aligned/$basename.sorted.bam"

	if [ ! -f "$1/../aligned/$basename.sorted.bam" ]; then
		    log  " "
    		log "Error:  Sorted Bam File not found!"
    		log " "
		else
    	    log  " successfull Sorting BAM  $basename"
		fi
	
	#remove duplicates


	# This would garuntee unique mapping, not sure if necessary # carefull changes affix to sorted.uniq.bam
	#log "Filtering uniquely mapped reads for $basename"
	#samtools view -b -q 30 "$1/../aligned/$basename.sorted.bam" > "$1/../aligned/$basename.sorted.uniq.bam"
	#log "End filtering uniquely mapped reads for $basename"


	
	log "Removing duplicates from $basename"
	#picard MarkDuplicates I="$1/../aligned/$basename.sorted.bam" O="$1/../aligned/$basename.sorted.rmdup.bam" M="$basename.marked_dup_metrics.txt" REMOVE_DUPLICATES=TRUE #old syntax
	#picard MarkDuplicates -I"$1/aligned/$basename.sorted.bam" -O"$1/aligned/$basename.sorted.rmdup.bam" -M"$basename.marked_dup_metrics.txt" -REMOVE_DUPLICATES TRUE #future syntax
	

	# need to specify the directory to picard
	java -jar /home/julie/anaconda3/envs/envCM/share/picard-3.3.0-0/picard.jar MarkDuplicates \
	-I "$1/../aligned/$basename.sorted.bam" \
    -O "$1/../aligned/$basename.sorted.rmdup.bam" \
    -M "$1/../aligned/$basename.marked_dup_metrics.txt" \
    --REMOVE_DUPLICATES true

	if [ ! -f "$1/../aligned/$basename.sorted.rmdup.bam" ]; then
		    log  " "
    		log "Error:  Unduplicated Bam File not found!"
    		log " "
		else
    	    log  " successfull Removing duplicates BAM  $basename"
		fi

	#indexing
	
	log "Indexing BAM for $basename"
	samtools index "$1/../aligned/$basename.sorted.rmdup.bam" -@20
	log " End Indexing BAM for $basename"
	#Making bigWig
	
	log "Switching to envDT..."
	conda activate envDT



	
	log "Making BigWig for $basename"
	# Normalization for RPKM 
	bamCoverage -b "$1/../aligned/$basename.sorted.rmdup.bam" -o "$1/../BigWig/$basename.bw"  --effectiveGenomeSize 2913022398 --normalizeUsing RPKM  --binSize 10 --extendReads 200 -p20
	
	# Normalization more for  Reads per million
	#bamCoverage -b "$1/../aligned/$CUR.sorted.rmdup.bam" -o "$1/../BigWig/$CUR.bw" --normalizeUsing CPM  --binSize 10 --extendReads 200 -p20
	
	if [ ! -f "$1/../BigWig/$basename.bw" ]; then
		    log  " "
    		log "Error:  BigWig File not found!"
    		log " "
		else
    	    log  " successfull making BigWig for  $basename"
		fi




	
done

log "Switching to envDT..."
	conda activate envDT
	
	log "END of Alignment"
	log2 "Active Conda Environment: $(conda info --envs | grep '*' | awk '{print $1}')"
    log2 "Conda Environment Details: $(conda list)"

	 log "END of script"
	
	