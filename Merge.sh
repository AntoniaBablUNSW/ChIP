#!/bin/sh

# Initialize conda (adjust the path to conda.sh as needed)
source /home/julie/anaconda3/etc/profile.d/conda.sh
conda activate envCM
# $1 = input directory (fastq)

mkdir -p "$1/../Log"
mkdir "$1/../Merged"
mkdir -p "$1/../BigWig"

log_file="$1/../Log/Log.txt"
stat_file="$1/../Log/Stats.txt"

log() {
	echo "$1"
    echo "Align: $(date '+%Y-%m-%d %H:%M:%S') - $1" >> "$log_file"
}
log2() {
	
    echo "Align: $(date '+%Y-%m-%d %H:%M:%S') - $1" >> "$log_file"
}

 log "Group by Prefix"

declare -A bam_groups

# Iterate over BAM files
for f in $1/../aligned/*.sorted.rmdup.bam; do
    # Extract the first two prefixes (assuming filenames are separated by underscores)
    prefix=$(basename "$f" | cut -d'_' -f1-2)

    # Append file to the corresponding group
    bam_groups["$prefix"]+="$f "
done

for prefix in "${!bam_groups[@]}"; do
    
    bam_list=(${bam_groups["$prefix"]})  # Convert space-separated string to an array
    num_files=${#bam_list[@]}  # Count number of BAM files

    if [ "$num_files" -eq 1 ]; then
        log "Skipping merge for $prefix (only one file: ${bam_list[0]})"
        continue  # Skip merging and move to the next prefix
    fi
    
    
    
    output_bam="$1/../Merged/${prefix}.merged.bam"


    log "Merging files for $prefix -> $output_bam"
    
    samtools merge -@ 4 "$output_bam" ${bam_groups["$prefix"]}

    if [ ! -f "$output_bam" ]; then
		    log  " "
    		log "Error: merged File  not found!"
    		log " "
		else
    	    log " Succsefull Merge for $prefix"
		fi
    

    log "Start sorting of $prefix"
    # Optional: Sort and index merged BAM
    samtools sort -@ 4 -o "${output_bam%.bam}.sorted.bam" "$output_bam"

    log "Start sindexing of $prefix"
    samtools index "${output_bam%.bam}.sorted.bam"


     if [ ! -f "$1/../Merged/${prefix}.merged.sorted.bam" ]; then
		    log  " "
    		log "Error: sorted merged File  not found!"
    		log " "
		else
    	    log " Succsefull Sorting for $prefix"

            # Delete unsorted merged BAM file
            rm "$output_bam"
            log "Deleted unsorted merged BAM: $output_bam"

		fi
    
    log "Finished merging $prefix"
done

log "All merging tasks completed."

    log "Switching to envDT..."
	conda activate envDT

	log2 "Active Conda Environment: $(conda info --envs | grep '*' | awk '{print $1}')"
    log2 "Conda Environment Details: $(conda list)"

for f  in $1/../Merged/*merged.sorted.bam ; do


    # Extract the basename of the BAM file (without path or extension)
    basename=$(basename "$f" .bam)  # Remove the .bam extension
    
    log "Making BigWig for $basename"
	# Normalization for RPKM 
	bamCoverage -b "$f" -o "$1/../BigWig/$basename.bw"  --effectiveGenomeSize 2913022398 --normalizeUsing RPKM  --binSize 10 --extendReads 200 -p20
	
	# Normalization more for  Reads per million
	#bamCoverage -b "$1/../aligned/$CUR.sorted.rmdup.bam" -o "$1/../BigWig/$CUR.bw" --normalizeUsing CPM  --binSize 10 --extendReads 200 -p20
	
	if [ ! -f "$1/../BigWig/$basename.bw" ]; then
		    log  " "
    		log "Error:  BigWig File not found!"
    		log " "
		else
    	    log  " Successfull making BigWig for  $f"
		fi

done


log "all merged files turned into BigWigs"
log "Script completed"
