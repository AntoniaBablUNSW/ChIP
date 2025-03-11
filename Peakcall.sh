#!/bin/bash

# Conda environment envPC

# Initialize conda (adjust the path to conda.sh as needed)
source /home/julie/anaconda3/etc/profile.d/conda.sh
conda activate envPC

# $1 = directory in aligned
# $2 = input for background files .bam IgG

mkdir -p "$1/../Log"
log_file="$1/../Log/Log.txt"
log() {
	echo "$1"
    echo "Peakcall: $(date '+%Y-%m-%d %H:%M:%S') - $1" >> "$log_file" 
}
log2() {
	
    echo "Peakcall: $(date '+%Y-%m-%d %H:%M:%S') - $1" >> "$log_file"
}
mkdir -p "$1/../Peaks"

# Log basic infos
log "Starting script Align"
log2 "User: $(whoami)"
log2 "Working Directory: $(pwd)"
log2 "System Info: $(uname -a)"
log2 "Active Conda Environment: $(conda info --envs | grep '*' | awk '{print $1}')"
log2 "Conda Environment Details: $(conda list)"

#LIST of Histone mods

HisMods=("H3" "H4" "K27me3" "K27ac" "K4me3" "k4me3" "k27ac")


for f in $1/../aligned/*.sorted.rmdup.bam; do
    
    
    
    basename="${f##*/}"  # Remove path
    basename="${basename%%.*}"  # Remove '.' suffix

 # Set variables
    CHIP_BAM="$1/../aligned/$basename.sorted.rmdup.bam"
    INPUT_BAM="$2"
    OUTPUT_DIR="$1/../Peaks"
    GENOME_SIZE="hs"  # Use "hs" for human, "mm" for mouse, or specify the number of base pairs
    NAME="$basename"



 # Check if basename starts with any histone modification
    is_histone_mod=0
    for mod in "${HisMods[@]}"; do
        if [[ "$basename" == "$mod"* ]]; then
            is_histone_mod=1
            break
        fi
    done

    if [[ $is_histone_mod -eq 1 ]]; then #Check if Histone_mod Was set to 1 

    log "Running MACS2 for broad peak calling on $basename..."
        # Run MACS2 for broad peak calling
    macs2 callpeak \
        -t "$CHIP_BAM" \
        -c "$INPUT_BAM" \
        -f BAMPE \
        -g "$GENOME_SIZE" \
        -n "$NAME" \
        --outdir "$OUTPUT_DIR" \
        -q 0.05 \
        --broad \
        -B \
        --keep-dup auto

        if [ -f "$OUTPUT_DIR/${NAME}_peaks.broadPeak" ]; then
            
         log " Broad peak calling successful for $NAME"
         
        else
            log " "
            echo "Error: Peak calling failed for $NAME."
            log " "
            fi




    
    else
 
    log "Running MACS2 for narrow peak calling on $basename..."

    # Run MACS2 for narrow peak calling
    macs2 callpeak \
        -t "$CHIP_BAM" \
        -c "$INPUT_BAM" \
        -f BAMPE \
        -g "$GENOME_SIZE" \
        -n "$NAME" \
        --outdir "$OUTPUT_DIR" \
        -q 0.05 \
        -B \
        --keep-dup auto\
        --call-summits

    
        if [ -f "$OUTPUT_DIR/${NAME}_peaks.narrowPeak" ]; then
            
         log "Narrow peak calling successful for $NAME"
         
        else
            log " "
            echo "Error: Peak calling failed for $NAME."
            log " "
            fi

    fi

done
log "End of Peak calling script"