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
log "Starting script Peakcallv2"
log2 "User: $(whoami)"
log2 "Working Directory: $(pwd)"
log2 "System Info: $(uname -a)"
log2 "Active Conda Environment: $(conda info --envs | grep '*' | awk '{print $1}')"
log2 "Conda Environment Details: $(conda list)"

#LIST of Histone mods
HisMods=("H3" "H4" "K27me3" "K27ac" "K4me3" "k4me3" "k27ac")

# Define prefix to q-value mapping
declare -A prefix_qvalue
prefix_qvalue=( 
    ["TAL1"]="0.05"
    ["RUNX1"]="0.01"
    ["GATA1"]="0.025"
    ["KLF1"]="0.05"
    ["BCL11A"]="0.01"
    # Add more prefix to q-value mappings as needed
)

# Default q-value if prefix not found in mapping
DEFAULT_QVALUE="0.05"

for f in $1/../aligned/*.sorted.rmdup.bam; do
    
    basename="${f##*/}"  # Remove path
    basename="${basename%%.*}"  # Remove '.' suffix

    # Extract the prefix (first part before underscore)
    prefix=$(echo "$basename" | cut -d'_' -f1)
    
    # Set q-value based on prefix
    if [[ -n "${prefix_qvalue[$prefix]}" ]]; then
        qvalue="${prefix_qvalue[$prefix]}"
        log "Using q-value of $qvalue for prefix $prefix"
    else
        qvalue="$DEFAULT_QVALUE"
        log "Using default q-value of $qvalue for prefix $prefix"
    fi

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

        log "Running MACS2 for broad peak calling on $basename with q-value $qvalue..."
        # Run MACS2 for broad peak calling
        macs2 callpeak \
            -t "$CHIP_BAM" \
            -c "$INPUT_BAM" \
            -f BAMPE \
            -g "$GENOME_SIZE" \
            -n "$NAME" \
            --outdir "$OUTPUT_DIR" \
            -q "$qvalue" \
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
        log "Running MACS2 for narrow peak calling on $basename with q-value $qvalue..."

        # Run MACS2 for narrow peak calling
        macs2 callpeak \
            -t "$CHIP_BAM" \
            -c "$INPUT_BAM" \
            -f BAMPE \
            -g "$GENOME_SIZE" \
            -n "$NAME" \
            --outdir "$OUTPUT_DIR" \
            -q "$qvalue" \
            -B \
            --keep-dup auto \
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