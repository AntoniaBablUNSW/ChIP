#!/bin/bash
# Doing the first QC and muliqc

#carefull here assumed ends of .fd.gz if need to change search and replace

# $1 Input directory to fastq

# Initialize conda (adjust the path to conda.sh as needed)
source /home/julie/anaconda3/etc/profile.d/conda.sh
conda activate envCM 
echo "Switching to envCM..."

shopt -s globstar

mkdir "$1/../Uncut"
mkdir "$1/../QCs"
mkdir "$1/../Log"

log_file="$1/../Log/Log.txt"
log() {
    echo  "$1"
    echo "Log: $(date '+%Y-%m-%d %H:%M:%S') - $1" >> "$log_file"
}

log2() {
    echo "Log: $(date '+%Y-%m-%d %H:%M:%S') - $1" >> "$log_file"
}



# Log basic infos
log "Starting script QC1"
log2 "User: $(whoami)"
log2 "Working Directory: $(pwd)"
log2 "System Info: $(uname -a)"
log2 "Active Conda Environment: $(conda info --envs | grep '*' | awk '{print $1}')"
log2 "Conda Environment Details: $(conda list)"


if test $# -lt 1  
then	
	echo "Need 1 arguments, input dir (most likley fastq folder)"
	exit 0
fi


if [ ! -d "$1" ]; then
    echo "Error: Directory $1 does not exist."
    exit 1
fi



for file in "$1"/*.fq.gz; do

 

    log "Starting fastQC1 of $file"
    fastqc -o $1/../QCs "$file"
    base=$(basename "$file" .fq.gz)
    mv "$1/../QCs/${base}_fastqc.zip" "$1/../QCs/QC1_${base}_fastqc.zip"
    mv "$1/../QCs/${base}_fastqc.html" "$1/../QCs/QC1_${base}_fastqc.html"
    if [ ! -f "$1/../QCs/QC1_${base}_fastqc.html"]; then
		    log  " "
    		log "Error:  QC1 File not found!"
    		log " "
		else
    	    log  " successfull making QC1 for  $base"
		fi


done

log "Starting multiQC1"
multiqc $1/../QCs -o $1/../QCs --filename QC1_multiqc
if [ ! -f "$1/../QCs/QC1_multiqc.html"]; then
		    log  " "
    		log "Error:  multiQC1 File not found!"
    		log " "
		else
    	    log  " successfull making MuliQC1"
		fi
log " End MultiQC1"


# CUTADAPT starts
log "Starting script cut adapt"
source /home/julie/anaconda3/etc/profile.d/conda.sh
conda activate envQC
echo "Switching to envQC..."



log2 "Active Conda Environment: $(conda info --envs | grep '*' | awk '{print $1}')"
log2 "Conda Environment Details: $(conda list)"

if test $# -lt 1  
then	
	echo "Need 1 arguments, input dir (most likley fastq folder)"
	exit 0
fi


if ! ls "$1" | grep -Eq '^(s_|p_)'; then

    log "Not all files sortet please check"

elif ls s_* &>/dev/null; then  # testing for single or paired reads, defined in MD5
   
    log "Single end read present. Please adjust Programm"

    # add here code for single end reads if necessary

else 
    
    log "All are pair end reads"


    for file in "$1"/*_1.*; do 
   



     if [[ -f "$file" && "$file" == "$1"/*_1.fq.gz ]]; then           # running cut adapt for _1 not clean samples
       
       # Extract the base name without _1 and the extension
        basename=$(basename "$file" "_1.fq.gz")
       log "Processing base name: $basename"
       
       log " ${basename}_1 is not clean"

        if [ -f "$1/${basename}_2.fq.gz" ]; then  # cheak if the partner has clean
            # Run cutadapt using the base name for paired-end reads

          log "${basename}_2 is not clean and Cutadapt started"  
            cutadapt -b CTGTCTCTTATACACATCT -B CTGTCTCTTATACACATCT \
         -o "$1/cut_${basename}_1.fq.gz" -p "$1/cut_${basename}_2.fq.gz" \
            "$1/${basename}_1.fq.gz" "$1/${basename}_2.fq.gz" -j 20
          log " End Cutadapt $basename"

            if [ ! -f "$1/cut_${basename}_1.fq.gz" ] || [ ! -f "$1/cut_${basename}_2.fq.gz" ]; then
		    log  " "
    		log "Error:  cut File not found!"
    		log " "
		    else
    	    log  " successfull cuting of $basename"
		    fi







        elif [ -f "$1/${basename}_2.clean.fq.gz" ]; then

        log " ${basename}_2 is clean and Cutadapt started"
            # Run cutadapt using the base name for paired-end reads
            cutadapt -b CTGTCTCTTATACACATCT -B CTGTCTCTTATACACATCT \
            -o "$1/cut_${basename}_1.fq.gz" -p "$1/cut_${basename}_2.clean.fq.gz" \
            "$1/${basename}_1.fq.gz" "$1/${basename}_2.clean.fq.gz" -j 20
            
            if [ ! -f "$1/cut_${basename}_1.fq.gz" ] || [ ! -f "$1/cut_${basename}_2.clean.fq.gz" ]; then
		    log  " "
    		log "Error:  cut File not found!"
    		log " "
		    else
    	    log  " successfull cuting of $basename"
		    fi

        else 
            log "No pair found for ${basename}"
        fi
    


     elif [[ -f "$file" && "$file" == "$1"/*_1.clean.fq.gz ]]; then

        # Extract the base name without _1 and the extension
        basename=$(basename "$file" "_1.clean.fq.gz")

       log "Processing base name: $basename"
       log "$ {basename}_1 is clean"
       if [ -f "$1/${basename}_2.fq.gz" ]; then  # cheak if the partner has clean
            # Run cutadapt using the base name for paired-end reads

          log " ${basename}_2 is not clean and Cutadapt started"  
            cutadapt -b CTGTCTCTTATACACATCT -B CTGTCTCTTATACACATCT \
         -o "$1/cut_${basename}_1.clean.fq.gz" -p "$1/cut_${basename}_2.fq.gz" \
            "$1/${basename}_1.clean.fq.gz" "$1/${basename}_2.fq.gz" -j 20
          
          if [ ! -f "$1/cut_${basename}_1.clean.fq.gz" ] || [ ! -f "$1/cut_${basename}_2.fq.gz" ]; then
		    log  " "
    		log "Error:  cut File not found!"
    		log " "
		    else
    	    log  " successfull cuting of $basename"
		    fi

        elif [ -f "$1/${basename}_2.clean.fq.gz" ]; then

        log "$ {basename}_2 is clean and Cutadapt started"
            # Run cutadapt using the base name for paired-end reads
            cutadapt -b CTGTCTCTTATACACATCT -B CTGTCTCTTATACACATCT \
            -o "$1/cut_${basename}_1.clean.fq.gz" -p "$1/cut_${basename}_2.clean.fq.gz" \
            "$1/${basename}_1.clean.fq.gz" "$1/${basename}_2.clean.fq.gz" -j 20
            
            if [ ! -f "$1/cut_${basename}_1.clean.fq.gz" ] || [ ! -f "$1/cut_${basename}_2.clean.fq.gz" ]; then
		    log  " "
    		log "Error:  cut File not found!"
    		log " "
		    else
    	    log  " successfull cuting of $basename"
		    fi

        else 
            log "No pair found for ${basename}"
        fi
     fi

    done
        
         # moves uncut files to folder Uncut in fastq
        for file in "$1"/*; do
        if [[ -f "$file" && "$file" != "$1/cut_"* ]]; then
        mv "$file" "$1/../Uncut"
        log "$file moved to uncut"
        fi
        done

    

       
fi

log "Script Cut adapt finished"


#Starting QC2

# Initialize conda (adjust the path to conda.sh as needed)
source /home/julie/anaconda3/etc/profile.d/conda.sh
echo "Switching to envCM..."
conda activate envCM 


log2 "Active Conda Environment: $(conda info --envs | grep '*' | awk '{print $1}')"
log2 "Conda Environment Details: $(conda list)"








# doing fast qc and safing Results in folder QCs with the prefix QC1
for file in $1/cut_* ; do
    log "Starting fastQC2 of $file"
    fastqc -o $1/../QCs "$file"
    base=$(basename "$file" .fq.gz)
    mv "$1/../QCs/${base}_fastqc.zip" "$1/../QCs/QC2_${base}_fastqc.zip"
    mv "$1/../QCs/${base}_fastqc.html" "$1/../QCs/QC2_${base}_fastqc.html"
    if [ ! -f "$1/../QCs/QC2_${base}_fastqc.html"]; then
		    log  " "
    		log "Error:  QC2 File not found!"
    		log " "
		else
    	    log  " successfull making QC2 for  $base"
		fi




done

log "Starting multiQC2"
multiqc $1/../QCs/QC2_* -o $1/../QCs --filename QC2_multiqc
log "Finished multiQC1"
if [ ! -f "$1/../QCs/QC2_multiqc.html"]; then
		    log  " "
    		log "Error:  multiQC1 File not found!"
    		log " "
		else
    	    log  " successfull making MuliQC1"
		fi


log "Script finished"
