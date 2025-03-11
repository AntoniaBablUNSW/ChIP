# Checking for md5 number and conncatanate 
# Carefull checking for files to concatanate 1. is used as identifier
# if syntax different please adjust to uniqu difference between read 1 and read 2

    


  # $1 = input directory (must be top level directory) (eg AB_Sep25) Data

mkdir "$1/../Log"

log_file="$1/../Log/Log.txt"
log() { echo $1
    echo "Md5 and Concat: $(date '+%Y-%m-%d %H:%M:%S') - $1" >> "$log_file"
}

# Log basic infos
log "Starting script MD5"
log "User: $(whoami)"
log "Working Directory: $(pwd)"
log "System Info: $(uname -a)"
log "Active Conda Environment: $(conda info --envs | grep '*' | awk '{print $1}')"
log "Conda Environment Details: $(conda list)"


 # Function tho identify singled and pair end read (carefull context dependent)

# Function to count .gz files and rename
rename_gz_files() {
    # Find all .gz files
    gz_files=(*.gz)

     # Count the number of .gz files
     count=${#gz_files[@]}

    # Check the count and rename the files accordingly
    if [ $count -eq 1 ]; then
     mv "${gz_files[0]}" "s_${gz_files[0]}"
         log"Renamed to s_${gz_files[0]}"
         log "single end reads"
    elif [ $count -eq 2 ]; then
        mv "${gz_files[0]}" "p_${gz_files[0]}"
        mv "${gz_files[1]}" "p_${gz_files[1]}"
        log "Renamed to p_${gz_files[0]} and p_${gz_files[1]}"
        log "Paired end reads"
    elif [ $count -gt 2 ]; then
        log "Error: More than 2 .gz files found. Please ensure there are only 1 or 2 files."
        log "More than 2 files"
    else
        log "No .gz files found."
    fi
    } 



if test $# -lt 1  
then	
	log "Need 1 arguments, input dir (most likley Data folder)"
	exit 0
fi




mkdir $1/../fastq  # folder to store md5 checked and conncenated files 

for folder in $1/*/; do   # opens every folder in Data subsuccessively
    if [ -d "$folder" ]; then
        cd "$folder"

        log "Analysing Folder $folder"

        for file in *.gz; do  # Goes throug every ziped file
            if [ -f "$file" ]; then
             
            log  "Ckecking $file" 

            MD5=$(md5sum "$file" | awk '{print $1}')  # extraacts MD5 out of file

                if grep -q "$MD5" MD5.txt; then                 # searches for MD5 in Md5.txt file
        
                    log  "MD5 correct for $file"
                else
                    
                    log  " "
                    log "ERROR: MD5 not found for $file" >&2
                    log  " "

                fi
            fi

            log  "End Ckecking $file" 

        done

       


        # concatanate files

        # files 1. 
        
        # Step 1: Find files containing '1.'
        files=($(ls *1.*.gz 2>/dev/null | sort))


        # Step 2: Check if there are more than 2 matching files
        if [ ${#files[@]} -ge 2 ]; then

            log "$folder nedds to concatanate Read 1"
        
            # Step 5: Get the first filename alphabetically
            first_file="${files[0]}"  # First file in sorted list

            # Step 3: Create the '1' directory if it doesn't exist
            mkdir -p 1

            # Step 4: Move the matching files into the '1' directory
            mv "${files[@]}" 1/

            # Step 5: Concatenate and compress the files inside '1'
            log "Starts to concatanate $first_file"
            zcat 1/* | gzip >  $first_file

            #    Step 7: Print success message
         log "Files moved to '1/' and concatenated in the original folder."
            log " End to concatanate $first_file"
            fi


        # files containing .2



            # Step 1: Find files containing '1.'
        files=($(ls *2.*.gz 2>/dev/null | sort))



        # Step 2: Check if there are more than 2 matching files
        if [ ${#files[@]} -ge 2 ]; then
        

            log "$folder nedds to concatanate Read 2"

            # Step 5: Get the first filename alphabetically
            first_file="${files[0]}"  # First file in sorted list

            # Step 3: Create the '1' directory if it doesn't exist
            mkdir -p 2

            # Step 4: Move the matching files into the '1' directory
            mv "${files[@]}" 2/

            # Step 5: Concatenate and compress the files inside '1'
            log "Starts to concatanate $first_file"
            zcat 2/* | gzip >  $first_file

            #    Step 7: Print success message
            log "Files moved to '2/' and concatenated in the original folder."
             log " End to concatanate $first_file"
        fi


    
        rename_gz_files


    


        #move all .gz files to fastq folder

         mv [sp]*.gz "$1/../fastq"

         log "files moved to fastq folder"
   
    fi

done

log "Script md5 check and concatanating finished"

