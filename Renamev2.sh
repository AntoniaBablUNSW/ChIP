#!/bin/bash

# Usage: ./rename_files.sh path_to_directory path_to_excel.csv
# $1 - Directory  fastq
# $2 - Excel CSV file containing the identifier and name

# 1. Check if there are two arguments
if [ $# -ne 2 ]; then
    echo "Usage: $0 <directory_with_files> <excel_csv_file>"
    exit 1
fi

log_file="$1/../Log/Log.txt"
log() { echo $1
    echo "Md5 and Concat: $(date '+%Y-%m-%d %H:%M:%S') - $1" >> "$log_file"
}

# 2. Loop through the files in the specified directory
for file in $1/../aligned/*.sorted.rmdup.bam; do

    #echo "Processing files in directory: $directory"
    # 3. Extract the identifier from the filename
    # Assuming the identifier is the first part before an underscore (change as necessary)
    identifier=$(basename "$file" | cut -d'_' -f3)
   # echo "Extracted identifier: $identifier"
    # 4. Search for the corresponding name in the CSV file using the identifier
    # Assuming the identifier is in column 1 and name is in column 4
    name=$(awk -F, -v id="$identifier" '$1 == id {print $5}' "$2")
    
    #echo "Found name: $name"

    # 5. Check if a name was found
     if [ -z "$name" ]; then
        log "No matching name found for identifier: $identifier"
    else

        #counter=1
        new_name="$1/../aligned/${name}_${identifier}"
        #while [ -e "$new_name" ]; do
        #    counter=$((counter + 1))
        #    new_name="$1/../aligned/${identifier}_${name}_$counter"
        #done

    
        
        mv "$file" "$new_name.sorted.rmdup.bam"
        log "Renamed Bam $file to $new_name"


       original_basename=$(basename "$file")

        # Now update the CSV file to include the original basename in column 8
        # First, we read the file and replace column 8 with the original basename
        awk -F, -v id="$identifier" -v orig_basename="$original_basename" 'BEGIN{OFS=","} 
        $1 == id {$7 = orig_basename} {print $0}' "$2" > "${2}.tmp" && mv "${2}.tmp" "$2"

        log "Updated CSV with original basename in column 8 for $identifier"


    fi
done

log "File renaming in aligned completed."

# 2. Loop through the files in the specified directory
for file in $1/../BigWig/*.bw; do

    #echo "Processing files in directory: $directory"
    # 3. Extract the identifier from the filename
    # Assuming the identifier is the first part before an underscore (change as necessary)
    identifier=$(basename "$file" | cut -d'_' -f3)
   # echo "Extracted identifier: $identifier"
    # 4. Search for the corresponding name in the CSV file using the identifier
    # Assuming the identifier is in column 1 and name is in column 4
    name=$(awk -F, -v id="$identifier" '$1 == id {print $5}' "$2")
    
    #echo "Found name: $name"

    # 5. Check if a name was found
     if [ -z "$name" ]; then
        log "No matching name found for identifier: $identifier"
    else

        #counter=1
        new_name="$1/../BigWig/${name}_${identifier}"
        #while [ -e "$new_name" ]; do
        #    counter=$((counter + 1))
        #    new_name="$1/../BigWig/${name}_${identifier}_${counter}"
        #done

    
        
        mv "$file" "$new_name.bw"
        log "Renamed Bigwig $file to $new_name"


    fi
done
