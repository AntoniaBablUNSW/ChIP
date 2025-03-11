#!/bin/sh

#CHeck permissions for everything and Programms exist and are named correctly

# $1 = input directory (DATA)
# $2 = ChIP overview Map
# $3 = Referenc genom (fasta)
# $4 =Input or IgG sample for reference in Peakcalling

if test $# -lt 3
then	
	log "Need 3 arguments, input dir to (DATA) and genome location (hg38) location (fasta file) and CVS ChIP overview file"
	log "Usage: $0 <input_directory> <reference_database>"
	exit 0
fi


source $1/../Programms/MD5.sh  $1

source $1/../Programms/QC1CutQC2.sh   $1/../fastq

source $1/../Programms/Align.sh   $1/../fastq $3

source $1/../Programms/Renamev2.sh   $1/../fastq $2

source $1/../Programms/Merge.sh   $1/../fastq

source $1/../Programms/Peakcall.sh   $1/../fastq $4

