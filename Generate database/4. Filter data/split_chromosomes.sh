#!/bin/bash

# split_chromosomes.sh by DT Laura Vo Ngoc
# Part of the pipeline used to generate the EDIR database
# 
#*******************************************************************************
#
# PURPOSE: Takes output from structureLengthResults.sh, splits them by chromosome
#
#
# Written to run from command line.
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# ./split_chromosomes.sh
#-----------------------------------
#
#*******************************************************************************

# Input and output file paths
infile="../Output/Genes/All"
outdir="../Output/Genes/All_split"

mkdir -p "$outdir"


for file in "$infile"/*; do
	echo "$file"
	
	# Get just the filename (without path or extension)
	filename=$(echo "$file" | cut -d'/' -f5 | cut -d'.' -f1)

	# Define output filename (path/filename)
	outfile="${outdir}/${filename}"

	# Split by chromosome and output to file
	awk -v outfile="$outfile"  'BEGIN{FS="\t";OFS="\t"}{print > outfile"_chr"$1".tsv"}' "$file"
done