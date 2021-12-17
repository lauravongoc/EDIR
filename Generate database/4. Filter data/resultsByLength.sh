#!/bin/bash

# resultsByLength.sh by DT Laura Vo Ngoc
# Part of the pipeline used to generate the EDIR database
# 
#*******************************************************************************
#
# PURPOSE: Takes output from geneAnnotate_features.R, and extracts the results
# 		per repeat length of 7 to 20bp.
#
#
# Written to run from command line, parallized in GNU parallel.
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# ./resultsByLength.sh
#-----------------------------------
#
#*******************************************************************************



# ***************************************
#               INITIALIZE
# ***************************************

# Output file path and define "dist" variable
export outprefix="../Output/Genes/All/"
export dist=1000

# Initialize parallel function parameters
PARALLEL="parallel -j -1 --bar "


# ***************************************
#               FUNCTION
# ***************************************
# Repeat length-specific output extracting function

lengths_extract() {
	# wc -l "$outfile"
	# Get gene name and format for exact grep
	length=$1

	# Initialize output and null_output files
	outfile="${outprefix}all_${length}bp_dist${dist}.tsv"
	touch "$outfile"
	chmod 775 "$outfile"

	# Get lines from directory files, and select results for repeat length "$length"
	## -v length="$length" -> define length variable in awk
	## $2 -> 2nd column (repeat_length)
	grep -r "$length" ../Output/expandSeq/annot | awk -F':' '{print $2}' | awk -v var="$length" 'int($2) == var' > "$outfile"
}
export -f lengths_extract


# ***************************************
#               MAIN SCRIPT
# ***************************************

# Combine 7bp results
awk FNR-1 ../Output/expandSeq/annot7/* > "${outprefix}all_7bp_dist${dist}.tsv"

# Run function in parallel
$PARALLEL lengths_extract ::: {8..20}