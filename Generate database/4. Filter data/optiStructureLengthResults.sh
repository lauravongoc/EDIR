#!/bin/bash

# optiStructureLengthResults.sh by DT Laura Vo Ngoc
# Part of the pipeline used to generate the EDIR database
# 
#*******************************************************************************
#
# PURPOSE: Takes output from geneAnnotate_features.R, and extracts the results
# 		for each gene given in the list allGenesByChrom.txt
#
#
# Written to run from command line, parallized in GNU parallel.
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# ./structureLengthResults.sh
#-----------------------------------
#
#*******************************************************************************

# ***************************************
#               INITIALIZE
# ***************************************

# Input file path
export genelist="./geneExtract/Input/allGenesByChrom.txt"

# Set permissions
chmod 775 "$genelist"
chmod 775 ../Output/expandSeq/*
chmod 775 ../Output/expandSeq/annotFin/*

# Check that input is correctly formatted
dos2unix "$genelist"

# Output file paths and define variable
export outprefix="../Output/Genes_fixed/"
export dist=1000

# Initialize parallel function parameters
PARALLEL="parallel -j -1 --resume --joblog ./geneExtract/structlogFIX.log --bar "


# ***************************************
#               FUNCTION
# ***************************************
# Gene-specific output extracting function
lengths_gene_extract() {

	chromosome=$(echo "$1" | cut -f1 -d" ")
	gene=$(echo "$1" | cut -f2 -d" ")
	
	for length in {7..20};
	do
		# echo "$length"
		# Initialize output and null_output files
		outfile="${outprefix}${length}bp_dist${dist}.tsv"
		null_outfile="${outprefix}${length}bp_dist${dist}_null.txt"
		touch "$outfile"
		touch "$null_outfile"
		chmod 775 "$outfile"
		chmod 775 "$null_outfile"
		
		infile="../Output/expandSeq/annotFin/${length}bp_dist${dist}_chr${chromosome}.tsv"
		results=$(grep -rw "$gene" "$infile")


		# Get lines from directory files, and select results for repeat length "$length"
		## -v length="$length" -> define length variable in awk
		## $2 -> 2nd column (repeat_length)
		# results=$(echo "$grepres"| awk -v var="$length" 'int($2) == var') # rows

		# If $results is not empty
		if [ ! -z "$results" ]; then
			# Get total count of instances per gene
			tot_struct=$(($(echo "$results" | wc -l)/2))

			# Concatenate and collapse all sequence with number of rows/sequence
			## cut -f # -> cut column #
			## sed command -> remove whitespace
			## tr -> replace space with :
			each_seq=$(echo "$results" | cut -f 3 | sort | uniq -c | sed "s/^[ \s]*//" | tr " " : | paste -sd "/" -)

			# Add column to $results with distance between occurences | Get avg distance
			avgdist=$(echo "$results" | awk ' {total+=$6; count++} END {printf(total/count"\n")}')

			# Append to output file
			echo "$gene"$'\t'"$tot_struct"$'\t'"$avgdist"$'\t'"$each_seq" | tee -ai "$outfile" > /dev/null
			# echo "$gene"$'\t'"$tot_struct"$'\t'"$avgdist"$'\t'"$each_seq" >> "$outfile"

		# If $results is empty (no results)
		else
			echo "$gene" >> "$null_outfile"
		fi
	done
}
export -f lengths_gene_extract



# ***************************************
#               MAIN SCRIPT
# ***************************************

# Run function in parallel
$PARALLEL -a "$genelist" lengths_gene_extract



#********************************************************************
###-***********************END MAIN SCRIPT***********************-###
#********************************************************************