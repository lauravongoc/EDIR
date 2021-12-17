#!/bin/bash

# ref_genome_loop.sh by DT Laura Vo Ngoc
# Part of the pipeline used to generate the EDIR database
# 
#*******************************************************************************
#
# PURPOSE: Unzip and format reference genome from Ensembl for downstream referencing
#
#
# Written to run from command line.
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# ./ref_genome_loop.sh
#-----------------------------------
#
#*******************************************************************************

echo "***Unzipping...***"
for filename in ../Genome_input/Homo_sapiens/*.fa.gz; 
do
	newdir="$(echo ${filename##*/})"
	new="$(echo "../Genome_input/Homo_sapiens/${newdir::len-3}")"
	echo "$new"
	gunzip -c -k "$filename" > "$new"
done

echo "***Removing \n characters...***"
for filename in ../Genome_input/Homo_sapiens/Homo_sapiens.GRCh38.dna.toplevel.fa.gz;
do
	# chrom="$(echo $filename | tr "." " " | awk '{ print $4 }')"
	# echo "$chrom"
	output="$(echo "../Genome_input/Homo_sapiens/Homo_sapiens.nonew.fa")"
	echo "$output"
	tr -d '\n' < "$filename" > "$output"
	awk '{printf "%s%s", (prev~/REF$/?RS:""), $0; prev=$0} END{print ""}' < "$filename" > "$output"
done

echo "***Ensuring FASTA format...***"
for file in ../Genome_input/Homo_sapiens/Homo_sapiens.nonew.fa;
do
	# chrom="$(echo $file | tr "." " " | awk '{ print $2}')"
	# echo "$chrom"
	finfasta="$(echo "../Genome_input/Homo_sapiens/Homo_sapiens.fa")"
	sed -e 's/REF/REF\n/g' "$file" > "$finfasta"
done


#******************************************************************
###-***********************END OF SCRIPT***********************-###
#******************************************************************