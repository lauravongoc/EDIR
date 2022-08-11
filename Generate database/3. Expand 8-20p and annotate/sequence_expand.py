#!/usr/bin/env python

import sys
import os
import pandas as pd
from Bio import SeqIO
import multiprocessing as mp


# Initialize directory containing input data
# indir = "../Output/fuzz/geneAnnotate/"
indir = "/theia/scratch/projects/brightcore/edir/geneAnnotate/"

# Initialize output directory
# output_file_name_prefix = "../Output/fuzz/expandSeq/"
output_file_name_prefix = "/theia/scratch/projects/brightcore/edir/expandSeq/"


# Initialize list of final colnames
cols = ['chromosome', 'repeat_length', 'repeat_seq', 'start', 'end', 'distance', 'ensembl_gene_id', 'hgnc_symbol', 'gene_range', 'ensembl_transcript_id', 'transcript_range', 'intron_exon', 'number', 'feature', 'stringdist']

# Import input file, reorder columns
def input_prep(
	infile):
	
	data_file = pd.read_csv(infile, sep="\t", dtype={"chromosome": "str"}) 
	print(data_file)

	# Put stringdist last
	data = data_file.iloc[:,[0,1,2,3,4,5,6,7,8,9,10,11,12,14,13]]

	return data


# Get 33-base window centered around the 7bp seed sequence (to expand to 20bp on either side)
def sequence_window(
	data):

	# Initialize list of 0s for results
	windows = [(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)]

	with open(source, "r") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			chrnum = record.id            # chromosome number
			fasta_seq = record.seq        # sequence

			# Extract data which match the chromosome number
			chrdata = data.loc[data['chromosome'] == chrnum]

			# Get 33-base window centered on 7bp (+13 left and right to get 20bp sequences)
			for index, row in chrdata.iterrows():
				if row['chromosome'] == str(chrnum):
					row['repeat_seq'] = str(fasta_seq[row['end']-22:row['start']+20])
					windows.append(row)
				else:
					break
	# Remove initial row of 0s
	windows.remove((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

	# Turn into a df
	windows = pd.DataFrame(windows, columns=cols)

	return windows


# Progressively expand one base at a time from the 7bp seed sequence
# Only retain if the total mismatch number is 1
def length_offset(
	data):
	
	# 7 bp start/end coordinates in window range
	ori_start = 14
	ori_end = 20
	
	# Initialize empty df for results
	expanded = pd.DataFrame(columns=cols)

	# Initialize data copies for forward/reverse expansion
	dataF = data.copy()
	dataR = data.copy()


	# For different length offsets (+1 for Python's non-inclusive range)
	for offset in range(2, 15):
	# for offset in range(2, 3):
		print("length: ", offset + 6)

		# FORWARD EXP + EXACT
		# Subset dataF with no mismatch
		subset = dataF.loc[dataF['stringdist'] == 0]

		# If subset is not empty
		if not (subset.empty):
			temp = pd.DataFrame(columns=cols)
			temp = temp.append(subset)

			# Set repeat seq to one base forward
			temp['repeat_seq'] = temp['repeat_seq'].str.slice(ori_start, ori_end + offset)   # Expanded sequence
			temp['end'] = temp['end'] + offset - 1                                           # New 'end' value
			temp['distance'] = temp['distance'] - offset + 1                                 # Adjust distance
			temp['add'] = temp['repeat_seq'].str.strip().str[-1]                             # Column with added base
			temp['repeat_length'] = 6 + offset                                               # Repeat length column

			# Update stringdist variable to account for added base
			# Store first "add" value (added base)
			base = temp.iloc[0,15]

			# Iterate over full table
			for index, row in temp.iterrows():
				# print(row['add'])
				# If odd row (because first row is the header)
				if (index % 2) == 0:
					# Store added base in variable
					base = row['add']
				# Even row
				else:
					# If added base is different
					if row['add'] != base:
						# Change stringdist in temp and dataF
						dataF.at[index-1, 'stringdist'] = 1
						dataF.at[index, 'stringdist'] = 1
						temp.at[index-1, 'stringdist'] = 1
						temp.at[index, 'stringdist'] = 1
						# print("mismatch introduced")

			# Append results to output df
			expanded = expanded.append(temp.iloc[:,:-1])

		# BACKWARDS EXP + EXACT
		subset = dataR.loc[dataR['stringdist'] == 0]

		if not (subset.empty):
			temp = pd.DataFrame(columns=cols)
			temp = temp.append(subset)

			# Set repeat seq to one base back
			temp['repeat_seq'] = temp['repeat_seq'].str.slice(ori_start-offset + 1, ori_end + 1)
			temp['start'] = temp['start'] - offset + 1
			temp['distance'] = temp['distance'] - offset + 1
			temp['add'] = temp['repeat_seq'].str.strip().str[0]
			temp['repeat_length'] = 6 + offset

			# Update stringdist variable to account for added base
			base = temp.iloc[0,15]
			for index, row in temp.iterrows():
				# print(row['add'])
				# If odd row
				if (index % 2) == 0:
					# Store added base in variable
					base = row['add']
				# Even row, if added base is different
				else:
					if row['add'] != base:
						# Change stringdist
						dataR.at[index-1, 'stringdist'] = 1
						dataR.at[index, 'stringdist'] = 1
						temp.at[index-1, 'stringdist'] = 1
						temp.at[index, 'stringdist'] = 1
						# print("mismatch introduced")

			# Append results to output df
			expanded = expanded.append(temp.iloc[:,:-1])


		# FORWARD EXP + MISMATCH
		# For subset with existing mismatch
		subset = dataF.loc[dataF['stringdist'] == 1]

		# If there are rows in subset
		if not (subset.empty):
			# temp variable to hold expansion
			temp = pd.DataFrame(columns=cols)
			temp = temp.append(subset)

			# Set repeat seq to one base forward
			temp['repeat_seq'] = temp['repeat_seq'].str.slice(ori_start, ori_end + offset)
			temp['end'] = temp['end'] + offset - 1
			temp['distance'] = temp['distance'] - offset + 1
			temp['add'] = temp['repeat_seq'].str.strip().str[-1]
			temp['repeat_length'] = 6 + offset

			for index, row in temp.iterrows():
				# If odd row
				if (index % 2) == 0:
					# Store added base in variable
					base = row['add']
				# Even row, if added base is different
				else:
					if row['add'] != base:
						# Change stringdist
						dataF.at[index-1, 'stringdist'] = 2
						dataF.at[index, 'stringdist'] = 2
						temp.at[index-1, 'stringdist'] = 2
						temp.at[index, 'stringdist'] = 2
						# print("extra mismatch")
						
			temp = temp[(temp['stringdist'] != 2)]

			# Append results to output df
			expanded = expanded.append(temp.iloc[:,:-1])


		# BACKWARDS EXP + MISMATCH
		# For subset with existing mismatch
		subset = dataR.loc[dataR['stringdist'] == 1]

		# If there are rows in subset
		if not (subset.empty):
			# temp variable to hold expansion
			temp = pd.DataFrame(columns=cols)
			temp = temp.append(subset)

			# Set repaet seq to one base back
			temp['repeat_seq'] = temp['repeat_seq'].str.slice(ori_start-offset + 1, ori_end + 1)
			temp['start'] = temp['start'] - offset + 1
			temp['distance'] = temp['distance'] - offset + 1
			temp['add'] = temp['repeat_seq'].str.strip().str[0]
			temp['repeat_length'] = 6 + offset

			# Update stringdist variable to account for added base
			base = temp.iloc[0,15]
			for index, row in temp.iterrows():
				# If odd row
				if (index % 2) == 0:
					# Store added base in variable
					base = row['add']
				# Even row, if added base is different
				else:
					if row['add'] != base:	
						# Change stringdist
						dataR.at[index-1, 'stringdist'] = 2
						dataR.at[index, 'stringdist'] = 2
						temp.at[index-1, 'stringdist'] = 2
						temp.at[index, 'stringdist'] = 2
						# print("extra mismatch")
						
			temp = temp[(temp['stringdist'] != 2)]


			# Append results to output df
			expanded = expanded.append(temp.iloc[:,:-1])

	# Remove row numbers
	expanded = expanded.reset_index(drop=True)

	# Return output df
	return expanded



# Expand initial data and export results
def repeats_length_expand(
	infile):

	sequence = infile.split("_")[0]
	print(sequence)

	infile = str(indir) + str(infile)	

	# Initialize empty df for results
	output = pd.DataFrame(columns=cols)

	# Reformat input
	indata = input_prep(infile)

	# Expand to 20bp
	output = length_offset(sequence_window(indata))
	
	# Append expanded output to initial 7bp
	output = indata.append(output)

	# Export filename
	outfile = "{prefix}expansion_results.{sequence}.{distance}bp.tsv".format(prefix=output_file_name_prefix, sequence=sequence,distance=distance)
	
	# Export to file
	output.to_csv(outfile, sep='\t',index = False)




if __name__ == "__main__" and '__file__' in globals():

	
	import argparse
	import time

	parser = argparse.ArgumentParser(prog='sequence_expansion.py',
		description="sequence_expansion.py takes \
		a short sequence represented as a string and a source FASTA-formatted \
		sequence (either a file or URL to a file) and makes an accounting of \
		the occurrences of that short sequence element (pattern) on both \
		strands of the main sequence. Exact matches only are counted as \
		occurrences. Matching is case-insensitive.")

	parser.add_argument("input_directory", help="file with content of the \
		directory containing the input files from the geneAnnotate.R output.", 
		metavar="INPUT_DIRECTORY")
	parser.add_argument("distance", help="Desired distance between results", 
		metavar="DISTANCE")


	# Check there is 2 arguments

	if len(sys.argv)==1:	#from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
		parser.print_help()
		sys.exit(1)
	args = parser.parse_args()

	# Input directory
	file = open(args.input_directory)

	# Read file containing list of input files
	infiles = file.read().splitlines()
	
	# Define variables
	source = "../Genome_input/Homo_sapiens/genome.fa" # Homo_sapiens GChR38
	distance = args.distance

	# obtain number of cores allocated to the job
	ncore = len(os.sched_getaffinity(0))

	# Making a pool object and set number of processes equal to number of cores
	pool = mp.Pool(processes=ncore)
	result = pool.map(repeats_length_expand, infiles)

