#!/usr/bin/env python

# sequence_expand.py by DT Laura Vo Ngoc
# Part of the pipeline used to generate the EDIR database
# 
#*******************************************************************************
# Written initially in Python 3.6
#
# PURPOSE: Takes list of output from geneAnnotate_optimized.R, a 
#		FASTA-formatted genome.fa reference genome file, and a given maximum 
#		distance, and returns the repeat sequences, incrementally expanded in 
#		length to a max of 20bp repeat length. Only identical sequences within 
#		the given distance from each other are retained.
#
#
# Written to run from command line, parallized in SCOOP.
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python3 -m scoop ./sequence_expand.py list_of_input_files.txt distance_number
#-----------------------------------
#
#*******************************************************************************


# Libraries and dependencies
import sys
import os
import re
import pandas as pd
from Bio import SeqIO
from scoop import shared

# Define output directory
output_file_name_prefix = "../Output/expandSeq/"



def input_prep(
	infile):
	'''
	Takes input file (output from geneAnnotate_optimized.R), reorders the columns, adds "repeat_length"
	column, and removes duplicates.
	'''

	data_file = pd.read_csv(infile, sep="\t")
	data = data_file.iloc[:,[4,0,1,2,6,7,8,9,10,11,12]]

	data.insert(loc=1, column='repeat_length', value=['' for i in range(data.shape[0])])

	data.drop_duplicates(keep="first")

	return data


def sequence_window(
	data):
	'''
	Obtains "window" of input repeat sequence + 13 bases on either sides by referencing source.
	'''

	source = shared.getConst("source")
	distance = shared.getConst("distance")

	windows = [(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)]

	with open(source, "r") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			chrnum = record.id
			fasta_seq = record.seq		

			chrdata = data.loc[data['chromosome'] == chrnum]
			for index, row in chrdata.iterrows():
				if row['chromosome'] == str(chrnum):
					start = row['start']
					stop = row['end']
					row['repeat_seq'] = str(fasta_seq[row['end']-22:row['start']+20])
					windows.append(row)
				else:
					break

	windows.remove((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

	windows = pd.DataFrame(windows, columns=['chromosome', 'repeat_length', 'repeat_seq', 'start', 'end', 'dist', 'ensembl_gene_id', 'hgnc_symbol', 'gene_range', 'ensembl_transcript_id', 'transcript_range', 'intron_exon', 'number'])

	return windows



def length_offset(
	data):
	'''
	Expands the original sequence within the given sequence window, incrementing the length from
	the original 7bp to a maximum of 20bp.
	'''

	ori_start = 14
	ori_end = 20
	
	expanded = pd.DataFrame(columns=['chromosome', 'repeat_length', 'repeat_seq', 'start', 'end', 'dist', 'ensembl_gene_id', 'hgnc_symbol', 'gene_range', 'ensembl_transcript_id', 'transcript_range', 'intron_exon', 'number'])

	for offset in range(2, 15):
		temp = pd.DataFrame(columns=['chromosome', 'repeat_length', 'repeat_seq', 'start', 'end', 'dist', 'ensembl_gene_id', 'hgnc_symbol', 'gene_range', 'ensembl_transcript_id', 'transcript_range', 'intron_exon', 'number'])
		temp = temp.append(data)
		
		temp['repeat_seq'] = temp['repeat_seq'].str.slice(ori_start, ori_end + offset)
		temp['end'] = temp['end'] + offset - 1

		temp_min = pd.DataFrame(columns=['chromosome', 'repeat_length', 'repeat_seq', 'start', 'end', 'dist', 'ensembl_gene_id', 'hgnc_symbol', 'gene_range', 'ensembl_transcript_id', 'transcript_range', 'intron_exon', 'number'])
		temp_min = temp_min.append(data)
		temp_min['repeat_seq'] = temp_min['repeat_seq'].str.slice(ori_start-offset + 1, ori_end + 1)
		temp_min['start'] = temp_min['start'] - offset + 1
		
		temp = temp.append(temp_min)
		temp['repeat_length'] = 6 + offset
		
		expanded = expanded.append(temp)

	expanded = expanded.reset_index(drop=True)

	return expanded		



def results_filter(
	data):
	'''
	Filters results to retain only those within the given distance in bases.
	'''

	distance = int(shared.getConst("distance"))

	output = [(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)]

	data['start'] = pd.to_numeric(data['start'])
	data['chromosome'] = pd.Categorical(data['chromosome'], 
		categories = ('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','X','Y'), 
		ordered=True)
	
	data.sort_values(by=['repeat_length', 'chromosome', 'repeat_seq', 'start'], inplace = True)
	data = data.reset_index(drop=True)
	
	prev_line = pd.DataFrame()

	for index, line in data.iterrows():
		# if prev_line in not empty:
		if prev_line.empty == False:
			dist = line['start'] - prev_line['end']
			line['dist'] = dist
			if line['chromosome'] == prev_line['chromosome']:
				if line['repeat_seq'] == prev_line['repeat_seq']:
					if dist <= distance:
						if output[-1][3] != prev_line['start']:
							if (index == 1) or (abs(prev_line['dist']) > distance):
								prev_line['dist'] = dist
							output.append(prev_line)
						output.append(line)
			prev_line = line
			prev_line['dist'] = int(prev_line['dist'])
		else:
			prev_line = line

	output.remove((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
	output = pd.DataFrame(output, columns=['chromosome', 'repeat_length', 'repeat_seq', 'start', 'end', 'dist', 'ensembl_gene_id', 'hgnc_symbol', 'gene_range', 'ensembl_transcript_id', 'transcript_range', 'intron_exon', 'number'])

	return output


def repeats_length_expand(
	infile):
	'''
	Main function of the script.
	'''

	infile = "../Output/geneAnnotate/" + str(infile)
	distance = shared.getConst("distance")
	sequence = infile.split("_")[2]
	print(sequence)
	output = results_filter(length_offset(sequence_window(input_prep(infile))))
	outfile = "{prefix}expansion_results.{sequence}.{distance}bp.tsv".format(prefix=output_file_name_prefix, sequence=sequence,distance=distance)
	output.to_csv(outfile, sep='\t',index = False)



if __name__ == "__main__" and '__file__' in globals():

	
	import argparse
	from scoop import futures
	import time

	parser = argparse.ArgumentParser(prog='sequence_expansion.py',
		description="Takes list of output from geneAnnotate_optimzed.R, a \
		FASTA-formatted genome.fa reference genome file, and a given maximum \
		distance, and returns the repeat sequences, incrementally expanded in \
		length to a max of 20bp repeat length. Only identical sequences within \
		the given distance from each other are retained.")
	parser.add_argument("input_directory", help="file with content of the \
		directory containing the input files from the geneAnnotate_optimzed.R \
		output.", 
		metavar="INPUT_DIRECTORY")
	parser.add_argument("distance", help="Desired distance between results", 
		metavar="DISTANCE")
  
	# Trigger error if no arguments are provided
	if len(sys.argv)==1:	#from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
		parser.print_help()
		sys.exit(1)

	# Parse arguments
	args = parser.parse_args()
	file = open(args.input_directory)
	infiles = file.read().splitlines()
	
	# Set shared constants (source genome.fa file, distance)
	shared.setConst(source="../Genome_input/Homo_sapiens/genome.fa")
	shared.setConst(distance=args.distance)

	# Run in parallel using SCOOP
	list(futures.map(repeats_length_expand, infiles))


#********************************************************************
###-***********************END MAIN SCRIPT***********************-###
#********************************************************************