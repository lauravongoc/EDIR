#!/usr/bin/env python


# Adapted by DT LAURA VO NGOC from original 
#   find_sequence_element_occurrences_in_sequence.py
#   author: Wayne Decatur
#   license: MIT
#   version: 0.1.0
#   source: https://github.com/fomightez/sequencework/tree/master/FindSequence




# find_sequence_element_occurrences_in_sequence.py
#
#*******************************************************************************
# Written initially in Python 3.
#
# PURPOSE: Takes a short sequence represented as a string and a FASTA-formatted
# sequence (from a file) and makes an accounting of the occurrences of that short
# sequence element on both strands of the main sequence. Exact matches only are 
# counted as occurrences. The script is case-insensitive, meaning the case of 
# either the element or FASTA sequence provided do not matter and matches will 
# be noted anyway. Only occurrences within a defined distance will be retained.
#
# As written, the script expects and only processes one FASTA-formatted 
# sequence. If your FASTA file has more than one sequence entry within it and 
# the one you want to have scanned is not the first, copy and paste it to a new 
# file and use that as the sequence to scan.
#
# Required arguments are a file containing a list of sequences, and an integer
# denoting maximun distance between occurrences.
# 
# Written to run from command line 
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python3 -m scoop find_seq_element_occurrences_in_seq.py sequences.txt 1000




#*******************************************************************************
##################################
#  USER ADJUSTABLE VALUES        #

##################################
#
output_file_name_prefix = "../Output/Homo_sapiens/"


limit_of_name = 17 # number of bases of the sequence element to limit to using 
# if the sequence element sequence is used to make the name for the output file


#
#*******************************************************************************
#**********************END USER ADJUSTABLE VARIABLES****************************













#*******************************************************************************
#*******************************************************************************
###DO NOT EDIT BELOW HERE - ENTER VALUES ABOVE###

import sys
import os
import re
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq # for reverse complement
from scoop import shared



###---------------------------HELPER FUNCTIONS---------------------------------###


def get_seq_element_representation(element):
    '''
    Takes `element` and returns a string representation for using
    in the output file name and the dataframe.
    '''
    if element:
        elem_id = element
    elif len(element) > limit_of_name:
        elem_id = element[:limit_of_name]+"..."
    else:
        elem_id = element
    return elem_id


def generate_output_file_name(element,id_of_seq_scanned,distance):
    '''
    Takes a file name as an argument and returns string for the name of the
    output file. The generated name is based on a prefix that can be adjusted
    under the section ' USER ADJUSTABLE VALUES ', plus the provided text in 
    the function call.
    If there is no `element` specified, the sequence element will be
    used, limited to the first number of bases specified in `limit_of_name`.
    Specific examples
    =================
    Calling function with
        ("elem1","chrmt")
    returns
        "seq_account_elem1_chrmt.tsv"
    Calling function with
        (None,"chrmt")
    returns
        "seq_account_GAATTC_chrmt.tsv"
    if `GAATTC` happened to be the provided sequence element.
    '''
    elem_id = element

    if elem_id.endswith("...") and (element == None):
        elem_id = elem_id[:-3] + "-" #add hyphen to indicate there is more 
    return "{prefix}{elem_id}_{seqid}.{distance}bp".format(
        prefix=output_file_name_prefix,elem_id=elem_id,
        seqid=id_of_seq_scanned,distance=distance)


def extract_id_of_seq_scanned(source):
    '''
    Take something like:
    https://downloads.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chrmt.fsa
    -or- 
    /local/directory1/directory2/chrmt.fa
    -or-
    chrmt.fa
    And return:
    chrmt
    '''
    if "/" in source:
        last_bit = source.split("/")[-1]
    else:
        last_bit = source
    if '.' in last_bit:
        main_part_of_name, file_extension = os.path.splitext(
        last_bit) #from http://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
        return main_part_of_name
    else:
        return last_bit


def get_seq_from_URL(url):
    '''
    takes a URL and gets the sequence
    '''
    try:
        from StringIO import StringIO
    except ImportError:
        from io import StringIO

    # Getting html
    try:
        # For Python 3.0 and later
        from urllib.request import urlopen
    except ImportError:
        # Fall back to Python 2's urllib2
        from urllib2 import urlopen
    html = urlopen(url)
    fasta_iterator = SeqIO.parse(StringIO(
        html.read().decode(encoding='UTF-8')), "fasta")

    record = next(fasta_iterator)
    return record.seq


def get_fasta_seq(source):
    '''
    Takes a source URL or filepath/ file name and gets sequence if it is in
    FASTA format.
    It won't return anything if what is provided isn't FASTA format
    '''

    # Read sequence, treating source as a filepath
    with open(source, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            # print(record.id) # for debugging
            return record.seq


def search_strand(pattern, sequence_to_scan, chromosome, distance, strand=1):
    '''
    Takes a sequence pattern (element) and find occurrences of that on the 
    provided, larger 5'-->3' sequence.
    Assumes strand is first unless provided.
    Tracks the start and end points of each occurrence, returning a list of
    that information where each element is a tuple of the start and end points
    along with the strand.
    Keeps sequences where positions are within "distance" from each other.
    '''

    # Initiate empty list of results and empty temporary variables
    occurrences = [(0, 0, 0, 0)]
    
    start1 = int(0)   # Temp variable that will hold previous match start position
    end1 = int(0)     # Temp variable that will hold previous match end position
    dist = int(distance)

    # For each part of sequence matching query, get start and end position
    for match in re.finditer(pattern.upper(), str(sequence_to_scan.upper())):   
        start_pos = match.start() + 1
        end_pos = match.end() + 1

        # If first match, initiate temp variable to its start position
        if start1 == 0:
            start1 = start_pos

        # If not first match, check if it is <= given distance from previous match
        elif start_pos - start1 <= dist:
            # If previous match is not already in results list, append it
            if start1 != occurrences[-1][1]:
                occurrences.append((chromosome, start1, end1, strand))
                if occurrences[-2][0] == 0:
                    occurrences.remove((0, 0, 0, 0))
            # Append current match to results list
            occurrences.append((chromosome, start_pos, end_pos, strand))
            
        # Set temp variables to current match
        start1 = start_pos
        end1 = end_pos
        
    return occurrences

###--------------------------END OF HELPER FUNCTIONS---------------------------###
###--------------------------END OF HELPER FUNCTIONS---------------------------###




#*******************************************************************************
###------------------------'main' function of script---------------------------##
def find_sequence_element_occurrences_in_sequence(
    element, return_dataframe = False):
    '''
    Main function of script. Scan a sequence (source) and report on occurrences 
    of a sub-sequence element in that sequence.
    Returns None
    Unless `return_dataframe = True`, and then it returns a dataframe of 
    accounting as well. That option being meant when using this script in a cell
    in a Jupyter notebook or importing it into another script.
    '''

    # Set shared SCOOP constants
    source = shared.getConst("source")
    distance = shared.getConst("distance")

    # Get the fasta_seq to scan
    fasta_seq = get_fasta_seq(source)
    

    assert fasta_seq, (
    "The provided source of the FASTA-formatted sequence seems invalid. Is it "
    "FASTA format?")

    #assert that element cannot be longer than fasta_seq
    assert len(element) < len(fasta_seq), (
    "the FASTA sequence has to be longer than the provided sequence element.\n"
    "The provided FASTA sequence, {0}, is {1} bases;\nthe provided sequence "
    "element '{2}' is {3} bases.".format(
        id_of_seq_scanned,len(fasta.seq),element,len(element)))

    totRes = 0
    with open(source, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            chrnum = record.id
            fasta_seq = record.seq

            occurrences = search_strand(element, fasta_seq, chrnum, distance) # first strand
            if occurrences:
                # make it into a dataframe since provides convenient options for 
                # handling
                df = pd.DataFrame(occurrences, columns=['chromosome','start_pos','end_pos','strand'])
                totRes+= len(df)

                # add some useful information to the dataframe
                df["element"] = element
                df = df[['element','chromosome','start_pos','end_pos','strand']] # 'seq.element'
                # column will show on far right otherwise

                # write to tab-delimited table file
                output_file_name = generate_output_file_name(
                    element,id_of_seq_scanned, distance)
                output_file_name = output_file_name + '.tsv'
               

                # If results file exists, append, else write to file
                resFile = os.path.isfile(output_file_name)
                if resFile:
                    with open(output_file_name, 'a') as output:
                        df.to_csv(output, sep='\t', index=False, header=False)
                else:
                    df.to_csv(output_file_name, sep='\t',index = False)
                    
                if return_dataframe:
                    print( "\n\nReturning a dataframe with the information "
                        "as well.")
                    return df
            else:
                # print( "\nNo occurrences of '{0}' "
                #         "found in the provided sequence.".format(element))
                if return_dataframe:
                    print( "\n\nNo data to return in a dataframe and so "
                        "returning `None`.")
                    return None


###--------------------------END OF MAIN FUNCTION----------------------------###
###--------------------------END OF MAIN FUNCTION----------------------------###







#*******************************************************************************
###------------------------'main' section of script---------------------------##
        

if __name__ == "__main__" and '__file__' in globals():
    """ This is executed when run from the command line """
       
    from scoop import futures

###-----------------for parsing command line arguments-----------------------###
    import argparse

    parser = argparse.ArgumentParser(prog='find_sequence_element_occurrences_in_sequence.py',
        description="find_sequence_element_occurrences_in_sequence.py takes \
        a short sequence represented as a string and a source FASTA-formatted \
        sequence (from a file) and makes an accounting of the occurrences of \
        that short sequence element (pattern) on both strands of the main \
        sequence. Exact matches only are counted as occurrences. Matching is \
        case-insensitive. Only occurrences within a defined distance will be \
        retained.")
    parser.add_argument("elements_file", help="File containing list of elements.", 
        metavar="ELEMENTS_FILE")
    parser.add_argument("distance", help="Desired distance between results", 
        metavar="DISTANCE")
    parser.add_argument('-id', '--element', action='store', type=str, 
        help="**OPTIONAL**Identifier \
        to reference sequence element. If none is provided, the sequence of \
        the user provided element will be used, limited to first {} bases if \
        exceeds that length.".format(limit_of_name))



    # Check that there are at least two arguments
    if len(sys.argv)==1:    #from http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    file = open(args.elements_file)
    elements = file.read().splitlines()


    # Set SCOOP shared variables
    shared.setConst(source="../Genome_input/Homo_sapiens/genome.fa") # Homo_sapiens GChR38
    shared.setConst(distance=args.distance)



    import time
    import resource

    start = time.time()
    
    # Parallelization function
    list(futures.map(find_sequence_element_occurrences_in_sequence, elements))
    end = time.time() 
    print("\n\nRunning time: {} sec".format(round(end-start, 2)))
    print("Memory used: {} kb".format(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))


#*******************************************************************************
###-***********************END MAIN PORTION OF SCRIPT***********************-###
#*******************************************************************************