#!/usr/bin/env python


# fuzzy_permutations.py
#
#*******************************************************************************
# Written initially in Python 3.
#
# PURPOSE: Takes a list of all possible permutations of 7bp sequences (16384 seqs)
# and retains only the first of sequences from which a list of sequences can be
# obtained using a maximum of 1 bp mismatch.
# 
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# python3 fuzz_permutations.py



def find_fuzzy(pattern, elements):
        
    # Find all 1-substitution fuzzy match
    for match in regex.findall("(" + str(pattern.upper()) + "){s<=1}", str(elements), overlapped=True):
        # Append current match to results list if it's not already there
        if "'" not in match:
            if match not in fuzzies:
                fuzzies.append(match)




if __name__ == "__main__" and '__file__' in globals():
    """ This is executed when run from the command line """
       

###-----------------for parsing command line arguments-----------------------###
    # import argparse
    import regex
 
    # Get file of 7bp permutations, and extract each line
    file = open("./seqFind/Input/all_sequences_length_7bp.txt")
    elements = [line.rstrip('\n') for line in file]

    # Possible fuzzies found
    fuzzies = []

    # "Unique" non-fuzzies
    nonfuzz = []

    for seq in elements:
        if seq not in fuzzies:
            nonfuzz.append(seq)
            find_fuzzy(seq, elements)

    # Print summary of results
    print("Fuzzies: ", len(fuzzies))
    print("Unique fuzz: ", len(list(set(fuzzies))))
    print(nonfuzz[0:20])
    print("Non fuzzies: ", len(nonfuzz))

    # Output results of minimum sequences to obtain full coverage
    with open('./seqFind/Input/fuzzy_sequences_7bp.txt', 'w') as f:
        for item in nonfuzz:
            f.write("%s\n" % item)