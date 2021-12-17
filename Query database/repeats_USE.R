#!/usr/bin/Rscript

# repeats_USE.R by DT Laura Vo Ngoc
# Part of the pipeline used to query the EDIR database
# 
#*******************************************************************************
# Written initially in R 3.5.0 Joy in Playing
#
# PURPOSE: Description of query function gene_lookup() and its parameters,
#       with examples of usage. 
#

# Written to run locally in RStudio.
#
#*******************************************************************************




# ********************************************
# ---------------- INITIALIZE ----------------
# ********************************************

# Set path of file location
path <- "C:/Users/rockp/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu18.04onWindows_79rhkp1fndgsc/LocalState/rootfs/home/laura/genRepeats/Genome-repeats/Repeats"

# Load required libraries and function
source(paste0(path, "/repeats_functions.R"))


# DESCRIPTION -------------------------------------------------------------

# ******************************
#---- FUNCTION: GENE LOOKUP ----
# ******************************

# Function to look up a specified gene in the dataset

# REQUIRED INPUT: 
#     gene ID (ENSEMBL or HGNC)
# OPTIONAL INPUT: 
#     length ---- repeat length (7-20), if not defined will output all results
#     mindist --- minimum distance, default at 0
#     maxdist --- maximum distance, default at 1000
#     summary --- whether to include summary in output, defaults to FALSE

# OUTPUT:
#     Prints to console:
#      - Gene ID (Ensembl ID / HGNC symbol)
#      - Gene length (in bp)
#      - Ensembl transcript ID
#      - Distance between repeats (default: 0-1000 bp)
#      - Overview of total results for the given repeat length

# If function output is assigned to a variable and if the "summary" 
# parameter is TRUE, the variable will contain two data frames:
#     var$summary -- same as the console-printed overview
#     var$results -- dataframe with all results of the query. This detailed 
#                    results is pairwise so each 2 rows is a structure.
# If the "summary" parameter is FALSE, variable will only contain detailed results

# All console outputs include the runtime



# ***************************************
# ---------------- USAGE ----------------
# ***************************************

# Only printing a single length result overview for given gene to console
gene_lookup("GAA", length = 7, mindist = 10, maxdist = 1000)
gene_lookup("GLA", mindist = 10, maxdist = 1000)
gene_lookup("GBA", length = 20, mindist = 10, maxdist = 1000)
gene_lookup("GAA", length = 20, mindist = 10, maxdist = 1000)
gene_lookup("gaa", length = 10, mindist = 10, maxdist = 1000)
gene_lookup("GAA", length = 7, mindist = 400, maxdist = 500)

# Only printing all length results overview for given gene to console
gene_lookup("GAA", mindist = 10, maxdist = 1000)
gene_lookup("SRY", mindist = 10, maxdist = 1000)
test <- gene_lookup("SRY", length = 9, mindist = 10, maxdist = 1000)

# Printing all length results overview for given gene to console
# Detailed results stored in variable
gaa7 <- gene_lookup("GAA", length = 7, mindist = 10, maxdist = 1000)
gaa20 <- gene_lookup("GAA", length = 20, mindist = 10, maxdist = 1000)
gaa10 <- gene_lookup("GAA", length = 10, mindist = 10, maxdist = 1000)

# Setting summary to TRUE will give tibble including the console summary and stored data
acan10 <- gene_lookup("ACAN", length = 10, mindist = 10, maxdist = 1000, summary = T)
head(acan10$summary)
head(acan10$results)