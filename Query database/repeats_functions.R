#!/usr/bin/Rscript

# repeats_functions.R by DT Laura Vo Ngoc
# Part of the pipeline used to query the EDIR database
# 
#*******************************************************************************
# Written initially in R 3.5.0 Joy in Playing
#
# PURPOSE: Functions and dependencies developed to query the EDIR database.
#
# Must be stored in same directory as INPUT folder, and repeats_USE.R script.
# 
# 
# Written to run locally in RStudio.
#
#*******************************************************************************



# ************************************************
# ---------------- LOAD LIBRAIRES ----------------
# ************************************************

if (!suppressMessages(suppressWarnings(require("dplyr")))) install.packages("dplyr")
if (!suppressMessages(suppressWarnings(require("tibble")))) install.packages("tibble")
if (!suppressMessages(suppressWarnings(require("tictoc")))) install.packages("tictoc")

# Load libraries
suppressMessages(suppressWarnings(library("dplyr")))               # v.0.8.3
suppressMessages(suppressWarnings(library("tibble")))              # v.2.1.3
suppressMessages(suppressWarnings(library("tictoc")))              # v.1.9

# Set path of file location
setwd(path)


# ********************************************
# ---------------- INITIALIZE ----------------
# ********************************************

# Set date variable
today <- format(Sys.Date(), format="%y%m%d")

# Set data name for output saving purposes
directory <- "./Input/"

# Load list genes with corresponding chromosome
genes_chr <- readRDS("./Input/genes_chromosome.rds")



# ******************************************************
# ---------------- GENE LOOKUP FUNCTION ----------------
# ******************************************************
# 
# Function to look up a specified gene in the dataset
# 
# REQUIRED INPUT: 
#     gene ID (ENSEMBL or HGNC)
# OPTIONAL INPUT: 
#     length ---- repeat length (7-20), if not defined will output all results
#     mindist --- minimum distance, default at 0
#     maxdist --- maximum distance, default at 1000
#     summary --- whether to include summary in output, defaults to FALSE
# 
# OUTPUT:
#     Prints to console:
#      - Gene ID (Ensembl ID / HGNC symbol)
#      - Gene length (in bp)
#      - Ensembl transcript ID
#      - Distance between repeats (default: 0-1000 bp)
#      - Overview of total results for the given repeat length
# 
# If function output is assigned to a variable and if the "summary" 
# parameter is TRUE, the variable will contain two data frames:
#     var$summary -- same as the console-printed overview
#     var$results -- dataframe with all results of the query. This detailed 
#                    results is pairwise so each 2 rows is a structure.
# If the "summary" parameter is FALSE, variable will only contain detailed results
# 
# All console outputs include the runtime
# 
# **********************************************************************************


gene_lookup <- function(gene, 
                        length = NA, 
                        mindist = 0, 
                        maxdist = 1000,
                        summary = F) {
    # USER INPUT CHECK 
    # Gene ID incorrect, or not Ensembl ID/HGNC symbol
    if (is.na(any(genes_chr==gene)) || !any(genes_chr==gene)) {
        stop("Incorrect gene ID: use ENSEMBL ID or HGNC symbol")
    }
    
    # Minimum distance >= maximumn distance
    if (mindist >= maxdist) {
        stop("mindist must be less than maxdist")
    }
    
    
    tic("Runtime")
    # Initialize: uppercase gene ID, get chromosome number
    gene <- toupper(gene)
    chromosome <- genes_chr$chromosome_name[which(genes_chr$ensembl_gene_id==gene | genes_chr$hgnc_symbol==gene)]
    
    # If length is defined
    if (!is.na(length)) {
        # Load file corresponding to defined length and gene's chromosome number
        files <- paste0("./Input/", list.files(path="./Input/", pattern=paste0("^", as.character(length), "bp_dist1000_chr", chromosome, ".rds")))
        numgenes <- readRDS(files)
    }
    
    # If length not defined
    else {
        # Load all files corresponding to gene's chromosome number
        files <- paste0("./Input/", list.files(path="./Input/", pattern=paste0("bp_dist1000_chr", chromosome, ".rds")))
        numgenes <- do.call(rbind,
                            lapply(files, readRDS))
    }
    
    # Filter only resuls of queried gene
    results <- numgenes[which(numgenes$ensembl_gene_id==gene | numgenes$hgnc_symbol==gene),]
    
    # Format numeric cols
    results$repeat_length <- as.numeric(results$repeat_length)
    results$distance <- as.numeric(results$distance)
    
    # Filter on mindist and maxdist
    results <- results[which(results$distance >= mindist & results$distance <= maxdist),]
    
    # Sort results by repeat length, chromosome number, repeat sequence, and start position
    results <- results[order(results$repeat_length, results$chromosome, results$repeat_seq, results$start),]
    
    # No results
    if(nrow(results)==0){
        # cat("Gene not in dataset\n\n")
        cat("NO RESULT FOUND\n\n")
        toc()
    } else {
        # summ df
        summ <- data.frame(repeat_length = sort(as.numeric(unique(results$repeat_length))), 
                              unique_seqs = NA)
        
        # For each repeat length in results: number unique repeat_seq, number unique instances, number structures, average distance
        for (i in summ$repeat_length) {
            temp <- results[which(results$repeat_length==i),]
            
            summ$unique_seqs[which(summ$repeat_length==i)] <- length(unique(temp$repeat_seq))
            summ$tot_instances[which(summ$repeat_length==i)] <- nrow(temp[!duplicated(temp[,c(1:5,7:13)]),c(1:5,7:13)])
            summ$tot_structures[which(summ$repeat_length==i)] <- nrow(temp)/2
            summ$avg_dist[which(summ$repeat_length==i)] <- mean(as.numeric(temp$distance))
        }
        
        # Results normalizd by gene length
        gene_length <- as.numeric(strsplit(results$gene_range[1],"-")[[1]][2]) - as.numeric(strsplit(results$gene_range[1],"-")[[1]][1])
        summ$norm_instances_bp <- summ$tot_instances/gene_length
        summ$norm_instances_Mb <- summ$tot_instances/gene_length*1e+06
        summ$norm_structures_bp <- summ$tot_structures/gene_length
        summ$norm_structures_Mb <- summ$tot_structures/gene_length*1e+06
        
        # Print summ to console
        if (is.na(length)) {
            cat("\tParameters\n") 
        } else {
            cat("\tParameters\n Repeat length:\t", length, " bp")
        }
        cat(" \n Gene:\t\t ", results$ensembl_gene_id[1], 
            " / ", results$hgnc_symbol[1],
            " \n Gene length:\t ", gene_length,
            " bp\n Transcript ID:\t ", results$ensembl_transcript_id[1],
            "\n Distance:\t ", mindist, "-", maxdist, " bp",
            "\n\n",
            sep="")
        
        # Print summ df
        print(summ)
        cat("\n")
        results <- results[order(as.numeric(rownames(results))),]
        toc()
        
        # Output list of summ and detailed results
        if (summary == F ){
             out <- results
        } else {
             summary <- summ
             out <- tibble::lst(summary, results)
        }
        invisible(out)
    }
}

