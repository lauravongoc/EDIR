#!/usr/bin/Rscript

# checkGene.R by DT Laura Vo Ngoc
# Part of the pipeline used to generate the EDIR database
# 
#*******************************************************************************
# Written initially in R 4.1.3 -- One Push-Up
#
# PURPOSE: Function that takes an input dataframe of chromosomal positions, 
#          and filters out rows not falling in the position ranges of known genes.
#
# Requires the input file "refgenomeAllGR.rds" 
# 
# 
# Written to receive and output in the Python script:
#      find_seq_element_occrrences_in_seq.py (line 306)
#
#*******************************************************************************



# Check that required packages are installed
if (!suppressMessages(suppressWarnings(require("BiocManager")))) install.packages("BiocManager")
if (!suppressMessages(suppressWarnings(require("GenomicRanges")))) BiocManager::install("GenomicRanges")

# Load libraries
suppressMessages(suppressWarnings(library("GenomicRanges")))         # v.1.11.18

# Load reference data
gr.intex <- readRDS("./Input/refgenomeAllGR.rds")


checkGene <- function(seqfind) {
    # FORMAT DATA FOR ANALYSIS
    # Change strand to * to omit filters on strand
    seqfind$strand <- "*"

    # Change colnames
    colnames(seqfind)[1] <- "repeat_seq"

    # Convert to GRanges object
    gr.seqfind <- makeGRangesFromDataFrame(seqfind,
                                           keep.extra.columns=TRUE,
                                           ignore.strand=FALSE,
                                           seqinfo=NULL,
                                           seqnames.field=c("seqnames", "seqname",
                                                            "chromosome", "chrom",
                                                            "chr", "chromosome_name",
                                                            "seqid"),
                                           start.field="start_pos",
                                           end.field=c("end_pos", "stop"),
                                           strand.field="strand",
                                           starts.in.df.are.0based=FALSE)

    annot.all <- subsetByOverlaps(gr.seqfind, gr.intex)

    # Convert to dataframe; remove "width" column, and reorder columns
    results <- as.data.frame(annot.all)
    results <- results[,-4]
    results <- results[,c(5,2:4,1)]
    colnames(results)[5] <- "chromosome"

    return(results)
}
