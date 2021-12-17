#!/usr/bin/Rscript

# geneAnnotate_features.R by DT Laura Vo Ngoc
# Part of the pipeline used to generate the EDIR database
# 
#*******************************************************************************
# Written initially in R 3.5.0 Joy in Playing
#
# PURPOSE: Takes input sequences and chromosomal positions, and annotates gene ID,
#       transcript, and intron/exon number. Pairs of identical repeat sequences are 
#       labeled with comparison features: in same intron or exon, spanning an intron 
#       and exon, or flanking an intron or exon. Then only pairs not in the same intron
#       are retained.
# 
# A summary output is also generated.
#

# Written to run from command line, parallized.
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# Rscript ./geneAnnotate_optimized.R input_directory
#-----------------------------------
#
#*******************************************************************************



# ***************************
# --- INIT FROM ARGUMENTS ---
# ***************************

args <- commandArgs(trailingOnly = T)

# Check if there is at least one argument: if not, return an error
if (length(args)==0) {
    stop("At least one argument must be supplied (input directory name).", call.=FALSE)
}


# **********************************
# --- LIBRARIES, REFERENCE DATA ----
# **********************************

# Check that required packages are installed
if (!suppressMessages(suppressWarnings(require("BiocManager")))) install.packages("BiocManager")
if (!suppressMessages(suppressWarnings(require("GenomicRanges")))) BiocManager::install("GenomicRanges")
if (!suppressMessages(suppressWarnings(require("dplyr")))) install.packages("dplyr")
if (!suppressMessages(suppressWarnings(require("doParallel")))) install.packages("doParallel")
if (!suppressMessages(suppressWarnings(require("foreach")))) install.packages("foreach")

# Load libraries
suppressMessages(suppressWarnings(library("GenomicRanges")))         # v.1.11.18
suppressMessages(suppressWarnings(library("dplyr")))                 # v.0.8.1
suppressMessages(suppressWarnings(library("doParallel")))            # v.1.0.15
suppressMessages(suppressWarnings(library("foreach")))               # v.1.4.7

# Load reference data
gr.intex <- readRDS("./geneAnnotate/Data_input/refgenomeAllGR.rds")

# mendeliome = list of genes in the Mendeliome
load("./geneAnnotate/Data_input/mendeliome_gene_list.RData")

mindist <- 0


# ********************************************
# ---------- LOAD DATA TO ANNOTATE ----------
# ********************************************

# genfiles <- list.files("../Output/seqFind/")
genfiles <- list.files(args)

# Set parallelization parameters
cl <- parallel::makeCluster(Sys.getenv("PBS_NP"))
registerDoParallel(cores = Sys.getenv("PBS_NP"))

para <- foreach (i = 1:length(genfiles)) %dopar% {
    worker <- Sys.getpid()
    
    # Load input
    genfile <- genfiles[i]
    seqfind <- as.data.frame(read.table(paste("../Output/seqFind/", genfile, sep=""), header=TRUE, stringsAsFactors=FALSE))
    
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
    
    
    # INITIALIZE OTHER VARIABLES
    # Set distance between instances
    maxdist <- as.integer(strsplit(strsplit(genfile,'[.]')[[1]][2], "bp")[[1]][1])
    # Get repeat seq
    repeat.seq <- strsplit(genfile,'[_]')[[1]][1]
    repeat.len <- nchar(repeat.seq)
    
    print(repeat.seq, quote = F)
    
    
    # ********************************************
    # ---------------- ANNOTATION ----------------
    # ********************************************
    
    # Find overlaps between seqFind output and reference genmome
    annot.all <- subsetByOverlaps(gr.seqfind, gr.intex)
    overlap.all <- findOverlaps(annot.all, gr.intex)
    
    # Initialize some column name groups
    info <- c("ensembl_gene_id", "hgnc_symbol")
    info.2 <- c("intron_exon", "number")
    
    # Add metadata columns to results: 
    # ensembl_gene_id, hgnc_symbol, concatenate gene start/end positions, ensembl_transcript_id, 
    # concatenate intron/exon start/end positions, whether the position refers to an intron or exon, 
    # and which intron/exon number corresponds to position. Adds a "feature" column for later filtering
    mcols(annot.all)[queryHits(overlap.all), info] <- mcols(gr.intex)[subjectHits(overlap.all), info]
    mcols(annot.all)[queryHits(overlap.all), "gene_range"] <- paste(mcols(gr.intex)[subjectHits(overlap.all), "start_position"], 
                                                                    "-", 
                                                                    mcols(gr.intex)[subjectHits(overlap.all), "end_position"], 
                                                                    sep="")
    mcols(annot.all)[queryHits(overlap.all), "ensembl_transcript_id"] <- mcols(gr.intex)[subjectHits(overlap.all), "ensembl_transcript_id"]
    mcols(annot.all)[queryHits(overlap.all), "transcript_range"] <- paste(start(gr.intex[subjectHits(overlap.all)]), 
                                                                          "-", 
                                                                          end(gr.intex[subjectHits(overlap.all)]), 
                                                                          sep="")
    
    mcols(annot.all)[queryHits(overlap.all), info.2] <- mcols(gr.intex)[subjectHits(overlap.all), info.2]
    
    # Convert to dataframe; remove "width" column, and reorder columns
    results <- as.data.frame(annot.all)
    results <- results[,-4]
    results <- results[,c(5,2:4,1,6:12)]
    colnames(results)[5] <- "chromosome"
    


    # ****************************************
    # ---------- FEATURES ANNOTATION ---------
    # ****************************************

    # Number code for features  
    #                       (1)                 (2)             (3)           (4)             (5)
    features <- c("spanning intron-exon", "same intron", "flanking exon", "same exon", "flanking intron")
    
    # Initialize output df
    out <- data.frame()
    
    for (j in 1:nrow(results)) {
        # print(i)
    # for (i in 1:10) {
        k <- j+1
        
        # Get row
        row <- results[j,]
        
        # While k within df; same repeat sequence, same gene
        while(k <= nrow(results) &
            results$repeat_seq[k] == row$repeat_seq &
            results$ensembl_gene_id[k] == row$ensembl_gene_id) {
            # print(k)
            row1 <- results[k,]
            
            # Define distance variable
            dist <- row1$start - row$end
            
            # If distance within range of interest
            if (dist >= mindist && dist <= maxdist) {
                
                # Define feature (see number code above)
                # If intron and exon
                feat <- ifelse(row$intron_exon != row1$intron_exon,
                               # TRUE: diff intron_exon
                               1,
                               # FALSE: same intron_exon
                               ifelse(row$intron_exon == "I",
                                      # TRUE: row intron
                                      ifelse(row$number == row1$number,
                                             # TRUE: same number
                                             2,
                                             # FALSE: diff number
                                             3),
                                      # FALSE: row exon
                                      ifelse(row$number == row1$number,
                                             # TRUE: same number
                                             4,
                                             # FALSE: diff number
                                             5)))
                
                # If feature is NOT "same intron" (2)
                if (feat != 2) {
                    # append row number, distance, and feature of both rows
                    out <- rbind(out, c(j, dist, feat), c(k, dist, feat))
                }
                
                rm(feat)
            }
            
            # Increment k
            k <- k+1
            
        } 
    }
    colnames(out) <- c("rownum", "dist", "feat")
    
    final.res <- results[out$rownum,]
    final.res$distance <- out$dist
    final.res$repeat_length <- 7
    final.res <- final.res[, c(5, 14, 1:3, 13, 6:12)]
    final.res$feature <- features[out$feat]
    rownames(final.res) <- NULL

    
    
    # *******************************************************
    # ----- GENERATE OUTPUT, MENDELIOME, SUMMARY TABLES -----
    # *******************************************************
    final.res <- data.frame(read.table(file="../Output/geneAnnotate/results_repeat_CGGGATC_dist_1000_bp.tsv", sep = "\t", header = T), stringsAsFactors = F)
    
    
    # Output final.res to results table
    write.table(final.res,
                paste("../Output/geneAnnotate/n-results_repeat", repeat.seq, "dist", maxdist, "bp.tsv", sep="_"), 
                sep="\t", row.names=F, col.names=T, quote=F)
    
    
    # Output results in Mendeliome list to a table
    res.mendeliome <- final.res[which(final.res$hgnc_symbol %in% mendeliome),]
  
    
    # Generate and append summary of all results to summary table
    # Set condition for whether to include colnames
    if (file.exists(paste("../Output/geneAnnotate/summary/", worker, "-summary_results.tsv", sep=""))) {
        headerrowsum <- FALSE
    } else {
        headerrowsum <- TRUE
    }
    # Generate and append
    summary <- data.frame(sequence=repeat.seq, 
                          abs_freq=nrow(final.res), 
                          genes=length(unique(final.res$hgnc_symbol)), 
                          mend_freq=nrow(res.mendeliome),
                          mendeliome=length(unique(res.mendeliome$hgnc_symbol)), 
                          fullGenelist=paste(unique(final.res$hgnc_symbol), collapse="/"),
                          mendeliomeGenes=paste(unique(res.mendeliome$hgnc_symbol), collapse="/"))
    write.table(summary, 
                paste("../Output/geneAnnotate/summary/", worker, "-summary_results.tsv", sep=""), 
                sep="\t", row.names=F, col.names=headerrowsum, quote=F, append=T)
}
parallel::stopCluster(cl)



#******************************************************************
###-***********************END OF SCRIPT***********************-###
#******************************************************************