#!/usr/bin/Rscript

# geneAnnotate_features.R by DT Laura Vo Ngoc
# Part of the pipeline used to generate the EDIR database
# 
#*******************************************************************************
# Written initially in R 3.5.0 Joy in Playing
#
# PURPOSE: Annotates pairs of repeat sequences as in same intron or exon,  
#       spanning an intron and exon, or flanking an intron or exon. Then only 
#       pairs not in the same intron are retained.
# 
# A summary output is also generated.
#

# Written to run from command line, parallized.
#
# TO RUN:
# Examples,
# Enter on the command line of your terminal, the line
#-----------------------------------
# Rscript ./geneAnnotate_features.R input_directory
#-----------------------------------
#
#*******************************************************************************



# ***************************
# --- INIT FROM ARGUMENTS ---
# ***************************

args <- commandArgs(trailingOnly = T)

# Check if there is at least one argument: if not, return an error
if (length(args) == 0) {
    stop("At least one argument must be supplied (input directory name).", call. = FALSE)
}

outdir <- args[2]
sumdir <- args[3]



# *************************************
# --- LIBRARIES AND REFERENCE DATA ----
# *************************************

# Check that required packages are installed
if (!suppressMessages(suppressWarnings(require("dplyr")))) install.packages("dplyr")
if (!suppressMessages(suppressWarnings(require("doParallel")))) install.packages("doParallel")
if (!suppressMessages(suppressWarnings(require("foreach")))) install.packages("foreach")

# Load libraries
suppressMessages(suppressWarnings(library("dplyr")))                 # v.0.8.1
suppressMessages(suppressWarnings(library("doParallel")))            # v.1.0.15
suppressMessages(suppressWarnings(library("foreach")))               # v.1.4.7


# mendeliome = list of genes in the Mendeliome
load("./geneAnnotate/Data_input/mendeliome_gene_list.RData")

mindist <- 0



# ********************************************
# ---------- LOAD DATA TO ANNOTATE ----------
# ********************************************

# Load data
genfiles <- list.files(args[1])

# Parallelization parameters
cl <- parallel::makeCluster(3)
registerDoParallel(cores = Sys.getenv("PBS_NP"))
# registerDoParallel(cores = 3)

para <- foreach (i = 1:length(genfiles)) %dopar% {
    worker <- Sys.getpid()
    
    # Load input
    genfile <- genfiles[i]
    results <- as.data.frame(read.table(paste(args[1], "/", genfile, sep = ""), header = TRUE, stringsAsFactors = FALSE))
    

    # INITIALIZE OTHER VARIABLES
    # Set distance between instances
    maxdist <- as.integer(strsplit(strsplit(genfile, "[.]")[[1]][3], "bp")[[1]][1])
    # Get repeat seq
    repeat.seq <- strsplit(genfile,'[.]')[[1]][2]
    print(repeat.seq, quote = F)
    
    # If input file not empty
    if (nrow(results) > 0) {
        
        # unique rows (excluding "distance" column)
        results <- unique(results[-6])
        results <- results[order(results$repeat_length, results$chromosome, results$repeat_seq, results$start),]
        


        # ****************************************
        # ---------- FEATURES ANNOTATION ---------
        # ****************************************

        # Number code for features  
        #                       (1)                 (2)             (3)           (4)             (5)
        features <- c("spanning intron-exon", "same intron", "flanking exon", "same exon", "flanking intron")
        
        # Initialize output df
        out <- data.frame()
        
        for (j in 1:nrow(results)) {
            k <- j+1
            
            # Get row
            row <- results[j,]
            
            # While k within df; same repeat sequence, same gene
            while(k < =  nrow(results) &
                  results$repeat_seq[k] = =  row$repeat_seq &
                  results$ensembl_gene_id[k] = =  row$ensembl_gene_id) {
                # print(k)
                row1 <- results[k,]
                
                # Define distance variable
                dist <- row1$start - row$end
                
                # If distance within range of interest
                if (dist > =  mindist && dist < =  maxdist) {
                    
                    # Define feature (see number code above)
                    # If intron and exon
                    feat <- ifelse(row$intron_exon ! =  row1$intron_exon,
                                   # TRUE: diff intron_exon
                                   1,
                                   # FALSE: same intron_exon
                                   ifelse(row$intron_exon = =  "I",
                                          # TRUE: row intron
                                          ifelse(row$number = =  row1$number,
                                                 # TRUE: same number
                                                 2,
                                                 # FALSE: diff number
                                                 3),
                                          # FALSE: row exon
                                          ifelse(row$number = =  row1$number,
                                                 # TRUE: same number
                                                 4,
                                                 # FALSE: diff number
                                                 5)))
                    
                    # If feature is NOT "same intron" (2)
                    if (feat ! =  2) {
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
        
        if (nrow(out) > 0) {
            print("yes")
        }
        

        final.res <- results[out$rownum,]
        final.res$distance <- out$dist
        final.res <- final.res[, c(1:5, 13, 6:12)]
        final.res$feature <- features[out$feat]
        rownames(final.res) <- NULL

    
    
        # *******************************************************
        # ----- GENERATE OUTPUT, MENDELIOME, SUMMARY TABLES -----
        # *******************************************************
    
        # Output final.res to results table
        write.table(final.res,
                    paste(outdir, "/expansion_repeat_", repeat.seq, "_dist", maxdist, "_bp.tsv", sep = ""),
                    sep = "\t", row.names = F, col.names = T, quote = F)
    
    
        # Output results in Mendeliome list to a table
        res.mendeliome <- final.res[which(final.res$hgnc_symbol %in% mendeliome),]
    
        # Generate and append summary of all results to summary table
        # Set condition for whether to include colnames
        if (file.exists(paste(sumdir, "/", worker, "-summary_expand_results.tsv", sep = ""))) {
            headerrowsum <- FALSE
        } else {
            headerrowsum <- TRUE
        }
        # Generate and append
        summary <- data.frame(repeat_length = NA,
                              sequence = names(table(final.res$repeat_seq)),
                              abs_freq = as.vector(table(final.res$repeat_seq)),
                              genes = NA,
                              mend_freq = NA,
                              mendeliome = NA,
                              fullGenelist = NA,
                              mendeliomeGenes = NA)
        
        # Populate summar table
        for (i in 1:nrow(summary)) {
            seq <- summary$sequence[i]
            summary$repeat_length[i] <- nchar(toString(seq))
            summary$genes[i] <- length(unique(final.res$hgnc_symbol[which(final.res$repeat_seq == seq)]))
            summary$fullGenelist[i] <- paste(unique(final.res$hgnc_symbol[which(final.res$repeat_seq == seq)]), collapse = "/")
            if (seq %in% res.mendeliome$repeat_seq) {
                summary$mend_freq[i] <- as.vector(table(res.mendeliome$repeat_seq)[which(names(table(res.mendeliome$repeat_seq)) == seq)])
                summary$mendeliome[i] <- length(unique(res.mendeliome$hgnc_symbol[which(res.mendeliome$repeat_seq == seq)]))
                summary$mendeliomeGenes[i] <- paste(unique(res.mendeliome$hgnc_symbol[which(res.mendeliome$repeat_seq == seq)]), collapse = "/")
            }
        }
        summary <- summary[order(summary$repeat_length, summary$sequence),]
    
        write.table(summary,
                    paste(sumdir, "/", worker, "-summary_expand_results.tsv", sep = ""),
                    sep = "\t", row.names = F, col.names = headerrowsum, quote = F, append = T)
    } else {
        
        # Add sequence to "null" list
        if (file.exists(paste(outdir, "/", worker, "-expand_results_NULL.tsv", sep = ""))) {
            headerrowsum <- FALSE
        } else {
            headerrowsum <- TRUE
        }
        
        write.table(repeat.seq,
                    paste(outdir, "/", worker, "-expand_results_NULL.tsv", sep = ""),
                    sep = "\t", row.names = F, col.names = headerrowsum, quote = F, append = T)
    }
}
parallel::stopCluster(cl)


#******************************************************************
###-***********************END OF SCRIPT***********************-###
#******************************************************************