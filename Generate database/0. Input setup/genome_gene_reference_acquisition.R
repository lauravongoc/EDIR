#!/usr/bin/Rscript

# genome_gene_reference_acquisition.R by DT Laura Vo Ngoc
# Part of the pipeline used to generate the EDIR database
# 
#*******************************************************************************
# Written initially in R 3.5.0 Joy in Playing
#
# PURPOSE: Using biomaRt to get a reference list of genes IDs, corresponding 
#       chromosomes, chromosomal positions, transcript IDs, exon start/end 
#       positions, and strand.
#

# Written to run locally from RStudio.
#
#*******************************************************************************


# ********************************************
# ---------------- INITIALIZE ----------------
# ********************************************

# Load libraries
library("biomaRt")          # v.2.39.2
library("dplyr")            # v.0.8.1
library("gtools")           # v.3.8.1
library("tictoc")           # v.1.0
library("GenomicRanges")    # v.1.34.0


# Set wd in find_seq_element_occurrences.py output folder
setwd("C:/Users/rockp/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu18.04onWindows_79rhkp1fndgsc/LocalState/rootfs/home/laura/genRepeats/Genome-repeats/")


# ********************************************
# ---- GENOME GENE INFORMATION AQUISITION ----
# ********************************************

# Get list of available datasets
ensembl <- useMart("ensembl")
listDatasets(ensembl)

# Define dataset to query
dataset <- "hsapiens_gene_ensembl"

# Query Ensembl for dataset
ensembl <- useMart("ensembl", dataset=dataset)


# Get all relevant gene info in each chr
head(getBM(c("chromosome_name","ensembl_gene_id" ,"hgnc_symbol", 'external_gene_name', "start_position", "end_position"),
           filters = c("chromosome_name"),
           values = list(1),
           mart = ensembl))

# Split chromosomes into 4 groups (importing all at once is impossible)
chromosomes1 <- c(noquote(1:5))
chromosomes2 <- c(noquote(6:10))
chromosomes3 <- c(noquote(11:15))
chromosomes4 <- c(noquote(16:22), "X", "Y")

# Get gene information from each chromosome group
chr.data1 <- getBM(c("chromosome_name","ensembl_gene_id","start_position","end_position",
                     "ensembl_transcript_id", "exon_chrom_start","exon_chrom_end","strand"),
                   filters = c("chromosome_name"),
                   values = list(chromosomes1),
                   mart = ensembl)
chr.data2 <- getBM(c("chromosome_name","ensembl_gene_id","start_position","end_position",
                     "ensembl_transcript_id","exon_chrom_start","exon_chrom_end","strand"),
                   filters = c("chromosome_name"),
                   values = list(chromosomes2),
                   mart = ensembl)
chr.data3 <- getBM(c("chromosome_name","ensembl_gene_id","start_position","end_position",
                     "ensembl_transcript_id","exon_chrom_start","exon_chrom_end","strand"),
                   filters = c("chromosome_name"),
                   values = list(chromosomes3),
                   mart = ensembl)
chr.data4 <- getBM(c("chromosome_name","ensembl_gene_id","start_position","end_position",
                     "ensembl_transcript_id","exon_chrom_start","exon_chrom_end","strand"),
                   filters = c("chromosome_name"),
                   values = list(chromosomes4),
                   mart = ensembl)

# Combine gene information into one df
chr.data <- rbind(chr.data1, chr.data2, chr.data3, chr.data4)


# Setup mart for gene symbols
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
              dataset = dataset, 
              host = "www.ensembl.org",
              ensemblRedirect = FALSE)


# Get gene symbols from Biomart
gene.symbols <- getBM(c("ensembl_gene_id", "hgnc_symbol"),
                      filters = c("ensembl_gene_id"),
                      values = list(chr.data$ensembl_gene_id),
                      mart = mart)

# Get gene symbols from Biomart -- Pan troglodytes
gene.symbols <- getBM(c("ensembl_gene_id", "external_gene_name"),
                      filters = c("ensembl_gene_id"),
                      values = list(chr.data$ensembl_gene_id),
                      mart = mart)

# Annotate genes table with symbols, reorder to put next to ensembl gene id
refInfo <- left_join(chr.data, gene.symbols)
refInfo <- refInfo[c(1,2,9,3:8)]

# Sort by start position, add an index, put index as first column
refInfo <- refInfo[order(refInfo$start_position),]
refInfo$index <- c(1:nrow(refInfo))
refInfo <- refInfo[,c(10,1:9)]

# Sort by chromosome number, then by gene start position, save as RData file
refInfo <- refInfo[with(refInfo, order(refInfo$chromosome_name, refInfo$start_position)),]

# Save df as RDS file
saveRDS(refInfo, file="./geneAnnotate/Data_input/refgenomeAllGR.rds")

refInfo <- readRDS("./geneAnnotate/Data_input/20210814_pan_troglodytes_biomart_gene_exons.rds")

single <- 0
df <- data.frame(gene=NA, num_transcr=NA)
df <- df[-1,]
for (gene in unique(refInfo$ensembl_gene_id)[1]) {
    transcr <- length(unique(refInfo$ensembl_transcript_id[which(refInfo$ensembl_gene_id==gene)]))
    if (transcr > 1) {
        df <- append(df, c(gene, transcr))
    } else {
        single <- single+1
    }
    
}

# ********************************************
# ----- GENE INFO SEPARATE BY CHROMOSOME -----
# ********************************************

# Initialize list of chromosome numbers/letters
nums <- c(1:22, "X", "Y")

# Loop through chromosome numbers
for (i in nums) {
    # Form chromosme name (chr#)
    chrom <- paste("chr", i, sep = "")
    # Initialize filename
    filename <- paste("./Data/", chrom, ".RData", sep="")
    
    # Select gene data for chromosome
    chrom_df <- refInfo[which(refInfo$chromosome_name==i),]
    
    # Retain only single row per gene, and remove transcript/exon info
    chrom_df <- chrom_df[!duplicated(chrom_df$ensembl_gene_id), 3:6]
    chrom_df <- chrom_df[which(chrom_df$hgnc_symbol!=""),]
    
    # Assigns df to appropriate variable name, and saves
    assign(chrom, chrom_df)
    save(list=chrom, file=filename)
}


# ********************************************
# ------- MENDELIOME GENE LIST EXTRACT -------
# ********************************************

# Read in text file as vector
# Mendeliome list obtained from BRIGHTcore
mendeliome <- scan(file="Mendeliome_gene_list.txt", what=character(), sep="\n")

# Save as RData file
save(mendeliome, file="./Data_input/mendeliome_gene_list.RData")




# *********************************
# ------- LOAD BIOMART DATA -------
# *********************************

# Load reference data (refInfo, chromosome-separated data,...)
for (data_file in list.files(path="./geneAnnotate/Data_input")) {
    load(paste("./geneAnnotate/Data_input/", data_file, sep=""))
}

# List chromosome-split gene df names
chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
                 "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
                 "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")

chrY$chromosome <- "Y"

refgenome <- rbind(chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10,
                   chr11, chr12, chr13, chr14, chr15, chr16, chr18, chr19, chr20,
                   chr21, chrX, chrY)

save(refgenome, file="./geneAnnotate/Data_input/refgenome_concat.RData")


# Extract just list of gene Ensembl IDs
gr.intex <- readRDS("./geneAnnotate/Data_input/refgenomeAllGR.rds")

refInfo <- as.data.frame(gr.intex)
data <- refInfo[which(refInfo$hgnc_symbol!=""),]
allgenes <- sort(refInfo$ensembl_gene_id)
allgenes <- unique(allgenes)
write.table(allgenes, file="./geneExtract/Input/EnsemblAllGenesSort.txt", col.names = F, row.names = F, quote = F)


#******************************************************************
###-***********************END OF SCRIPT***********************-###
#******************************************************************