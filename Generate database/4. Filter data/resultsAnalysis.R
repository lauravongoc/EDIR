#!/usr/bin/Rscript

# resultsAnalysis.R by DT Laura Vo Ngoc
# Part of the pipeline used to generate the EDIR database
# 
#*******************************************************************************
# Written initially in R 3.5.0 Joy in Playing
#
# PURPOSE: Takes final output of the pipeline (resultsByLength.sh and 
#       optiStructureLengthResults.sh), and generates plots of the data.
# 
#
# 
# Written to run locally in RStudio.
##
#*******************************************************************************


# **********************************************
# ---- INIT LIBRARIES AND WORKING DIRECTORY ----
# *********************************************

library("ggplot2")
library("reshape2")

setwd("C:/Users/rockp/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu18.04onWindows_79rhkp1fndgsc/LocalState/rootfs/home/laura/genRepeats/Genome-repeats")


# **************
# ---- INIT ----
# **************

# Set date variable
today <- format(Sys.Date(), format="%y%m%d")

# Set data name for output saving purposes
directory <- "../Output/Analysis_fixed/"
parameters <- "dist1000"
filename <- paste(today, parameters, sep="_")
dist <- 1000

# List of chromosomes
chr_label <- c(1:22, "X", "Y")

# Genes, transcripts, intron/exon information
reference <- data.frame(readRDS("./geneAnnotate/Data_input/refgenomeAllGR.rds"), stringsAsFactors = F)

# Mendeliome genes list
mendeliome <- read.table("./geneAnnotate/Data_input/Mendeliome_gene_list.txt", stringsAsFactors = F)
mendeliome <- mendeliome$V1


# *********************************
# ---- INIT GENEEXTRACT OUTPUT ----
# *********************************

# List files in input directory
files <- list.files(path="../Output/Genes_fixed/", pattern="1000.tsv")

for (file in files) {
    print(file)

    # Init: filename, repeat length, import table
    replen <- as.numeric(strsplit(file, "bp")[[1]][1])
    genes <- as.data.frame(read.table(paste("../Output/Genes_fixed/", file, sep=""), sep="\t", stringsAsFactors = F))
    colnames(genes) <- c("ensembl_gene_id", "tot_structures", "avg_dist", "sequences")
    genes$repeat_length <- replen
    numgenes <- genes[, c(5, 1:3)]
    rm(genes)
    
    # Label with gene info: ensembl gene ID, HGNC smbol, start/end, gene length
    numgenes$chromosome <- reference$seqname[match(numgenes$ensembl_gene_id, reference$ensembl_gene_id)]
    numgenes$hgnc_symbol <- reference$hgnc_symbol[match(numgenes$ensembl_gene_id, reference$ensembl_gene_id)]
    numgenes$start_position <- reference$start_position[match(numgenes$ensembl_gene_id, reference$ensembl_gene_id)]
    numgenes$end_position <- reference$end_position[match(numgenes$ensembl_gene_id, reference$ensembl_gene_id)]
    numgenes$gene_length <- numgenes$end_position-numgenes$start_position
    numgenes <- numgenes[with(numgenes, order(chromosome, start_position)),]
    
    # Normalize number of structures by gene length
    numgenes$norm_struct_bp <- numgenes$tot_structures/numgenes$gene_length
    numgenes$norm_struct_Mb <- numgenes$tot_structures/numgenes$gene_length*1e+06
    
    # Remove duplicated rows; output
    numgenes <- numgenes[!is.na(numgenes$tot_structures),]
    numgenes <- unique(numgenes)
    saveRDS(numgenes, file=paste0(directory, replen, "bp_", parameters, "_genome.rds"))
    
    # Order by number of structures
    ord.genes <- numgenes[order(-numgenes$tot_structures),]
    ord.genes <- unique(ord.genes)
    write.csv(ord.genes, file=paste(directory, replen, "bp_", parameters, "_ord_genes.csv", sep=""), row.names = F, quote = F)
    
    # Mendeliome genes only
    mendgenes <- numgenes[which(numgenes$hgnc_symbol %in% mendeliome),]
    mendgenes <- unique(mendgenes)
    saveRDS(mendgenes, file=paste0(directory, replen, "bp_", parameters, "_mendeliome.rds"))

    # Order Mendeliome genes by number of structures
    ord.mendeliome <- mendgenes[order(-mendgenes$tot_structures),]
    ord.mendeliome <- unique(ord.mendeliome)
    write.csv(ord.mendeliome, file=paste(directory, replen, "bp_", parameters, "_ord_mendeliome.csv", sep=""), row.names = F, quote = F)
    
}
rm(file, files, replen, genes, numgenes, mendgenes, ord.genes, ord.mendeliome)


# ********************************
# ---- Load by lengths output ----
# ********************************

# FOR FUNCTION INPUT

# files <- list.files(path="../Output/Genes/", pattern="^all_")
files <- list.files(path="../Output/Genes/All_split")
for (file in files) {
    print(file)
    replen <- as.numeric(strsplit(strsplit(file, "_") [[1]][2], "bp")[[1]][1])
    genes <- as.data.frame(read.table(paste("../Output/Genes/All_split/", file, sep=""), sep="\t", stringsAsFactors = F))
    colnames(genes) <- c("chromosome", "repeat_length", "repeat_seq", "start", "end", "distance", "ensembl_gene_id", 
                         "hgnc_symbol", "gene_range", "ensembl_transcript_id", "transcript_range", "in_ex", 
                         "number", "feature")
    genes$intron_exon <- paste0(genes$in_ex, genes$number)
    genes <- genes[,c(1:11, 15, 14)]
    
    namet <- strsplit(file, "\\.")[[1]][1]
    name <- substr(namet, 5, nchar(namet))
    
    saveRDS(genes, file=paste0("../Output/repeats/", name, ".rds"))
}
rm(file, files, replen, genes, numgenes, mendgenes, ord.genes, ord.mendeliome)


# ********************************
# ---- REPEAT LENGTH OVERVIEW----
# ********************************

# ---- ~ Instances summary ----
summary <- as.data.frame(read.table(file="../Output/expandSeq/summary_results_7-20bp_1000bp.tsv", 
                                    sep="\t", header = T, stringsAsFactors = F))
summary.num <- summary[,1:6]

files <- paste0("../Output/Analysis_fixed/", list.files(path="../Output/Analysis_fixed/", 
                                                        pattern=paste0("^", length,"bp_dist", dist)))
numgenes <- readRDS(files[1])
mendgenes <- readRDS(files[2])

repeatLen <- data.frame(length=7:20, tot_instances=NA, mend_instances=NA, tot_structures=NA, mend_structures=NA)

files <-  list.files(path="./Repeats/Input/", pattern="dist1000_chr")

for (file in files) {
    data <- readRDS(paste0("./Repeats/Input/", file))
    filename <- 
    write.table(data, file=paste0("../Output/repeats_fix/", strsplit(file, "\\.")[[1]][1], ".tsv"), 
        sep="\t", quote=F, row.names = F, col.names = F)
}


for (length in 8) {
    files <- paste0("./Repeats/Input/", list.files(path="./Repeats/Input/", 
                                                   pattern=paste0("^",length,"bp_dist", dist, "_chr")))
    files <- files[c(1,12,16:22,2:11,13:15,23:24)]
    numgenes <- do.call(rbind,
                        lapply(files, readRDS))

    # Tot instances
    repeatLen[length-6,2] <- nrow(unique(numgenes[,c(1:5, 7:12)]))


    mendgenes <- numgenes[which(numgenes$hgnc_symbol %in% mendeliome),]
    rm(numgenes)
    mendgenes <- unique(mendgenes[,c(1:5, 7:12)])
    
    # Mend instances
    repeatLen[length-6,3] <- nrow(unique(mendgenes))

    files <- paste0("../Output/Analysis_fixed/", list.files(path="../Output/Analysis_fixed/",
                                                            pattern=paste0("^", length,"bp_dist", dist)))
    numgenes <- readRDS(files[1])
    mendgenes <- readRDS(files[2])

    # Tot structures
    repeatLen[length-6,4] <- sum(numgenes$tot_structures)

    # Mend structures
    repeatLen[length-6,5] <- sum(mendgenes$tot_structures)
}


for (i in 7:20) {
    repeatLen[i-6, 2] <- sum(summary.num$abs_freq[which(summary.num$repeat_length==i)], na.rm=T)
    repeatLen[i-6, 3] <- sum(summary.num$mend_freq[which(summary.num$repeat_length==i)], na.rm=T)
}


write.table(repeatLen, file="../Output/summary_structures_by_length_fixed.txt", sep="\t", col.names = T, row.names = F, quote = F)
repeatLen <- data.frame(read.table("../Output/summary_structures_by_length_fixed.txt", sep="\t", header = T), stringsAsFactors = F)

# ---- ~ Structures/gene ----
# ---- ~~ Genome ----
# 35584 = max number of genes from 7bp dataset
tot.struct <- data.frame(order=c(1:35584), struct_7bp=NA, struct_8bp=NA, struct_9bp=NA, struct_10bp=NA, struct_11bp=NA, struct_12bp=NA, 
                       struct_13bp=NA, struct_14bp=NA, struct_15bp=NA, struct_16bp=NA, struct_17bp=NA, struct_18bp=NA, struct_19bp=NA, struct_20bp=NA)

for (length in 7:20) {
    files <- paste0(directory, list.files(path="../Output/Analysis_fixed/", pattern=paste0("^", as.character(length), "bp_")))
    ord.genes <- as.data.frame(read.csv(files[3], header = T, stringsAsFactors = F))
    
    tot.struct[which(tot.struct$order<=nrow(ord.genes)),length-5] <- ord.genes$tot_structures
    
}

tot.struct.melt <- melt(tot.struct, id.vars="order")
colnames(tot.struct.melt) <- c("order", "repeat_length", "structures")

struct.ymax <- 500*ceiling(max(tot.struct[,c(2:15)], na.rm=T)/500)

# ---- ~~ Norm genome ----
# 35584 = max number of genes from 7bp dataset
norm.struct <- data.frame(order=c(1:35584), struct_7bp=NA, struct_8bp=NA, struct_9bp=NA, struct_10bp=NA, struct_11bp=NA, struct_12bp=NA, 
                       struct_13bp=NA, struct_14bp=NA, struct_15bp=NA, struct_16bp=NA, struct_17bp=NA, struct_18bp=NA, struct_19bp=NA, struct_20bp=NA)

for (length in 7:20) {
    files <- paste0(directory, list.files(path="../Output/Analysis_fixed/", pattern=paste0("^", as.character(length), "bp_")))
    ord.genes <- as.data.frame(read.csv(files[3], header = T, stringsAsFactors = F))
    ord.normgenes <- ord.genes[order(-ord.genes$norm_struct_Mb),]
    # ord.normmend <- ord.mendeliome[order(-ord.mendeliome$norm_struct_Mb),]
    norm.struct[which(norm.struct$order<=nrow(ord.normgenes)),length-5] <- ord.normgenes$norm_struct_Mb
    
}

norm.struct.melt <- melt(tot.struct, id.vars="order")
colnames(norm.struct.melt) <- c("order", "repeat_length", "norm_struct")


# ---- ~~ Mendeliome ----

# 35584 = max number of genes from 7bp dataset
mend.struct <- data.frame(order=c(1:3874), struct_7bp=NA, struct_8bp=NA, struct_9bp=NA, struct_10bp=NA, struct_11bp=NA, struct_12bp=NA, 
                        struct_13bp=NA, struct_14bp=NA, struct_15bp=NA, struct_16bp=NA, struct_17bp=NA, struct_18bp=NA, struct_19bp=NA, struct_20bp=NA)

for (length in 7:20) {
    files <- paste0(directory, list.files(path="../Output/Analysis_fixed/", pattern=paste0("^", as.character(length), "bp_")))
    ord.mendeliome <- as.data.frame(read.csv(files[4], header = T, stringsAsFactors = F))
    
    mend.struct[which(mend.struct$order<=nrow(ord.mendeliome)),length-5] <- ord.mendeliome$tot_structures
    
}

mend.struct.melt <- melt(mend.struct, id.vars="order")
colnames(mend.struct.melt) <- c("order", "repeat_length", "structures")


# ---- ~~ Norm Mendeliome ----
# 35584 = max number of genes from 7bp dataset
normmend.struct <- data.frame(order=c(1:3874), struct_7bp=NA, struct_8bp=NA, struct_9bp=NA, struct_10bp=NA, struct_11bp=NA, struct_12bp=NA, 
                        struct_13bp=NA, struct_14bp=NA, struct_15bp=NA, struct_16bp=NA, struct_17bp=NA, struct_18bp=NA, struct_19bp=NA, struct_20bp=NA)

for (length in 7:20) {
    files <- paste0(directory, list.files(path="../Output/Analysis_fixed/", pattern=paste0("^", as.character(length), "bp_")))
    ord.mendeliome <- as.data.frame(read.csv(files[4], header = T, stringsAsFactors = F))
    ord.normmend <- ord.mendeliome[order(-ord.mendeliome$norm_struct_Mb),]
    
    normmend.struct[which(normmend.struct$order<=nrow(ord.normmend)),length-5] <- ord.normmend$norm_struct_Mb
    
}

normmend.struct.melt <- melt(tot.struct, id.vars="order")
colnames(normmend.struct.melt) <- c("order", "repeat_length", "norm_struct")


# ---- ~ Plot ----
# Set y-axis maxima
struct.ymax <- 500*ceiling(max(tot.struct[,c(2:15)], na.rm=T)/500)
mend.ymax <- 500*ceiling(max(mend.struct[,c(2:15)], na.rm=T)/500)
norm.ymax <- 500*ceiling(max(norm.struct[,c(2:15)], na.rm=T)/500)
normmend.ymax <- 500*ceiling(max(normmend.struct[,c(2:15)], na.rm=T)/500)
structures.ymax <- 500*ceiling(max(repeatLen$tot_structures, na.rm=T)/500)
structmend.ymax <- 500*ceiling(max(repeatLen$mend_structures, na.rm=T)/500)



pdf(paste("./Figures/", filename, "_structures_per_length_log.pdf", sep=""), h=8, w=15)
# Total structures/repeat length
print(ggplot(repeatLen, aes(x = length, y=tot_structures)) +
          geom_bar(stat="identity", fill="#0084ff") +
          ggtitle(paste0("Structures per repeat length (dist=", dist, "bp)")) +
          xlab("Repeat length (bp)") +
          ylab("Structures (log)") +
          # coord_cartesian(xlim=c(6.5,20.5), ylim=c(0,structures.ymax), expand=FALSE) +
          scale_x_continuous(breaks=c(7:20), expand=c(0,0)) +
          scale_y_log10(expand=c(0,0)) + 
          theme_classic() +
          theme(text = element_text(size=20)))
dev.off()

# Total Mendeliome structures/repeat length
print(ggplot(repeatLen, aes(x = length, y=mend_structures)) +
          geom_bar(stat="identity", fill="forestgreen") +
          ggtitle(paste0("Mendeliome structures per repeat length (dist=", dist, "bp)")) +
          xlab("Repeat length (bp)") +
          ylab("Structures") +
          coord_cartesian(xlim=c(6.5,20.5), ylim=c(0,structmend.ymax), expand=FALSE) +
          scale_x_continuous(breaks=c(7:20)) +
          theme_bw() +
          theme(text = element_text(size=15)))
# print(ggplot(struct, aes(x = length, y=tot_structs)) +
#           geom_bar(stat="identity") +
#           ggtitle(paste0("Structures per repeat length (dist=", dist, "bp)")) +
#           xlab("Repeat length (bp)") +
#           ylab("Structures") +
#           coord_cartesian(xlim=c(6.5,20.5), ylim=c(0,structures.ymax), expand=FALSE) +
#           scale_y_continuous(limits=c(0, structures.ymax), breaks=seq(0, structures.ymax, by=1e06)) +
#           scale_x_continuous(breaks=c(7:20)) +
#           theme_bw() +
#           theme(text = element_text(size=15)))
# Line plot total structures/gene
print(ggplot(na.omit(tot.struct.melt), aes(x = order, y = as.numeric(structures))) +
          geom_line(aes(color=repeat_length), size=1) +
          ggtitle(paste0("Structures per gene (dist ", dist, "bp)")) +
          xlab("Genes") +
          ylab("structures") +
          scale_color_discrete(name = "Repeat length",
                               labels = c("7bp", "8bp", "9bp", "10bp", "11bp", "12bp", "13bp", "14bp", 
                                          "15bp", "16bp", "17bp", "18bp", "19bp", "20bp")) +
          # coord_cartesian(ylim=c(0.11, struct.ymax), expand=FALSE) +
          scale_x_continuous(expand = c(0, 0)) +
          scale_y_log10() + 
          theme_bw() +
          theme(text = element_text(size=15),
                axis.text.x = element_blank(),
                axis.ticks = element_blank()))

# Line plot normalized structures/gene
print(ggplot(na.omit(norm.struct.melt), aes(x = order, y = norm_struct)) +
          geom_line(aes(color=repeat_length), size=1) +
          ggtitle(paste0("Normalized structures per gene (dist ", dist, "bp)")) +
          xlab("Genes") +
          ylab("structures per Mb") +
          # coord_cartesian(ylim=c(0, norm.ymax), expand=FALSE) +
          scale_color_discrete(name = "Repeat length",
                               labels = c("7bp", "8bp", "9bp", "10bp", "11bp", "12bp", "13bp", "14bp", 
                                          "15bp", "16bp", "17bp", "18bp", "19bp", "20bp")) +
          scale_x_continuous(expand = c(0, 0)) +
          scale_y_log10() +
          theme_bw() +
          theme(text = element_text(size=15),
                axis.text.x = element_blank(),
                axis.ticks = element_blank()))

# Line plot total structures/Mendeliome gene
print(ggplot(na.omit(mend.struct.melt), aes(x = order, y = structures)) +
          geom_line(aes(color=repeat_length), size=1) +
          ggtitle(paste0("Structures per Mendeliome gene (dist ", dist, "bp)")) +
          xlab("Genes") +
          ylab("structures") +
          # coord_cartesian(ylim=c(0, mend.ymax), expand=FALSE) +
          scale_color_discrete(name = "Repeat length",
                               labels = c("7bp", "8bp", "9bp", "10bp", "11bp", "12bp", "13bp", "14bp", 
                                          "15bp", "16bp", "17bp", "18bp", "19bp", "20bp")) +
          scale_x_continuous(expand = c(0, 0)) +
          scale_y_log10() +
          theme_bw() +
          theme(text = element_text(size=15),
                axis.text.x = element_blank(),
                axis.ticks = element_blank()))

# Line plot normalized structures/Mendeliome gene
print(ggplot(na.omit(normmend.struct.melt), aes(x = order, y = norm_struct)) +
          geom_line(aes(color=repeat_length), size=1) +
          ggtitle(paste0("Normalized structures per Mendeliome gene (dist ", dist, "bp)")) +
          xlab("Genes") +
          ylab("structures per Mb") +
          # coord_cartesian(ylim=c(0, normmend.ymax), expand=FALSE) +
          scale_color_discrete(name = "Repeat length",
                               labels = c("7bp", "8bp", "9bp", "10bp", "11bp", "12bp", "13bp", "14bp", 
                                          "15bp", "16bp", "17bp", "18bp", "19bp", "20bp")) +
          scale_x_continuous(expand = c(0, 0)) +
          scale_y_log10() +
          theme_bw() +
          theme(text = element_text(size=15),
                axis.text.x = element_blank(),
                axis.ticks = element_blank()))

dev.off()



# ---- LENGTHS LOOP ----
# Input = repeat length, output = plots
for (length in 7:20) {
    # for (length in 8) {
    files <- paste0(directory, list.files(path="../Output/Analysis_fixed/", pattern=paste0("^", as.character(length), "bp_")))
    numgenes <- readRDS(files[1])
    mendgenes <- readRDS(files[2])
    ord.genes <- as.data.frame(read.csv(files[3], header = T, stringsAsFactors = F))
    dist.order <- ord.genes[!is.na(ord.genes$avg_dist),]
    dist.order <- dist.order[order(-dist.order$avg_dist),]
    dist.order$hgnc_symbol <- reorder(dist.order$hgnc_symbol, dist.order$avg_dist)
    
    ord.mendeliome <- as.data.frame(read.csv(files[4], header = T, stringsAsFactors = F))
    mend.dist.order <- ord.mendeliome[!is.na(ord.mendeliome$avg_dist),]
    mend.dist.order <- mend.dist.order[order(-mend.dist.order$avg_dist),]
    
    ord.normgenes <- ord.genes[order(-ord.genes$norm_struct_Mb),]
    ord.normmend <- ord.mendeliome[order(-ord.mendeliome$norm_struct_Mb),]
    
    filename <- paste0(today, "_", parameters, "_", length, "bp")
    
    # ---- ~ Overview/gene plots ----
    struct.ymax <- 500*ceiling(max(ord.genes$tot_structures)/500)
    mend.ymax <- 500*ceiling(max(ord.mendeliome$tot_structures)/500)
    norm.ymax <- 500*ceiling(max(ord.normgenes$norm_struct_Mb)/500)
    normmend.ymax <- 500*ceiling(max(ord.normmend$norm_struct_Mb)/500)
    
    # Plotting 
    pdf(paste("./Figures/", filename, "_gene.pdf", sep=""), h=8, w=14)
    # Structures 
    print(ggplot(ord.genes, aes(x = reorder(ensembl_gene_id, -tot_structures), y = tot_structures)) +
              geom_bar(stat="identity") +
              ggtitle(paste0("Structures per gene (dist ", dist, "bp, repeat ", length, "bp)")) +
              xlab("Genes") +
              ylab("Structures") +
              coord_cartesian(ylim=c(0, struct.ymax), expand=FALSE) +
              theme_classic() +
              theme(text = element_text(size=15),
                    axis.text.x = element_blank(),
                    axis.ticks = element_blank()))
    
    # Structures top 50
    print(ggplot(ord.genes[1:50,], aes(x = reorder(hgnc_symbol, -tot_structures), y = tot_structures)) +
              geom_bar(stat="identity") +
              ggtitle(paste0("Structures per gene (top 50; dist ", dist, "bp, repeat ", length, "bp)")) +
              xlab("") +
              ylab("Structures") +
              coord_cartesian(ylim=c(0, struct.ymax), expand=FALSE) +
              theme_bw() +
              theme(text = element_text(size=15),
                    axis.text.x = element_text(angle = 90, hjust=0.95, vjust=0.2),
                    axis.ticks = element_blank()))
    
    # Normalized (by Mb)
    print(ggplot(ord.normgenes, aes(x = reorder(ensembl_gene_id, -norm_struct_Mb), y = norm_struct_Mb)) +
              geom_bar(stat="identity", fill="slategray") +
              ggtitle(paste0("Normalized structures per gene (dist ", dist, "bp, repeat ", length, "bp)")) +
              xlab("Genes") +
              ylab("Structures per Mb") +
              coord_cartesian(ylim=c(0, norm.ymax), expand=FALSE) +
              scale_y_continuous(labels = function(x) format(x, scientific=T)) +
              theme_classic() +
              theme(text = element_text(size=15),
                    axis.text.x = element_blank(),
                    axis.ticks = element_blank()))
    
    # Normalized (by Mb) top 50
    print(ggplot(ord.normgenes[1:50,], aes(x = reorder(hgnc_symbol, -norm_struct_Mb), y = norm_struct_Mb)) +
              geom_bar(stat="identity", fill="slategray") +
              ggtitle(paste0("Normalized structures per gene (top 50; dist ", dist, "bp, repeat ", length, "bp)")) +
              xlab("") +
              ylab("Structures per Mb") +
              coord_cartesian(ylim=c(0, norm.ymax), expand=FALSE) +
              scale_y_continuous(labels = function(x) format(x, scientific=T)) +
              theme_bw() +
              theme(text = element_text(size=15),
                    axis.text.x = element_text(angle = 90, hjust=0.95, vjust=0.2),
                    axis.ticks = element_blank()))
    
    # Avg distance
    print(ggplot(dist.order, aes(x = reorder(ensembl_gene_id, -avg_dist), y = avg_dist)) +
              geom_bar(stat="identity") +
              ggtitle(paste0("Average distance per gene (dist ", dist, "bp, repeat ", length, "bp)")) +
              xlab("Genes") +
              ylab("Average distance") +
              coord_cartesian(ylim=c(0, dist), expand=FALSE) +
              scale_y_continuous(breaks = seq(0, dist, by = 50)) +
              theme_classic() +
              # theme(axis.text.x = element_text(angle = 90, hjust=0.95, vjust=0.2),
              theme(text = element_text(size=15),
                    axis.text.x = element_blank(),
                    axis.ticks = element_blank()))
    
    # Avg distance top 50
    print(ggplot(dist.order[1:50,], aes(x = reorder(hgnc_symbol, -avg_dist), y = avg_dist)) +
              geom_bar(stat="identity") +
              ggtitle(paste0("Average distance per gene (top 50; dist ", dist, "bp, repeat ", length, "bp)")) +
              xlab("") +
              ylab("Average distance") +
              coord_cartesian(ylim=c(0, dist), expand=FALSE) +
              theme_bw() +
              # theme(axis.text.x = element_text(angle = 90, hjust=0.95, vjust=0.2),
              theme(text = element_text(size=15),
                    axis.text.x = element_text(angle = 90, hjust=0.95, vjust=0.2),
                    axis.ticks = element_blank()))
    
    # ---- ~ Mendeliome ----
    # Structures Mendeliome
    print(ggplot(ord.mendeliome, aes(x = reorder(ensembl_gene_id, -tot_structures), y = tot_structures)) +
              geom_bar(stat="identity", fill="forestgreen") +
              ggtitle(paste0("Structures per Mendeliome gene (dist ", dist, "bp, repeat ", length, "bp)")) +
              xlab("Genes") +
              ylab("Structures") +
              coord_cartesian(ylim=c(0, mend.ymax), expand=FALSE) +
              theme_classic() +
              theme(text = element_text(size=15),
                    axis.text.x = element_blank(),
                    axis.ticks = element_blank()))
    
    # Structures Mendeliome top 50
    print(ggplot(ord.mendeliome[1:50,], aes(x = reorder(hgnc_symbol, -tot_structures), y = tot_structures)) +
              # print(ggplot(ord.mendeliome[1:50,], aes(x = hgnc_symbol, y = tot_structures)) +
              geom_bar(stat="identity", fill="forestgreen") +
              ggtitle(paste0("Structures per Mendeliome gene (top 50; dist ", dist, "bp, repeat ", length, "bp)")) +
              xlab("") +
              ylab("Structures") +
              coord_cartesian(ylim=c(0, mend.ymax), expand=FALSE) +
              theme_bw() +
              theme(text = element_text(size=15),
                    axis.text.x = element_text(angle = 90, hjust=0.95, vjust=0.2),
                    axis.ticks = element_blank()))
    
    # Normalized (by Mb) Mendeliome
    print(ggplot(ord.normmend, aes(x = reorder(ensembl_gene_id, -norm_struct_Mb), y = norm_struct_Mb)) +
              geom_bar(stat="identity", fill="seagreen") +
              ggtitle(paste0("Normalized structures per Mendeliome gene (dist ", dist, "bp, repeat ", length, "bp)")) +
              xlab("Genes") +
              ylab("Structures per Mb") +
              coord_cartesian(ylim=c(0, normmend.ymax), expand=FALSE) +
              scale_y_continuous(labels = function(x) format(x, scientific=T)) +
              theme_classic() +
              theme(text = element_text(size=15),
                    axis.text.x = element_blank(),
                    axis.ticks = element_blank()))
    
    # Normalized (by Mb) Mendeliome top 50
    print(ggplot(ord.normmend[1:50,], aes(x = reorder(hgnc_symbol, -norm_struct_Mb), y = norm_struct_Mb)) +
              geom_bar(stat="identity", fill="seagreen") +
              ggtitle(paste0("Normalized structures per Mendeliome gene (top 50; dist ", dist, "bp, repeat ", length, "bp)")) +
              xlab("") +
              ylab("Structures per Mb") +
              coord_cartesian(ylim=c(0, normmend.ymax), expand=FALSE) +
              scale_y_continuous(labels = function(x) format(x, scientific=T)) +
              theme_bw() +
              theme(text = element_text(size=15),
                    axis.text.x = element_text(angle = 90, hjust=0.95, vjust=0.2),
                    axis.ticks = element_blank()))
    
    # Avg distance Mendeliome
    print(ggplot(mend.dist.order, aes(x = reorder(ensembl_gene_id, -avg_dist), y = avg_dist)) +
              geom_bar(stat="identity", fill="forestgreen") +
              ggtitle(paste0("Average distance per Mendeliome gene (dist ", dist, "bp, repeat ", length, "bp)")) +
              xlab("Genes") +
              ylab("Average distance") +
              coord_cartesian(ylim=c(0, dist), expand=FALSE) +
              theme_classic() +
              # theme(axis.text.x = element_text(angle = 90, hjust=0.95, vjust=0.2),
              theme(text = element_text(size=15),
                    axis.text.x = element_blank(),
                    axis.ticks = element_blank()))
    
    # Avg distance Mendeliome top 50
    print(ggplot(mend.dist.order[1:50,], aes(x = reorder(hgnc_symbol, -avg_dist), y = avg_dist)) +
              geom_bar(stat="identity", fill="forestgreen") +
              ggtitle(paste0("Average distance per Mendeliome gene (top 50; dist ", dist, "bp, repeat ", length, "bp)")) +
              xlab("") +
              ylab("Average distance") +
              coord_cartesian(ylim=c(0, dist), expand=FALSE) +
              theme_bw() +
              # theme(axis.text.x = element_text(angle = 90, hjust=0.95, vjust=0.2),
              theme(text = element_text(size=15),
                    axis.text.x = element_text(angle = 90, hjust=0.95, vjust=0.2),
                    axis.ticks = element_blank()))
    
    dev.off()
}



# ---- MINDIST INCREMENTS ----

struct <- data.frame(read.table(file="../Output/summary_structures_by_length_mindists_incrementsFIN.txt", sep="\t", header = T, stringsAsFactors = F))
# struct <- struct[3:nrow(struct),]
for (i in 1:ncol(struct)){
    struct[i] <- as.numeric(unlist(struct[i]))
}
colnames(struct) <- c("length", seq(0, 900, 50))


# ---- ~ Barplot per length ----
pdf(paste("./Figures/", filename, "_structures_per_length_mindists_increments.pdf", sep=""), h=8, w=14)
for (i in 2:ncol(struct)) {
    mindist <- as.numeric(colnames(struct)[i])
    print(ggplot(struct, aes(x = length, y=struct[,i])) +
         geom_bar(stat="identity") +
         ggtitle(paste0("Structures per repeat length (dist = ", mindist, "-", dist, " bp)")) + 
         xlab("Repeat length (bp)") +
         ylab("Structures") +
         coord_cartesian(xlim=c(6.5,20.5), ylim=c(0,500*ceiling(max(struct[,i], na.rm=T)/500)), expand=FALSE) +
         scale_y_continuous(limits = c(0, 500*ceiling(max(struct[,i], na.rm=T)/500)),
                            breaks = seq(0, 500*ceiling(max(struct[,i], na.rm=T)/500), by=1e06)) + 
         scale_x_continuous(breaks=c(7:20)) +
         theme_bw() +
         theme(text = element_text(size=15)))
}
dev.off()

pdf(paste("./Figures/", filename, "_structures_per_length_small.pdf", sep=""), h=8, w=8)
for (i in 2) {
    mindist <- as.numeric(colnames(struct)[i])
    print(ggplot(struct, aes(x = length, y=struct[,i])) +
              geom_bar(stat="identity") +
              # ggtitle(paste0("Structures per repeat length (dist = ", mindist, "-", dist, " bp)")) + 
              xlab("Repeat length (bp)") +
              ylab("Structures") +
              coord_cartesian(xlim=c(6.5,20.5), ylim=c(0,500*ceiling(max(as.numeric(struct[,i]), na.rm=T)/500)), expand=FALSE) +
              # scale_y_log10() +
              # scale_y_continuous(limits = c(0, 500*ceiling(max(struct[,i], na.rm=T)/500))
                                 # breaks = seq(0, 500*ceiling(max(struct[,i], na.rm=T)/500), by=1e06)) + 
              scale_x_continuous(breaks=c(7:20)) +
              theme_bw() +
              theme(text = element_text(size=30)))
}
dev.off()


# ~ Line plot per mindist -------------------------------------------------

struct.melt <- melt(struct, id = "length", value.name = "structures")
colnames(struct.melt)[2] <- "mindist"
struct.melt$length <- paste0(struct.melt$length, "bp")
struct.melt$length <- factor(struct.melt$length, levels=c("7bp", "8bp", "9bp", "10bp", "11bp", "12bp", "13bp", "14bp",
                                                          "15bp", "16bp", "17bp", "18bp", "19bp", "20bp"))

pdf(paste("./Figures/", filename, "_structures_per_length_mindists_increments_overview.pdf", sep=""), h=8, w=15)
# Overview per repeat length and mindist LINEAR
print(ggplot(struct.melt, aes(x = mindist, y = structures, group = length)) +
          geom_line(aes(color=length), size=1) +
          geom_point(aes(color=length)) +
          ggtitle(paste0("Structures per minimum distance (dist ", dist, "bp)")) +
          xlab("Min distance (bp)") +
          ylab("Structures") +
          scale_y_continuous(limits = c(0, 500*ceiling(max(struct.melt$structures, na.rm=T)/500))) +
          scale_color_discrete(name = "Repeat length") +
          theme_bw() +
          theme(text = element_text(size=15)))
# Overview per repeat length and mindist LOG
print(ggplot(struct.melt, aes(x = mindist, y = structures, group = length)) +
          geom_line(aes(color=length), size=1) +
          geom_point(aes(color=length)) +
          ggtitle(paste0("Structures per minimum distance (dist ", dist, "bp)")) +
          xlab("Min distance (bp)") +
          ylab("Structures") +
          scale_y_log10() +
          scale_color_discrete(name = "Repeat\nlength") +
          theme_bw() +
          theme(text = element_text(size=15)))

# Per mindist, separate plot per repeat length
for (i in struct.melt$length[1:14]) {
    data <- struct.melt[which(struct.melt$length==i),]
    print(ggplot(data, aes(x = mindist, y = structures)) +
              # geom_line(aes(color=length), size=1, group=1) +
              geom_line(size=1, group=1, linetype="dotted") +
              geom_point() +
              ggtitle(paste0("Structures per minimum distance (repeat ", i, ", dist ", dist, "bp)")) +
              xlab("Min distance (bp)") +
              ylab("Structures") +
              coord_cartesian(ylim=c(0, 500*ceiling(max(data$structures, na.rm=T)/500))) +
              theme_bw() +
              theme(text = element_text(size=15)))
}
dev.off()





# BY CHROMOSOMES ----------------------------------------------------------

# Number of genes per chromosomes
# Source: https://www.ncbi.nlm.nih.gov/books/NBK22266/
geneschr <- c(3000, 2500, 1900, 1600, 1700, 1900, 1800, 1400, 1400, 1400, 
              2000, 1600, 800, 1200, 1200, 1300, 1600, 600, 1700, 900, 
              400, 800, 1400, 200)

genes.chr <- data.frame(chromosome=chr_label, genes=geneschr, 
                        seven=NA, eight=NA, nine=NA, ten=NA, eleven=NA, twelve=NA,
                        thirteen=NA, fourteen=NA, fifteen=NA, sixteen=NA, seventeen=NA,
                        eighteen=NA, nineteen=NA, twenty=NA)

struct.chr <- data.frame(chromosome=chr_label, genes=geneschr, 
                         seven=NA, eight=NA, nine=NA, ten=NA, eleven=NA, twelve=NA,
                         thirteen=NA, fourteen=NA, fifteen=NA, sixteen=NA, seventeen=NA,
                         eighteen=NA, nineteen=NA, twenty=NA)


for (length in 7:20) {
    
    files <- paste0("./Repeats/Input/", list.files(path="./Repeats/Input/", 
                                                   pattern=paste0("^",length,"bp_dist", dist, "_chr")))
    files <- files[c(1,12,16:22,2:11,13:15,23:24)]
    
    for (i in 1:24) {
        data <- readRDS(files[i])
        genes.chr[i, length-4] <- length(unique(data$ensembl_gene_id))
        struct.chr[i, length-4] <- nrow(data)/2
    }
    
}

write.table(genes.chr, file="../Output/genes_per_chromosomes.txt", sep="\t", 
            col.names = T, row.names = F, quote = F)

write.table(struct.chr, file="../Output/structures_per_chromosomes.txt", 
            sep="\t", col.names = T, row.names = F, quote = F)


# ---- ~ Genes per chr ----
genes.chr$chromosome <- factor(genes.chr$chromosome, levels = chr_label)

# ---- ~~ Plot ----
lengths <- c(7:20)

pdf(paste("./Figures/", filename, "_genes_per_chromosome.pdf", sep=""), h=8, w=14)
for (i in 3:ncol(genes.chr)) {
    maxlim <- 500*ceiling(max(genes.chr[,i], na.rm=T)/500)
    if (maxlim <- 3000) {
        maxlim <- 3100
    }
    
    print(ggplot(genes.chr, aes(x = chromosome, y=genes.chr[,i])) +
              geom_bar(stat="identity", fill="gray") +
              geom_line(aes(x = as.numeric(chromosome), y = genes)) +
              ggtitle(paste0("Genes per chromosome (dist = ", dist, "bp, repeat = ", lengths[i-2], "bp)")) + 
              xlab("") +
              ylab("Number of genes") +
              coord_cartesian(xlim=c(0.5,24.5), ylim=c(0,maxlim), expand=FALSE) +
              # scale_y_continuous(limits = c(0, 500*ceiling(max(genes.chr[,i], na.rm=T)/500)),
              #                    breaks = seq(0, 500*ceiling(max(genes.chr[,i], na.rm=T)/500), by=500)) +
              # scale_x_continuous(breaks=c(7:20)) +
              theme_bw() +
              theme(text = element_text(size=15)))
}
dev.off()


# ---- ~ Structures per chr ----
# Source: https://www.ncbi.nlm.nih.gov/books/NBK22266/
struct.chr <- data.frame(read.table(file="../Output/structures_per_chromosomes.txt", header = T), stringsAsFactors = F)
struct.chr$chromosome <- factor(struct.chr$chromosome, levels = chr_label)

# ---- ~~ Plot ----

lengths <- c(7:20)

pdf(paste("./Figures/", filename, "_structures_per_chromosome.pdf", sep=""), h=8, w=14)
for (i in 3:ncol(struct.chr)) {
    maxlim <- 500*ceiling(max(struct.chr[,i], na.rm=T)/500)
    
    if (maxlim < 3000) {
        maxlim <- 3100
    }
    
    print(ggplot(struct.chr, aes(x = chromosome, y=struct.chr[,i])) +
              geom_bar(stat="identity", fill="#0084ff") +
              # geom_line(aes(x = as.numeric(chromosome), y = genes)) +
              ggtitle(paste0("Structures per chromosome (dist = ", dist, "bp, repeat = ", lengths[i-2], "bp)")) + 
              xlab("") +
              ylab("Number of structures") +
              coord_cartesian(xlim=c(0.5,24.5), ylim=c(0,maxlim), expand=FALSE) +
              # scale_y_continuous(limits = c(0, 500*ceiling(max(struct.chr[,i], na.rm=T)/500)),
              #                    breaks = seq(0, 500*ceiling(max(struct.chr[,i], na.rm=T)/500), by=500)) +
              # scale_x_continuous(breaks=c(7:20)) +
              theme_classic() +
              theme(text = element_text(size=15)))
}
dev.off()


# ---- MISC ----
# Normalized
ord.normgenes <- ord.genes[order(-ord.genes$norm_inst_bp),]
norm.ymax <- ceiling(max(ord.normgenes$norm_inst_bp))

pdf(paste("./Figures/", filename, "_instances_normgene.pdf", sep=""), h=8, w=15)
ggplot(ord.normgenes, aes(x = reorder(hgnc_symbol, -norm_inst_bp), y = norm_inst_bp)) +
    geom_bar(stat="identity", fill="steelblue") +
    ggtitle(paste0("Normalized instances per gene (dist ", dist, "bp, repeat ", length, "bp)")) + 
    ggtitle("Normalized instances per gene") + 
    xlab("") +
    ylab("Instances normalized by gene length") +
    coord_cartesian(ylim=c(0,norm.ymax), expand=FALSE) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank())
dev.off()

pdf(paste("./Figures/", filename, "_instances_normgene_50.pdf", sep=""), h=8, w=15)
ggplot(ord.normgenes[1:50,], aes(x = reorder(hgnc_symbol, -norm_inst_bp), y = norm_inst_bp)) +
    geom_bar(stat="identity", fill="steelblue") +
    ggtitle(paste0("Normalized instances per gene (top 50; dist ", dist, "bp, repeat ", length, "bp)")) + 
    xlab("") +
    ylab("Instances normalized by gene length") +
    coord_cartesian(ylim=c(0,norm.ymax), expand=FALSE) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust=0.95, vjust=0.2),
          axis.ticks = element_blank())
dev.off()

# Mendeliome
pdf(paste("./Figures/", filename, "_instances_gene_mendeliome.pdf", sep=""), h=8, w=15)
ggplot(ord.mendeliome, aes(x = reorder(hgnc_symbol, -tot_instances), y = tot_instances)) +
    geom_bar(stat="identity", fill="forestgreen") +
    ggtitle(paste0("Instances per Mendeliome gene (dist ", dist, "bp, repeat ", length, "bp)")) + 
    ggtitle("Instances per Mendeliome gene") + 
    xlab("") +
    ylab("Instances") +
    coord_cartesian(ylim=c(0,plot.ymax), expand=FALSE) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank())
dev.off()

pdf(paste("./Figures/", filename, "_instances_gene_mendeliome_50.pdf", sep=""), h=8, w=15)
ggplot(ord.mendeliome[1:50,], aes(x = reorder(hgnc_symbol, -tot_instances), y = tot_instances)) +
    geom_bar(stat="identity", fill="forestgreen") +
    ggtitle(paste0("Instances per Mendeliome gene (top 50; dist ", dist, "bp, repeat ", length, "bp)")) + 
    xlab("") +
    ylab("Instances") +
    coord_cartesian(ylim=c(0,plot.ymax), expand=FALSE) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust=0.95, vjust=0.2),
          axis.ticks = element_blank())
dev.off()


pdf(paste(directory, "Plots/", filename, "_avgdist_mendeliome_50.pdf", sep=""), h=8, w=15)
ggplot(mend.dist.order[1:50,], aes(x = reorder(hgnc_symbol, -avg_dist), y = avg_dist)) +
    geom_bar(stat="identity", fill="forestgreen") +
    ggtitle("Average distance per Mendeliome gene (top 50)") + 
    xlab("") +
    ylab("Average distance") +
    coord_cartesian(ylim=c(0,800), expand=FALSE) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust=0.95, vjust=0.2),
          axis.ticks = element_blank())
dev.off()

write.csv(ord.mendeliome, file=paste(directory, "Output/", filename, "_ord_mendeliome_genes.csv", sep=""), row.names = F, quote = F)