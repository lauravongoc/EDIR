

#' Look Up a Gene in EDIR Dataset
#'
#' This function searches for a specified gene in the EDIR dataset. A gene name
#' is a required parameter.
#'
#' Summary of results printed to console includes gene name, gene length (bp),
#' Ensembl transcript ID, queried distance between repeats (default: 0-1000 bp),
#' and an overview of total results for the given repeat length. Console outputs
#' include runtime.
#'
#' @param gene The gene name (ENSEMBL ID or HGNC symbol)
#' @param length Repeat sequence length, must be between 7 and 20. Defaults to
#' NA. If NA, results will include all available lengths in dataset for queried
#' gene.
#' @param mindist Minimum spacer distance between repeats. Defaults to 0.
#' @param maxdist Maximum spacer distance between repeats. Defaults to 1000.
#' @param summary Logical value indicating whether to store summary. Defaults
#' to FALSE.
#' @return A data.frame of the results from the EDIR database. If
#' `summary = TRUE`, returns a tibble containing the summary
#' (`$summary`), and query results (`$results`).
#' @param mismatch Logical value indicating whether to allow 1 mismatch in
#' sequences. Defaults to TRUE.
#' @param path String containing path to directory holding downloaded dataset
#' files. Defaults to NA. If not provided (`path` = NA),
#' `gene_lookup()` will use subset of data provided as example.
#' @examples
#' ## With given repeat length,
#' gene_lookup("GAA", length = 7, mindist = 10, maxdist = 1000,
#'             mismatch = TRUE)
#'
#' ## Without specified repeat length
#' gene_lookup("GAA", mindist = 0, maxdist = 1000, mismatch = TRUE)
#'
#' ## To access query results, store in variable
#' output <- gene_lookup("GAA", length = 7, mindist = 10, maxdist = 1000,
#'                         mismatch = FALSE)
#' head(output)
#'
#' ## With summary = TRUE
#' output <- gene_lookup("GAA", length = 10, mindist = 10, maxdist = 1000,
#'                         summary = TRUE,
#'                         mismatch = TRUE)
#' output$summary
#' head(output$results)
#' @export gene_lookup
#'
#' @import readr
#' @import tibble
#' @import tictoc
#' @import utils
#' @importFrom stats setNames
#'


gene_lookup <- function(gene,
                        length = NA,
                        mindist = 0,
                        maxdist = 1000,
                        summary = FALSE,
                        mismatch = TRUE,
                        path = NA) {
    #Initialize: uppercase gene ID, get chromosome number
    gene <- toupper(gene)

    # USER INPUT CHECK
    # Gene ID incorrect, or not Ensembl ID/HGNC symbol
    if (is.na(any(gene_chr == gene)) || !any(gene_chr == gene)) {
        stop("Incorrect gene ID: use ENSEMBL ID or HGNC symbol")
    }

    # Minimum distance >= maximumn distance
    if (mindist >= maxdist) {
        stop("mindist must be less than maxdist")
    }

    tic("Runtime")


    # Initialize filter and coltypes for readr
    gen_filter <- function(gene) {
        function(x, pos) {
            # subset(x, hgnc_symbol %in% gene | ensembl_gene_id %in% gene)
            subset(x, hgnc_symbol %in% gene |
                        ensembl_gene_id %in% gene)
        }
    }
    col_types <- cols(
        chromosome = col_double(),
        repeat_length = col_double(),
        repeat_seq = col_character(),
        start = col_double(),
        end = col_double(),
        repeat_seq2 = col_character(),
        start2 = col_double(),
        end2 = col_double(),
        distance = col_double(),
        ensembl_gene_id = col_character(),
        hgnc_symbol = col_character(),
        gene_range = col_character(),
        ensembl_transcript_id = col_character(),
        transcript_range = col_character(),
        intron_exon = col_character(),
        intron_exon2 = col_character(),
        feature = col_character(),
        mismatch = col_double()
    )

    # Get chromosome number
    chromosome <-
        gene_chr$chromosome_name[which(gene_chr$ensembl_gene_id == gene |
                                            gene_chr$hgnc_symbol == gene)]



    # If a path to downloaded data is given
    if (!is.na(path)) {
        # Check if last character of path is "/", if not, append "/"
        path <- ifelse(substring(path, nchar(path)) == "/",
                        path,
                        paste0(path, "/"))

        # If length is defined
        if (!is.na(length)) {
            if (length < 7 | length > 20) {
                stop("Length must be between 7 and 20")
            }

            # Define file path
            file <- paste0(path,
                            as.character(length),
                            "bp_dist1000_chr",
                            chromosome,
                            ".txt")

            # Read in file by chunks and select only rows corresponding to gene
            numgenes <- as.data.frame(
                readr::read_delim_chunked(
                    file,
                    readr::DataFrameCallback$new(gen_filter(gene)),
                    chunk_size = 10000,
                    delim = c("\t"),
                    col_types = col_types
                ),
                stringsAsFactors = FALSE
            )
            rm(file)
        }

        # If length not defined
        else {
            # Get all files of gene's chromosome number
            files <-
                paste0(path, list.files(
                    path = path,
                    pattern = paste0("bp_dist1000_chr",
                                    chromosome, ".txt")
                ))

            # Initialize empty df for relevant input
            numgenes <- data.frame()

            for (i in 7:20) {
                file <- paste0(path,
                                as.character(i),
                                "bp_dist1000_chr",
                                chromosome,
                                ".txt")

                temp <- as.data.frame(
                    readr::read_delim_chunked(
                        file,
                        readr::DataFrameCallback$new(gen_filter(gene)),
                        chunk_size = 10000,
                        delim = c("\t"),
                        col_types = col_types
                    ),
                    stringsAsFactors = FALSE
                )

                # Stop loop if current length has no results
                if (nrow(temp) == 0) {
                    break
                }

                numgenes <- rbind(numgenes, temp)
                rm(temp)
            }
            rm(files, file)

        }


        # If no given path to downloaded data, load included data subset
    } else {
        if (!(gene %in% c("GAA", "ENSG00000171298"))) {
            stop("Use gene GAA / ENSG00000171298, or download full dataset and
                specify path")
        }

        # If length is defined
        if (!is.na(length)) {
            if (length < 7 | length > 20) {
                stop("Length must be between 7 and 20")
            }

            # Load file corresponding to defined length and gene's chromosome
            numgenes <- get(paste0(as.character(length),
                                    "bp_dist1000_chr17"))

            # Filter only resuls of queried gene
            numgenes <-
                numgenes[which(numgenes$ensembl_gene_id == gene |
                                numgenes$hgnc_symbol == gene),]

        }
        # If length not defined
        else {
            # Load all internal files, corresponding to gene's chromosome number
            numgenes <- data.frame()

            for (i in 7:20) {
                input <- get(paste0(as.character(i),
                                    "bp_dist1000_chr17"))
                temp <- input[which(input$ensembl_gene_id == gene |
                                        input$hgnc_symbol == gene),]
                numgenes <- rbind(numgenes, temp)
            }
            rm(i, input, temp)
        }
    }
    # Clear memory
    gc(verbose = FALSE)

    # Define mismatch checking variable (1 for T, 0 for F)
    checkMismatch <- ifelse(mismatch == TRUE, 1, 0)

    # Filter by mismatch value, either none, or both 1 and 0
    results <- numgenes[which(numgenes$mismatch <= checkMismatch),]
    rm(numgenes)

    # Format numeric cols
    results$repeat_length <- as.numeric(results$repeat_length)
    results$distance <- as.numeric(results$distance)

    # Filter on mindist and maxdist
    results <- results[which(results$distance >= mindist &
                                results$distance <= maxdist),]

    # Sort results by repeat length, chromosome number, repeat sequence,
    # and start position
    results <-
        results[order(results$repeat_length,
                        results$chromosome,
                        results$repeat_seq,
                        results$start),]

    # No results
    if (nrow(results) == 0) {
        # cat("Gene not in dataset\n\n")
        cat("NO RESULT FOUND\n\n")
        toc()
    } else {
        # summary df
        summ <- data.frame(repeat_length =
                                sort(as.numeric(unique(
                                    results$repeat_length
                                ))),
                            unique_seqs = NA)

        # For each repeat length in results: number unique repeat_seq, number
        # unique instances, number structures, average distance
        for (i in summ$repeat_length) {
            temp <- results[which(results$repeat_length == i),]

            # Temp variable with just the sequences and start/ends
            inst <- rbind(temp[, c(3:5)],
                            setNames(temp[, c(6:8)],
                                    c("repeat_seq", "start", "end")))

            summ$unique_seqs[which(summ$repeat_length == i)] <-
                length(unique(inst$repeat_seq))
            summ$tot_instances[which(summ$repeat_length == i)] <-
                nrow(inst[!duplicated(inst),])
            summ$tot_structures[which(summ$repeat_length == i)] <-
                nrow(temp)
            summ$avg_dist[which(summ$repeat_length == i)] <-
                mean(as.numeric(temp$distance))
        }

        # Results normalizd by gene length
        gene_length <-
            as.numeric(strsplit(results$gene_range[1], "-")[[1]][2]) -
            as.numeric(strsplit(results$gene_range[1], "-")[[1]][1])
        summ$norm_instances_bp <- summ$tot_instances / gene_length
        summ$norm_instances_Mb <-
            summ$tot_instances / gene_length * 1e+06
        summ$norm_structures_bp <- summ$tot_structures / gene_length
        summ$norm_structures_Mb <-
            summ$tot_structures / gene_length * 1e+06

        # Print summ to console
        if (is.na(length)) {
            cat("\tParameters\n")
        } else {
            cat("\tParameters\n Repeat length:  ", length, " bp")
        }
        cat(
            " \n Gene:            ",
            results$ensembl_gene_id[1],
            " / ",
            results$hgnc_symbol[1],
            " \n Gene length:     ",
            gene_length,
            " bp\n Transcript ID:   ",
            results$ensembl_transcript_id[1],
            "\n Distance:        ",
            mindist,
            "-",
            maxdist,
            " bp",
            "\n Mismatch:        ",
            mismatch,
            "\n\n",
            sep = ""
        )

        # Print summ df
        print(summ)
        cat("\n")
        results <- results[order(as.numeric(rownames(results))),]
        toc()

        # Output list of summ and detailed results
        if (summary == FALSE) {
            out <- results
        } else {
            summary <- summ
            out <- tibble::lst(summary, results)
        }
        invisible(out)
    }
}
