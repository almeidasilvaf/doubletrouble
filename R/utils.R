
#' Get a list of anchor pairs for each species
#'
#' @param blast_list A list of data frames containing BLAST tabular output
#' for intraspecies comparisons.
#' Each list element corresponds to the BLAST output for a given species,
#' and names of list elements must match the names of list elements in
#' `annotation`. BLASTp, DIAMOND or simular programs must be run on processed
#' sequence data as returned by \code{process_input()}.
#' @param annotation A processed GRangesList or CompressedGRangesList object as
#' returned by \code{syntenet::process_input()}.
#' @param evalue Numeric scalar indicating the E-value threshold. 
#' Default: 1e-10.
#' @param anchors Numeric indicating the minimum required number of genes
#' to call a syntenic block, as in \code{syntenet::infer_syntenet}. 
#' Default: 5.
#' @param max_gaps Numeric indicating the number of upstream and downstream
#' genes to search for anchors, as in \code{syntenet::infer_syntenet}. 
#' Default: 25.
#' @param collinearity_dir Character indicating the path to the directory
#' where .collinearity files will be stored. If NULL, files will
#' be stored in a subdirectory of \code{tempdir()}. Default: NULL.
#' 
#' @return A list of data frames representing intraspecies anchor pairs.
#' @importFrom syntenet intraspecies_synteny
#' @export
#' @rdname get_anchors_list 
#' @examples 
#' data(diamond_intra)
#' data(yeast_annot)
#' data(yeast_seq)
#' blast_list <- diamond_intra
#' 
#' # Get processed annotation for S. cerevisiae
#' annotation <- syntenet::process_input(yeast_seq, yeast_annot)$annotation
#' 
#' # Get list of intraspecies anchor pairs
#' anchorpairs <- get_anchors_list(blast_list, annotation)
get_anchors_list <- function(
        blast_list = NULL, annotation = NULL,
        evalue = 1e-10, anchors = 5, max_gaps = 25,
        collinearity_dir = NULL
) {
    
    # Create output directory
    intradir <- collinearity_dir
    if(is.null(intradir)) {
        daytime <- format(Sys.time(), "%d_%b_%Y_%Hh%M")
        intradir <- file.path(tempdir(), paste0("intra_", daytime))
    }

    # Filter DIAMOND list by e-value
    fblast <- lapply(blast_list, function(x) return(x[x$evalue <= evalue, ]))
    
    # Get .collinearity files for intragenome comparisons
    col_files <- syntenet::intraspecies_synteny(
        blast_intra = fblast, 
        annotation = annotation,
        intra_dir = intradir, 
        anchors = anchors, 
        max_gaps = max_gaps
    )
    
    # Parse files
    anchors <- lapply(col_files, syntenet::parse_collinearity)
    names(anchors) <- gsub("\\.collinearity", "", basename(col_files))
    
    return(anchors)
}


#' Parse .collinearity files into a data frame of syntenic blocks
#'
#' @param collinearity_paths Character vector of paths to .collinearity files.
#'
#' @return A 4-column data frame with the variables:
#' \describe{
#'   \item{block}{Syntenic block}
#'   \item{anchor1}{Anchor pair 1}
#'   \item{anchor2}{Anchor pair 2}
#' }
#'
#' @importFrom utils read.table
#' @noRd
collinearity2blocks <- function(collinearity_paths = NULL) {
    
    fname <- gsub("\\.collinearity", "", basename(collinearity_paths))
    names(collinearity_paths) <- fname
    
    blocks <- lapply(seq_along(collinearity_paths), function(x) {
        lines <- readLines(collinearity_paths[x])
        nlines <- length(lines[!startsWith(lines, "#")])
        
        df <- NULL
        if(nlines > 0) {
            df <- read.table(
                collinearity_paths[x], sep = "\t", comment.char = "#"
            )
            df <- df[, c(1, 2, 3)]
            names(df)[c(2, 3)] <- c("anchor1", "anchor2")
            
            # Get syntenic block IDs
            df$V1 <- gsub(":", "", df$V1)
            block_ids <- strsplit(df$V1, "-")
            blocks <- lapply(block_ids, function(x) return(as.numeric(x[1])))
            
            # Add syntenic block IDs to data frame
            df$block <- unlist(blocks)
            df <- df[, c("block", "anchor1", "anchor2")]
        }
        return(df)
    })
    blocks <- Reduce(rbind, blocks)
    return(blocks)
}


#' Get a data frame of intron counts per gene
#' 
#' @param txdb A `txdb` object with transcript annotations. See details below 
#' for examples on how to create `txdb` objects from different kinds of input.
#' 
#' 
#' @return A data frame with intron counts per gene, with variables:
#' \describe{
#'   \item{gene}{Character with gene IDs.}
#'   \item{introns}{Numeric with number of introns per gene.}
#' }
#' 
#' @details
#' The family of functions \code{makeTxDbFrom*} from 
#' the \strong{GenomicFeatures} package can be used to create `txdb` objects
#' from a variety of input data types. You can create `txdb` objects
#' from e.g., `GRanges` objects (\code{makeTxDbFromGRanges()}),
#' GFF files (\code{makeTxDbFromGFF()}), 
#' an Ensembl database (\code{makeTxDbFromEnsembl}), and
#' a Biomart database (\code{makeTxDbFromBiomart}).
#' 
#' @rdname get_intron_counts
#' @export
#' @importFrom GenomicFeatures intronsByTranscript
#' @importFrom AnnotationDbi select
#'
#' @examples
#' data(yeast_annot)
#' 
#' # Create txdb object from GRanges
#' txdb <- GenomicFeatures::makeTxDbFromGRanges(yeast_annot[[1]])
#'
#' # Get intron counts
#' intron_counts <- get_intron_counts(txdb)
get_intron_counts <- function(txdb) {
    
    # Get a data frame with the number of introns per transcript
    introns_by_tx <- intronsByTranscript(txdb, use.names = TRUE)
    
    intron_counts <- data.frame(
        tx = names(introns_by_tx),
        introns = lengths(introns_by_tx)
    )
    
    # Create a data frame of transcript-to-gene mapping
    suppressMessages({
        tx2gene <- AnnotationDbi::select(
            txdb, 
            keys = unique(intron_counts$tx),
            columns = "GENEID", 
            keytype = "TXNAME"
        )
    })
    names(tx2gene) <- c("tx", "gene")
    
    # Create a data frame of intron counts per gene
    intron_counts_gene <- merge(intron_counts, tx2gene)[, c("gene", "introns")]
    intron_counts_gene <- intron_counts_gene[order(-intron_counts_gene$introns), ]
    intron_counts_gene <- intron_counts_gene[!duplicated(intron_counts_gene$gene), ]
    rownames(intron_counts_gene) <- NULL
    
    return(intron_counts_gene)
}



#' Find line intersect between pairs of Gaussian mixtures
#'
#' This function finds x-axis coordinate of n-1 intersections between lines 
#' of n Gaussian mixtures. Thus, it will find 1 intersection for Ks distros
#' with 2 peaks, 2 intersections for distros with 2 peaks, and so on.
#' 
#' @param peaks A list with elements \strong{mean}, \strong{sd}, 
#' \strong{lambda}, and \strong{ks}, as returned by the 
#' function \code{fins_ks_peaks()}.
#'
#' @return A numeric scalar or vector with the x-axis coordinates of the 
#' intersections.
#' @importFrom ggplot2 ggplot_build
#' @noRd
#' @rdname find_intersect_mixtures
#' @examples
#' data(fungi_kaks)
#' scerevisiae_kaks <- fungi_kaks$saccharomyces_cerevisiae
#' ks <- scerevisiae_kaks$Ks
#' 
#' # Find 2 peaks in Ks distribution
#' peaks <- find_ks_peaks(ks, npeaks = 2)
#'
#' # Get intersects
#' inter <- find_intersect_mixtures(peaks)
find_intersect_mixtures <- function(peaks) {
    
    p <- plot_ks_peaks(peaks)
    npeaks <- length(peaks$mean)
    if(npeaks == 1) {
        stop("Cannot find intersect of peaks with only 1 peak.")
    }
    
    # Create list of density line indices to iterate through
    iteration_list <- list(
        c(2,3), c(3,4), c(4,5), c(5,6), c(6,7), c(7,8), c(8,9)
    )
    iteration_list <- iteration_list[seq_len(npeaks-1)]
    
    # Get intersection between density line i and density line i+1    
    ints <- unlist(lapply(iteration_list, function(x) {
        l1 <- x[1]
        l2 <- x[2]
        line_df <- data.frame(
            x = ggplot_build(p)$data[[l1]]$x,
            line1 = ggplot_build(p)$data[[l1]]$y,
            line2 = ggplot_build(p)$data[[l2]]$y
        )
        # Get minimal distance between lines along y axis
        line_df$delta <- line_df$line1 - line_df$line2
        
        # Get x value for minimal delta y
        int <- line_df$x[which(diff(sign(diff((abs(line_df$delta))))) == 2)+1]
        return(int)
    }))
    return(ints)
}


#' Get a duplicate count matrix for each genome
#'
#' @param duplicate_list A list of data frames with the duplicated genes or
#' gene pairs and their modes of duplication as returned 
#' by \code{classify_gene_pairs()} or \code{classify_genes()}.
#' @param shape Character specifying the shape of the output data frame.
#' One of "long" (data frame in the long shape, in the tidyverse sense),
#' or "wide" (data frame in the wide shape, in the tidyverse sense).
#' Default: "long".
#' 
#' @return If \strong{shape = "wide"}, a count matrix containing the 
#' frequency of duplicated genes (or gene pairs) by mode for each species, 
#' with species in rows and duplication modes in columns.
#' If \strong{shape = "long"}, a data frame in long format with the following
#' variables:
#' \describe{
#'   \item{type}{Factor, type of duplication.}
#'   \item{n}{Numeric, number of duplicates.}
#'   \item{species}{Character, species name}
#' }
#' 
#' @export
#' @rdname duplicates2counts
#' @examples
#' data(fungi_kaks)
#' 
#' # Get unique duplicates
#' duplicate_list <- classify_genes(fungi_kaks)
#' 
#' # Get count table
#' counts <- duplicates2counts(duplicate_list)
duplicates2counts <- function(duplicate_list, shape = "long") {
    
    # Get factor levels for variable `type`
    tlevels <- lapply(duplicate_list, function(x) return(levels(x$type)))
    tlevels <- tlevels[[names(sort(lengths(tlevels), decreasing = TRUE)[1])]]
    
    counts <- Reduce(rbind, lapply(seq_along(duplicate_list), function(x) {
        
        species <- names(duplicate_list)[x]
        
        dup_table <- duplicate_list[[x]]
        dup_table$type <- factor(dup_table$type, levels = tlevels)
        
        if(shape == "long") {
            final_dups <- as.data.frame(table(dup_table$type))
            names(final_dups) <- c("type", "n")
            final_dups$species <- species
        } else if(shape == "wide") {
            final_dups <- t(as.matrix(table(dup_table$type)))
            final_dups <- cbind(species, as.data.frame(final_dups))
        } else {
            stop("Argument 'format' must be one of 'long' or 'wide'.")
        }
        
        return(final_dups)
    }))
    
    return(counts)
}

