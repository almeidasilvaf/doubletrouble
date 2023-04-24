
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
get_anchors_list <- function(blast_list = NULL, annotation = NULL,
                             evalue = 1e-10, anchors = 5, max_gaps = 25) {
    
    # Create output directory
    daytime <- format(Sys.time(), "%d_%b_%Y_%Hh%M")
    intradir <- file.path(tempdir(), paste0("intra_", daytime))
    
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

#' Get gene pairs derived from whole-genome and small-scale duplications
#'
#' @param anchor_pairs A 2-column data frame with anchor pairs in columns 1
#' and 2.
#' @param duplicate_pairs A 2-column data frame with all duplicate pairs. This
#' is equivalent to the first 2 columns of the tabular output of BLAST-like
#' programs.
#'
#' @return A 3-column data frame with the variables:
#' \describe{
#'   \item{dup1}{Duplicated gene 1}
#'   \item{dup2}{Duplicated gene 2}
#'   \item{type}{Duplication type, which can be either 
#'               "WGD" (whole-genome duplication) or 
#'               "SSD" (small-scale duplication).}
#' }
#' @rdname get_wgd_pairs
#' @export
#' @examples
#' data(diamond_intra)
#' data(yeast_annot)
#' data(yeast_seq)
#' blast_list <- diamond_intra
#' 
#' # Get processed annotation for S. cerevisiae
#' annotation <- syntenet::process_input(yeast_seq, yeast_annot)$annotation[1]
#' 
#' # Get list of intraspecies anchor pairs
#' anchor_pairs <- get_anchors_list(blast_list, annotation)
#' anchor_pairs <- anchor_pairs[[1]][, c(1, 2)]
#' 
#' # Get duplicate pairs from DIAMOND output
#' duplicates <- diamond_intra[[1]][, c(1, 2)]
#' dups <- get_wgd_pairs(anchor_pairs, duplicates)
get_wgd_pairs <- function(anchor_pairs = NULL, duplicate_pairs = NULL) {
    
    p <- duplicate_pairs
    anchorp <- anchor_pairs
    names(p) <- c("dup1", "dup2")
    if(is.null(anchorp)) {
        duplicates <- p
        duplicates$type <- "SSD"
    } else {
        names(anchorp) <- c("anchor1", "anchor2")
        
        # Look for anchor pairs in duplicate pairs - vector-based approach
        p_vector <- paste0(p$dup1, p$dup2)
        anchor_vector <- c(paste0(anchorp$anchor1, anchorp$anchor2),
                           paste0(anchorp$anchor2, anchorp$anchor1))
        check <- which(p_vector %in% anchor_vector)
        
        # Create 2 data frames: WGD- and SSD-derived duplicates
        wgd <- p[check, ]
        wgd$type <- "WGD"
        ssd <- p[-check, ]
        ssd$type <- "SSD"
        
        # Combine the two data frames into one
        duplicates <- rbind(wgd, ssd)
    }
    return(duplicates)
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

#' Get transposed duplicate pairs
#'
#' @param pairs A 2-column data frame with duplicated gene 1 and 2 in
#' columns 1 and 2, respectively.
#' @param blast_inter A list of data frames of length 1 
#' containing BLAST tabular output for the comparison between target
#' species and outgroup. Names of list elements must match the names of 
#' list elements in `annotation`. BLASTp, DIAMOND or simular programs must 
#' be run on processed sequence data as returned by \code{process_input()}.
#' @param annotation A processed GRangesList or CompressedGRangesList object as
#' returned by \code{syntenet::process_input()}.
#'
#' @return A 3-column data frame with the following variables:
#' \describe{
#'   \item{dup1}{Duplicated gene 1}
#'   \item{dup2}{Duplicated gene 2}
#'   \item{type}{Duplication type, which can be either 
#'               "TRD" (transposed duplication) or 
#'               "DD" (dispersed duplication).}
#' }
#' @importFrom syntenet interspecies_synteny
#' @export
#' @rdname get_transposed
#' @examples 
#' data(diamond_inter)
#' data(diamond_intra)
#' data(yeast_seq)
#' data(yeast_annot)
#' blast_inter <- diamond_inter
#' 
#' # Get processed annotation
#' pdata <- syntenet::process_input(yeast_seq, yeast_annot)
#' annotation <- pdata$annotation
#' 
#' # Get duplicated pairs
#' annot <- pdata$annotation["Scerevisiae"] 
#' pairs_all <- classify_gene_pairs(diamond_intra, annot)
#' pairs <- pairs_all$Scerevisiae[pairs_all$Scerevisiae$type == "DD", c(1, 2)]
#'
#' trd <- get_transposed(pairs, blast_inter, annotation)
get_transposed <- function(pairs, blast_inter, annotation) {
    
    if(length(blast_inter) > 1) {
        stop("The list in `blast_inter` must have only 1 element.")
    }
    names(pairs)[c(1, 2)] <- c("dup1", "dup2")
    
    # Detect synteny between target species and outgroup
    target <- unlist(strsplit(names(blast_inter), "_"))[1]
    outgroup <- unlist(strsplit(names(blast_inter), "_"))[2]
    
    syn <- syntenet::interspecies_synteny(
        blast_inter,
        annotation = annotation[c(target, outgroup)]
    )
    
    # Read and parse interspecies synteny results
    parsed_syn <- collinearity2blocks(syn)[, c("anchor2", "block")]
    parsed_syn <- parsed_syn[!duplicated(parsed_syn$anchor2), ]
    final <- pairs
    
    if(!is.null(parsed_syn)) {
        # Find TRD-derived pairs (syntenic in outgroup, but not in target)
        pairs_ancestral <- merge(
            pairs[, c(1, 2)], parsed_syn, by.x = "dup1", by.y = "anchor2",
            all.x = TRUE
        )
        names(pairs_ancestral)[3] <- "block1"
        
        pairs_ancestral2 <- merge(
            pairs_ancestral, parsed_syn, by.x = "dup2", by.y = "anchor2",
            all.x = TRUE
        )
        names(pairs_ancestral2)[4] <- "block2"
        
        final <- pairs_ancestral2
        final$type <- ifelse(final$block1 == final$block2, "TRD", "DD")
        final$type[is.na(final$type)] <- "DD"
        final <- final[, c("dup1", "dup2", "type")]
    }
    return(final)
}


#' Classify small-scale duplication-derived gene pairs into subcategories
#' 
#' SSD-derived gene pairs are classified into tandem, proximal, and dispersed
#' duplicates (TD, PD, and DD, respectively).
#' 
#' @param ssd_pairs A 2-column data frame with SSD-derived gene pairs.
#' This data frame can be obtained by filtering the output of
#' \code{get_wgd_pairs()} to keep only rows where type == "SSD".
#' @param annotation_granges A processed GRanges object as in each element 
#' of the list returned by \code{syntenet::process_input()}.
#' @param annotation A processed GRangesList or CompressedGRangesList 
#' object as returned by \code{syntenet::process_input()}, which must
#' contain the gene ranges for all species.
#' @param proximal_max Numeric scalar with the maximum distance (in number
#' of genes) between two genes to consider them as proximal duplicates.
#' Default: 10.
#' @param blast_inter A list of data frames containing the tabular output
#' of interspecies BLAST/DIAMOND searches, as returned by \code{run_diamond()}.
#' Each element must contain the pairwise comparison between a target species
#' and its outgroup, which will be used to identify duplicated genes
#' derived from transpositions (TRD). If this parameter is NULL, 
#' this function will not identify TRD genes.
#'
#' @return A 3-column data frame with the variables:
#' \describe{
#'   \item{dup1}{Duplicated gene 1}
#'   \item{dup2}{Duplicated gene 2}
#'   \item{type}{Duplication type, which can be
#'               "TD" (tandem duplication), 
#'               "PD" (proximal duplication), 
#'               "TRD" (transposed duplication), and
#'               "DD" (dispersed duplication).}
#' }
#' @rdname classify_ssd_pairs
#' @export
#' @examples
#' data(diamond_intra)
#' data(diamond_inter)
#' data(yeast_annot)
#' data(yeast_seq)
#' blast_list <- diamond_intra
#' blast_inter <- diamond_inter
#' 
#' # Get processed annotation for S. cerevisiae
#' pdata <- annotation <- syntenet::process_input(yeast_seq, yeast_annot)
#' annotation <- pdata$annotation[1]
#' 
#' # Get list of intraspecies anchor pairs
#' anchor_pairs <- get_anchors_list(blast_list, annotation)
#' anchor_pairs <- anchor_pairs[[1]][, c(1, 2)]
#' 
#' # Get duplicate pairs from DIAMOND output and classify them
#' duplicates <- diamond_intra[[1]][, c(1, 2)]
#' dups <- get_wgd_pairs(anchor_pairs, duplicates)
#' ssd_pairs <- dups[dups$type == "SSD", ]
#' 
#' # Get GRanges
#' annotation_granges <- pdata$annotation[["Scerevisiae"]]
#'
#' # Get annotation list
#' annotation <- pdata$annotation
#'
#' # Classify SSD-derived gene pairs
#' ssd_classes <- classify_ssd_pairs(
#'     ssd_pairs, annotation_granges, annotation, blast_inter = blast_inter
#' )
classify_ssd_pairs <- function(ssd_pairs = NULL, annotation_granges = NULL,
                               annotation = NULL,
                               proximal_max = 10, blast_inter = NULL) {
    
    annot <- as.data.frame(annotation_granges)
    annot <- annot[, c("seqnames", "gene", "start", "end")]
    ssd_pairs <- ssd_pairs[, c(1, 2)] # Just in case df has >2 columns
    
    # Add chromosome number and order in the chromosome for each gene pair
    annot <- annot[order(annot$seqnames, annot$start), ]
    annot_bychr <- split(annot, annot$seqnames)
    annot_order <- Reduce(rbind, lapply(annot_bychr, function(x) {
        x$order <- seq_len(nrow(x))
        return(x[, c("seqnames", "gene", "order")])
    }))
    
    ssd_pos <- merge(ssd_pairs, annot_order, by.x = "dup1", by.y = "gene")
    names(ssd_pos)[c(3, 4)] <- c("chr_dup1", "order_dup1")
    ssd_pos <- merge(ssd_pos, annot_order, sort = FALSE,
                     by.x = "dup2", by.y = "gene")[, c(2, 1, 3, 4, 5, 6)]
    names(ssd_pos)[c(5, 6)] <- c("chr_dup2", "order_dup2")
    
    td <- NULL
    pd <- NULL
    others <- NULL
    #----1) Find TD-derived gene pairs--------------------------------------
    same_chr <- ssd_pos[ssd_pos$chr_dup1 == ssd_pos$chr_dup2, ]
    if(nrow(same_chr) != 0) {
        same_chr$dist <- abs(same_chr$order_dup1 - same_chr$order_dup2)
        td <- same_chr[same_chr$dist == 1, c(1, 2)]
        if(nrow(td) != 0) {
            td$type <- "TD"
        }
        
        #----2) Find PD-derived gene pairs--------------------------------------
        pd <- same_chr[same_chr$dist > 1 & same_chr$dist <= proximal_max, c(1, 2)]
        if(nrow(pd) != 0) {
            pd$type <- "PD"
        }
    }

    #----3) Find TRD-derived and DD-derived gene pairs--------------------------
    others <- rbind(
        ssd_pos[ssd_pos$chr_dup1 != ssd_pos$chr_dup2, c(1, 2)], # different chroms
        same_chr[same_chr$dist > proximal_max, c(1, 2)] # same chr, too distant
    )
    if(nrow(others) != 0) {
        others$type <- "DD"
        if(!is.null(blast_inter)) {
            others <- get_transposed(others[, c(1, 2)], blast_inter, annotation)
        }
    }

    ssd_dups <- rbind(td, pd, others)
    return(ssd_dups)
}







