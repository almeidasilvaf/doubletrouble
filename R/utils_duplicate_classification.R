
#' Get a list of anchor pairs for each species
#'
#' @param annotation A processed GRangesList or CompressedGRangesList object as
#' returned by \code{syntenet::process_input()}.
#' @param blast_list A list of data frames containing BLAST tabular output
#' for intraspecies comparisons.
#' Each list element corresponds to the BLAST output for a given species,
#' and names of list elements must match the names of list elements in
#' `annotation`. BLASTp, DIAMOND or simular programs must be run on processed
#' sequence data as returned by \code{process_input()}.
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
#' data(sce_diamond)
#' data(sce_annotation)
#' blast_list <- sce_diamond
#' annotation <- sce_annotation
#' anchorpairs <- get_anchors_list(blast_list, annotation)
get_anchors_list <- function(blast_list = NULL, annotation = NULL,
                             evalue = 1e-10, anchors = 5, max_gaps = 25) {
    
    # Filter by e-value and convert GRanges to data frame
    fblast <- lapply(blast_list, function(x) return(x[x$evalue <= evalue, ]))
    fannotation <- lapply(annotation, function(x) {
        return(as.data.frame(x)[, c("seqnames", "gene", "start", "end")])
    })
    
    # Detect synteny
    daytime <- format(Sys.time(), "%d_%b_%Y_%Hh%M")
    intradir <- file.path(tempdir(), paste0("intra_", daytime))
    anchorp <- lapply(seq_along(fannotation), function(x) {
        sp <- names(fannotation)[x]
        blast <- fblast[grep(sp, names(fblast))] # Do names match? Get it.
        annot <- fannotation[x]
        anch <- syntenet::intraspecies_synteny(
            blast, intradir, annot, anchors = anchors, max_gaps = max_gaps
        )
        return(syntenet::parse_collinearity(anch))
    })
    names(anchorp) <- names(fannotation)
    return(anchorp)
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
#' data(sce_anchors)
#' data(sce_duplicates)
#' dups <- get_wgd_pairs(sce_anchors, sce_duplicates)
get_wgd_pairs <- function(anchor_pairs = NULL, duplicate_pairs = NULL) {
    
    p <- duplicate_pairs
    anchorp <- anchor_pairs
    names(p) <- c("dup1", "dup2")
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
    return(duplicates)
}

#' Classify small-scale duplication-derived gene pairs into subcategories
#' 
#' SSD-derived gene pairs are classified into tandem, proximal, and dispersed
#' duplicates (TD, PD, and DD, respectively).
#' 
#' @param ssd_pairs A 2-column data frame with SSD-derived gene pairs.
#' This data frame can be obtained by filtering the output of
#' \code{get_wgd_pairs()} to keep only rows where type == "SSD".
#' @param annotation A processed GRanges object as in each element of the list
#' returned by \code{syntenet::process_input()}.
#' @param proximal_max Numeric scalar with the maximum distance (in number
#' of genes) between two genes to consider them as proximal duplicates.
#' Default: 10.
#'
#' @return A 3-column data frame with the variables:
#' \describe{
#'   \item{dup1}{Duplicated gene 1}
#'   \item{dup2}{Duplicated gene 2}
#'   \item{type}{Duplication type, which can be either 
#'               "TD" (tandem duplication), "PD" (proximal duplication), and
#'               "DD" (dispersed duplication).}
#' }
#' @rdname classify_ssd_pairs
#' @export
#' @examples
#' data(sce_annotation)
#' annotation <- sce_annotation[[1]]
#' data(sce_anchors)
#' data(sce_duplicates)
#' # Get SSD-derived gene pairs
#' dups <- get_wgd_pairs(sce_anchors, sce_duplicates)
#' ssd_pairs <- dups[dups$type == "SSD", 1:2]
#'
#' # Classify SSD-derived gene pairs
#' ssd_classes <- classify_ssd_pairs(ssd_pairs, annotation)
classify_ssd_pairs <- function(ssd_pairs = NULL, annotation = NULL,
                               proximal_max = 10) {
    
    annot <- as.data.frame(annotation)[, c("seqnames", "gene", "start", "end")]
    ssd_pairs <- ssd_pairs[, 1:2] # Just in case df has >2 columns
    
    # Add chromosome number and order in the chromosome for each gene pair
    annot <- annot[order(annot$seqnames, annot$start), ]
    annot_bychr <- split(annot, annot$seqnames)
    annot_order <- Reduce(rbind, lapply(annot_bychr, function(x) {
        x$order <- seq_len(nrow(x))
        return(x[, c("seqnames", "gene", "order")])
    }))
    
    ssd_pos <- merge(ssd_pairs, annot_order, by.x = "dup1", by.y = "gene")
    names(ssd_pos)[3:4] <- c("chr_dup1", "order_dup1")
    ssd_pos <- merge(ssd_pos, annot_order, sort = FALSE,
                     by.x = "dup2", by.y = "gene")[, c(2, 1, 3:6)]
    names(ssd_pos)[5:6] <- c("chr_dup2", "order_dup2")
    
    #----1) Find TD-derived gene pairs--------------------------------------
    same_chr <- ssd_pos[ssd_pos$chr_dup1 == ssd_pos$chr_dup2, ]
    same_chr$dist <- abs(same_chr$order_dup1 - same_chr$order_dup2)
    
    td <- same_chr[same_chr$dist == 1, 1:2]
    td$type <- "TD"
    
    #----2) Find PD-derived gene pairs--------------------------------------
    pd <- same_chr[same_chr$dist > 1 & same_chr$dist <= proximal_max, 1:2]
    pd$type <- "PD"
    
    #----3) Find DD-derived gene pairs--------------------------------------
    dd1 <- ssd_pos[ssd_pos$chr_dup1 != ssd_pos$chr_dup2, 1:2]
    dd2 <- same_chr[same_chr$dist > proximal_max, 1:2]
    dd <- rbind(dd1, dd2)
    dd$type <- "DD"
    
    ssd_dups <- rbind(td, pd, dd)
    return(ssd_dups)
}







