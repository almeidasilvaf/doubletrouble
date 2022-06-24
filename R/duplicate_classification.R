

#' Classify duplicate gene pairs based on their modes of duplication
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
#' @param binary Logical indicating whether to perform a binary classification
#' (i.e., whole-genome and small-scale duplication) or not. If FALSE,
#' small-scale duplications are subdivided into tandem, proximal, and dispersed
#' duplications. Default: FALSE.
#' @param proximal_max Numeric scalar with the maximum distance (in number
#' of genes) between two genes to consider them as proximal duplicates.
#' Default: 10.
#'  
#' @return A list of 3-column data frames of duplicated gene pairs 
#' (columns 1 and 2), and their modes of duplication (column 3).
#' @export
#' @importFrom GenomicRanges GRangesList
#' @rdname classify_gene_pairs
#' @examples 
#' data(sce_diamond)
#' data(sce_annotation)
#' blast_list <- sce_diamond
#' annotation <- sce_annotation
#' duplicates <- classify_gene_pairs(blast_list, annotation)
classify_gene_pairs <- function(blast_list = NULL, annotation = NULL,
                                evalue = 1e-10, anchors = 5, max_gaps = 25,
                                binary = FALSE, proximal_max = 10) {
    
    anchorp <- get_anchors_list(
        blast_list, annotation, evalue, anchors, max_gaps
    )
    
    # Get duplicate pairs and filter duplicate entries
    pairs <- lapply(blast_list, function(x) {
        fpair <- x[x$evalue <= evalue, 1:2]
        fpair <- fpair[fpair[, 1] != fpair[, 2], ]
        fpair <- fpair[!duplicated(t(apply(fpair, 1, sort))),]
        names(fpair) <- c("dup1", "dup2")
        return(fpair)
    })

    dups <- lapply(seq_along(anchorp), function(x) {
        # Find WGD-derived gene pairs
        sp <- names(anchorp)[x]
        p <- pairs[[grep(sp, names(pairs))]]
        dups <- get_wgd_pairs(anchorp[[x]], p)

        if(!binary) {
            ssd <- dups[dups$type == "SSD", ]
            wgd <- dups[dups$type == "WGD", ]
            
            # Classify SSD-derived gene pairs
            annot <- annotation[[sp]]
            ssd_classes <- classify_ssd_pairs(ssd, annot)
            
            dups <- rbind(wgd, ssd_classes)
        }
        rownames(dups) <- NULL
        return(dups)
    })
    names(dups) <- names(anchorp)
    return(dups)
}


#' Classify genes into unique modes of duplication
#'
#' @param gene_pairs_list List of classified gene pairs as returned 
#' by \code{classify_gene_pairs()}.
#' 
#' @return A list of 2-column data frames with variables \strong{gene} 
#' and \strong{type} representing gene ID and duplication type, respectively.
#'
#' @rdname classify_genes
#' @export
#' @examples
#' data(sce_diamond)
#' data(sce_annotation)
#' blast_list <- sce_diamond
#' annotation <- sce_annotation
#' gene_pairs_list <- classify_gene_pairs(blast_list, annotation)
#' class_genes <- classify_genes(gene_pairs_list)
classify_genes <- function(gene_pairs_list = NULL) {
    
    # Classify genes following the order of priority: WGD > TD > PD > DD
    class_genes <- lapply(gene_pairs_list, function(x) {

        pairs_by_type <- split(x, x$type)
        gene_type <- Reduce(rbind, lapply(pairs_by_type, function(x) {
            genes <- unique(c(x$dup1, x$dup2))
            genes_df <- data.frame(gene = genes, type = x$type[1])
            genes_df <- genes_df[!duplicated(genes_df$gene), ]
            return(genes_df)
        }))
        
        # Reorder 'type' variable based on
        ref <- c("WGD", "TD", "PD", "DD")
        if(length(unique(gene_type$type)) == 2) {
            ref <- c("WGD", "SSD")
        }
        gene_type <- gene_type[order(match(gene_type$type, ref)), ]
        gene_type <- gene_type[!duplicated(gene_type$gene), ]
    })
    return(class_genes)
}









