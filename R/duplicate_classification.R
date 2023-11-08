

#' Classify duplicate gene pairs based on their modes of duplication
#'
#' @param annotation A processed GRangesList or CompressedGRangesList object as
#' returned by \code{syntenet::process_input()}.
#' @param blast_list A list of data frames containing BLAST tabular output
#' for intraspecies comparisons.
#' Each list element corresponds to the BLAST output for a given species,
#' and names of list elements must match the names of list elements in
#' \strong{annotation}. BLASTp, DIAMOND or simular programs must be run 
#' on processed sequence data as returned by \code{process_input()}.
#' @param scheme Character indicating which classification scheme to use.
#' One of "binary", "standard", "extended", or "full". See details below
#' for information on what each scheme means. Default: "standard".
#' @param blast_inter (Only valid if \code{scheme == "extended" or "full"}).
#' A list of data frames containing BLAST tabular output 
#' for the comparison between target species and outgroups. 
#' Names of list elements must match the names of 
#' list elements in `annotation`. BLASTp, DIAMOND or simular programs must 
#' be run on processed sequence data as returned by \code{process_input()}.
#' @param intron_counts (Only valid if \code{scheme == "full"}). 
#' A list of 2-column data frames with the number of
#' introns per gene as returned by \code{get_intron_counts()}. Names
#' of list elements must match names of \strong{annotation}.
#' @param evalue Numeric scalar indicating the E-value threshold. 
#' Default: 1e-10.
#' @param anchors Numeric indicating the minimum required number of genes
#' to call a syntenic block, as in \code{syntenet::infer_syntenet}. 
#' Default: 5.
#' @param max_gaps Numeric indicating the number of upstream and downstream
#' genes to search for anchors, as in \code{syntenet::infer_syntenet}. 
#' Default: 25.
#' @param proximal_max Numeric scalar with the maximum distance (in number
#' of genes) between two genes to consider them as proximal duplicates.
#' Default: 10.
#' @param collinearity_dir Character indicating the path to the directory
#' where .collinearity files will be stored. If NULL, files will
#' be stored in a subdirectory of \code{tempdir()}. Default: NULL.
#'  
#' @return A list of 3-column data frames of duplicated gene pairs 
#' (columns 1 and 2), and their modes of duplication (column 3).
#' 
#' @details
#' The classification schemes increase in complexity (number of classes)
#' in the order 'binary', 'standard', 'extended', and 'full'.
#' 
#' For classification scheme "binary", duplicates are classified into
#' one of 'SD' (segmental duplications) or 'SSD' (small-scale duplications).
#' 
#' For classification scheme "standard" (default), duplicates are
#' classified into 'SD' (segmental duplication), 'TD' (tandem duplication),
#' 'PD' (proximal duplication), and 'DD' (dispersed duplication).
#' 
#' For classification scheme "extended", duplicates are classified into
#' 'SD' (segmental duplication), 'TD' (tandem duplication), 
#' 'PD' (proximal duplication), 'TRD' (transposon-derived duplication), 
#' and 'DD' (dispersed duplication).
#' 
#' Finally, for classification scheme "full", duplicates are classified into
#' 'SD' (segmental duplication), 'TD' (tandem duplication), 
#' 'PD' (proximal duplication), 'rTRD' (retrotransposon-derived duplication), 
#' 'dTRD' (DNA transposon-derived duplication), and 
#' 'DD' (dispersed duplication).
#' 
#' @export
#' @rdname classify_gene_pairs
#' @examples 
#' # Load example data
#' data(diamond_intra)
#' data(diamond_inter)
#' data(yeast_annot)
#' data(yeast_seq)
#' 
#' # Get processed annotation data
#' annotation <- syntenet::process_input(yeast_seq, yeast_annot)$annotation
#' 
#' # Get list of intron counts
#' txdb_list <- lapply(yeast_annot, GenomicFeatures::makeTxDbFromGRanges)
#' intron_counts <- lapply(txdb_list, get_intron_counts)
#' 
#' # Classify duplicates - full scheme
#' dup_class <- classify_gene_pairs(
#'     annotation = annotation, 
#'     blast_list = diamond_intra, 
#'     scheme = "full",
#'     blast_inter = diamond_inter, 
#'     intron_counts = intron_counts
#' )
#' 
#' # Check number of gene pairs per class
#' table(dup_class$Scerevisiae$type)
#' 
classify_gene_pairs <- function(
        annotation = NULL, blast_list = NULL, scheme = "standard",
        blast_inter = NULL, intron_counts,
        evalue = 1e-10, anchors = 5, max_gaps = 25, proximal_max = 10,
        collinearity_dir = NULL
) {
    
    anchorp <- get_anchors_list(
        blast_list, annotation, evalue, anchors, max_gaps, collinearity_dir
    )
    
    # Get duplicate pairs and filter duplicate entries
    pairs <- lapply(blast_list, function(x) {
        fpair <- x[x$evalue <= evalue, c(1, 2)]
        fpair <- fpair[fpair[, 1] != fpair[, 2], ]
        fpair <- fpair[!duplicated(t(apply(fpair, 1, sort))), ]
        names(fpair) <- c("dup1", "dup2")
        return(fpair)
    })
    
    
    dup_list <- lapply(seq_along(anchorp), function(x) {
        # 1) Get segmental duplicates
        sp <- names(anchorp)[x]
        p <- pairs[[grep(paste0(sp, "$"), names(pairs))]]
        
        dups <- get_segmental(anchorp[[x]], p)
        if(scheme == "binary") {
            dups$type <- gsub("DD", "SSD", dups$type)
            dups$type <- factor(dups$type, levels = c("SD", "SSD"))
        } else {
            # 2) Get tandem and proximal duplicates
            dups <- get_tandem_proximal(
                dups, annotation_granges = annotation[[sp]], 
                proximal_max = proximal_max
            )
            
            if(scheme %in% c("extended", "full")) {
                # 3) Get transposed duplicates
                binter <- blast_inter[startsWith(names(blast_inter), paste0(sp, "_"))]
                dups <- get_transposed(
                    dups, binter, annotation, evalue = evalue,
                    anchors = anchors, max_gaps = max_gaps,
                    collinearity_dir = collinearity_dir
                )
                
                if(scheme == "full") {
                    # 4) Get TRD classes (rTRD and dTRD)
                    dups <- get_transposed_classes(dups, intron_counts[[sp]])
                }
            }
        }
        
        return(dups)
    })
    names(dup_list) <- names(anchorp)
    
    return(dup_list)
}


#' Classify genes into unique modes of duplication
#'
#' @param gene_pairs_list List of classified gene pairs as returned 
#' by \code{classify_gene_pairs()}.
#' 
#' @return A list of 2-column data frames with variables \strong{gene} 
#' and \strong{type} representing gene ID and duplication type, respectively.
#' 
#' @details
#' If a gene is present in pairs with different duplication modes, the gene
#' is classified into a unique mode of duplication following the order
#' of priority indicated in the levels of the factor \strong{type}.
#' 
#' For scheme "binary", the order is SD > SSD.
#' For scheme "standard", the order is SD > TD > PD > DD.
#' For scheme "extended", the order is SD > TD > PD > TRD > DD.
#' For scheme "full", the order is SD > TD > PD > rTRD > dTRD > DD.
#'
#' @rdname classify_genes
#' @export
#' @importFrom GenomicRanges GRangesList
#' @examples
#' data(scerevisiae_kaks)
#' 
#' cols <- c("dup1", "dup2", "type")
#' gene_pairs_list <- list(Scerevisiae = scerevisiae_kaks[, cols])
#' 
#' class_genes <- classify_genes(gene_pairs_list)
classify_genes <- function(gene_pairs_list = NULL) {
    
    # Classify genes into unique modes used factor levels as priority order
    class_genes <- lapply(gene_pairs_list, function(x) {

        pairs_by_type <- split(x, x$type)
        gene_type <- Reduce(rbind, lapply(pairs_by_type, function(x) {
            genes <- unique(c(x$dup1, x$dup2))
            genes_df <- data.frame(gene = genes, type = x$type[1])
            genes_df <- genes_df[!duplicated(genes_df$gene), ]
            return(genes_df)
        }))

        # For genes assigned to multiple classes, keep the first (level order)        
        ref <- levels(x$type)
        gene_type <- gene_type[order(match(gene_type$type, ref)), ]
        gene_type <- gene_type[!duplicated(gene_type$gene), ]
    })
    
    return(class_genes)
}









