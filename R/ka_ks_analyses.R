

#' Calculate Ka, Ks, and Ka/Ks from duplicate gene pairs
#'
#' @param gene_pairs_list List of data frames containing duplicated gene pairs
#' as returned by \code{classify_gene_pairs()}.
#' @param cds List of DNAStringSet objects containing the coding sequences 
#' of each gene.
#' @param model Character scalar indicating which codon model to use.
#' Possible values are "Li", "NG86", "NG", "LWL", "LPB", "MLWL", "MLPB", "GY", 
#' "YN", "MYN", "MS", "MA", "GNG", "GLWL", "GLPB", "GMLWL", "GMLPB", "GYN", 
#' and "GMYN". Default: "MYN".
#' @param bp_param BiocParallel back-end to be used. 
#' Default: `BiocParallel::SerialParam()`.
#' 
#' @return A list of data frames containing gene pairs and their Ka, Ks,
#' and Ka/Ks values.
#' @importFrom MSA2dist dnastring2kaks
#' @importFrom Biostrings width
#' @importFrom BiocParallel SerialParam bplapply
#' @export
#' @rdname pairs2kaks
#' @examples 
#' data(diamond_intra)
#' data(diamond_inter)
#' data(yeast_annot)
#' data(yeast_seq)
#' data(cds_scerevisiae)
#' blast_list <- diamond_intra
#' blast_inter <- diamond_inter
#' 
#' pdata <- syntenet::process_input(yeast_seq, yeast_annot)
#' annot <- pdata$annotation["Scerevisiae"]
#' 
#' # Binary classification scheme
#' gene_pairs_list <- classify_gene_pairs(annot, blast_list)
#' gene_pairs_list <- list(
#'     Scerevisiae = gene_pairs_list[[1]][seq(1, 5, by = 1), ]
#' )
#' 
#' cds <- list(Scerevisiae = cds_scerevisiae)
#' 
#' kaks <- pairs2kaks(gene_pairs_list, cds)
pairs2kaks <- function(
        gene_pairs_list, cds, model = "MYN", 
        bp_param = BiocParallel::SerialParam()
) {
    
    kaks_list <- lapply(seq_along(gene_pairs_list), function(x) {
        
        # Get pairs for species x
        species <- names(gene_pairs_list)[x]
        pairs <- gene_pairs_list[[x]]
        names(pairs)[c(1, 2)] <- c("dup1", "dup2")
        pairs$dup1 <- gsub("^[a-zA-Z]{2,5}_", "", pairs$dup1)
        pairs$dup2 <- gsub("^[a-zA-Z]{2,5}_", "", pairs$dup2)
        
        # Remove CDS that are not multiple of 3
        fcds <- cds[[species]]
        m3 <- Biostrings::width(fcds) %% 3
        remove <- which(m3 != 0)
        if(length(remove) != 0) {
            message(
                "For species ", species, ", the lengths of ", length(remove), 
                " CDS are not multiples of 3. Removing them..."
            )
            pairs <- pairs[!pairs$dup1 %in% names(fcds)[remove], ]
            pairs <- pairs[!pairs$dup2 %in% names(fcds)[remove], ]
            fcds <- fcds[-remove]
        }
        
        # Calculate Ka, Ks, and Ka/Ks for each gene pair
        seq_list <- lapply(seq_len(nrow(pairs)), function(y) {
            return(fcds[as.character(pairs[y, c(1, 2)])])
        })
        kaks <- BiocParallel::bplapply(seq_list, function(z, model) {
            
            rates <- MSA2dist::dnastring2kaks(
                z, model = model, isMSA = FALSE, verbose = FALSE
            )
            rates <- data.frame(
                dup1 = rates$seq1,
                dup2 = rates$seq2,
                Ka = as.numeric(ifelse(rates$Ka == "NA", NA, rates$Ka)),
                Ks = as.numeric(ifelse(rates$Ks == "NA", NA, rates$Ks)),
                Ka_Ks = as.numeric(ifelse(rates[["Ka/Ks"]] == "NA", NA, rates[["Ka/Ks"]]))
            )
            
            return(rates)
        }, BPPARAM = bp_param, model = model)
        kaks <- Reduce(rbind, kaks)
        
        if("type" %in% names(pairs)) {
            kaks$type <- pairs$type
        }
        
        return(kaks)
    })
    names(kaks_list) <- names(gene_pairs_list)
    
    return(kaks_list)
}



#' Find peaks in a Ks distribution with Gaussian Mixture Models
#'
#' @param ks A numeric vector of Ks values.
#' @param npeaks Numeric scalar indicating the number of peaks in 
#' the Ks distribution. If you don't know how many peaks there are, 
#' you can include a range of values, and the number of peaks that produces
#' the lowest BIC (Bayesian Information Criterion) will be selected as the
#' optimal. Default: 2.
#' @param min_ks Numeric scalar with the minimum Ks value. Removing
#' very small Ks values is generally used to avoid the incorporation of allelic 
#' and/or splice variants and to prevent the fitting of a component to infinity.
#' Default: 0.01.
#' @param max_ks Numeric scalar indicating the maximum Ks value. Removing
#' very large Ks values is usually performed to account for Ks saturation.
#' Default: 4.
#' @param verbose Logical indicating if messages should be printed on screen.
#' Default: FALSE.
#' 
#' @return A list with the following elements:
#' \describe{
#'   \item{mean}{Numeric with the estimated means.}
#'   \item{sd}{Numeric with the estimated standard deviations.}
#'   \item{lambda}{Numeric with the estimated mixture weights.}
#'   \item{ks}{Numeric vector of filtered Ks distribution based on
#'             arguments passed to min_ks and max_ks.}
#' }
#' @importFrom mclust densityMclust
#' @export
#' @rdname find_ks_peaks
#' @examples 
#' data(fungi_kaks)
#' scerevisiae_kaks <- fungi_kaks$saccharomyces_cerevisiae
#' ks <- scerevisiae_kaks$Ks
#' 
#' # Find 2 peaks in Ks distribution
#' peaks <- find_ks_peaks(ks, npeaks = 2)
#' 
#' # From 2 to 4 peaks, verbose = TRUE to show BIC values
#' peaks <- find_ks_peaks(ks, npeaks = c(2, 3, 4), verbose = TRUE)
find_ks_peaks <- function(ks, npeaks = 2, min_ks = 0.01, max_ks = 4,
                          verbose = FALSE) {
    
    # Data preprocessing
    ks <- ks[!is.na(ks)]
    fks <- ks[ks >= min_ks]
    fks <- fks[fks <= max_ks]
    
    # Find peaks
    peaks <- mclust::densityMclust(
        fks, G = npeaks, verbose = FALSE, plot = FALSE
    )
    
    if(verbose & length(npeaks) > 1) {
        message("Optimal number of peaks: ", peaks$G)
        print(peaks$BIC)
    }
    
    # Create result list
    peak_list <- list(
        mean = peaks$parameters$mean, 
        sd = sqrt(peaks$parameters$variance$sigmasq), 
        lambda = peaks$parameters$pro,
        ks = as.numeric(peaks$data[,1])
    )
    return(peak_list)
}





#' Split gene pairs based on their Ks peaks
#' 
#' The purpose of this function is to classify gene pairs by age when there
#' are 2+ Ks peaks. This way, newer gene pairs are found within a 
#' certain number of standard deviations from the highest peak, 
#' and older genes are found close within smaller peaks.
#'
#' @param ks_df A 3-column data frame with gene pairs in columns 1 and 2,
#' and Ks values for the gene pair in column 3.
#' @param peaks A list with mean, standard deviation, and amplitude of Ks
#' peaks as generated by \code{find_ks_peaks}.
#' @param nsd Numeric with the number of standard deviations to consider
#' for each peak.
#' @param binwidth Numeric scalar with binwidth for the histogram.
#' Default: 0.05.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{pairs}{A 4-column data frame with the variables 
#'                \strong{dup1} (character), \strong{dup2} (character), 
#'                \strong{ks} (numeric), and \strong{peak} (numeric),
#'                representing duplicate gene pair, Ks values, and peak ID,
#'                respectively.}
#'   \item{plot}{A ggplot object with Ks peaks as returned by 
#'               \code{plot_ks_peaks}, but with dashed red lines indicating
#'               boundaries for each peak.}
#' }
#' 
#' @importFrom ggplot2 geom_vline
#' @export
#' @rdname split_pairs_by_peak
#' @examples
#' data(fungi_kaks)
#' scerevisiae_kaks <- fungi_kaks$saccharomyces_cerevisiae
#'
#' # Create a data frame of duplicate pairs and Ks values
#' ks_df <- scerevisiae_kaks[, c("dup1", "dup2", "Ks")]
#'
#' # Create list of peaks
#' peaks <- find_ks_peaks(ks_df$Ks, npeaks = 2)
#' 
#' # Split pairs
#' spairs <- split_pairs_by_peak(ks_df, peaks) 
split_pairs_by_peak <- function(ks_df, peaks, nsd = 2, binwidth = 0.05) {
    
    names(ks_df) <- c("dup1", "dup2", "ks")
    npeaks <- length(peaks$mean)
    
    # Filter Ks data frame as done in find_ks_peaks()
    max_ks <- max(peaks$ks)
    min_ks <- min(peaks$ks)
    ks_df <- ks_df[!is.na(ks_df$ks), ]
    ks_df <- ks_df[ks_df$ks >= min_ks & ks_df$ks <= max_ks, ]
    
    # Get minimum, intersection points, and maximum
    min_boun <- peaks$mean[1] - nsd * peaks$sd[1]
    if(min_boun < 0) { min_boun <- 0 }
    max_boun <- peaks$mean[npeaks] + nsd * peaks$sd[npeaks]
    if(max_boun > max(ks_df$ks)) { max_boun <- max(ks_df$ks) }
    
    if(npeaks == 1) {
        cutpoints <- c(min_boun, max_boun)
    } else {
        inter <- find_intersect_mixtures(peaks)
        cutpoints <- c(min_boun, inter, max_boun)
    }

    # Plot histogram with cutpoints in "brown2" dashed lines
    p <- plot_ks_peaks(peaks, binwidth = binwidth)
    for(i in seq_along(cutpoints)) {
        p <- p + geom_vline(xintercept = cutpoints[i],
                            linetype = "dashed", color = "brown2")
    }
    
    # Create list of intervals
    int_list <- lapply(seq_len(length(cutpoints)-1), function(x) {
        return(c(cutpoints[x], cutpoints[x] + 1))
    })
    
    # Create list of data frames for each interval
    split_pairs <- Reduce(rbind, lapply(seq_along(int_list), function(x) {
        ivec <- int_list[[x]]
        pairs <- ks_df[ks_df$ks >= ivec[1] & ks_df$ks < ivec[2], ]
        pairs$peak <- x
        return(pairs)
    }))
    
    result_list <- list(pairs = split_pairs, plot = p)
}


