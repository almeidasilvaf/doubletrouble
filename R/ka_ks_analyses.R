

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
#' @param threads Numeric scalar indicating the number of threads to use.
#' Default: 1.
#' 
#' @return A list of data frames containing gene pairs and their Ka, Ks,
#' and Ka/Ks values.
#' @importFrom MSA2dist dnastring2kaks
#' @importFrom Biostrings width
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
#' gene_pairs_list <- classify_gene_pairs(blast_list, annot, binary = TRUE)
#' gene_pairs_list <- list(
#'     Scerevisiae = gene_pairs_list[[1]][seq(1, 5, by = 1), ]
#' )
#' 
#' cds <- list(Scerevisiae = cds_scerevisiae)
#' 
#' kaks <- pairs2kaks(gene_pairs_list, cds)
pairs2kaks <- function(gene_pairs_list, cds, model = "MYN", threads = 1) {
    
    kaks_list <- lapply(seq_along(gene_pairs_list), function(x) {
        species <- names(gene_pairs_list)[x]
        pairs <- gene_pairs_list[[x]]
        names(pairs)[c(1, 2)] <- c("dup1", "dup2")
        
        # Remove species ID from gene IDs in gene pairs
        pairs$dup1 <- gsub("[a-zA-Z]{2,5}_", "", pairs$dup1)
        pairs$dup2 <- gsub("[a-zA-Z]{2,5}_", "", pairs$dup2)
        
        kaks <- Reduce(rbind, lapply(seq_len(nrow(pairs)), function(y) {
            genes <- as.character(pairs[y, c(1, 2)])
            cds_genes <- cds[[species]][genes]
            
            multiple_3 <- Biostrings::width(cds_genes) %% 3
            if(sum(multiple_3) > 0) {
                vals <- NULL
                pgenes <- paste0(genes, collapse = ", ")
                message("CDS length is not a multiple of 3 for pair ", pgenes)
            } else {
                vals <- MSA2dist::dnastring2kaks(
                    cds_genes, model = "MYN", isMSA = FALSE, threads = threads
                )
                vals <- as.data.frame(vals)
                vals <- vals[, c("seq1", "seq2", "Ka", "Ks", "Ka/Ks")]
                names(vals) <- c("dup1", "dup2", "Ka", "Ks", "Ka_Ks")
                rownames(vals) <- NULL
                if("type" %in% names(pairs)) { vals$type <- pairs[y, "type"] }
            }
            return(vals)
        }))
        kaks$Ka <- gsub("NA", NA, kaks$Ka)
        kaks$Ka <- as.numeric(kaks$Ka)
        kaks$Ks <- gsub("NA", NA, kaks$Ks)
        kaks$Ks <- as.numeric(kaks$Ks)
        kaks$Ka_Ks <- gsub("NA", NA, kaks$Ka_Ks)
        kaks$Ka_Ks <- as.numeric(kaks$Ka_Ks)
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
#' data(scerevisiae_kaks)
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


#' Plot histogram of Ks distribution with peaks
#'
#' @param peaks A list with elements \strong{mean}, \strong{sd}, 
#' \strong{lambda}, and \strong{ks}, as returned by the 
#' function \code{fins_ks_peaks()}.
#' @param binwidth Numeric scalar with binwidth for the histogram.
#' Default: 0.05.
#'
#' @return A ggplot object with a histogram and lines for each Ks peak.
#'
#' @importFrom ggplot2 ggplot aes_ geom_histogram ggplot stat_function
#' labs theme_bw
#' @importFrom stats dnorm
#' @rdname plot_ks_peaks
#' @export
#' @examples 
#' data(scerevisiae_kaks)
#' ks <- scerevisiae_kaks$Ks
#' 
#' # Find 2 peaks in Ks distribution
#' peaks <- find_ks_peaks(ks, npeaks = 2)
#'
#' # Plot
#' plot_ks_peaks(peaks, binwidth = 0.05)
plot_ks_peaks <- function(peaks = NULL, binwidth = 0.05) {
    
    ks_df <- data.frame(ks = peaks$ks)
    
    # Define color palette
    pal <- c(
        "#6A6599FF", "#79AF97FF", "#B24745FF", "#00A1D5FF", 
        "#DF8F44FF", "#374E55FF", "#F39B7FFF", "#3C5488FF"
    )
    
    pal <- as.list(rev(pal[seq_along(peaks$mean)]))
    
    # Plot 
    p <- ggplot(ks_df, aes_(x = ~ks)) +
        geom_histogram(binwidth = binwidth, color = "black", fill = "grey80") +
        mapply(function(mean, sd, lambda, n, binwidth, color) {
            stat_function(geom = "line", fun = function(x) {
                (dnorm(x, mean = mean, sd = sd)) * n * binwidth * lambda
            }, 
            color = color, size = 1.5)
        }, mean = peaks$mean, sd = peaks$sd, lambda = peaks$lambda,
        n = length(ks_df$ks), binwidth = binwidth,
        color = pal) +
        theme_bw() +
        labs(title = "Ks distribution with peaks", y = "Frequency",
             x = "Ks values")
    
    return(p)
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
#' data(scerevisiae_kaks)
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
#' data(scerevisiae_kaks)
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


