

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
#' @importFrom ggplot2 ggplot aes geom_histogram ggplot stat_function
#' labs theme_bw
#' @importFrom stats dnorm
#' @importFrom rlang .data
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
    p <- ggplot(ks_df, aes(x = .data$ks)) +
        geom_histogram(binwidth = binwidth, color = "black", fill = "grey80") +
        mapply(function(mean, sd, lambda, n, binwidth, color) {
            stat_function(geom = "line", fun = function(x) {
                (dnorm(x, mean = mean, sd = sd)) * n * binwidth * lambda
            }, 
            color = color, linewidth = 1.5)
        }, mean = peaks$mean, sd = peaks$sd, lambda = peaks$lambda,
        n = length(ks_df$ks), binwidth = binwidth,
        color = pal) +
        theme_bw() +
        labs(title = "Ks distribution with peaks", y = "Frequency",
             x = "Ks values")
    
    return(p)
}