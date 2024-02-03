
#' Create a named vector with a color palette for duplication modes
#' 
#' @return A named character vector with colors for each duplication mode.
#' @noRd
#' 
dup_palette <- function() {
    
    pal <- c(
        All = "gray20",
        SD = "#EFC000FF",
        TD = "#CD534CFF",
        PD = "#79AF97FF",
        TRD = "#7AA6DCFF",
        rTRD = "#7AA6DCFF",
        dTRD = "#003C67FF",
        DD = "#6A6599FF",
        SSD = "#7AA6DCFF"
    )
    
    return(pal)
}


#' Plot frequency of duplicates per mode for each species
#'
#' @param dup_counts A data frame in long format with the number of
#' duplicates per mode for each species, as returned by 
#' the function \code{duplicates2counts}.
#' @param plot_type Character indicating how to plot frequencies. One of
#' 'facet' (facets for each level of the variable \strong{type}),
#' 'stack' (levels of the variable \strong{type} as stacked bars), or
#' 'stack_percent' (levels of the variable \strong{type} as stacked bars,
#' with x-axis representing relative frequencies). Default: 'facet'.
#' @param remove_zero Logical indicating whether or not to remove rows
#' with zero values. Default: TRUE.
#' 
#' @return A ggplot object.
#' 
#' @importFrom ggplot2 ggplot aes geom_bar facet_wrap theme_bw theme labs
#' scale_fill_manual element_blank
#' @importFrom rlang .data
#' @export
#' @rdname plot_duplicate_freqs
#' @examples
#' data(fungi_kaks)
#' 
#' # Get unique duplicates
#' duplicate_list <- classify_genes(fungi_kaks)
#' 
#' # Get count table
#' dup_counts <- duplicates2counts(duplicate_list)
#' 
#' # Plot counts
#' plot_duplicate_freqs(dup_counts, plot_type = "stack_percent")
plot_duplicate_freqs <- function(
        dup_counts, plot_type = "facet", remove_zero = TRUE
) {
    
    # Define palette
    pal <- dup_palette()
    
    # Remove zeros
    if(remove_zero) { dup_counts <- dup_counts[dup_counts$n != 0, ] }
    
    if(plot_type == "facet") {
        p <- ggplot(dup_counts, aes(x = .data$n, y = .data$species)) +
            geom_bar(
                aes(fill = .data$type), stat = "identity", color = "grey20",
                show.legend = FALSE
            ) +
            facet_wrap("type", nrow = 1, scales = "free_x") +
            labs(y = "", x = "Absolute frequency")
        
    } else if(plot_type == "stack") {
        p <- ggplot(
            dup_counts, aes(x = .data$n, y = .data$species, fill = .data$type)
        ) +
            geom_bar(color = "gray20", position = "stack", stat = "identity") +
            labs(fill = "Type", y = "", x = "Absolute frequency")
        
    } else if(plot_type == "stack_percent") {
        p <- ggplot(
            dup_counts, aes(x = .data$n, y = .data$species, fill = .data$type)
        ) +
            geom_bar(position = "fill", stat = "identity", color = "gray20") +
            labs(fill = "Type", y = "", x = "Relative frequency")
        
    } else {
        stop("Input to argument 'plot_type' must be one of 'facet', 'stack', or 'stack_percent'.")
    }
    
    p <- p + 
        scale_fill_manual(values = pal) +
        theme_bw() +
        theme(panel.grid = element_blank())

    return(p)
}


#' Plot distribution of synonymous substitution rates (Ks)
#' 
#' @param ks_df A data frame with Ks values for each gene pair
#' as returned by \code{pairs2kaks()}.
#' @param min_ks Numeric indicating the minimum Ks value to keep. 
#' Default: 0.01.
#' @param max_ks Numeric indicating the maximum Ks value to keep.
#' Default: 2.
#' @param bytype Logical indicating whether or not to plot the distribution
#' by type of duplication (requires a column named `type`).
#' @param type_levels (Only valid if \strong{bytype} is not NULL) Character
#' indicating which levels of the variable specified in 
#' parameter \strong{group_by} should be kept. By default, all levels are kept.
#' @param plot_type Character indicating the type of plot to create. 
#' If \strong{bytype = TRUE}, possible types are "histogram" or "violin".
#' If \strong{bytype = FALSE}, possible types are "histogram", "density",
#' or "density_histogram". Default: "histogram".
#' @param binwidth (Only valid if \strong{plot_type = "histogram"}) 
#' Numeric indicating the bin width. Default: 0.03.
#' 
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_density geom_histogram facet_grid
#' theme theme_bw geom_violin geom_boxplot scale_x_continuous vars
#' scale_y_continuous after_stat
#' @importFrom stats density
#' @export
#' @rdname plot_ks_distro
#' @examples
#' data(fungi_kaks)
#' ks_df <- fungi_kaks$saccharomyces_cerevisiae
#' 
#' # Plot distro
#' plot_ks_distro(ks_df, bytype = TRUE)
plot_ks_distro <- function(
        ks_df, min_ks = 0.01, max_ks = 2,
        bytype = FALSE, type_levels = NULL,
        plot_type = "histogram",
        binwidth = 0.03
) {
    
    pal <- dup_palette()
    
    # Filter Ks values
    filt_ks <- ks_df[ks_df$Ks >= min_ks & ks_df$Ks <= max_ks, ]
    filt_ks <- filt_ks[!is.na(filt_ks$Ks), ]
    
    if(bytype) {
        # Add level for all combined
        ks_all <- filt_ks
        ks_all$type <- factor("All", levels = "All")
        filt_ks <- rbind(ks_all, filt_ks)
        
        # Keep only desired levels (optional)
        if(!is.null(type_levels)) {
            filt_ks <- filt_ks[filt_ks$type %in% type_levels, ]
            filt_ks$type <- droplevels(filt_ks$type)
        }
        
        # Plot
        if(plot_type == "histogram") {
            p <- ggplot(filt_ks, aes(x = .data$Ks)) + 
                geom_histogram(
                    aes(fill = .data$type),
                    color = "gray30", binwidth = binwidth, show.legend = FALSE
                ) +
                scale_fill_manual(values = pal) +
                facet_grid(rows = vars(.data$type), scales = "free_y") +
                labs(y = "Count")
        } else if(plot_type == "violin") {
            p <- ggplot(filt_ks, aes(x = .data$Ks, y = .data$type)) +
                geom_violin(aes(fill = .data$type), show.legend = FALSE) +
                scale_fill_manual(values = pal) +
                labs(y = "Type")
        } else {
            stop("When plotting by groups, `plot_type` must be either 'histogram' or 'violin'.")
        }
        p <- p + labs(title = "Ks distribution for gene pairs by mode")
    } else {
        
        if(plot_type == "histogram") {
            p <- ggplot(filt_ks, aes(x = .data$Ks)) + 
                geom_histogram(
                    fill = "#9196ca", color = "#3e57a7", binwidth = binwidth
                ) +
                labs(y = "Count")
        } else if(plot_type == "density") {
            p <- ggplot(filt_ks, aes(x = .data$Ks)) +
                geom_density(
                    fill = "#9196ca", color = "#3e57a7"
                ) +
                labs(y = "Density")
        } else if(plot_type == "density_histogram") {
            p <- ggplot(filt_ks, aes(x = .data$Ks)) + 
                geom_histogram(
                    aes(y = after_stat(density)), alpha = 0.5,
                    fill = "#9196ca", color = "#3e57a7", binwidth = binwidth
                ) +
                geom_density(color = "gray30", linewidth = 1) +
                labs(y = "Density")
        } else {
            stop("Without groups, `plot_type` must be one of 'histogram', 'density', or 'density_histogram'.")
        }
        p <- p + scale_y_continuous(expand = c(1e-2, 1e-2)) +
            labs(title = "Ks distribution for gene pairs")
    }
    
    # Polish plot
    p <- p +
        theme_bw() +
        theme(panel.grid = element_blank()) +
        scale_x_continuous(expand = c(1e-2, 1e-2)) +
        labs(x = expression(K[s]))
    
    return(p)
}


#' Plot distributions of substitution rates (Ka, Ks, or Ka/Ks) per species
#'
#' @param kaks_list A list of data frames with substitution rates per gene
#' pair in each species as returned by \code{pairs2kaks()}.
#' @param rate_column Character indicating the name of the column to plot.
#' Default: "Ks".
#' @param bytype Logical indicating whether or not to show distributions by
#' type of duplication. Default: FALSE.
#' @param range Numeric vector of length 2 indicating the minimum and maximum
#' values to plot. Default: \code{c(0, 2)}.
#' @param fill Character with color to use for the fill aesthetic. Ignored
#' if \strong{bytype = TRUE}. Default: "deepskyblue3".
#' @param color Character with color to use for the color aesthetic. Ignored
#' if \strong{bytype = FALSE}. Default: "deepskyblue4".
#' 
#' @return A ggplot object.
#'
#' @details
#' Data will be plotted using the species order of the list. To change the
#' order of the species to plot, reorder the list elements 
#' in \strong{kaks_list}.
#' 
#' @importFrom ggplot2 geom_violin scale_fill_manual facet_wrap
#' @export
#' @rdname plot_rates_by_species
#' @examples
#' data(fungi_kaks)
#'
#' # Plot rates
#' plot_rates_by_species(fungi_kaks, rate_column = "Ka_Ks") 
plot_rates_by_species <- function(
        kaks_list, rate_column = "Ks", bytype = FALSE, range = c(0, 2),
        fill = "deepskyblue3", color = "deepskyblue4"
) {
    
    xl <- switch(
        rate_column,
        "Ks" = expression(K[s]),
        "Ka" = expression(K[a]),
        "Ka_Ks" = expression(K[a] / K[s]),
        stop("Input to 'rate_column' must be one of 'Ka', 'Ks', or 'Ka_Ks'.")
    )
    
    # From list to data frame
    rate_df <- Reduce(rbind, lapply(seq_along(kaks_list), function(x) {
        
        df <- kaks_list[[x]]
        df$species <- names(kaks_list)[x]
        
        return(df)
    }))
    rate_df <- rate_df[!is.na(rate_df[[rate_column]]), ]
    rate_df <- rate_df[rate_df[[rate_column]] >= range[1] & 
                           rate_df[[rate_column]] <= range[2], ]
    rate_df$species <- factor(rate_df$species, levels = names(kaks_list))
    
    # Create plot
    if(bytype) {
        p <- ggplot(rate_df, aes(x = .data[[rate_column]], y = .data$species)) +
            geom_violin(aes(fill = .data$type), show.legend = FALSE) +
            scale_fill_manual(values = dup_palette()) +
            facet_wrap("type", nrow = 1) 
        
    } else {
        p <- ggplot(rate_df, aes(x = .data[[rate_column]], y = .data$species)) +
            geom_violin(fill = fill, color = color)
    }
    
    # Polish the plot
    p <- p +
        theme_bw() +
        labs(x = xl, y = "") +
        theme(panel.grid = element_blank())
    
    return(p)
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
#' @importFrom ggplot2 ggplot aes geom_histogram ggplot stat_function
#' labs theme_bw
#' @importFrom stats dnorm
#' @importFrom rlang .data
#' @rdname plot_ks_peaks
#' @export
#' @examples 
#' data(fungi_kaks)
#' scerevisiae_kaks <- fungi_kaks$saccharomyces_cerevisiae
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
        theme(panel.grid = element_blank()) +
        labs(title = "Ks distribution with peaks", y = "Frequency",
             x = "Ks values")
    
    return(p)
}