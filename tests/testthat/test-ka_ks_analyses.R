
#----Load data------------------------------------------------------------------
data(gma_dups_kaks)
ks <- gma_dups_kaks$Ks
ks <- ks[!is.na(ks)]
 
# Remove Ks values > 1 for testing purposes
ks <- ks[ks <= 1]


#----Start tests----------------------------------------------------------------
test_that("find_peak_number() returns a non-negative integer", {
    
    # Using everything small for testing purposes
    bootstraps <- 1
    npeaks <- find_peak_number(ks[ks < 0.4], bootstraps, max_components = 1)
    
    expect_equal(class(npeaks), "numeric")
    expect_equal(npeaks, round(npeaks))
    expect_true(npeaks > 0)
})

test_that("find_ks_peaks() returns a list of mean, sd, and amplitudes", {
    
    peaks <- find_ks_peaks(ks, npeaks = 2)
    
    expect_equal(class(peaks), "list")
    expect_equal(names(peaks), c("mean", "sd", "lambda"))
    expect_equal(length(peaks$mean), 2)
})

test_that("plot_ks_peaks() returns a ggplot object", {
    
    peaks <- list(
        mean = c(0.118717925754829, 0.534196999662316), 
        sd = c(0.054568151633283, 0.227909257694474), 
        lambda = c(0.49001230165837, 0.509987698341629)
    )
    
    peaks_plot <- plot_ks_peaks(ks, peaks, binwidth = 0.05)
    expect_true("ggplot" %in% class(peaks_plot))
        
})

test_that("find_intersect_mixtures() returns a numeric scalar", {
    
    peaks <- find_ks_peaks(ks, npeaks = 2)
    inters <- find_intersect_mixtures(ks, peaks)
    
    expect_equal(class(inter), "numeric")
    expect_equal(length(inter), 1)
})

test_that("split_pairs_by_peak() returns a list", {
    
    # Create a data frame of duplicate pairs and Ks values
    ks_df <- gma_dups_kaks[!is.na(gma_dups_kaks$Ks), c("dup1", "dup2", "Ks")]
     
    # Remove Ks values >1 for testing purposes
    ks_df <- ks_df[ks_df$Ks <= 1, ]
    
    # Create list of peaks
    peaks <- find_ks_peaks(ks_df$Ks, npeaks = 2)
    spairs <- split_pairs_by_peak(ks_df, peaks) 
    
    expect_equal(class(spairs), "list")
    expect_equal(names(spairs), c("pairs", "plot"))
    expect_equal(class(spairs$pairs), "data.frame")
    expect_equal(ncol(spairs$pairs), 4)
    expect_true("ggplot" %in% class(spairs$plot))
    
    
})
