
#----Load data------------------------------------------------------------------
data(scerevisiae_kaks)
data(diamond_intra)
data(diamond_inter)
data(yeast_annot)
data(yeast_seq)
data(cds_scerevisiae)
blast_list <- diamond_intra
blast_inter <- diamond_inter
 
pdata <- syntenet::process_input(yeast_seq, yeast_annot)
annot <- pdata$annotation["Scerevisiae"]

gene_pairs_list <- classify_gene_pairs(annot, blast_list)
gene_pairs_list <- list(Scerevisiae = gene_pairs_list[[1]][1:2, ])

cds <- list(Scerevisiae = cds_scerevisiae)

ks <- scerevisiae_kaks$Ks

## Simulate gene with CDS that is not a multiple of 3
cds2 <- cds
cds2$Scerevisiae$YGR032W <- Biostrings::subseq(
    cds2$Scerevisiae$YGR032W, 1, length(cds2$Scerevisiae$YGR032W)-1
)

#----Start tests----------------------------------------------------------------
test_that("pairs2kaks() returns a data frame with Ka, Ks, and Ka/Ks", {
    
    kaks <- pairs2kaks(gene_pairs_list, cds)
    kaks2 <- pairs2kaks(gene_pairs_list, cds2)
    
    expect_equal(class(kaks), "list")
    expect_equal(class(kaks[[1]]), "data.frame")
    expect_equal(nrow(kaks[[1]]), 2)
    expect_true("Ks" %in% names(kaks[[1]]))
    expect_true("Ka" %in% names(kaks[[1]]))
    expect_true("Ka_Ks" %in% names(kaks[[1]]))
    expect_equal(class(kaks2), "list")
})


test_that("find_ks_peaks() returns a list of mean, sd, amplitudes, ks vals", {
    
    expect_message(
        peaks <- find_ks_peaks(ks, npeaks = 1:2, verbose = TRUE)
    )

    expect_equal(class(peaks), "list")
    expect_equal(names(peaks), c("mean", "sd", "lambda", "ks"))
    expect_equal(length(peaks$mean), 2)
})


test_that("plot_ks_peaks() returns a ggplot object", {
    
    peaks <- list(
        mean = c(0.118717925754829, 0.534196999662316), 
        sd = c(0.054568151633283, 0.227909257694474), 
        lambda = c(0.49001230165837, 0.509987698341629)
    )
    
    peaks_plot <- plot_ks_peaks(peaks, binwidth = 0.05)
    expect_true("ggplot" %in% class(peaks_plot))
        
})

test_that("find_intersect_mixtures() returns a numeric scalar", {
    
    peaks <- find_ks_peaks(ks, npeaks = 2)
    peak1 <- find_ks_peaks(ks, npeaks = 1)
    inters <- find_intersect_mixtures(peaks)
    
    expect_error(
        find_intersect_mixtures(peak1)
    )
    
    expect_equal(class(inters), "numeric")
    expect_equal(length(inters), 1)
})

test_that("split_pairs_by_peak() returns a list", {
    
    ks_df <- scerevisiae_kaks[, c("dup1", "dup2", "Ks")]
    peaks <- find_ks_peaks(ks_df$Ks, npeaks = 2)
    peak1 <- find_ks_peaks(ks, npeaks = 1)
    
    spairs <- split_pairs_by_peak(ks_df, peaks) 
    spairs1 <- split_pairs_by_peak(ks_df, peak1) 
    
    expect_equal(class(spairs1), "list")
    expect_equal(class(spairs), "list")
    expect_equal(names(spairs), c("pairs", "plot"))
    expect_equal(class(spairs$pairs), "data.frame")
    expect_equal(ncol(spairs$pairs), 4)
    expect_true("ggplot" %in% class(spairs$plot))
    
})
