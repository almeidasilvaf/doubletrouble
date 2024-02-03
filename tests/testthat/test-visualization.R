
# Load data ----
data(fungi_kaks)

# Start tests ----
test_that("duplicates2counts() returns counts", {
    
    duplicate_list <- classify_genes(fungi_kaks)
    
    # Get count table
    d1 <- duplicates2counts(duplicate_list, shape = "wide")
    d2 <- duplicates2counts(duplicate_list, shape = "long")
    
    expect_true(is.data.frame(d1))
    expect_true(is.data.frame(d2))
    
    expect_error(duplicates2counts(duplicate_list, shape = "error"))
})


test_that("plot_duplicate_freqs() returns a ggplot object", {
    
    duplicate_list <- classify_genes(fungi_kaks)
    
    # Get count table
    dup_counts <- duplicates2counts(duplicate_list)
    
    # Plot
    p1 <- plot_duplicate_freqs(dup_counts, plot_type = "facet")
    p2 <- plot_duplicate_freqs(dup_counts, plot_type = "stack")
    p3 <- plot_duplicate_freqs(dup_counts, plot_type = "stack_percent")
    
    expect_true("ggplot" %in% class(p1))
    expect_true("ggplot" %in% class(p2))
    expect_true("ggplot" %in% class(p3))
    
    expect_error(plot_duplicate_freqs(dup_counts, plot_type = "error"))
})


test_that("plot_ks_distro() returns a ggplot object", {
    
    df <- fungi_kaks$saccharomyces_cerevisiae
    
    p1 <- plot_ks_distro(df, bytype = TRUE)
    p2 <- plot_ks_distro(
        df, bytype = TRUE, plot_type = "violin", type_levels = c("SD", "All")
    )
    p3 <- plot_ks_distro(df, bytype = FALSE, plot_type = "histogram")
    p4 <- plot_ks_distro(df, bytype = FALSE, plot_type = "density")
    p5 <- plot_ks_distro(df, bytype = FALSE, plot_type = "density_histogram")
    
    expect_true("ggplot" %in% class(p1))
    expect_true("ggplot" %in% class(p2))
    expect_true("ggplot" %in% class(p3))
    expect_true("ggplot" %in% class(p4))
    expect_true("ggplot" %in% class(p5))
    
    expect_error(plot_ks_distro(df, bytype = TRUE, plot_type = "error"))
    expect_error(plot_ks_distro(df, bytype = FALSE, plot_type = "error"))
})


test_that("plot_rates_by_species() returns a ggplot object", {
    
    p1 <- plot_rates_by_species(fungi_kaks, rate_column = "Ka_Ks")
    p2 <- plot_rates_by_species(fungi_kaks, rate_column = "Ks", bytype = TRUE)
    p3 <- plot_rates_by_species(fungi_kaks, rate_column = "Ka")
    
    expect_true("ggplot" %in% class(p1))
    expect_true("ggplot" %in% class(p2))
    expect_true("ggplot" %in% class(p3))

    expect_error(plot_rates_by_species(fungi_kaks, rate_column = "Kw"))
})



