
#----Load data------------------------------------------------------------------
data(sce_diamond)
data(sce_annotation)
data(sce_anchors)
data(sce_duplicates)

annotation_ranges <- sce_annotation[[1]]
dups <- get_wgd_pairs(sce_anchors, sce_duplicates)
ssd_pairs <- dups[dups$type == "SSD", 1:2]
blast_list <- sce_diamond
annotation <- sce_annotation

#----Start tests----------------------------------------------------------------
test_that("get_anchors_list() returns a list of anchor pairs", {
    anchorpairs <- get_anchors_list(blast_list, annotation)
    
    expect_equal(class(anchorpairs), "list")
    expect_equal(class(anchorpairs[[1]]), "data.frame")
    expect_equal(ncol(anchorpairs[[1]]), 2)
})

test_that("get_wgd_pairs() returns a data frame with WGD and SSD pairs", {
    dups <- get_wgd_pairs(sce_anchors, sce_duplicates)
    
    expect_equal(class(dups), "data.frame")
    expect_equal(ncol(dups), 3)
    expect_equal(names(dups), c("dup1", "dup2", "type"))
    expect_equal(length(unique(dups$type)), 2)
})

test_that("classify_ssd_pairs() returns a data frame of SSD subclasses", {
    ssd_classes <- classify_ssd_pairs(ssd_pairs, annotation_ranges)
    
    expect_equal(class(ssd_classes), "data.frame")
    expect_equal(ncol(ssd_classes), 3)
    expect_equal(names(ssd_classes), c("dup1", "dup2", "type"))
    expect_true(length(unique(ssd_classes$type)) > 2)
})

test_that("classify_gene_pairs() returns a data frame", {
    duplicates <- classify_gene_pairs(blast_list, annotation)
    
    expect_equal(class(duplicates), "list")
    expect_equal(class(duplicates[[1]]), "data.frame")
    expect_equal(ncol(duplicates[[1]]), 3)
    expect_equal(names(duplicates[[1]]), c("dup1", "dup2", "type"))
    expect_true(length(unique(duplicates[[1]]$type)) > 2)
})

test_that("classify_genes() returns a unique duplication mode for each gene", {
    gene_pairs_list <- classify_gene_pairs(blast_list, annotation)
    class_genes <- classify_genes(gene_pairs_list)
    
    n_dup_genes <- length(unique(c(gene_pairs_list[[1]]$dup1,
                                   gene_pairs_list[[1]]$dup2)))
    
    expect_equal(class(class_genes), "list")
    expect_equal(class(class_genes[[1]]), "data.frame")
    expect_equal(ncol(class_genes[[1]]), 2)
    expect_equal(nrow(class_genes[[1]]), n_dup_genes)
})


