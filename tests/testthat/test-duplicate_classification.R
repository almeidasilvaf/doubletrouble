
#----Load data------------------------------------------------------------------
data(diamond_intra)
data(scerevisiae_kaks)
data(diamond_inter)
data(yeast_annot)
data(yeast_seq)
blast_list <- diamond_intra
blast_inter <- diamond_inter

pdata <- syntenet::process_input(yeast_seq, yeast_annot)
annotation <- pdata$annotation
annotation_granges <- pdata$annotation[["Scerevisiae"]]

## Get anchor pairs
all <- scerevisiae_kaks
all$dup1 <- paste0("Sce_", all$dup1)
all$dup2 <- paste0("Sce_", all$dup2)
anchor_pairs <- all[all$type == "WGD", 1:2]
ssd <- all[all$type != "WGD", 1:2]

## Get duplicate pairs from DIAMOND output
duplicates <- diamond_intra[[1]][, 1:2]


#----Start tests----------------------------------------------------------------
test_that("get_anchors_list() returns a list of anchor pairs", {
    anchorpairs <- get_anchors_list(blast_list, annotation)
    
    expect_equal(class(anchorpairs), "list")
    expect_equal(class(anchorpairs[[1]]), "data.frame")
    expect_equal(ncol(anchorpairs[[1]]), 2)
})

test_that("get_wgd_pairs() returns a data frame with WGD and SSD pairs", {
    dups <- get_wgd_pairs(anchor_pairs, duplicates)
    
    expect_equal(class(dups), "data.frame")
    expect_equal(ncol(dups), 3)
    expect_equal(names(dups), c("dup1", "dup2", "type"))
    expect_equal(length(unique(dups$type)), 2)
})

test_that("get_transposed() finds transposon-derived pairs", {
    trd <- get_transposed(ssd, blast_inter, annotation)
    
    expect_equal(class(trd), "data.frame")
    expect_equal(ncol(trd), 3)
    expect_equal(names(trd), c("dup1", "dup2", "type"))
})

test_that("classify_ssd_pairs() returns a data frame of SSD subclasses", {
    dups <- get_wgd_pairs(anchor_pairs, duplicates)
    ssd_pairs <- dups[dups$type == "SSD", ]
    
    ssd_classes <- classify_ssd_pairs(
        ssd_pairs, annotation_granges, annotation, blast_inter = blast_inter
    )
    
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


