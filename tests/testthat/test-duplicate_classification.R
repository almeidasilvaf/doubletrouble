
#----Load data------------------------------------------------------------------
data(diamond_intra)
data(fungi_kaks)
data(diamond_inter)
data(yeast_annot)
data(yeast_seq)
blast_list <- diamond_intra
blast_inter <- diamond_inter

scerevisiae_kaks <- fungi_kaks$saccharomyces_cerevisiae

pdata <- syntenet::process_input(yeast_seq, yeast_annot)
annotation <- pdata$annotation
annotation_granges <- pdata$annotation[["Scerevisiae"]]

## Get anchor pairs
all <- scerevisiae_kaks
all$dup1 <- paste0("Sce_", all$dup1)
all$dup2 <- paste0("Sce_", all$dup2)
anchor_pairs <- all[all$type == "SD", 1:2]
ssd <- all[all$type != "SD", 1:2]

## Get duplicate pairs from DIAMOND output
duplicates <- diamond_intra[[1]][, 1:2]


txdb <- GenomicFeatures::makeTxDbFromGRanges(yeast_annot[[1]])
intron_counts <- get_intron_counts(txdb)

ic_list <- list(Scerevisiae = intron_counts)

#----Start tests----------------------------------------------------------------
test_that("get_* functions returned classified duplicates", {
    
    # 1) get_anchors_list()
    anchorpairs <- get_anchors_list(blast_list, annotation)
    
    expect_equal(class(anchorpairs), "list")
    expect_equal(class(anchorpairs[[1]]), "data.frame")
    expect_equal(ncol(anchorpairs[[1]]), 2)
    
    # 2) get_segmental()
    dups <- get_segmental(anchor_pairs, duplicates)
    dups2 <- get_segmental(NULL, duplicates)
    
    expect_equal(class(dups2), "data.frame")
    expect_equal(ncol(dups), 3)
    expect_equal(names(dups), c("dup1", "dup2", "type"))
    expect_equal(length(unique(dups$type)), 2)
    
    # 3) get_tandem_proximal()
    td_pd <- get_tandem_proximal(dups, annotation_granges)
    
    expect_equal(class(td_pd), "data.frame")
    expect_equal(ncol(td_pd), 3)
    expect_equal(names(td_pd), c("dup1", "dup2", "type"))
    expect_equal(length(unique(td_pd$type)), 4)
    
    # 4) get_transposed
    trd <- get_transposed(td_pd, blast_inter, annotation)
    
    binter2 <- list(e1 = blast_inter[[1]], e2 = blast_inter[[1]])
    expect_error(get_transposed(td_pd, binter2, annotation))
    
    expect_equal(class(trd), "data.frame")
    expect_equal(ncol(trd), 3)
    expect_equal(names(trd), c("dup1", "dup2", "type"))
    expect_equal(length(unique(trd$type)), 5)
    
    # 5) get_transposed_classes
    trdc <- get_transposed_classes(trd, intron_counts)
    
    expect_equal(class(trdc), "data.frame")
    expect_equal(ncol(trdc), 3)
    expect_equal(names(trdc), c("dup1", "dup2", "type"))
    expect_equal(length(unique(trdc$type)), 6)
})


test_that("classify_gene_pairs() and classify_genes() return a data frame", {
    
    # 1) classify_gene_pairs()
    dup_full <- classify_gene_pairs(
        annotation = annotation,
        blast_list = diamond_intra,
        scheme = "full",
        blast_inter = diamond_inter,
        intron_counts = ic_list
    )
    
    dup_binary <- classify_gene_pairs(
        annotation = annotation,
        blast_list = diamond_intra,
        scheme = "binary"
    )
    
    expect_equal(class(dup_full), "list")
    expect_equal(class(dup_full[[1]]), "data.frame")
    expect_equal(ncol(dup_full[[1]]), 3)
    
    expect_equal(class(dup_binary), "list")
    expect_equal(class(dup_binary[[1]]), "data.frame")
    expect_equal(ncol(dup_binary[[1]]), 3)
    
    # 2) classify_genes()
    dup_genes <- classify_genes(dup_full)
    
    expect_equal(class(dup_genes[[1]]), "data.frame")
    expect_equal(class(dup_genes), "list")
    expect_equal(ncol(dup_genes[[1]]), 2)
})

