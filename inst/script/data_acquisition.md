Data acquisition
================

# Data in data/

## sce_annotation.rda

The organism will be the budding yeast *Saccharomyces cerevisiae strain
S288C*, and the data will be obtained from pico-PLAZA 3.0. Here, we will
process the proteomes and annotation with `syntenet`.

``` r
library(syntenet)
library(Biostrings)

#----1) Get and clean annotation data-------------------------------------------
## Read annotation
annotation <- rtracklayer::import(
    "ftp://ftp.psb.ugent.be/pub/plaza/plaza_pico_03/GFF/sac/annotation.selected_transcript.exon_features.sac.gff3.gz"
)

## Remove non-coding genes
coding <- as.data.frame(readr::read_delim(
    "ftp://ftp.psb.ugent.be/pub/plaza/plaza_pico_03/Annotation/annotation.selected_transcript.sac.csv.gz", skip = 8, delim = ";", show_col_types = FALSE
))
coding <- coding[coding$type == "coding", 1]

## Filter GRanges (include only gene ranges and relevant metadata)
annotation <- annotation[annotation$type == "gene" & 
                             annotation$ID %in% coding]
annotation$score <- NULL
annotation$phase <- NULL
annotation$SGD <- NULL
annotation$pid <- NULL
annotation$Uniprot <- NULL
annotation$Name <- NULL
annotation$Parent <- NULL
annotation$old_gi <- NULL

sce_annotation <- GenomicRanges::GRangesList(Sce = annotation)

#----2) Get and clean proteome data---------------------------------------------
proteome <- readAAStringSet(
    "ftp://ftp.psb.ugent.be/pub/plaza/plaza_pico_03/Fasta/proteome.selected_transcript.sac.fasta.gz"
)
names(proteome) <- gsub(".* \\| ", "", names(proteome))
sce_proteome <- list(Sce = proteome)


#----3) Process data------------------------------------------------------------
check_input(sce_proteome, sce_annotation)
sce_pdata <- process_input(sce_proteome, sce_annotation)
sce_annotation <- sce_pdata$annotation

# Save data
usethis::use_data(sce_annotation, compress = "xz", overwrite = TRUE)
```

## sce_diamond.rda

Here, we will create a list containing a data frame of DIAMOND tabular
output. As above, we will get the data for *S. cerevisiae* from
pico-PLAZA.

``` r
# Run DIAMOND
sce_diamond <- run_diamond(sce_pdata$seq, ... = "--sensitive")

usethis::use_data(sce_diamond, compress = "xz", overwrite = TRUE)
```

## sce_duplicates.rda

This is a 2-column data frame with all duplicate gene pairs (from
DIAMOND output).

``` r
# Filter by e-value, and remove duplicate and redundant entries
data(sce_diamond)
sce_duplicates <- lapply(sce_diamond, function(x) {
    fpair <- x[x$evalue <= 1e-10, 1:2]
    fpair <- fpair[fpair[, 1] != fpair[, 2], ]
    fpair <- fpair[!duplicated(t(apply(fpair, 1, sort))),]
    names(fpair) <- c("dup1", "dup2")
    return(fpair)
})
sce_duplicates <- sce_duplicates[[1]]

usethis::use_data(sce_duplicates, compress = "xz", overwrite = TRUE)
```

## sce_anchors.rda

This is a 2-column data frame with anchor pairs for S. cerevisiae.

``` r
data(sce_diamond)
data(sce_annotation)
blast_list <- sce_diamond
annotation <- sce_annotation
sce_anchors <- get_anchors_list(blast_list, annotation)[[1]]

usethis::use_data(sce_anchors, compress = "xz", overwrite = TRUE)
```

## gma_dups_kaks.rda

``` r
urls <- c(
    WGD = "https://raw.githubusercontent.com/almeidasilvaf/GmPR1/main/data/duplicated_genes_kaks/wgd_kaks.txt",
    TD = "https://raw.githubusercontent.com/almeidasilvaf/GmPR1/main/data/duplicated_genes_kaks/td_kaks.txt",
    PD = "https://raw.githubusercontent.com/almeidasilvaf/GmPR1/main/data/duplicated_genes_kaks/pd_kaks.txt",
    TRD = "https://raw.githubusercontent.com/almeidasilvaf/GmPR1/main/data/duplicated_genes_kaks/trd_kaks.txt",
    DD = "https://raw.githubusercontent.com/almeidasilvaf/GmPR1/main/data/duplicated_genes_kaks/dd_kaks.txt"
)

gma_dups_kaks <- Reduce(rbind, lapply(seq_along(urls), function(x) {
    pairs <- read.csv(urls[x], header = TRUE, sep = "\t", skip = 1)[, c(1:5)]
    names(pairs) <- c("dup1", "dup2", "Ka", "Ks", "Ka_Ks")
    mode <- names(urls)[x]
    pairs$mode <- mode
    return(pairs)
}))

usethis::use_data(gma_dups_kaks, compress = "xz", overwrite = TRUE)
```
