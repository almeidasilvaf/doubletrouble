Data acquisition
================

# Data in data/

Here, we will use genome data for two yeast species:

- *Saccharomyces cerevisiae*
- *Candida glabrata*

Data will be obtained from Ensembl Fungi.

First of all, let’s obtain a list of only protein-coding genes for each
species.

``` r
library(tidyverse)

# Get character vector of protein coding gene IDs
## S. cerevisiae
scerevisiae_coding <- as.data.frame(read_delim(
    "ftp://ftp.psb.ugent.be/pub/plaza/plaza_pico_03/Annotation/annotation.selected_transcript.sac.csv.gz", skip = 8, delim = ";", show_col_types = FALSE
))
scerevisiae_coding <- scerevisiae_coding[scerevisiae_coding$type == "coding", 1]

## C. glabrata
cglabrata_coding <- rtracklayer::import(
    "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-54/fungi/gff3/candida_glabrata/Candida_glabrata.GCA000002545v2.54.gff3.gz"
)
cglabrata_coding <- cglabrata_coding[cglabrata_coding$type == "gene", ]
cglabrata_coding <- cglabrata_coding[cglabrata_coding$biotype == "protein_coding", ]
```

## yeast_annot.rda

The object `yeast_annot` is a `GRangesList` object with elements
*Scerevisiae* and *Cglabrata*. Only ranges for protein-coding genes are
included.

``` r
library(rtracklayer)

# Get gene ranges
scerevisiae_annot <- import(
    "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-54/fungi/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.54.gff3.gz"
)
cglabrata_annot <- import(
    "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-54/fungi/gff3/candida_glabrata/Candida_glabrata.GCA000002545v2.54.gff3.gz"
)

# Filter GRanges (include only protein-coding genes and relevant metadata)
## Combine GRanges objects in a list
yeast_annot <- list(
    Scerevisiae = scerevisiae_annot,
    Cglabrata = cglabrata_annot
)

## Filter data
yeast_annot <- lapply(yeast_annot, function(x) {
    
    # Get ranges for coding genes, and use them to extract exons, mRNA, etc.
    gene_ranges <- x[x$biotype == "protein_coding" & x$type == "gene"]
    cranges <- subsetByOverlaps(x, gene_ranges)
    
    # Remove exons for TEs (to avoid warnings when building TxDb)
    te_tx <- cranges[cranges$type == "transposable_element", ]$transcript_id
    if(length(te_tx) > 0) {
        te_exonid <- paste0(rep(te_tx, each = 9), paste0("-E", 1:9))
        cranges <- cranges[-which(cranges$Name %in% te_exonid)]
    }

    # Remove unnecessary columns (for package size issues)
    cols <- c(
        "type", "phase", "ID", "Parent", "Name", 
        "gene_id", "transcript_id", "exon_id", "protein_id"
    )
    cranges <- cranges[, cols]
    
    return(cranges)
})

yeast_annot <- GenomicRanges::GRangesList(yeast_annot)

# Save data
usethis::use_data(yeast_annot, compress = "xz", overwrite = TRUE)
```

## yeast_seq.rda

The object `yeast_seq` is a list of `AAStringSet` objects with elements
*Scerevisiae* and *Cglabrata*. Only translated sequences for primary
transcripts (protein-coding only) are included.

``` r
library(Biostrings)

# Define small function to keep only longest isoform
ensembl_longest_isoform <- function(proteome = NULL) {

    pnames <- gsub(".*gene:", "", names(proteome))
    pnames <- gsub(" .*", "", pnames)

    names(proteome) <- pnames
    proteome <- proteome[order(Biostrings::width(proteome), decreasing = TRUE),]
    proteome <- proteome[!duplicated(names(proteome)), ]
    return(proteome)
}

# Get proteome data
scerevisiae_proteome <- readAAStringSet(
    "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-54/fungi/fasta/saccharomyces_cerevisiae/pep/Saccharomyces_cerevisiae.R64-1-1.pep.all.fa.gz"
) |> ensembl_longest_isoform()

cglabrata_proteome <- readAAStringSet(
    "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-54/fungi/fasta/candida_glabrata/pep/Candida_glabrata.GCA000002545v2.pep.all.fa.gz"
) |> ensembl_longest_isoform()

# Remove non-coding genes
scerevisiae_proteome <- scerevisiae_proteome[names(scerevisiae_proteome) %in%
                                                 scerevisiae_annot$gene_id, ]

cglabrata_proteome <- cglabrata_proteome[names(cglabrata_proteome) %in% 
                                             cglabrata_annot$gene_id]

# Store AAStringSet objects in a list
yeast_seq <- list(
    Scerevisiae = scerevisiae_proteome,
    Cglabrata = cglabrata_proteome
)

# Save object
usethis::use_data(yeast_seq, compress = "xz", overwrite = TRUE)
```

## diamond_intra.rda and diamond_inter.rda

The object `diamond_intra` is a list of DIAMOND data frames for
intraspecies comparisons of *S. cerevisiae*, while `diamond_inter`
contains the DIAMOND output of a comparison between *S. cerevisiae* and
*C. glabrata*.

``` r
# Load and process data
data(yeast_seq)
data(yeast_annot)

pdata <- process_input(yeast_seq, yeast_annot)

# Intraspecies DIAMOND
diamond_intra <- run_diamond(
    seq = pdata$seq["Scerevisiae"],
    compare = "intraspecies", 
    outdir = file.path(tempdir(), "diamond_intra_data"),
    ... = "--sensitive"
)

# Interspecies DIAMOND
comparisons <- data.frame(
    species = "Scerevisiae",
    outgroup = "Cglabrata"
)

diamond_inter <- run_diamond(
    seq = pdata$seq,
    compare = comparisons,
    outdir = file.path(tempdir(), "diamond_inter_data"),
    ... = "--sensitive"
)

# Save data
usethis::use_data(diamond_intra, compress = "xz", overwrite = TRUE)
usethis::use_data(diamond_inter, compress = "xz", overwrite = TRUE)
```

## cds_scerevisiae.rda

This is a `DNAStringSet object` containing the CDS of duplicated genes
in the S. cerevisiae genome.

``` r
library(Biostrings)

# Get duplicated genes
data(scerevisiae_kaks)
c_full <- scerevisiae_kaks[, c("dup1", "dup2", "type")]

dup_genes <- unique(c(c_full$dup1, c_full$dup2))
dup_genes <- gsub(".*_", "", dup_genes)

dup_sd <- c_full[c_full$type == "SD", ]
dup_sd <- unique(c(dup_sd$dup1, dup_sd$dup2))
dup_sd <- gsub(".*_", "", dup_sd)

# Get CDS and keep only longest isoform
cds_scerevisiae_full <- readDNAStringSet(
    "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-54/fungi/fasta/saccharomyces_cerevisiae/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz"
) |> ensembl_longest_isoform()

# Keep only duplicated genes
cds_scerevisiae <- cds_scerevisiae_full[names(cds_scerevisiae_full) %in% dup_wgd]

usethis::use_data(cds_scerevisiae, compress = "xz", overwrite = TRUE)
```

## scerevisiae_kaks.rda

This object is a data frame of duplicate pairs and their Ks values.

``` r
# Get all duplicated gene pairs
library(Biostrings)

data(yeast_seq)
data(yeast_annot)
data(diamond_intra)
data(diamond_inter)
pdata <- syntenet::process_input(yeast_seq, yeast_annot)

# Classify genes into the extended scheme
c_extended <- classify_gene_pairs(
    blast_list = diamond_intra,
    annotation = pdata$annotation,
    scheme = "extended",
    blast_inter = diamond_inter
)

# Get CDS
cds <- list(Scerevisiae = cds_scerevisiae_all)

# Calculate Ks values
scerevisiae_kaks_list <- pairs2kaks(c_extended, cds)
scerevisiae_kaks <- scerevisiae_kaks_list$Scerevisiae

usethis::use_data(scerevisiae_kaks, compress = "xz", overwrite = TRUE)
```

## gmax_ks.rda

This object is a 3-column data frame of duplicate pairs and their Ks
values for *Glycine max* (soybean).

``` r
# Get data
annot <- rtracklayer::import(
    "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/gtf/glycine_max/Glycine_max.Glycine_max_v2.1.53.gtf.gz"
)
annot <- list(Gmax = annot)

seq <- Biostrings::readAAStringSet(
    "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/fasta/glycine_max/pep/Glycine_max.Glycine_max_v2.1.pep.all.fa.gz"
) |> ensembl_longest_isoform()
seq <- list(Gmax = seq)

cds <- Biostrings::readDNAStringSet(
    "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/fasta/glycine_max/cds/Glycine_max.Glycine_max_v2.1.cds.all.fa.gz"
) |> ensembl_longest_isoform()

# Process data
pdata <- syntenet::process_input(seq, annot)

# Intraspecies comparison
diamond_intra <- run_diamond(
    seq = pdata$seq["Gmax"],
    compare = "intraspecies", 
    outdir = file.path(tempdir(), "diamond_intra_data"),
    ... = "--sensitive"
)

# Binary classification
c_binary <- classify_gene_pairs(
    blast_list = diamond_intra,
    annotation = pdata$annotation,
    binary = TRUE
)

cds <- list(Gmax = cds)

# Calculate Ks values
gmax_kaks_list <- pairs2kaks(c_binary, cds)
gmax_ks <- gmax_kaks_list$Gmax
gmax_ks <- gmax_ks[, c("dup1", "dup2", "Ks")]

gmax_ks <- gmax_ks[gmax_ks$Ks <= 2, ]
gmax_ks <- gmax_ks[!is.na(gmax_ks$Ks), ]

usethis::use_data(gmax_ks, compress = "xz", overwrite = TRUE)
```
