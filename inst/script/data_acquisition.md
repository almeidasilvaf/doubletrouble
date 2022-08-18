Data acquisition
================

# Data in data/

Here, we will use genome data for two yeast species:

-   *Saccharomyces cerevisiae*
-   *Candida glabrata*

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
*Scerevisiae* and *Cglabrata*. Only ranges for primary transcripts
(protein-coding only) are included.

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
## S. cerevisiae
scerevisiae_annot <- scerevisiae_annot[scerevisiae_annot$type == "gene" & 
                             scerevisiae_annot$biotype == "protein_coding"]
scerevisiae_annot <- scerevisiae_annot[, c("type", "gene_id")]

## C. glabrata
cglabrata_annot <- cglabrata_annot[cglabrata_annot$type == "gene" & 
                             cglabrata_annot$biotype == "protein_coding"]
cglabrata_annot <- cglabrata_annot[, c("type", "gene_id")]

# Combine GRanges objects in a GRangesList
yeast_annot <- GenomicRanges::GRangesList(
    Scerevisiae = scerevisiae_annot,
    Cglabrata = cglabrata_annot
)

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
