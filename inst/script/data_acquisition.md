Data acquisition
================

# Data in data/

Here, we will use genome data for two yeast species:

-   *Saccharomyces cerevisiae*
-   *Schizosaccahromyces pombe*

Data will be obtained from Pico-PLAZA 3.0.

First of all, let’s obtain a list of only protein-coding genes for each
species.

``` r
# Get character vector of protein coding gene IDs
## S. cerevisiae
scerevisiae_coding <- as.data.frame(readr::read_delim(
    "ftp://ftp.psb.ugent.be/pub/plaza/plaza_pico_03/Annotation/annotation.selected_transcript.sac.csv.gz", skip = 8, delim = ";", show_col_types = FALSE
))
scerevisiae_coding <- scerevisiae_coding[scerevisiae_coding$type == "coding", 1]

## S. pombe
spombe_coding <- as.data.frame(readr::read_delim(
    "ftp://ftp.psb.ugent.be/pub/plaza/plaza_pico_03/Annotation/annotation.selected_transcript.scp.csv.gz", skip = 8, delim = ";", show_col_types = FALSE
))
spombe_coding <- spombe_coding[spombe_coding$type == "coding", 1]
```

## yeast_seq.rda

The object `yeast_seq` is a list of `AAStringSet` objects with elements
*Scerevisiae* and *Spombe*. Only translated sequences for primary
transcripts (protein-coding only) are included.

``` r
library(Biostrings)

# Get proteome data
scerevisiae_proteome <- readAAStringSet(
    "ftp://ftp.psb.ugent.be/pub/plaza/plaza_pico_03/Fasta/proteome.selected_transcript.sac.fasta.gz"
)
spombe_proteome <- readAAStringSet(
    "ftp://ftp.psb.ugent.be/pub/plaza/plaza_pico_03/Fasta/proteome.selected_transcript.scp.fasta.gz"
)
names(scerevisiae_proteome) <- gsub(".* \\| ", "", names(scerevisiae_proteome))
names(spombe_proteome) <- gsub(".* \\| ", "", names(spombe_proteome))

# Remove non-coding genes
scerevisiae_proteome <- scerevisiae_proteome[names(scerevisiae_proteome) %in%
                                                 scerevisiae_coding, ]
spombe_proteome <- spombe_proteome[names(spombe_proteome) %in% 
                                       spombe_coding]

# Store AAStringSet objects in a list
yeast_seq <- list(
    Scerevisiae = scerevisiae_proteome,
    Spombe = spombe_proteome
)

# Save object
usethis::use_data(yeast_seq, compress = "xz", overwrite = TRUE)
```

## yeast_annot.rda

The object `yeast_annot` is a `GRangesList` object with elements
*Scerevisiae* and *Spombe*. Only ranges for primary transcripts
(protein-coding only) are included.

``` r
library(rtracklayer)

# Get gene ranges
scerevisiae_annot <- import(
    "ftp://ftp.psb.ugent.be/pub/plaza/plaza_pico_03/GFF/sac/annotation.selected_transcript.exon_features.sac.gff3.gz"
)
spombe_annot <- import(
    "ftp://ftp.psb.ugent.be/pub/plaza/plaza_pico_03/GFF/scp/annotation.selected_transcript.exon_features.scp.gff3.gz"
)

# Filter GRanges (include only protein-coding genes and relevant metadata)
## S. cerevisiae
scerevisiae_annot <- scerevisiae_annot[scerevisiae_annot$type == "gene" & 
                             scerevisiae_annot$ID %in% scerevisiae_coding]
scerevisiae_annot$source <- NULL
scerevisiae_annot$score <- NULL
scerevisiae_annot$phase <- NULL
scerevisiae_annot$SGD <- NULL
scerevisiae_annot$pid <- NULL
scerevisiae_annot$Uniprot <- NULL
scerevisiae_annot$Name <- NULL
scerevisiae_annot$ID <- NULL
scerevisiae_annot$Parent <- NULL
scerevisiae_annot$old_gi <- NULL

## S. pombe
spombe_annot <- spombe_annot[spombe_annot$type == "gene" & 
                             spombe_annot$ID %in% spombe_coding]
spombe_annot$source <- NULL
spombe_annot$score <- NULL
spombe_annot$phase <- NULL
spombe_annot$SGD <- NULL
spombe_annot$pid <- NULL
spombe_annot$Uniprot <- NULL
spombe_annot$Name <- NULL
spombe_annot$ID <- NULL
spombe_annot$Parent <- NULL
spombe_annot$old_gi <- NULL


# Combine GRanges objects in a GRangesList
yeast_annot <- GenomicRanges::GRangesList(
    Scerevisiae = scerevisiae_annot,
    Spombe = spombe_annot
)

# Save data
usethis::use_data(yeast_annot, compress = "xz", overwrite = TRUE)
```

## diamond_intra.rda and diamond_inter.rda

The object `diamond_intra` is a list of DIAMOND data frames for
intraspecies comparisons of *S. cerevisiae*, while `diamond_inter`
contains the DIAMOND output of a comparison between *S. cerevisiae* and
*S. pombe*.

``` r
# Load and process data
data(yeast_seq)
data(yeast_annot)

pdata <- process_input(yeast_seq, yeast_annot)

# Intraspecies DIAMOND
diamond_intra <- run_diamond(
    seq = pdata$seq["Scerevisiae"],
    compare = "intraspecies", 
    outdir = file.path(tempdir(), "diamond_intra_example"),
    ... = "--sensitive"
)

# Interspecies DIAMOND
comparisons <- data.frame(
    species = "Scerevisiae",
    outgroup = "Spombe"
)

diamond_inter <- run_diamond(
    seq = pdata$seq,
    compare = comparisons,
    outdir = file.path(tempdir(), "diamond_inter_example"),
    ... = "--sensitive"
)

# Save data
usethis::use_data(diamond_intra, compress = "xz", overwrite = TRUE)
usethis::use_data(diamond_inter, compress = "xz", overwrite = TRUE)
```
