
#' Protein sequences of the yeast species S. cerevisiae and C. glabrata
#'
#' Data obtained from Ensembl Fungi. Only translated sequences of primary
#' transcripts were included. 
#' 
#' @name yeast_seq
#' @format A list of AAStringSet objects with the elements
#' \strong{Scerevisiae} and \strong{Cglabrata}.
#' @examples
#' data(yeast_seq)
#' @usage data(yeast_seq)
"yeast_seq"


#' Genome annotation of the yeast species S. cerevisiae and C. glabrata
#'
#' Data obtained from Ensembl Fungi. Only annotation data for primary
#' transcripts were included.
#' 
#' @name yeast_annot
#' @format A CompressedGRangesList containing 
#' the elements \strong{Scerevisiae} and \strong{Cglabrata}.
#' @examples
#' data(yeast_annot)
#' @usage data(yeast_annot)
"yeast_annot"


#' Intraspecies DIAMOND output for S. cerevisiae
#'
#' List obtained with \code{run_diamond()}.
#' 
#' @name diamond_intra
#' @format A list of data frames (length 1) containing the whole paranome of
#' S. cerevisiae resulting from intragenome similarity searches.
#' @examples 
#' data(diamond_intra)
#' @usage data(diamond_intra)
"diamond_intra"


#' Interspecies DIAMOND output for yeast species
#' 
#' This list contains a similarity search of S. cerevisiae against
#' C. glabrata, and it was obtained with \code{run_diamond()}.
#' 
#' @name diamond_inter
#' @format A list of data frames (length 1) containing the output of a 
#' DIAMOND search of S. cerevisiae against C. glabrata (outgroup).
#' @examples 
#' data(diamond_inter)
#' @usage data(diamond_inter)
"diamond_inter"


#' Coding sequences (CDS) of S. cerevisiae
#' 
#' Data were obtained from Ensembl Fungi, and only CDS of primary transcripts
#' were included.
#' 
#' @name cds_scerevisiae 
#' @format A DNAStringSet object with CDS of S. cerevisiae.
#' @examples 
#' data(cds_scerevisiae)
#' @usage data(cds_scerevisiae)
"cds_scerevisiae"


#' Duplicate pairs and Ka, Ks, and Ka/Ks values for S. cerevisiae
#'
#' This data set was obtained with \code{classify_gene_pairs()} followed
#' by \code{pairs2kaks()}.
#' 
#' @name scerevisiae_kaks
#' @format A data frame with the following variables:
#' \describe{
#'   \item{dup1}{Character, duplicated gene 1.}
#'   \item{dup2}{Character, duplicated gene 2.}
#'   \item{Ka}{Numeric, Ka values.}
#'   \item{Ks}{Numeric, Ks values.}
#'   \item{Ka_Ks}{Numeric, Ka/Ks values.}
#'   \item{type}{Character, mode of duplication}
#' }
#' @examples 
#' data(scerevisiae_kaks)
#' @usage data(scerevisiae_kaks)
"scerevisiae_kaks"


#' Duplicate pairs and Ks values for Glycine max
#'
#' This data set was obtained with \code{classify_gene_pairs()} followed
#' by \code{pairs2kaks()}.
#' 
#' @name gmax_ks
#' @format A data frame with the following variables:
#' \describe{
#'   \item{dup1}{Character, duplicated gene 1.}
#'   \item{dup2}{Character, duplicated gene 2.}
#'   \item{Ks}{Numeric, Ks values.}
#' }
#' @examples 
#' data(gmax_ks)
#' @usage data(gmax_ks)
"gmax_ks"

