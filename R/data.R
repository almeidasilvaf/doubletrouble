
#' Processed genome annotation for Saccharomyces cerevisiae
#'
#' Data obtained from Pico-PLAZA 3.0. Only annotation data for primary
#' transcripts were included. The data was processed with 
#' syntenet::process_input().
#' 
#' @name sce_annotation
#' @format A CompressedGRangesList containing 
#' the element `Sce`.
#' @references
#' Van Bel, M., Silvestri, F., Weitz, E. M., Kreft, L., Botzki, A.,
#' Coppens, F., & Vandepoele, K. (2021). PLAZA 5.0: extending the scope
#' and power of comparative and functional genomics in plants.
#' Nucleic acids research.
#' @examples
#' data(sce_annotation)
#' @usage data(sce_annotation)
"sce_annotation"


#' List of data frames containing BLAST-like tabular output
#'
#' Only the DIAMOND output for S. cerevisiae against itself was included.
#' Here, we used DIAMOND in sensitive mode.
#' 
#' @name sce_diamond
#' @format A list of data frames (length 1) containing the whole paranome of
#' S. cerevisiae resulting from intragenome similarity searches.
#' @examples 
#' data(sce_diamond)
#' @usage data(sce_diamond)
"sce_diamond"


#' All duplicate gene pairs for S. cerevisiae
#'
#' Duplicate gene pairs were obtained from the output of DIAMOND in sensitive
#' mode with an E-value filter of 1e-10.
#' 
#' @name sce_duplicates
#' @format A 2-column data frame with each gene of the duplicate pair in each
#' column. Columns are named "dup1" and "dup2", respectively.
#' @examples 
#' data(sce_duplicates)
#' @usage data(sce_duplicates)
"sce_duplicates"


#' Anchor pairs in the S. cerevisiae genome
#'
#' Anchor pairs were obtained with \code{get_anchors_list()}.
#'
#' @name sce_anchors
#' @format A 2-column data frame with each gene of the anchor pair in each 
#' column. Columns are named "Anchor1" and "Anchor2".
#' @examples 
#' data(sce_anchors)
#' @usage data(sce_anchors)
"sce_anchors"


#' Soybean duplicated gene pairs and their Ka, Ks, and Ka/Ks values
#'
#' @name gma_dups_kaks
#' @format A data frame with the following columns:
#' \describe{
#'   \item{dup1}{Character, duplicated gene 1.}
#'   \item{dup2}{Character, duplicated gene 2.}
#'   \item{Ka}{Numeric, Ka values.}
#'   \item{Ks}{Numeric, Ks values.}
#'   \item{Ka_Ks}{Numeric, Ka/Ks values.}
#'   \item{mode}{Mode of duplication.}
#' }
#' @references 
#' Almeida-Silva, F., Moharana, K. C., Machado, F. B., & 
#' Venancio, T. M. (2020). Exploring the complexity of soybean (Glycine max) 
#' transcriptional regulation using global gene co-expression networks. 
#' Planta, 252(6), 1-12.
#' @examples 
#' data(gma_dups_kaks)
#' @usage data(gma_dups_kaks)
"gma_dups_kaks"

