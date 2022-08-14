
#' Protein sequences of the yeast species S. cerevisiae and S. pombe
#'
#' Data obtained from Pico-PLAZA 3.0. Only translated sequences of primary
#' transcripts were included. 
#' 
#' @name yeast_seq
#' @format A list of AAStringSet objects with the elements
#' \strong{Scerevisiae} and \strong{Spombe}.
#' @references
#' Van Bel, M., Silvestri, F., Weitz, E. M., Kreft, L., Botzki, A.,
#' Coppens, F., & Vandepoele, K. (2021). PLAZA 5.0: extending the scope
#' and power of comparative and functional genomics in plants.
#' Nucleic acids research.
#' @examples
#' data(yeast_seq)
#' @usage data(yeast_seq)
"yeast_seq"


#' Genome annotation of the yeast species S. cerevisiae and S. pombe
#'
#' Data obtained from Pico-PLAZA 3.0. Only annotation data for primary
#' transcripts were included.
#' 
#' @name yeast_annot
#' @format A CompressedGRangesList containing 
#' the elements \strong{Scerevisiae} and \strong{Spombe}.
#' @references
#' Van Bel, M., Silvestri, F., Weitz, E. M., Kreft, L., Botzki, A.,
#' Coppens, F., & Vandepoele, K. (2021). PLAZA 5.0: extending the scope
#' and power of comparative and functional genomics in plants.
#' Nucleic acids research.
#' @examples
#' data(yeast_annot)
#' @usage data(yeast_annot)
"yeast_annot"


#' List of data frames containing intraspecies DIAMOND output for S. cerevisiae
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


#' List of data frames containing interspecies DIAMOND output for yeast species
#' 
#' This list contains a similarity search of S. cerevisiae against
#' S. pombe, and it was obtained with \code{run_diamond()}.
#' 
#' @name diamond_inter
#' @format A list of data frames (length 1) containing the output of a 
#' DIAMOND search of S. cerevisiae against S. pombe (outgroup).
#' @examples 
#' data(diamond_inter)
#' @usage data(diamond_inter)
"diamond_inter"

