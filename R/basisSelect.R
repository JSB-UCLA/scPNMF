#' @title Perform basis selection
#' 
#' @description Wrapper of basis annotation and data-driven basis selection.
#'
#' @param W The weight matrix output by PNMF
#' @param S The score matrix output by PNMF
#' @param X The input gene by cell logcount matrix
#' @param toTest Whether to select bases by Pearson correlation w/ cell library size and test of multimodality
#' @param cor_thres Pearson correlation w/ cell library size cutoff. Default is 0.7.
#' @param pval_thres Adjusted p-value cutoff on test of multimodality. Default is 0.01.
#' @param return_fig Whether to print scatter plot of score vector against cell library size and distribution of score vectors for each basis.
#' @param adj_method P-value correction method. Default is "BH".
#' @param mc.cores The number of cores to use for function mclapply().
#' @param ncol Columns for facets in plots.
#' @param toAnnotate Whether to perform Gene Ontology (GO) enrichment analysis on each basis
#' @param dim_use The bases (columns) to be used in the selected weight matrix, \code{NULL} value uses all bases
#'
#'
#' @return The selected weight matrix
#' @export 
#'
#' 
#' 
#' 
basisSelect <- function(W, S = NULL, X = NULL, 
                        toTest = TRUE, cor_thres = 0.7, pval_thres = 0.01, 
                        return_fig = FALSE, adj_method = "BH", mc.cores = 10, ncol = 4,
                        toAnnotate = FALSE, dim_use = NULL
) {
  
  dim_selected <- 1:(dim(W)[2])
  if (toAnnotate & !is.null(dim_use)) {
    dim_selected <- dim_use
  }
  
  if (toTest) {
    stopifnot(!is.null(S) & !is.null(X))
    resTest <- basisTest(S = S, X = X, return_fig = return_fig, adj_method = adj_method, mc.cores = mc.cores, ncol = ncol)
    dim_Test <- which(resTest$cor_value <= cor_thres | resTest$adj_pvalue <= pval_thres)
    dim_selected <- intersect(dim_selected, dim_Test)
    if (return_fig) {
      print(resTest$p_scatter)
      print(resTest$p_hist)
    }
  }
  
  W_s <- W[, dim_selected]
  
  return(W_s)
}