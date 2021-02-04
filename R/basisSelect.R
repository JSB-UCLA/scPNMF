#' Basis selection function
#'
#' @param W 
#' @param S 
#' @param X 
#' @param toTest 
#' @param cor_thres 
#' @param pval_thres 
#' @param return_fig 
#' @param adj_method 
#' @param mc.cores 
#' @param nrow 
#' @param toAnnotate 
#' @param dim_use 
#'
#' @return The selected weight matrix
#' @export 
#'
#' 
#' 
#' 
basisSelect <- function(W, S = NULL, X = NULL, 
                        toTest = TRUE, cor_thres = 0.7, pval_thres = 0.01, 
                        return_fig = FALSE, adj_method = "BH", mc.cores = 10, nrow = 4,
                        toAnnotate = FALSE, dim_use = NULL
) {
  
  dim_selected <- 1:(dim(W)[2])
  if (toAnnotate & !is.null(dim_use)) {
    dim_selected <- dim_use
  }
  
  if (toTest) {
    stopifnot(!is.null(S) & !is.null(X))
    resTest <- basisTest(S = S, X = X, return_fig = return_fig, adj_method = adj_method, mc.cores = mc.cores, nrow = nrow)
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