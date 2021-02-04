.MTruncation <- function(W_s, M = 100) {
  stopifnot(!is.null(rownames(W_s)))
  stopifnot(M <= dim(W_s)[1])
  
  row_max <- apply(W_s, 1, max)
  row_sum <- apply(W_s, 1, sum)
  wM <- sort(row_max, decreasing = TRUE)[M]
  orderM <- order(row_max, row_sum, decreasing = TRUE)
  idxM <- sort(orderM[1:M])
  
  InfoGene <- rownames(W_s)[idxM]
  W_sM <- W_s[idxM, ]
  W_sM[W_sM < wM] <- 0
  
  return(list(InfoGene = InfoGene, W_sM = W_sM))
}


#' Get informative genes
#'
#' @param W_s 
#' @param M 
#'
#' @return A vector containing M top informative genes
#' @export
#'
#'
#'
getInfoGene <- function(W, M = 100, by_basis = FALSE, dim_use = NULL) {
  W_s <- W
  if (!is.null(dim_use)) {
    W_s <- W[, dim_use]
  }
  ResMTrun <- MTruncation(W_s = W_s, M = M)
  W_sM <- ResMTrun$W_sM
  
  if (by_basis) {
    return(lapply(1:dim(W_sM)[2], function(k) {
      ResMTrun$InfoGene[W_sM[,k] > 0]
      }))
  } else {
    return(ResMTrun$InfoGene)
  }
}

#' Get data projection
#'
#' @param X 
#' @param W_s 
#' @param M 
#'
#' @return A K0 * cells matrix representing low-dimensional data projection
#' @export
#'
#'
getProjection <- function(X, W_sM) {
  stopifnot(!is.null(rownames(X)))
  stopifnot(!is.null(rownames(W_sM)))
  
  X_M <- X[rownames(X) %in% rownames(W_sM), ]
  S_M <- t(W_sM) %*% X_M
  
  return(S_M)
}