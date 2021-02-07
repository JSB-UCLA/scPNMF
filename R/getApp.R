#' Get informative genes
#'
#' @param W The weight matrix output by PNMF
#' @param M The user-defined informative gene number
#' @param by_basis Return informative genes by basis or not 
#' @param return_trunW Return the truncated weight matrix or not 
#' @param dim_use The bases (columns) to be used in the selected weight matrix, \code{NULL} value uses all bases
#'
#' @return A vector containing M top informative genes
#' @export
#'
#'
#'
getInfoGene <- function(W, M = 100, by_basis = FALSE, return_trunW = FALSE, dim_use = NULL) {
  W_s <- W
  if (!is.null(dim_use)) {
    W_s <- W[, dim_use]
  }
  ResMTrun <- .MTruncation(W_s = W_s, M = M)
  W_sM <- ResMTrun$W_sM
  output_W_sM <- NULL
  if (return_trunW) {
    output_W_sM <- W_sM
  }
  
  if (by_basis) {
    output_InfoGene <- lapply(1:dim(W_sM)[2], function(k) {
      ResMTrun$InfoGene[W_sM[,k] > 0]
    })
    names(output_InfoGene) <- colnames(W)
    return(list(InfoGene = output_InfoGene, trunW = output_W_sM))
  } else {
    return(list(InfoGene = ResMTrun$InfoGene, trunW = output_W_sM))
  }
}

#' Get data projection
#'
#' @param X The input gene by cell logcount matrix
#' @param W_sM The truncated, selected weight matrix
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

