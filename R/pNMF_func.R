#' @title Fitting PNMF Models
#'
#' @description 
#' Fast Projective Nonnegative Matrix Factorization Realizatiton based on Euclidean Distance / KL Divergence / Discriminant pNMF.
#'
#' @param X Input data matrix, where rows represent features (genes), columns represent samples (cells).
#' @param K Specification of the factorization rank (number of low dimension).
#' @param tol A threshold below which would be considered converged. Default is 1e-3.
#' @param maxIter Number of max iteration times. Default is 500.
#' @param verboseN A boolean value indicating whether to print number of iterations.
#' @param zerotol A threshold on basis loadings below which would be considered zero. Default is 1e-10.
#' @param method A character string indicating which method to be used. One of \code{"EucDist"}, \code{"KL"}, or \code{"DPNMF"}.
#' @param label A character vector indicating the cluster type for each cell. Required only when \code{method = "DPNMF"}.
#' @param mu A numerical value which controls the penalty term. Larger \code{mu} represents haivier penalization of class distances in \code{DPNMF}. Default is 1.
#' @param lambda A numerical value which controls the magnituide of within class distances. Larger \code{lambda} represents larger proportion of within class distances in the total penalty term. Default is 0.01.
#' @param seed Random seed of the initialization.
#' 
#' @details 
#' Given a data matrix (rows as features and columns as samples), this function
#' computes the Projective Nonnegative Matrix Factorization (PNMF). Based on different objective functions,
#' the choices are Euclidean distance (\code{"EucDist"}), KL divergence (\code{"KL"}) (Yang, Zhirong, and Erkki Oja. 2010), or Discriminant PNMF (\code{"DPNMF"}) (Guan, Naiyang, et al. 2013). 
#' \code{"EucDist"} is supposed to be the most common one;
#' \code{"KL"} is similar to KL-NMF (Poisson-NMF), and may work better for count data; 
#' \code{"DPNMF"} requires the predefined labels.
#' 
#' The fitting result of PNMF shares characteristics of both PCA and NMF. The model returns
#' a \code{basis} matrix, which is similar to the loading matrix in PCA. However, notice that
#' unlike in PCA the first PC always represents the largest variation, each basis vectors are equal in PNMF.
#' 
#' @return A list with components:
#' \describe{
#'   \item{\code{Weight}}{The basis of model fit (\eqn{W}).}
#'   \item{\code{Score}}{The mapped scores (\eqn{W^TX}), which is the dimension reduced result.}
#' }
#' 
#' @importFrom irlba irlba
#' @export
#' 
#' @author Kexin Li, \email{aileenlikexin@@outlook.com}
#' @author Dongyuan Song, \email{dongyuansong@@g.ucla.edu}
#'
#' @references
#' \itemize{
#' \item Yang, Z., & Oja, E. (2010). Linear and nonlinear projective nonnegative matrix factorization. IEEE Transactions on Neural Networks, 21(5), 734-749.
#' \item Guan, N., Zhang, X., Luo, Z., Tao, D., & Yang, X. (2013). Discriminant projective non-negative matrix factorization. PloS one, 8(12), e83291.
#' \item \url{https://github.com/richardbeare/pNMF}
#' }
#' 
#'
#' @examples 
#' data(zheng4)
#' X <- SummarizedExperiment::assay(zheng4, "logcounts")
#' res_pnmf <- scPNMF::PNMFfun(X = X,
#'                            K = 3,
#'                            method="EucDist", 
#'                            tol=1e-4, 
#'                            maxIter=100,
#'                            verboseN = FALSE)
PNMFfun <- function(X, 
                    K=10, 
                    tol=1e-3, 
                    maxIter=500, 
                    verboseN=FALSE, 
                    zerotol=1e-10, 
                    method="EucDist", 
                    label=NULL, 
                    mu=1, 
                    lambda=0.01, 
                    seed=123) {
  if (is.null(rownames(X))) {
    stop("Gene names missing!")
  }
  if (is.null(colnames(X))) {
    stop("Cell names missing!")
  }
  
  #nmfmod <- NMF::nmf(X, rank)
  set.seed(seed)
  Init <- irlba(X, nv = K)
  Winit <- Init$u
  Winit <- abs(Winit)
  
  if (method == "EucDist") {
    W <- PNMF_EucDistC(X, Winit, tol, maxIter, verboseN, zerotol) 
    W <- W/norm(W, "2")
    ld <- t(X) %*% W
  }
  else if (method == "KL") {
    W <- PNMF_KLC(X, Winit, tol, maxIter, verboseN, zerotol) 
    W <- W/norm(W, "2")
    ld <- t(X) %*% W
  }
  else if (method == "DPNMF") {
    if (length(label) != dim(X)[2]) {
      stop("Cluster labels must have same length as number of cells.")
    }
    cluvec <- as.factor(label)
    cluvec.num <- as.numeric(cluvec)
    cluvec.ord <- order(cluvec.num)
    Xord <- X[, cluvec.ord]
    clunum <- as.integer(table(cluvec))
    
    W <- DPNMFC(X, Winit, tol, maxIter, verboseN, zerotol, Xord, clunum, mu, lambda) 
    W <- W/norm(W, "2")
    ld <- t(X) %*% W
  }
  
  rownames(W) <- rownames(X)
  rownames(ld) <- colnames(X) 
  colnames(W) <- paste0("Basis", 1:K)
  colnames(ld) <- paste0("Basis", 1:K)
  
  return(list(Weight=W, Score=ld))
}







