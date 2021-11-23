#' Title \eqn{K} selection based on *dev.ortho*
#' 
#' @description 
#' A low rank \eqn{K} selection method based on *dev.ortho*, the normalized difference between \eqn{\mathbf{W}^T\mathbf{W}} and \eqn{\mathbf{I}}. It runs the PNMF algorithm on every candidate \eqn{K} provided by the user, and find the elbow point \eqn{K} as the suggested one. 
#' 
#' However, this process can be time-consuming as it needs to run the PNMF algorithm on every candidate \eqn{K}. Empirically we recommend the users to skip this section, and directly use \eqn{K=20} in function [PNMFfun()], as we found it reaches stability for most scRNA-seq data.
#' 
#'
#' @param X Input data matrix, where rows represent features (genes), columns represent samples (cells). Rownames (gene names) and colnames (cell names) are required.
#' @param K.seq An integer vector of candidate \eqn{K}'s. All elements need to be distinct integers, and the total number needs to be no fewer than 5.
#' @param ncores The number of cores to use.
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
#' @return A list with components:
#' \describe{
#'   \item{\code{BestK}}{The selected elbow point \eqn{K}.}
#'   \item{\code{BestRes}}{The list containing the \code{Weight} and \code{Score} returned by the PNMF algorithm corresponding to \code{BestK}.}
#'   \item{\code{plot.dev.ortho}}{\code{dev.ortho} plot.}
#'   \item{\code{AllRes}}{A list containing all PNMF outputs for all candidate \eqn{K}'s in \code{K.seq}.}
#' }
#' @importFrom parallel mclapply
#' @importFrom irlba irlba
#' @importFrom akmedoids elbow_point
#' @import ggplot2
#' @export
#' 
#' @author Kexin Li, \email{aileenlikexin@@outlook.com}
#'
#' @examples
K_selection <- function(X, K.seq = 10L:30L, ncores = 1, 
                        tol = 1e-3, 
                        maxIter = 500, 
                        verboseN = FALSE, 
                        zerotol = 1e-10, 
                        method = "EucDist", 
                        label = NULL, 
                        mu = 1, 
                        lambda = 0.01, 
                        seed = 123)  {
  
  dtmat <- as.matrix(X)
  if (is.null(rownames(dtmat))) {
    stop("Gene names missing!")
  }
  if (is.null(colnames(dtmat))) {
    stop("Cell names missing!")
  }
  if (!all(K.seq == round(K.seq))) {
    K.seq <- round(K.seq)
    warning("K.seq are not all integers!")
  }
  if (length(unique(K.seq)) <= 4) {
    stop("More candidate K's are needed to determine the elbow point!")
  }
  
  
  ### run all K case
  res_list <- parallel::mclapply(X = K.seq, function(k, dtmat, tol, maxIter, verboseN, zerotol, method, label, mu, lambda, seed) {
    tryCatch({
      res_pnmfEuc <- scPNMF::PNMFfun(X = dtmat, K = k,
                             tol = tol, maxIter = maxIter, verboseN = verboseN, zerotol = zerotol,
                             method = method, label = label, mu = mu, lambda = lambda, seed = seed)
      res_pnmfEuc
    },
    error = function(e) {
      list(Weight = NULL, Score = NULL)
    })
  }, dtmat = dtmat, tol = tol, maxIter = maxIter, verboseN = verboseN, zerotol = zerotol,
  method = method, label = label, mu = mu, lambda = lambda, seed = seed, 
  mc.cores = ncores)
  # a list with length = length(K.seq); every element = list(Weight, Score).
  
  
  ### Calculate dev.ortho for every K
  dev.ortho.list <- parallel::mclapply(1:length(K.seq), function(idx, K.seq, res_list){
    k1 <- K.seq[idx]
    k2 <- dim(res_list[[idx]]$Weight)[2]
    if (k2 != k1) {
      return(c(k1, NA))
    } else {
      bct <- crossprod(res_list[[idx]]$Weight)
      diag(bct) <- 1
      norm_bct <- norm(bct - diag(k1))
      norm_bct_k <- norm_bct/(k1^2 - k1)
      return(c(k1,norm_bct_k))
    }
  }, K.seq = K.seq, res_list = res_list, mc.cores = ncores)
  dev.ortho.df <- data.frame(matrix(unlist(dev.ortho.list), ncol=2, byrow=TRUE))
  colnames(dev.ortho.df) <- c("K", "dev.ortho")
  
  ### Find the elbow point of K
  # elbowpoints FROM akmedoids
  dev.ortho.df.naomit <- na.omit(dev.ortho.df)
  ep <- akmedoids::elbow_point(dev.ortho.df.naomit[, 1], dev.ortho.df.naomit[, 2])$x
  ep <- round(ep) # the elbow point K value
  ptsize <- rep(1, nrow(dev.ortho.df.naomit))
  ptsize[dev.ortho.df.naomit$K == ep] <- 3
  dev.ortho.df.naomit$ptsize <- ptsize
  
  p.ortho <- ggplot(dev.ortho.df.naomit, aes(x=K, y=dev.ortho)) + geom_line(color = "#FC4E07") + geom_point(aes(size=ptsize), color = "#FC4E07") + 
    ggtitle("Deviation from Orthogonality") + theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5), legend.position="none") 
  
  ### Output all results
  besti <- which(K.seq == ep)
  if (length(besti) == 1) {
    best.list <- list(Weight = res_list[[besti]]$Weight, Score = res_list[[besti]]$Score)
  } else {
    best.list <- scPNMF::PNMFfun(X = dtmat, K = ep,
                                 tol = tol, maxIter = maxIter, verboseN = verboseN, zerotol = zerotol,
                                 method = method, label = label, mu = mu, lambda = lambda, seed = seed)
  }
  names(res_list) <- sapply(K.seq, function(k) {paste0("K = ", k)})
  
  
  
  ### Return the results
  return(list(BestK = ep, BestRes = best.list, plot.dev.ortho = p.ortho, AllRes = res_list))
  
  
}





