context("Run scPNMF")
library(scPNMF)

test_that("scPNMF works", {
  data(zheng4)
  
  X <- zheng4@assays@data@listData$logcounts[1:200, ]
  
  res_pnmf <- scPNMF::PNMFfun(X = X,
                              K = 10,
                              method="EucDist", tol=1e-4, maxIter=1000,
                              verboseN = TRUE)
  W <- res_pnmf$Weight
  S <- res_pnmf$Score
  
  #rownames(W) <- rownames(zheng4)[1:1000]
  #colnames(W) <- paste0("Basis", 1:10)
  
  expect_equal(is.matrix(W), TRUE)
  
  
  res_annotation <- scPNMF::basisAnnotate(W = W, 
                                          dim_use = 1:10,
                                          id_type = "ENSEMBL")
  #print(res_annotation)
  #expect_equal(class(res_annotation), "compareClusterResult")
  
  
  res_test <- scPNMF::basisTest(S,
                                X,
                                return_fig = TRUE,
                                mc.cores = 1)
  expect_equal(length(res_test), 3)
  
  
  res_select <- scPNMF::basisSelect(W,
                                    S,
                                    X,mc.cores = 1)
  expect_equal(is.matrix(res_select), TRUE)             
  
  
  res_gene <- scPNMF::getInfoGene(W = res_select)
  expect_equal(is.vector(res_gene$InfoGene), TRUE)
  
})