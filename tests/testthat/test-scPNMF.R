context("Run scPNMF")
library(scPNMF)

test_that("scPNMF works", {
  data(zheng4)
  X <- zheng4@assays@data@listData$logcounts
  
  res_pnmf <- scPNMF::PNMFfun(X = X,
                              rank = 10,
                              method = "EucDist") 
  
})
