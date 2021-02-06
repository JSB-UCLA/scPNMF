#' @title Annotate the bases of scPNMF
#' @description Perform Gene Ontology (GO) enrichment analysis on each basis to provide biological interpretation
#' 
#' @param W The weight matrix from scPNMF.
#' @param DE_prop The proportion of genes to be treated as "DE" genes for GO analysis.
#' @param dim_use Which bases to annotate.
#' @param id_type The gene id types. For example, official gene symbol is \code{SYMBOL}. See details in package \code{clusterProfiler}.
#' @param ont Which gene ontology to use. One of "BP", "MF" or "CC". The default is "BP".
#' @param OrgDb OrgDb to use.
#' @param pvalueCutoff P-value cutoff on on enrichment tests to report. Default is 0.1.
#' @param qvalueCutoff Q-value cutoff on on enrichment tests to report. Default is 0.1.
#' @param minGSSize Minimal size of genes annotated by Ontology term for testing. Default is 20.
#' @param maxGSSize Maximal size of genes annotated by Ontology term for testing. Default is 100.
#' @param simp If simplifying the GO terms. Default is FALSE.
#' @param return_fig If returning the visualization plots. Default is FALSE.
#' @param cat_num How many GO terms will be shown for basis comparison plot.
#' @param word_num How many GO terms will be used for word cloud plot.
#' 
#' @details This function annotates the weight matrix \eqn{W}. Each basis of \eqn{W} represents a functional gene cluster. Users can use this annotation to decide which bases should be kept for further analysis.
#' See details about the GO analysis in \code{clusterProfiler}.
#' 
#' @importFrom magrittr %>%
#' @importFrom methods is
#' @importFrom tibble as_tibble
#' @import ggplot2
#'  
#' 
#' @export 
#' @author Dongyuan Song


basisAnnotate <- function(W,
                          DE_prop = 0.1,
                          dim_use = NULL,
                          id_type = "SYMBOL",
                          ont = "BP",
                          OrgDb = "org.Hs.eg.db",
                          pvalueCutoff = 0.1,
                          qvalueCutoff = 0.1,
                          minGSSize    = 20,
                          maxGSSize    = 100,
                          simp = FALSE,
                          return_fig = FALSE,
                          cat_num = 10, 
                          word_num = 50
) {
  
  Cluster <- ID <- Description <- word <- n <- NULL
  
  if(is.null(rownames(W))) {
    stop("Input matrix must contain row names!")
  }
  
  gene_num <- dim(W)[1]
  basis_num <- dim(W)[2]
  
  colnames(W) <- paste0("Basis", 1:basis_num)
  
  if(is.null(dim_use)) dim_use <- 1:(dim(W)[2])
  
  
  gene_list <- getInfoGene(W = W,
                           M = round(DE_prop*gene_num),
                           by_basis = TRUE,
                           dim_use = dim_use)
  
  name_vec <- colnames(W)[dim_use]
  names(gene_list) <- name_vec
  
  
  pnmf_comp <- clusterProfiler::compareCluster(geneClusters = gene_list, 
                                               keyType = id_type, 
                                               fun = "enrichGO", 
                                               OrgDb = OrgDb, 
                                               ont = ont, 
                                               pvalueCutoff = pvalueCutoff, 
                                               qvalueCutoff = qvalueCutoff, 
                                               universe = rownames(W), 
                                               pAdjustMethod = "BH",
                                               minGSSize = minGSSize,
                                               maxGSSize = maxGSSize,
                                               readable = TRUE)
  
  if(simp) pnmf_comp <- clusterProfiler::simplify(pnmf_comp, cutoff=0.7, by="p.adjust", select_fun=min)
  
  
  if(return_fig) {
    p_comp <- enrichplot::dotplot(pnmf_comp, showCategory = cat_num, font.size = 8)
    
    p_wordcloud <- lapply(name_vec, function(x) {
      dat <- pnmf_comp@compareClusterResult %>% dplyr::filter(Cluster == x)
      
      num_use <- min(word_num, dim(dat)[1])
      
      dat_p <- dat %>% dplyr::select(ID, Description)
      dat_p <- dat_p[1:num_use,]
      
      p <-  dat_p %>% tidytext::unnest_tokens(word, Description, token = "ngrams", n = 2) %>%
        #anti_join(stop_words) %>%
        dplyr::count(word, sort = TRUE) %>%
        dplyr::filter(n >= 2) %>%
        dplyr::filter(!sapply(word, function(x) {
          (grepl("^of", x) | grepl("of$", x)  | grepl("^and", x) | grepl("and$", x) | grepl("^to", x) | grepl("to$", x) | grepl("^via", x) | grepl("via$", x)   )
        })) %>%
        with(ggwordcloud::ggwordcloud(word, n, random.order = FALSE, scale = c(4, 1)))
      p
    })
    
    names(p_wordcloud) <- name_vec
    
    return(list(p_comp = p_comp, p_wordcloud = p_wordcloud, go_res = pnmf_comp))
  }
  else return(pnmf_comp)
}



