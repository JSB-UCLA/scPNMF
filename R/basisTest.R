#' @title Data-driven basis selection
#' @description Perform data-driven basis selection to remove useless bases and keep informative bases. It contains two approaches:
#' (1) examine bases by correlations with cell library sizes; (2) examine bases by multimodality tests.
#' @param S The score matrix of scPNMF.
#' @param X The expression matrix of scPNMF.
#' @param return_fig If returning the visualization plots. Default is FALSE.
#' @param adj_method Multiple test correction method.
#' @param mc.cores The number of cores to use.
#' @param ncol Columns for facets in plots.
#' 
#' @importFrom magrittr %>%
#' @importFrom methods is
#' @import ggplot2
#' @importFrom parallel mclapply
#' @import tibble
#'  
#' 
#' @export
#' @author Dongyuan Song



basisTest <- function(S,
                      X,
                      return_fig = FALSE,
                      adj_method = "BH",
                      mc.cores = 10,
                      ncol = 4
) {
  
  basis <- Score <- NULL
  cell_lib_size <- colSums(X)
  
  cell_num <- dim(S)[1]
  basis_num <- dim(S)[2]
  colnames(S) <- paste0("Basis", 1:basis_num)
  
  dat_coef_list <- as.list(data.frame(S))
  
  cor_value <- simplify2array(lapply(dat_coef_list, function(x) {stats::cor(x, cell_lib_size)}))
  
  
  multimode_ACR_pvalue <- simplify2array(mclapply(dat_coef_list, function(x) {
    multimode::modetest(x, mod0 = 1, method ="ACR", B=999)$p.value}, mc.cores = mc.cores))
  
  adj_pvalue <- stats::p.adjust(multimode_ACR_pvalue, method = adj_method)
  
  if(return_fig) {
    dat <- S %>% as_tibble() %>% dplyr::mutate(cell_lib_size = cell_lib_size)
    dat_long <- dat %>% tidyr::pivot_longer(dplyr::starts_with("Basis"), names_to = "basis", values_to = "Score") %>% dplyr::mutate(basis = factor(basis, levels = colnames(S)))
    
    p_scatter <- dat_long %>% ggplot(aes(y = Score, x = cell_lib_size)) + geom_point(alpha = 0.5, size = 0.2) + facet_wrap(~basis,  scales = "free", ncol = ncol) + xlab("Cell lib size") + ggtitle("Correlations with cell library sizes") +theme_bw() + theme(legend.position = "bottom", aspect.ratio = 1, axis.text.x = element_text(angle = 45, hjust = 1)) + guides(color = guide_legend(override.aes = list(size=5)))
    
    p_hist <- dat_long %>% ggplot(aes(x = Score)) + geom_histogram(color = "black", fill = "white") + facet_wrap(~basis, scales = "free", ncol = ncol) + ggtitle("Multimodality") +theme_bw() + theme(legend.position = "bottom", aspect.ratio = 1)+ guides(color = guide_legend(override.aes = list(size=5))) + ylab("Count")
    
    p <- gridExtra::grid.arrange(p_scatter, p_hist, nrow = 2)
    
    return(list(cor_value = cor_value, adj_pvalue = adj_pvalue, fig = p))
  }
  else return(list(cor_value = cor_value, adj_pvalue = adj_pvalue, fig = NULL))
  
  
}











