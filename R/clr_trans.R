#' The centered log-ratio transformation
#'
#' @param eBay.comps the estimated microbial compositions by empirical Bayes.
#' @return Returns the clr transfromation of microbial compositions.
#' @examples
#' library(phyloseq)
#' otu_tab <- t(combo.phyloseq.obj@otu_table@.Data)
#' tree <- phy_tree(combo.phyloseq.obj)
#' eBay.comps <- eBay_comps(otu_tab,prior = "Dirichlet", tree=NULL)
#' eBay.comps.clr <- clr_trans(eBay.comps)
#' @export

### centered log-ratio (clr)###
clr_trans <- function(eBay.comps){
  comps_clr <- apply(eBay.comps, 1, function(x){log(x) - mean(log(x))}) %>% t
  rownames(comps_clr) <- rownames(eBay.comps)
  return(comps_clr)
}

