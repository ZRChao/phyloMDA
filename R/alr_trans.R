#' The additive log-ratio transformation
#'
#' @param eBay.comps the estimated microbial compositions by empirical Bayes.
#' @return Returns the alr transfromation of microbial compositions.
#' @examples
#' library(phyloseq)
#' otu_tab <- t(combo.phyloseq.obj@otu_table@.Data)
#' tree <- phy_tree(combo.phyloseq.obj)
#' eBay.comps <- eBay_comps(otu_tab,prior = "Dirichlet", tree=NULL)
#' eBay.comps.alr <- alr_trans(eBay.comps)
#' @export

### additive log-ratio (alr) ###
alr_trans <- function(eBay.comps){
  comps_alr <- apply(eBay.comps, 1, function(x){log(x[-length(x)]) - log(x[length(x)])}) %>% t
  return(comps_alr)
  }


