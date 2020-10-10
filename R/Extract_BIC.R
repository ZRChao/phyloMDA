#' Calculates the Bayesian information criterion for a list of MGLM objects.
#'
#' This function extracts BIC from a list of MGLM objects.
#'
#' @param x A list of MGLM objects. 'MGLMdtmFit', 'MGLMdtmReg', 'MGLMdtmSparseReg', or 'MGLMdtmTune'.
#' @return Returns the BIC of a MGLMdtm object (a list of MGLM objects).
#' @examples
#' library(phyloseq)
#' otu.tab <- t(otu_table(combo.phyloseq.obj))
#' tree <- phy_tree(combo.phyloseq.obj)
#' fit <- MGLMdtmFit(otu.tab, tree)
#' Extract_BIC(fit)
#' @export
#' @import phyloseq MGLM

Extract_BIC <- function(x) {

  ic <- sapply(x, function(x) {
    if(is.null( attributes(x)$select) ) c(x@BIC)
    else  c( x@select@BIC)
  })

  ic[is.infinite(ic)] <- NA
   bic <- sum(ic, na.rm = T)
   return(bic)
}
