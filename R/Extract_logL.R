#' Calculates the loglikelihood for a list of MGLM objects.
#'
#' This function extracts loglikelihood from a list of MGLM objects.
#'
#' @param x A list of MGLM objects. 'MGLMdtmFit', 'MGLMdtmReg', 'MGLMdtmSparseReg', or 'MGLMdtmTune'.
#' @return Returns the loglikelihood of a MGLMdtm object (a list of MGLM objects).
#' @examples
#' library(phyloseq)
#' otu.tab <- t(otu_table(combo.phyloseq.obj))
#' tree <- phy_tree(combo.phyloseq.obj)
#' fit <- MGLMdtmFit(otu.tab, tree)
#' Extract_logL(fit)
#' @export
#' @import phyloseq MGLM

Extract_logL <- function(x) {

  ic <- sapply(x, function(x) {
    if(is.null( attributes(x)$select) ) c(x@logL)
    else  c( x@select@logL)
  })

  ic[is.infinite(ic)] <- NA
  logL <- sum(ic, na.rm = T)
  return(logL)
}
