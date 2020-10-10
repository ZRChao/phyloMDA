#' Fit the Dirichlet-tree multinomial (DTM) regression.
#'
#' @param otu.tab a data frame or matrix containing the count data. Rows of the matrix represent observations and columns are the taxa.
#' @param X a data frame or matrix containing the covariates. Rows of the matrix represent observations.
#' @param tree the phylogenetic tree.
#' @return Returns a list of MGLM objects ('MGLMreg') corresponding to internal nodes of the tree.
#' @references Wang, T., Zhao, H.: A Dirichlet-tree multinomial regression model for associating dietary nutrients with gut microorganisms. Biometrics 73(3), 792–801 (2017).
#' @examples
#' library(phyloseq)
#' otu.tab <- t(otu_table(combo.phyloseq.obj))
#' tree <- phy_tree(combo.phyloseq.obj)
#' sodium <- sample_data(combo.phyloseq.obj)$sodium
#' fit <- MGLMdtmReg(otu.tab, sodium, tree)
#' Extract_logL(fit)
#' @export
#' @import phyloseq MGLM



MGLMdtmReg <- function(otu.tab, X, tree) {
  ytree <- YtreeFun(otu.tab, tree)
  fit <- lapply(ytree, function(y)  MGLMreg(y~1+X, dist='DM'))
  fit
}
