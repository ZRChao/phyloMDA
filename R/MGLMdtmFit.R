#' Fit the Dirichlet-tree multinomial (DTM) distribution.
#'
#' @param otu.tab a data frame or matrix containing the count data. Rows of the matrix represent observations and columns are the taxa.
#' @param tree the phylogenetic tree.
#' @return Returns a list of MGLM objects ('MGLMfit') corresponding to internal nodes of the tree.
#' @references Wang, T., Zhao, H.: A Dirichlet-tree multinomial regression model for associating dietary nutrients with gut microorganisms. Biometrics 73(3), 792–801 (2017).
#' @examples
#' library(phyloseq)
#' otu.tab <- t(otu_table(combo.phyloseq.obj))
#' tree <- phy_tree(combo.phyloseq.obj)
#' fit <- MGLMdtmFit(otu.tab, tree)
#' Extract_logL(fit)
#' @export
#' @import phyloseq MGLM


MGLMdtmFit <- function(otu.tab, tree) {
  ytree <- YtreeFun(otu.tab, tree)
  fit <- lapply(ytree, function(x) MGLMfit(x, dist = 'DM'))
  fit
}
