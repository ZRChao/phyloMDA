#' Fit the Tree-guided Automatic Sub-composition Selection Operator (TASSO) model.
#'
#' @param y the univariate response.
#' @param X a data frame or matrix containing the compositional data. Rows of the matrix represent observations and columns are the taxa.
#' @param tree when tree=NULL, it reduces to constrained lasso.
#' @return Returns a genlasso object.
#' @references Wang, T., Zhao, H.: Structured subcomposition selection in regression and its application to microbiome data analysis. The Annals of Applied Statistics 11(2), 771–791 (2017).
#' @examples
#' library(phyloseq)
#' otu.tab <- t(otu_table(combo.phyloseq.obj))
#' X <- eBay_comps(otu.tab, prior = "Dirichlet")
#' tree <- phy_tree(combo.phyloseq.obj)
#' metadata <- sample_data(combo.phyloseq.obj)
#' y <- metadata$bmi
#' fit <- TASSO(y, X, tree)
#' IC <- ICgenlasso(fit, IC='AIC')
#' IC <- ICgenlasso(fit, IC='BIC')
#' @export
#' @import phyloseq MGLM
#' @import caper genlasso


TASSO <- function(y, X, tree=NULL) {
  # type, default constrained lasso (genlasso)
  # construct the penalty matrix
  K <- ncol(X)
  D1 <- diag(K) - 1 / K
  rownames(D1) <- 1 : K

  if(!is.null(tree)) {
    tree.edge <- tree$edge
    N.node <- tree$Nnode

    node.tips <- clade.members.list(tree, tips = F, tip.labels = F)
    M <- matrix(0, N.node, K)
    rownames(M) <- names(node.tips)
    root <- NULL
    for (j in 1 : N.node) {
      node.tips.j <- node.tips[[j]]
      if (length(node.tips.j) == K) root <- j
      M[j, ] <- -length(node.tips.j) / K
      M[j, node.tips.j] <- M[j, node.tips.j] + 1
    }
    D1 <- rbind(D1, M[-root, ])
  }

  alpha <- 0.0001
  centered_Y <- c(y-mean(y), rep(0, K))
  centered_X <- rbind(scale(clr_trans(X), T, F), alpha*diag(K)) #
  fit <- genlasso(centered_Y, centered_X, D1)
  fit
}
