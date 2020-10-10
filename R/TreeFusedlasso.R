#' Fit the tree-guided fused lasso model.
#'
#' @param y the univariate response.
#' @param X a data frame or matrix containing the compositional data. Rows of the matrix represent observations and columns are the taxa.
#' @param tree the phylogenetic tree.
#' @param type two different penalty construct methods : 1 or 2.
#' @return Returns a genlasso object.
#' @references Wang, T., Zhao, H.: Constructing predictive microbial signatures at multiple taxonomic levels. Journal of the American Statistical Association 112(519), 1022–1031 (2017).
#' @examples
#' library(phyloseq)
#' otu.tab <- t(otu_table(combo.phyloseq.obj))
#' X <- eBay_comps(otu.tab, prior = "Dirichlet")
#' tree <- phy_tree(combo.phyloseq.obj)
#' metadata <- sample_data(combo.phyloseq.obj)
#' y <- metadata$bmi
#' fit <- TreeFusedlasso(y, X, tree)
#' IC <- ICgenlasso(fit, IC='AIC')
#' IC <- ICgenlasso(fit, IC='BIC')
#' @export
#' @import phyloseq MGLM
#' @import genlasso ape caper

TreeFusedlasso <- function(y, X, tree, type = 1) {
  K <- length(tree$tip.label)
  alpha <- 0.0001
  centered_Y <- c(y-mean(y), rep(0, K))
  centered_X <- rbind(scale(X, T, F), alpha*diag(K))

  node.tips <- clade.members.list(tree, tips = T, tip.labels = F)
  N.node <- tree$Nnode
  tree.edges <- tree$edge
  dist.mat <- dist.nodes(tree)

  if(type == 1) { ## tree-guided fused lasso penalty 1
    weights <- rep(NA, N.node)
    D.mat <- matrix(0, N.node, K)
    for (i in 1 : N.node) {
      i.nodes <- tree.edges[tree.edges[, 1] == (K + i), 2]
      weights[i] <- dist.mat[i.nodes[1], i.nodes[2]]
      tips.i.1 <- node.tips[[i.nodes[1]]]
      tips.i.2 <- node.tips[[i.nodes[2]]]
      D.mat[i, tips.i.1] <-  1 / length(tips.i.1)
      D.mat[i, tips.i.2] <- -1 / length(tips.i.2)
    }
    D.mat.w <- D.mat * weights^(-1)
  }
  if(type == 2) { ## tree-guided fused lasso penalty 2
    group.size <- NULL
    for (g in 1 : length(node.tips)) {
      node.tips.g <- node.tips[[g]]
      group.size <- c(group.size, length(node.tips.g))
    }
    level <- unique(sort(group.size))
    weights <- rep(NA, tree$Nnode + K)
    D.mat <- D.temp <- diag(1, tree$Nnode + K, K)
    for (i in level[-1]) {
      i.nodes <- which(group.size == i)
      for (j in i.nodes) {
        i.nodes.j <- tree.edges[tree.edges[, 1] == j, 2]
        weights[j] <- dist.mat[i.nodes.j[1], i.nodes.j[2]]
        D.mat[j, ] <- D.temp[i.nodes.j[1], ] - D.temp[i.nodes.j[2], ]
        D.temp[j, ] <-  abs(D.mat[j, ]) / 2
      }
    }
    weights <- weights[-(1 : K)]
    D.mat <- D.mat[-(1 : K), ]
    D.mat.w <- D.mat * weights^(-2)
  }

  fit <- genlasso(centered_Y, centered_X, D.mat.w)
  fit
}
