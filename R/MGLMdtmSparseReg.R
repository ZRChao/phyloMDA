#' Fit the Dirichlet-tree multinomial (DTM) sparse regression.
#'
#' @param otu.tab a data frame or matrix containing the count data. Rows of the matrix represent observations and columns are the taxa.
#' @param X a data frame or matrix containing the covariates. Rows of the matrix represent observations.
#' @param tree the phylogenetic tree.
#' @param penalty penalty type for the regularization term. Can be chosen from 'sweep', 'group', or 'nuclear'. See Details for the description of each penalty type of MGLM.
#' @param lambda penalty parameter.
#' @return Returns a list of MGLM objects ('MGLMsparsereg') corresponding to internal nodes of the tree.
#' @references Wang, T., Zhao, H.: A Dirichlet-tree multinomial regression model for associating dietary nutrients with gut microorganisms. Biometrics 73(3), 792–801 (2017).
#' @examples
#' library(phyloseq)
#' library(magrittr)
#' otu.tab <- t(otu_table(combo.phyloseq.obj))
#' metadata <- sample_data(combo.phyloseq.obj)
#' tree <- phy_tree(combo.phyloseq.obj)
#' X <- as.matrix(metadata)[,18:37] %>% apply(., 2, as.numeric)
#' fit <- MGLMdtmSparseReg(otu.tab, X, tree, penalty = 'sweep', lambda = Inf)
#' Extract_logL(fit)
#' @export
#' @import phyloseq MGLM


MGLMdtmSparseReg <- function(otu.tab, X, tree, penalty = 'sweep', lambda = Inf) {
  ytree <- YtreeFun(otu.tab, tree)
  fit <- list()
  for(i in 1: tree$Nnode) {
    y <- ytree[[i]]
    t1 <- try(f1 <- MGLMsparsereg(formula = y~1+X, dist = 'DM', lambda = lambda,penalty = penalty, penidx=c(F, rep(T, ncol(X)))), silent = T)
    if(class(t1) == 'try-error') {
      t2 <- try(f1 <- MGLMsparsereg(formula = y~1+X, dist = 'MN', lambda = lambda, penalty = penalty, penidx=c(F, rep(T, ncol(X)))), silent = T)
      if(class(t2) == 'try-error') f1 <-  MGLMreg(y~1, dist = 'DM')

    }
    fit[[i]] <- f1
  }
  fit
}
