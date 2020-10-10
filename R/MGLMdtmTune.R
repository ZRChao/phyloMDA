#' Finds the tuning parameter value that yields the smallest BIC for the Dirichlet-tree multinomial (DTM) sparse regression.
#'
#' @param otu.tab a data frame or matrix containing the count data. Rows of the matrix represent observations and columns are the taxa.
#' @param X a data frame or matrix containing the covariates. Rows of the matrix represent observations.
#' @param tree the phylogenetic tree.
#' @param penalty penalty type for the regularization term. Can be chosen from 'sweep', 'group', or 'nuclear'. See Details for the description of each penalty type of MGLM.
#' @param ngridpt an optional numeric variable specifying the number of grid points to tune. see MGLMtune.
#' @return Returns a list of MGLM objects ('MGLMtune') corresponding to internal nodes of the tree.
#' @references Wang, T., Zhao, H.: A Dirichlet-tree multinomial regression model for associating dietary nutrients with gut microorganisms. Biometrics 73(3), 792–801 (2017).
#' @examples
#' library(phyloseq)
#' library(magrittr)
#' otu.tab <- t(otu_table(combo.phyloseq.obj))
#' metadata <- sample_data(combo.phyloseq.obj)
#' tree <- phy_tree(combo.phyloseq.obj)
#' X <- as.matrix(metadata)[,18:37] %>% apply(., 2, as.numeric)
#' # fit <- MGLMdtmTune(otu.tab, X, tree, penalty = "sweep")
#' # Extract_logL(fit)
#' @export
#' @import phyloseq MGLM
#' @import foreach



MGLMdtmTune <- function(otu.tab, X, tree, penalty = 'sweep', ngridpt = 20) {
  ytree <- YtreeFun(otu.tab, tree)

  MGLMdtmTune_sub <- function(i){
    y <- ytree[[i]]
    t1 <- try(f1 <- MGLMtune(formula = y~1+X, dist = 'DM', penalty = penalty, ngridpt = ngridpt, penidx=c(F, rep(T, ncol(X)))), silent = T)
    if(class(t1) == 'try-error') {
      t2 <- try(f1 <- MGLMtune(formula = y~1+X, dist = 'MN', penalty = penalty, ngridpt = ngridpt, penidx=c(F, rep(T, ncol(X)))), silent = T)
      if(class(t2) == 'try-error')
        t2 <- try(f1 <- MGLMtune(formula = y~1+X, dist = 'DM', penalty = penalty, ngridpt = 1, penidx=c(F, rep(T, ncol(X)))), silent = T)
    }
    f1
  }

  fit <-  foreach(i=1:tree$Nnode)%dopar% MGLMdtmTune_sub(i)

  return(fit)
}

