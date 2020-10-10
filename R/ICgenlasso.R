#' Calculates the Bayesian information criterion or Akaike's information criterion for a genlasso object.
#'
#' This function extracts BIC or AIC from a genlasso object.
#'
#' @param fit a genlasso object.
#' @param IC IC = BIC for Bayesian information criterion and IC = AIC for Akaike's information criterion.
#' @return Returns the minimal AIC or BIC and the corresponding lambda.
#' @examples
#' library(phyloseq)
#' otu.tab <- t(otu_table(combo.phyloseq.obj))
#' X <- eBay_comps(otu.tab, prior = "Dirichlet")
#' tree <- phy_tree(combo.phyloseq.obj)
#' metadata <- sample_data(combo.phyloseq.obj)
#' y <- metadata$bmi
#' fit <- TASSO(y, X, tree)
#' IC <- ICgenlasso(fit, IC="AIC")
#' IC <- ICgenlasso(fit, IC="BIC")
#' @export
#' @import phyloseq MGLM
#' @import genlasso


ICgenlasso <- function(fit, IC = 'BIC') {
  lambda <- fit$lambda
  res <- summary(fit)
  n <- length(fit$y)
  rss <- res[, 3]
  df <- res[, 1] - (n - dim(fit$X)[2])
  if (IC == 'AIC'){
    ic <- log(rss) + 2/n * df
    ic_res <- cbind(lambda, ic)
    colnames(ic_res) <- c('lambda', 'AIC')
  }
  if (IC == 'BIC'){
    ic <- log(rss) + log(n)/n*df
    ic_res <- cbind(lambda, ic)
    colnames(ic_res) <- c('lambda', 'BIC')
    }
  return(data.frame(ic_res))
}
