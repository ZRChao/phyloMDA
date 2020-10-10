# phyloMDA

An R implementation of statistical methods from “Statistical Methods for Analyzing Tree Structured Microbiome Data” by Tao Wang and Hongyu Zhao.

To install this package on your machine, first install the devtools package
```r
install.packages("devtools")  
```
and then use
```r
library(devtools)
install_github("liudoubletian/phyloMDA")  
library(phyloMDA)  
```
Here, we show a brief example.

Load Example Data


```r
#phyloseq: A tool to import, store, analyze, and display phylogenetic sequencing data
library(phyloseq); packageVersion("phyloseq")

#phyloMDA: An R package for phylogeny-aware microbiome data analysis
library(phyloMDA); packageVersion("phyloMDA")
(phyloseq.obj <- combo.phyloseq.obj)

#Plot the phylogenetic tree
plot_tree(phyloseq.obj, "treeonly", nodeplotblank, label.tips="taxa_names")
tree <- phy_tree(phyloseq.obj)

#Heatmap of microbial counts
plot_heatmap(phyloseq.obj, taxa.order=taxa_names(phyloseq.obj))
otu_tab <- t(phyloseq.obj@otu_table@.Data)

#Metadata
metadata <- sample_data(phyloseq.obj)
```

Multinomial-logit regression

```r
# MGLM: A package for multivariate response GLMs
library(MGLM); packageVersion("MGLM")

fit_mn <- MGLMfit(data=otu_tab, dist="MN")
fit_mn@logL # MN loglikelihood

sodium <- metadata$sodium
reg_mn <- MGLMreg(otu_tab~1+sodium, dist="MN")
reg_mn@logL # simple MN regression loglikelihood
```


Dirichlet-multinomial regression
```r
fit_dm <- MGLMfit(data=otu_tab, dist="DM")
fit_dm@logL # DM loglikelihood

reg_dm <- MGLMreg(otu_tab~1+sodium, dist="DM")
reg_dm@logL # simple DM regression loglikelihood

#MGLMsparsereg and MGLMtune fit sparse regression
Nutrs <- metadata[, 18:37] # first 20 nutrients
Nutrs <- as.matrix(data.frame(Nutrs))
idx <- c(F, rep(T, dim(Nutrs)[2]))

sreg_dm <- MGLMsparsereg(otu_tab~1+Nutrs, dist="DM", lambda=Inf, penalty="sweep", penidx=idx)
sreg_dm@logL

sreg_dm_tune <- MGLMtune(otu_tab~1+Nutrs, dist="DM", penalty="sweep", penidx=idx)
sreg_dm_tune@select@logL

```

Dirichlet-tree multinomial regression
```r
fit_dtm <- MGLMdtmFit(otu_tab, tree)
Extract_logL(fit_dtm) # DTM loglikelihood

reg_dtm <- MGLMdtmReg(otu_tab, sodium, tree)
Extract_logL(reg_dtm) # DTM regression loglikelihood

sreg_dtm <- MGLMdtmSparseReg(otu_tab, Nutrs, tree, lambda=Inf, penalty="sweep")
Extract_logL(sreg_dtm)

sreg_dtm_tune <- MGLMdtmTune(otu_tab, Nutrs, tree, penalty="sweep")
Extract_logL(sreg_dtm_tune)

```

Empirical Bayes normalization
```r
eBay.comps <- eBay_comps(otu_tab, prior="Dirichlet")

eBay.tree.comps <- eBay_comps(otu_tab, prior="Dirichlet-tree", tree=tree)

```

Log-ratio transformations
```r
eBay.comps.alr <- alr_trans(eBay.comps)
eBay.comps.clr <- clr_trans(eBay.comps)
```


Constrained lasso and log-ratio lasso
```r
library(logratiolasso)
packageVersion("logratiolasso")

y <- metadata$bmi
x <- log(eBay.comps) # log of estimated compositions
centered_y <- y - mean(y)
centered_x <- scale(x, center=T, scale=F)

#constrained lasso
classo <- glmnet.constr(centered_x, centered_y)
set.seed(10)
cv_constr_lasso <- cv.glmnet.constr(classo, centered_x, centered_y)

# two-stage log-ratio lasso
set.seed(10)
cv_ts_lasso <- cv_two_stage(centered_x, centered_y, k_max = 7)

```

TASSO
```r

fit_tasso <- TASSO(y, eBay.comps, tree)
fit_tasso

fit_classo <- TASSO(y, eBay.comps, tree = NULL) # constrained lasso
fit_classo

```

Tree-guided fused lasso
```
fit_tflasso1 <- TreeFusedlasso(y, eBay.comps, tree)
fit_tflasso1

fit_tflasso2 <- TreeFusedlasso(y, eBay.comps, tree, type = 2)
fit_tflasso2

```

