---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# MultiPhenoAssoc

<!-- badges: start -->
<!-- badges: end -->

The goal of MultiPhenoAssoc is to ...

## Installation

You can install the development version of MultiPhenoAssoc like so:
```{r}
library(devtools)
install_github("HungChingChang/MultiPhenoAssoc")
```

## Example

Multivariate phenotype-gene association analysis
```{r}
library(MultiPhenoAssoc)
data(exampleData)
Expr <- exampleData$Gene.expr
Pheno <- exampleData$Pheno
Confounder <- exampleData$Confounder
Pheno.type <- c("continuous", "count", "binary", "survival")
AWFisher.MultiPheno(expr = Expr,
                    pheno = Pheno,
                    confounder = Confounder,
                    Pheno.type = Pheno.type,
                    method = "AFp",
                    num.bootstrap = 10,
                    ncores = 5)
```

prepare distance matrix for tight clustering
```{r}
load("AFp_result.Rdata")
load("Distance.matrix.AFp.Rdata")
AFp.BH <- p.adjust(AFp.result$AFp.pvalue, "BH")
sig.gene.index <- which(AFp.BH < 0.01)
Distance.matrix <- 1 - Distance.matrix.AFp[sig.gene.index, sig.gene.index]
colnames(Distance.matrix) <- rownames(Distance.matrix) <- rep("",dim(Distance.matrix)[1])
```

tight clustering
```{r}
res <- tight.clust(Distance.matrix, 3, 10, random.seed=12345)
save(res, file="tight_clustering.Rdata")
set1.index <- which(res$cluster==1)
set2.index <- which(res$cluster==2)
set3.index <- which(res$cluster==3)
tight.cluster.gene <- c(set1.index, set2.index, set3.index)
colors <- rep(as.character(1:3),
              c(length(set1.index), length(set2.index), length(set3.index)))
tight.cluster.distance.matrix <- Distance.matrix[tight.cluster.gene, tight.cluster.gene]

library(gplots)
heatmap.2(1-tight.cluster.distance.matrix,
          trace = "none", col="greenred", Rowv = NULL, dendrogram = "none",
          Colv=NULL, RowSideColors = colors, ColSideColors = colors,
          density.info = "none", key.title = NA, key.xlab = "Co-membership", key.ylab = NA,
          keysize=1, key.par = list(cex=1))
legend("top", legend = paste0("C",1:3),
       col = as.character(1:3),
       cex=1, lty=1, lwd=3, horiz = T)
```

