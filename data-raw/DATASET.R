# make toy example data from MetaBric
library(dplyr)
setwd("~/OneDrive/Pitts/side_project/GENES/MultiPhenoAssoc/")
load("../AWFISHER real data/MetaBric_gene.Rdata")
load("../AWFISHER real data/Complete_METABRIC_Clinical_Features_Data.rbin")
load("../AWFISHER real data/Complete_METABRIC_Clinical_Survival_Data_OS.rbin")
set.seed(12345)
# selection criteria: grade = 2 or 3
grade23.index <- which(Complete_METABRIC_Clinical_Features_Data$grade != "1")
Data.Gene <- t(MetaBric_gene[,grade23.index])
coef_variation <- apply(Data.Gene, 2, function(x) sd(x)/mean(x))
Data.Gene <- Data.Gene[,order(coef_variation, decreasing = T)[1:3000]] # 3k most variable genes
survival.outcome <- Complete_METABRIC_Clinical_Survival_Data_OS
Data.Clinical <- data.frame(Complete_METABRIC_Clinical_Features_Data[grade23.index,],
                            OS = survival.outcome[grade23.index])
Pheno.variable<-c("size","lymph_nodes_positive","grade","OS")
confounder.variable<-c("age_at_diagnosis")
Data.confounder <- Data.Clinical %>%  select(all_of(confounder.variable))
Pheno <- Data.Clinical %>% select(all_of(Pheno.variable))
Pheno$grade <- as.numeric(Pheno$grade == 3) # change grade variable to be binary
head(Pheno)

missing.confounder <- which(apply(Data.confounder,1,function(x){sum(is.na(x))})!=0)
missing.Pheno <- which(apply(Pheno,1,function(x){sum(is.na(x))})!=0)
missing.expression <- which(apply(Data.Gene,1,function(x){sum(is.na(x))})!=0)
index_delete <- Reduce(union, list(missing.confounder, missing.Pheno, missing.expression))
Data.Gene <- Data.Gene[-index_delete,]
Data.Clinical <- Data.Clinical[-index_delete,]
Data.confounder <- as.matrix(Data.confounder[-index_delete,])
Pheno <- Pheno[-index_delete,]
Gene.expr <- Data.Gene
Confounder <- cbind(Data.confounder,
                    rnorm(nrow(Data.confounder), 0, 1))
dim(Pheno);dim(Gene.expr);dim(Confounder)  ## 1708 completed samples
data.type <- c("continuous", "count", "binary", "survival")

# subsample 50 samples from grade2 and 50 samples from grade3
grade2 <- sample(rownames(Pheno[Pheno$grade == 0,]), 50)
grade3 <- sample(rownames(Pheno[Pheno$grade == 1,]), 50)
subsample <- which(rownames(Pheno)%in% c(grade2, grade3))
Pheno <- Pheno[subsample,]
#Pheno$size <- abs(rnorm(length(Pheno$size), 20, 10))
Gene.expr <- Gene.expr[subsample,]
Confounder <- Confounder[subsample,]
dim(Pheno);dim(Gene.expr);dim(Confounder) # 200 subsamples with 5000 genes
exampleData <- list(Pheno = Pheno,
                    Gene.expr = Gene.expr,
                    Confounder = Confounder)
usethis::use_data(exampleData, overwrite = TRUE)
