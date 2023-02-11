# p-value of association between multiple phenotype and gene expression
Getpvalue <- function(X,Y,Z, data.type){
  pvalue<-matrix(0,nrow=ncol(X),ncol=ncol(Y))
  Coef<-matrix(0,nrow=ncol(X),ncol=ncol(Y))
  t.stat<-matrix(0,nrow=ncol(X),ncol=ncol(Y))
  # show ProgressBar
  pb <- txtProgressBar(min=0, max=ncol(X)*ncol(Y), style=3)
  pbcount <- 0
  setTxtProgressBar(pb, pbcount)
  for(i in 1:ncol(Y)){
    for(j in 1:ncol(X)){
      if(is.null(Z)){
        data.temp<-cbind(data.frame(Y=Y[,i],X=X[,j]))
      }else{
        data.temp<-data.frame(Y=Y[,i],X=X[,j],Z)
      }
      if(data.type[i] == "continuous"){
        model<-lm(Y~.,data=data.temp)
        Coef[j,i]<-model$coefficients[2]
        model1<-summary(model)
        t.stat[j,i]<-model1$coefficients[2,3]
        pvalue[j,i]<-model1$coefficients[2,4]
      }
      if(data.type[i] == "count"){
        model <- glm(Y ~ ., data=data.temp,  family = "poisson")
        Coef[j,i]<-model$coefficients[2]
        model1<-summary(model)
        t.stat[j,i]<-model1$coefficients[2,3]
        pvalue[j,i]<-model1$coefficients[2,4]
      }
      if(data.type[i] == "binary"){
        model <- glm(Y ~ ., data=data.temp, family = "binomial")
        Coef[j,i]<-model$coefficients[2]
        model1<-summary(model)
        t.stat[j,i]<-model1$coefficients[2,3]
        pvalue[j,i]<-model1$coefficients[2,4]
      }
      if(data.type[i] == "survival"){
        model <- coxph(Y ~ ., data=data.temp)
        Coef[j,i]<-model$coefficients[1]
        model1<-summary(model)
        t.stat[j,i]<-model1$coefficients[1,4]
        pvalue[j,i]<-model1$coefficients[1,5]
      }
      pbcount <- pbcount+1
      setTxtProgressBar(pb, pbcount)
    }
  }
  res<-list(pvalue=pvalue,t.stat=t.stat,coef=Coef)
  return(res)
}

# all possible combinations of multiple phenotype
GetSortMatrix_exact <- function(n){
  W<-do.call(expand.grid, rep(list(c(0, 1)), n))[-1,]
  return(W)
}

# permute gene expression data
Permutation <- function(X){
  n<-dim(X)[1]
  index.permute<-sample(1:n,replace = F)
  X<-X[index.permute,]

  return(X)
}

# permutation procedure to control type I error
GetPermutationPvalue <- function(Data.Gene,Y.residual,num.perm,ncores=1){
  G=ncol(Data.Gene)
  K<-ncol(Y.residual)
  Pvalue.permutation<-NULL
  Tstat.permutation<-NULL
  pvalue<-matrix(0,nrow=ncol(Data.Gene),ncol=ncol(Y.residual))
  t.stat<-matrix(0,nrow=ncol(Data.Gene),ncol=ncol(Y.residual))
  wcs<-pbmclapply(1:num.perm,function(t){
    set.seed(t)
    print(t)
    Data.Gene.perm<-Permutation(Data.Gene)
    for(i in 1:ncol(Y.residual)){
      res <- regression(x=Data.Gene.perm,y=Y.residual[,i])
      pvalue[,i] <- res[,2]
      t.stat[,i] <- res[,1]
    }
    return(list(pvalue = pvalue, t.stat = t.stat))
  },mc.cores = ncores)

  for(i in 1:num.perm){
    Pvalue.permutation[[i]]<-wcs[[i]]$pvalue
    Tstat.permutation[[i]]<-wcs[[i]]$t.stat
  }
  res<-list(Pvalue=Pvalue.permutation,Tstat=Tstat.permutation)
  return(res)
}

# AFp method to aggregate p-values from each dependent phenotype
AFp <- function(pvalue,pvalue.perm){
  W<-GetSortMatrix_exact(ncol(pvalue))
  W<-t(W)
  TEMP1<-(-2*log(pvalue.perm))%*%W
  num.nonzero<-apply(W,2,function(x){sum(x>0)})
  TStat<-(-2*log(pvalue)%*%W)

  # show ProgressBar
  pb <- txtProgressBar(min=0, max=2*ncol(TStat)+nrow(pvalue), style=3)
  pbcount <- 0
  setTxtProgressBar(pb, pbcount)

  AFp<-matrix(NA,nrow=nrow(TStat),ncol=ncol(TStat))
  for(j in 1:ncol(AFp)){
    score.perm<-TEMP1[,j]
    score<-TStat[,j]
    index<-order(c(score,score.perm),decreasing = F)
    index1<-match(1:length(score),index)
    for(i in 1:nrow(AFp)){
      large.all<-(nrow(TEMP1)+nrow(TStat)-index1[i])
      large.sample<-sum(index1>index1[i])
      AFp[i,j]<-(large.all-large.sample)/nrow(TEMP1)
      if(AFp[i,j]==0){{AFp[i,j]<-1/nrow(TEMP1)}}
    }
    pbcount <- pbcount+1
    setTxtProgressBar(pb, pbcount)
  }
  AFp.Stat<-apply(AFp,1,min)#get the test statistic
  final.AFp<-matrix(0,nrow=nrow(AFp),ncol=ncol(pvalue))
  for(i in 1:nrow(final.AFp)){
    index<-which(AFp[i,]==AFp.Stat[i])
    num1<-num.nonzero[index]
    Stat<-TStat[i,index]
    index.num<-which.max(num1)
    if(length(index.num)>1){index.num<-index.num[which.max(Stat[index.num])]}
    index1<-index[index.num]
    final.AFp[i,]<-W[,index1]
    index<-which(final.AFp[i,]==1)
    final.AFp[i,index]<-ifelse(pvalue[i,index]>0.1,0,1)
    pbcount <- pbcount+1
    setTxtProgressBar(pb, pbcount)
  }

  #----------------For permuted data
  AFp.perm<-matrix(NA,nrow=nrow(TEMP1),ncol=ncol(TEMP1))
  for(j in 1:ncol(AFp.perm)){
    score.perm<-TEMP1[,j]
    index<-match(score.perm,sort(score.perm,decreasing = F))
    AFp.perm[,j]<-(nrow(TEMP1)-index+1)/nrow(TEMP1)
    pbcount <- pbcount+1
    setTxtProgressBar(pb, pbcount)
  }
  AFp.Stat.perm<-apply(AFp.perm,1,min)
  pvalue.res<-unlist(lapply(AFp.Stat,function(x){sum(AFp.Stat.perm<x)/nrow(pvalue.perm)}))
  pvalue.res<-unlist(lapply(pvalue.res,function(x){
    if(x==0){
      x<-1/nrow(pvalue.perm)
    }
    return(x)
    }))
  res<-list(weight=final.AFp, pvalue=pvalue.res)
  return(res)
}

# AFz method to aggregate p-values from each dependent phenotype
AFz <- function(pvalue,pvalue.perm){
  W<-GetSortMatrix_exact(ncol(pvalue))
  W<-t(W)
  TEMP1<-(-2*log(pvalue.perm))%*%W
  TStat<-(-2*log(pvalue)%*%W)

  # show ProgressBar
  pb <- txtProgressBar(min=0, max=2*ncol(TStat)+nrow(pvalue), style=3)
  pbcount <- 0
  setTxtProgressBar(pb, pbcount)

  AFz<-matrix(NA,nrow=nrow(TStat),ncol=ncol(TStat))
  for(j in 1:ncol(AFz)){
    AFz[,j]<-abs((TStat[,j]-mean(TEMP1[,j]))/sd(TEMP1[,j]))
    pbcount <- pbcount+1
    setTxtProgressBar(pb, pbcount)
  }
  AFz.Stat<-apply(AFz,1,max)#get the test statistic
  final.AFz<-matrix(NA,nrow=nrow(TStat),ncol=ncol(pvalue))
  for(i in 1:nrow(AFz)){
    final.AFz[i,]<-W[,which.max(AFz[i,])]
    index<-which(final.AFz[i,]==1)
    final.AFz[i,index]<-ifelse(pvalue[i,index]>0.1,0,1)
    pbcount <- pbcount+1
    setTxtProgressBar(pb, pbcount)
  }

  #----------------For permuted data
  AFz.perm<-matrix(NA,nrow=nrow(TEMP1),ncol=ncol(TEMP1))
  for(j in 1:ncol(AFz.perm)){
    AFz.perm[,j]<-abs((TEMP1[,j]-mean(TEMP1[,j]))/sd(TEMP1[,j]))
    pbcount <- pbcount+1
    setTxtProgressBar(pb, pbcount)
  }

  AFz.Stat.perm<-apply(AFz.perm,1,max)
  pvalue.res<-unlist(lapply(AFz.Stat,function(x){sum(AFz.Stat.perm>x)/nrow(pvalue.perm)}))
  pvalue.res<-unlist(lapply(pvalue.res,function(x){
    if(x==0){
      x<-1/nrow(pvalue.perm)
    }
    return(x)
  }))
  res<-list(weight=final.AFz, pvalue=pvalue.res)
  return(res)
}
