#problem1
AStest=function(y,X,C){ #y represent phenotype, X represent genotype, C represent covariates.
  y=y[,-1] #remove taxa name
  C=C[,-1] #remove taxa name of covariates
  G=X[,-1] #these code come from lecture of GLM
  n=nrow(G)
  m=ncol(G)
  P=matrix(NA,1,m)
  for (i in 1:m){
    x=G[,i]
    if(max(x)==min(x)){
      p=1}else{
        X1=cbind(1, C[,1:2],x)
        LHS=t(X1)%*%as.matrix(X1)
        C1=solve(LHS)
        RHS=t(X1)%*%y
        b=C1%*%RHS
        yb=as.matrix(X1)%*%b
        e=y-yb
        n=length(y)
        ve=sum(e^2)/(n-1)
        vt=C1*ve
        t=b/sqrt(diag(vt))
        p=2*(1-pt(abs(t),n-2))
      } #end of testing variation
    P[i]=p[length(p)]
  } #end of looping for markers
  return(P)
}
myP<-AStest(y=myY,X=myGD,C=myC) #test

#problem2
GWAS=function(y,X,C){ #y represent phenotype, X represent genotype, C represent covariates.
  y=y[,-1] #remove taxa name
  C=C[,-1] #remove taxa name of covariates
  G=X[,-1]
  ###remove PC that are in linear dependent to the covariates
  PCA=prcomp(G)
  r<-cor(PCA$x,C)
  index1=r==1
  r[index1]=NA
  r_remain=na.omit(r) #I don't know how to remove rows with value of 1, so I change 1 to NA, than remove them.
  tPC=t(PCA$x)
  keep_tPC<-tPC[rownames(tPC) %in% rownames(r_remain),]
  keep_PC<-t(keep_tPC)
  #GWAS related content
  n=nrow(G)
  m=ncol(G)
  P=matrix(NA,1,m)
  for (i in 1:m){
    x=G[,i]
    if(max(x)==min(x)){
      p=1}else{
        X1=cbind(1,keep_PC[,1:3],x) #only set PC1-3 as cofactor,right?
        LHS=t(X1)%*%as.matrix(X1)
        C1=solve(LHS)
        RHS=as.matrix(t(X1))%*%y
        b=C1%*%RHS
        yb=as.matrix(X1)%*%b
        e=y-yb
        n=length(y)
        ve=sum(e^2)/(n-1)
        vt=C1*ve
        t=b/sqrt(diag(vt))
        p=2*(1-pt(abs(t),n-2))
      } #end of testing variation
    P[i]=p[length(p)]
  } #end of looping for markers
  return(P)
}
myGWAS<-GWAS(y=myY,X=myGD,C=myC) #test
