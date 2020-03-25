# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

hello <- function() {
  print("Hello, world! fan!")
}
#problem1
myGD=read.table(file="http://zzlab.net/GAPIT/data/mdp_numeric.txt",head=T)
myGM=read.table(file="http://zzlab.net/GAPIT/data/mdp_SNP_information.txt",head=T)
myY=read.table(file="http://zzlab.net/GAPIT/data/CROP545_Phenotype.txt",head=T)
myC=read.table(file="http://zzlab.net/GAPIT/data/CROP545_Covariates.txt",head=T)
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
#problem3

#problem4
setwd("/Users/fanzhang/Desktop/teach_PPT/homework2")
genotype<-data.frame(read.table("GAPIT.Genotype.Numerical.txt", header = T, stringsAsFactors = F, sep = "\t"))
phenotype<-read.table("220pheno.txt",header = T)
myPC<-prcomp(genotype[,-1])
myC<-myPC$x[,100:220] #set PC100 to PC220 as cofactor
myGWAS<-GWAS(y=phenotype,X=genotype,C=myC)
#problem5
myGD=read.table(file="http://zzlab.net/GAPIT/data/mdp_numeric.txt",head=T)
myGM=read.table(file="http://zzlab.net/GAPIT/data/mdp_SNP_information.txt",head=T)
source("http://zzlab.net/StaGen/2020/R/G2P.R")
source("http://zzlab.net/StaGen/2020/R/GWASbyCor.R")
X=myGD[,-1]
index1to5=myGM[,2]<6
X1to5 = X[,index1to5]
set.seed(9916)
mySim=G2P(X= X1to5,h2=.75,alpha=1,NQTN=10,distribution="norm")
p_cor= GWASbyCor(X=X,y=mySim$y)
color.vector <- rep(c("deepskyblue","orange","forestgreen","indianred3"),10)
m=nrow(myGM)
par(mfrow=c(2,1), oma =c(0,1,0,1), mar = c(4,2,3,2)+0.1)
plot(t(-log10(p_cor))~seq(1:m),col=color.vector[myGM[,2]],main="GWASbyCor")
abline(v=mySim$QTN.position, lty = 2, lwd=2, col = "black")
##myGWAS
myC=read.table(file="http://zzlab.net/GAPIT/data/CROP545_Covariates.txt",head=T)
y=cbind(myGD[,1],mySim$y)
X=myGD
myPC<-prcomp(myGD[,-1])
myC<-myPC$x[,100:281] #set PC100 to PC220 as cofactor
myP<-GWAS(y=y,X=X,C=myC)
plot(t(-log10(myP))~seq(1:m),col=color.vector[myGM[,2]])
abline(v=mySim$QTN.position, lty = 2, lwd=2, col = "black")

