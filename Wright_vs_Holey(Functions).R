library(rethinking);library(mvtnorm);library(coda);library(loo);library(MASS);library(gdata);library(Rcpp);library(gtools)

#### cppFunction ####
cppFunction('NumericVector SeqInVecOpt(NumericVector myVector, NumericVector mySequence) {

            int vecSize = myVector.size();
            int seqSize = mySequence.size();
            NumericVector comparison(seqSize);
            NumericVector res(vecSize);
            int foundCounter = 0;
            
            for (int i = 0; i < vecSize; i++ ) {
            
            if (myVector[i] == mySequence[0]) {
            for (int j = 0; j < seqSize; j++ ) {
            comparison[j] = mySequence[j] == myVector[i + j];
            }
            
            if (sum(comparison) == seqSize) {
            for (int j = 0; j < seqSize; j++ ) {
            res[foundCounter] = i + j + 1;
            foundCounter++;
            }
            }
            }
            }
            
            IntegerVector idx = seq(0, (foundCounter-1));
            return res[idx];
            }')

#### Functions ####
cor2cov<-function(vars,cormat){   
  sdMat<-diag(sqrt(vars))
  corMat<-cormat
  mat<-sdMat %*% corMat %*% t(sdMat)
  return(mat)
}

pop.gen<-function(N=100,n.traits=10,h2=0.5,Emu=0,Evar=1){
  G<-matrix(runif(n.traits*N,-1,1),nrow=N)
  E<-matrix(rnorm(n.traits*N,Emu,Evar),nrow=N)
  init.pop<-cbind((h2*G)+((1-h2)*E),G,E)
}

fitness.AL.Lande1980<-function(z,optimum,omega){
  z<-matrix(z,ncol=1)
  if(dim(z)[1]!=dim(omega)[1]){
    stop("z and omega are for different numbers of traits")
  }
  if(dim(z)[1]!=dim(optimum)[1]){
    stop("z and theta are for different numbers of traits")
  }
  exp(-.5*(t(z-optimum)%*%solve(omega)%*%(z-optimum)))
}

fitness.HL<-function(mat,pattern,HLS.fit){
  out=SeqInVecOpt(mat,pattern)
  which.one=(matrix(out,ncol=n.traits,byrow=T)[,1]-1)/n.traits
  this.one=na.omit(which.one[is.wholenumber(which.one)=="TRUE"]+1)
  HLS.fit[this.one]
}

trunc.select<-function(indivs,prop=0.2,decreasing=TRUE){
  fit.indivs<-indivs[order(indivs[,1],decreasing=decreasing),]
  n.selected<-round(prop*length(indivs[,1]),0)
  Selected<-fit.indivs[c(1:n.selected),]
  return(Selected)
}

make.binary<-function(x){ifelse(x<0,0,1)}

is.wholenumber<-function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol