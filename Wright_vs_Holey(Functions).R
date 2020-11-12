library(mvtnorm);library(coda);library(loo);library(MASS);library(gdata);library(Rcpp);library(gtools)

#### cppFunction ####
#This function matches phenotype sequence (e.g. 0011010001 to corresponding fitness)
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

#LKJ onion method for generating random correlation matrices,
#third-hand from McElreath's "rethinking" package
rlkjcorr <- function (n, K, eta = 1) {

  stopifnot(is.numeric(K), K >= 2, K == as.integer(K))
  stopifnot(eta > 0)
  #if (K == 1) return(matrix(1, 1, 1))

  f <- function() {
    alpha <- eta + (K - 2)/2
    r12 <- 2 * rbeta(1, alpha, alpha) - 1
    R <- matrix(0, K, K) # upper triangular Cholesky factor until return()
    R[1,1] <- 1
    R[1,2] <- r12
    R[2,2] <- sqrt(1 - r12^2)
    if(K > 2) for (m in 2:(K - 1)) {
      alpha <- alpha - 0.5
      y <- rbeta(1, m / 2, alpha)

      # Draw uniformally on a hypersphere
      z <- rnorm(m, 0, 1)
      z <- z / sqrt(crossprod(z)[1])

      R[1:m,m+1] <- sqrt(y) * z
      R[m+1,m+1] <- sqrt(1 - y)
    }
    return(crossprod(R))
  }
  R <- replicate( n , f() )
  if ( dim(R)[3]==1 ) {
    R <- R[,,1]
  } else {
    # need to move 3rd dimension to front, so conforms to array structure that Stan uses
    R <- aperm(R,c(3,1,2))
  }
  return(R)
}

#convert correlation matrix to covariance matrix
cor2cov<-function(vars,cormat){   
  sdMat<-diag(sqrt(vars))
  corMat<-cormat
  mat<-sdMat %*% corMat %*% t(sdMat)
  return(mat)
}

#Generate matrix of individual phenotypes--rows correspond to 
#individuals, 1:k columns are P, k+1:2k = G,  2k+1:3k = E.
pop.gen<-function(N=100,n.traits=10,h2=0.5,Emu=0,Evar=1){
  G<-matrix(runif(n.traits*N,-1,1),nrow=N)
  E<-matrix(rnorm(n.traits*N,Emu,Evar),nrow=N)
  init.pop<-cbind((h2*G)+((1-h2)*E),G,E)
}

#Calculate fitness based on a Gaussian landscape given  specified optima (optimum)
#and specified correlation matrix (omega)
fitness.AL.Lande1980<-function(z,optimum,omega){
  z<-matrix(z,ncol=1)
  if(dim(z)[1]!=dim(omega)[1]){
    stop("z and omega different numbers of traits")
  }
  if(dim(z)[1]!=dim(optimum)[1]){
    stop("z and theta are for different numbers of traits")
  }
  exp(-.5*(t(z-optimum)%*%solve(omega)%*%(z-optimum)))
}

#Calculate fitness [0, 1] based on a holey landscape 
#uses the C++ function above to match the phenotypes to fitness
fitness.HL<-function(mat,pattern,HLS.fit){
  out=SeqInVecOpt(mat,pattern)
  which.one=(matrix(out,ncol=n.traits,byrow=T)[,1]-1)/n.traits
  this.one=na.omit(which.one[is.wholenumber(which.one)=="TRUE"]+1)
  HLS.fit[this.one]
}

#sorts individuals according to the first column, takes the proportion ("prop")
#with highest fitness
trunc.select<-function(indivs,prop=0.2,decreasing=TRUE){
  fit.indivs<-indivs[order(indivs[,1],decreasing=decreasing),]
  n.selected<-round(prop*length(indivs[,1]),0)
  Selected<-fit.indivs[c(1:n.selected),]
  return(Selected)
}

make.binary<-function(x){ifelse(x<0,0,1)}

is.wholenumber<-function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
