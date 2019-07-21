library(CorShrink)
library(glasso)
library(corpcor)
library(Matrix)
library(psych)


hub_sim = function(N, P, block){
  mat <- 0.3*diag(1,block) + 0.7*rep(1,block) %*% t(rep(1, block))
  Sigma <-   bdiag(mat, mat, mat, mat, mat, mat, mat, mat, mat, mat)
  corSigma <- cov2cor(Sigma)
  data <- MASS::mvrnorm(N,rep(0,P),corSigma)
  return(data)
}

banded_prec_sim = function(N, P){
  diags <- list()
  diags[[1]] <- rep(1, P)
  diags[[2]] <- rep(-0.5, P)
  Kinv <- bandSparse(P, k = -(0:1), diag = diags, symm = TRUE)
  K <- solve(Kinv)
  corSigma <- cov2cor(K)
  data <- MASS::mvrnorm(N,rep(0,P),corSigma)
  return(data)
}

DM_toeplitz = function(n,P){
  library("MASS")
  index1=sort(sample(seq(1:n),(n/2)))
  index2=seq(1:n)[-index1]
  Sigmatp=function(P){
    a=array(0,dim=c(P,P))
    for(i in 1:P){
      for(j in 1:P){
        a[i,j]=max(1-0.1*(abs(i-j)),0)
      }
    }
    return(a)
  }
  Sigma = Sigmatp(P)
  data = mvrnorm(n,rep(0,P),Sigma)
  Xtest = data[index2,]
  Xtrain = data[index1,]
  Omega = solve(Sigma)
  return(list(Xtrain = Xtrain, Xtest = Xtest, Sigma = Sigma))
}


toeplitz_sim = function(N, P){
  ll <- DM_toeplitz(n=N, P=P)
  data <- rbind(ll$Xtrain, ll$Xtest)
  Sigma <- ll$Sigma
  corSigma <- cov2cor(Sigma)
  return(data)
}


nsamples = 100
nfeatures = 100
res1 = hub_sim(nsamples, nfeatures, 10)
res2 = banded_prec_sim(nsamples, nfeatures)
res3 = toeplitz_sim(nsamples, nfeatures)



