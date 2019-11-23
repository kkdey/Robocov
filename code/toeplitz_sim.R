options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
N <- as.numeric(toString(args[1]))
P <- as.numeric(toString(args[2]))
prop_missing <- as.numeric(toString(args[3]))



library(glasso)
library(corpcor)
library(Matrix)
library(psych)
library(CVXR)
library(Robocov)
library(CorShrink)

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
  ll = list("dat" = data, "cor" = corSigma)
  return(ll)
}


corSigma = as.matrix(toeplitz_sim(N, P)$cor)

nloglik = function(data, cormat){
  llik = 0
  for(m in 1:nrow(data)){
    idx = which(!is.na(data[m,]))
    if(length(idx) > 2){
      llik = llik + emdbook::dmvnorm(data[m, idx], rep(0, length(idx)), cormat[idx, idx], log = T)
    }
  }
  return(-llik)
}

angle_norm = function(S, Sigma){
  dist = 1 - (tr(as.matrix(cov2cor(S)%*%cov2cor(Sigma))))/(norm(cov2cor(S), type = "F")* norm(Sigma, type = "F"))
  return(dist)
}


NUM_SIM=10
ll <- vector(mode="list", length=NUM_SIM)

for(nsim in 1:NUM_SIM){
  data = toeplitz_sim(N, P)$dat

  #######################   Turn some of the entries to NA   ###################################

  data_missing = t(apply(data, 1, function(x){
    if(prop_missing > 0){
      rand = sample(1:length(x), floor(prop_missing*length(x)), replace = F)
      y = x
      y[rand] = NA
      return(y)
    }else{
      return(x)
    }}))



  standard_cor = cor(data_missing, use = "pairwise.complete.obs")
  cov_sample_ML <-  CorShrinkData(data_missing, sd_boot = FALSE,
                                  ash.control = list())
  corshrink_cor = cov2cor(cov_sample_ML$cor)
  robocov_box_cor = Robocov_box(data_with_missing = data_missing)

  alpha_vec = c(1e-02, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100)
  nlik_list = rep(0, length(alpha_vec))
  for(m in 1:length(alpha_vec)){
    temp_cor = Robocov_box_slack(data_missing, alpha=alpha_vec[m])
    nlik_list[m] = nloglik(data_missing, temp_cor)
    cat("Passed alpha value: ", alpha_vec[m], "\n")
  }
  final_alpha = alpha_vec[which.min(nlik_list)]
  robocov_box_slack_cor =  Robocov_box_slack(data_missing, alpha=final_alpha)


  alpha_vec = c(1e-02, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100)
  nlik_list = rep(0, length(alpha_vec))
  for(m in 1:length(alpha_vec)){
    temp_cor = Robocov_local(data_missing, alpha=alpha_vec[m])
    nlik_list[m] = nloglik(data_missing, as.matrix(nearPD(temp_cor)$mat))
    cat("Passed alpha value: ", alpha_vec[m], "\n")
  }
  final_alpha = alpha_vec[which.min(nlik_list)]
  robocov_local_cor =  Robocov_local(data_missing, alpha=final_alpha)




  df = rbind(c(norm(robocov_box_cor - corSigma, type = "F"),
               norm(robocov_local_cor - corSigma, type = "F"),
               norm(robocov_box_slack_cor - corSigma, type = "F"),
               norm(standard_cor - corSigma, type = "F"),
               norm(corshrink_cor - corSigma, type = "F")),
             c(
               norm(robocov_box_cor - corSigma, type = "2"),
               norm(robocov_local_cor - corSigma, type = "2"),
               norm(robocov_box_slack_cor - corSigma, type = "2"),
               norm(standard_cor - corSigma, type = "2"),
               norm(corshrink_cor - corSigma, type = "2")
             ),
             c(
               norm(robocov_box_cor - corSigma, type = "M"),
               norm(robocov_local_cor - corSigma, type = "M"),
               norm(robocov_box_slack_cor - corSigma, type = "M"),
               norm(standard_cor - corSigma, type = "M"),
               norm(corshrink_cor - corSigma, type = "M")
             ),
             c(
               norm(robocov_box_cor - corSigma, type = "I"),
               norm(robocov_local_cor - corSigma, type = "I"),
               norm(robocov_box_slack_cor - corSigma, type = "I"),
               norm(standard_cor - corSigma, type = "I"),
               norm(corshrink_cor - corSigma, type = "I")
             ),
             c(
               angle_norm(robocov_box_cor, corSigma),
               angle_norm(robocov_local_cor, corSigma),
               angle_norm(robocov_box_slack_cor, corSigma),
               angle_norm(standard_cor, corSigma),
               angle_norm(corshrink_cor, corSigma)
             ))
  df = t(df[,c(4, 5, 1, 3, 2)])

  rownames(df) = c("Standard", "CorShrink",
                   "Robocov_box", "Robocov_box_slack_cor", "Robocov_local")
  colnames(df) = c("Frobenius", "Operator", "Max.Modulus", "Infinity", "Angle")

  ll[[nsim]] = df
  cat("We are at simulation trial:", nsim, "\n")
}

save(df, file = paste0("/n/groups/price/kushal/Robocov/output/robocov_sim_toeplitz_n_",
                       N, "_p_", P, "_prop_", prop_missing, ".rda"))

