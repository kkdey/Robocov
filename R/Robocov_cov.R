#' @title Robocov covariance estimation using box constraints on Fisher Z-scores
#' @description A robust estimation of covariance matrix for data with missing entries using box
#' constraint on the difference between the population covariance matrix and pairwise sample covariance
#' matrix.
#'
#' @param data_with_missing Samples by features data matrix. May contain missing entries
#'        (NA) values.
#' @param loss Specify if we minimize L-1 (`lasso`), L-2 (`ridge`) or elastic-net (`elasticnet`)
#'             loss functions.
#'
#' @examples
#' data("sample_by_feature_data")
#' out = Robocov_cov(sample_by_feature_data)
#' corrplot::corrplot(as.matrix(cov2cor(out)), diag = FALSE,
#'         col = colorRampPalette(c("blue", "white", "red"))(200),
#'         tl.pos = "td", tl.cex = 0.4, tl.col = "black",
#'         rect.col = "white",na.label.col = "white",
#'         method = "color", type = "upper")
#'
#' @keywords box shrinkage, correlation
#' @import CVXR
#' @importFrom corrplot corrplot
#' @importFrom stats cor sd cov2cor
#' @export

Robocov_cov <- function(data_with_missing,
                        loss = c("lasso", "ridge", "elasticnet")){

  if(length(loss) == 3){
    loss = "lasso"   ## L1 norm penalty as default to induce sparsity
  }
  ##################  Building matrix of common samples for pairwise comparisons  ####################

  ##  B: binary matrix B_{N x P} similar to X_{N x P}
  ##  B^{T}B (P x P ) matrix with each entry equal to n_{ij}
  ##  n_{ii} = 0 y construction

  binary_indicator = matrix(1, nrow(data_with_missing), ncol(data_with_missing))
  binary_indicator[is.na(data_with_missing)]= 0
  common_samples = t(binary_indicator)%*%binary_indicator
  diag(common_samples) = 0

  #################  Pairwise correlations computation  ###############################

  ## C = cov(X) pairwise sample cov. matrix:

  pairwise_cov = cov(data_with_missing, use = "pairwise.complete.obs")
  pairwise_cov[is.na(pairwise_cov)] = 0
  sdvec = sqrt(diag(pairwise_cov))
  pairwise_cor = cov2cor(pairwise_cov)
  pairwise_cor[pairwise_cor > 0.95] = 0.95
  pairwise_cor[pairwise_cor < -0.95] = -0.95
  diag(pairwise_cor) = 1
  common_samples[common_samples <= 2] = 2
  pairwise_cov = diag(sdvec) %*% pairwise_cor %*% diag(sdvec)

  ################# Computing sample Fisher Z scores   ###########################

  ## Compute Z = 0.5 log ((1+r)/(1-r))

  pairwise_zscores = apply(pairwise_cor, c(1,2), function(x) return (0.5*log((1+x)/(1-x))))
  diag(pairwise_zscores) = 0

  ###############  Bounds on the correlations  ####################################

  ## Compute the C_{ij} constant upper bound

  bound1 = 12*exp(2*pairwise_zscores)/((exp(2*pairwise_zscores) + 1)^2)
  zscores_sd = sqrt(1/(common_samples - 1) + 2/(common_samples - 1)^2)
  bound2 = bound1*zscores_sd
  bound3 = zscores_sd^2*2*sqrt(3)
  overall_bound = bound2 + bound3
  constrained_overall_bound = apply(overall_bound, c(1,2), function(x) return(min(2,x)))
  diag(constrained_overall_bound) = 0
  constrained_overall_bound2 = diag(sdvec) %*% constrained_overall_bound %*% diag(sdvec)

  ###############  Convex optimization  ######################
  #library(CVXR)

  Sigma <- Semidef(dim(pairwise_cor)[1])
  if(loss == "lasso"){
    obj <- Minimize(p_norm(Sigma, 1))
  }else if (loss == "ridge"){
    obj <- Minimize(p_norm(Sigma, 2))
  }else if (loss == "elasticnet"){
    obj <- Minimize(0.5*p_norm(Sigma, 2) + 0.5*p_norm(Sigma,1))
  }else{
    stop("loss must be one of lasso, ridge or elasticnet.")
  }

  constraints = list(Sigma <= pairwise_cov + constrained_overall_bound2,
                     Sigma >= pairwise_cov - constrained_overall_bound2)
  prob <- Problem(obj, constraints)
  result <- solve(prob)
  Sigma_hat = as.matrix(result$getValue(Sigma))
  return(Sigma_hat)
}
