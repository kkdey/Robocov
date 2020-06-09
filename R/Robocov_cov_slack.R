#' @title  Robocov covariance estimation using box constraints on Fisher Z-scores with slack variables
#' @description A robust estimation of covariance matrix for data with missing entries using box
#' constraint on the difference between the population covariance matrix and pairwise sample covariance
#' matrix, with added L-1 penalty on the slack variables. This is a more flexible version of the
#' Robocov_cov function.
#'
#' @param data_with_missing Samples by features data matrix. May contain missing entries (NA) values.
#' @param alpha The tuning parameter for the L-1 shrinkage of the slack variables.
#' @param loss Specify if we minimize L-1 (`lasso`), L-2 (`ridge`) or elastic-net (`elasticnet`)
#'             loss functions.
#' @examples
#' data("sample_by_feature_data")
#' out = Robocov_cov_slack(sample_by_feature_data, alpha = 1)
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

Robocov_cov_slack <- function(data_with_missing,
                              alpha = 1,
                              loss = c("lasso", "ridge", "elasticnet")){

  ##################  Building matrix of common samples for pairwise comparisons  ####################

  binary_indicator = matrix(1, nrow(data_with_missing), ncol(data_with_missing))
  binary_indicator[is.na(data_with_missing)]= 0
  common_samples = t(binary_indicator)%*%binary_indicator
  diag(common_samples) = 0

  #################  Pairwise correlations computation  ###############################

  pairwise_cov = cov(data_with_missing, use = "pairwise.complete.obs")
  sdvec = sqrt(diag(pairwise_cov))
  pairwise_cor = cov2cor(pairwise_cov)
  pairwise_cor[is.na(pairwise_cor)] = 0
  pairwise_cor[pairwise_cor > 0.95] = 0.95
  pairwise_cor[pairwise_cor < -0.95] = -0.95
  diag(pairwise_cor) = 1
  common_samples[common_samples <= 2] = 2
  pairwise_cov = diag(sdvec) %*% pairwise_cor %*% diag(sdvec)

  ################# Computing sample Fisher Z scores   ###########################

  pairwise_zscores = apply(pairwise_cor, c(1,2), function(x) return (0.5*log((1+x)/(1-x))))
  diag(pairwise_zscores) = 0

  ###############  Bounds on the correlations  ####################################

  bound1 = 4*3.3*exp(2*pairwise_zscores)/((exp(2*pairwise_zscores) + 1)^2)
  zscores_sd = sqrt(1/(common_samples - 1) + 2/(common_samples - 1)^2)
  bound2 = bound1*zscores_sd
  bound3 = zscores_sd^2*4.2
  overall_bound = bound2 + bound3
  constrained_overall_bound = apply(overall_bound, c(1,2), function(x) return(min(2,x)))
  diag(constrained_overall_bound) = 0
  constrained_overall_bound2 = diag(sdvec) %*% constrained_overall_bound %*% diag(sdvec)


  ###############  Convex optimization  ######################

  #library(CVXR)

  Sigma <- Semidef(dim(pairwise_cor)[1])
  slack <- Variable(dim(pairwise_cor)[1], dim(pairwise_cor)[1])
  if(loss == "lasso"){
    obj <- Minimize(p_norm(Sigma, 1))
  }else if (loss == "ridge"){
    obj <- Minimize(p_norm(Sigma, 2))
  }else if (loss == "elasticnet"){
    obj <- Minimize(0.5*p_norm(Sigma, 2) + 0.5*p_norm(Sigma,1))
  }else{
    stop("loss must be one of lasso, ridge or elasticnet.")
  }


  constraints = list(diag(slack) == 0,
                     Sigma <= pairwise_cov + constrained_overall_bound2 + slack,
                     Sigma >= pairwise_cov - constrained_overall_bound2 - slack)
  prob <- Problem(obj, constraints)
  result <- solve(prob)
  Sigma_hat = as.matrix(result$getValue(Sigma))

  if(!is.null(colnames(data_with_missing))){
    rownames(Sigma_hat) = colnames(data_with_missing)
    colnames(Sigma_hat) = colnames(data_with_missing)
  }

  return(Sigma_hat)
}
