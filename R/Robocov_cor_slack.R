#' @title  Robocov correlation estimation using box constraints on Fisher Z-scores with slack variables
#' @description A robust estimation of correlation matrix for data with missing entries using box
#' constraint on the difference between the population correlation matrix and pairwise sample correlation
#' matrix, with added L-1 penalty on the slack variables. This is a more flexible version of the
#' Robocov_cor function.
#'
#' @param data_with_missing Samples by features data matrix. May contain missing entries (NA) values.
#' @param alpha The tuning parameter for the L-1 shrinkage of the slack variables.
#' @param loss Specify if we minimize L-1 (`lasso`), L-2 (`ridge`) or elastic-net (`elasticnet`)
#'             loss functions.
#' @examples
#' data("sample_by_feature_data")
#' out = Robocov_cor_slack(sample_by_feature_data, alpha = 1)
#' corrplot::corrplot(as.matrix(out), diag = FALSE,
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

Robocov_cor_slack <- function(data_with_missing,
                              alpha = 1,
                              loss = c("lasso", "ridge", "elasticnet")){

  ##################  Building matrix of common samples for pairwise comparisons  ####################

  binary_indicator = matrix(1, nrow(data_with_missing), ncol(data_with_missing))
  binary_indicator[is.na(data_with_missing)]= 0
  common_samples = t(binary_indicator)%*%binary_indicator
  diag(common_samples) = 0

  #################  Pairwise correlations computation  ###############################

  pairwise_cor = cor(data_with_missing, use = "pairwise.complete.obs")
  pairwise_cor[is.na(pairwise_cor)] = 0
  pairwise_cor[pairwise_cor > 0.95] = 0.95
  pairwise_cor[pairwise_cor < -0.95] = -0.95
  diag(pairwise_cor) = 1
  common_samples[common_samples <= 2] = 2

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

  ###############  Convex optimization  ######################

  #library(CVXR)
  R <- Semidef(dim(pairwise_cor)[1])
  slack <- Variable(dim(pairwise_cor)[1], dim(pairwise_cor)[1])

  if(loss == "lasso"){
    obj <- Minimize(p_norm(R, 1) + alpha*p_norm(slack, 1))
  }else if (loss == "ridge"){
    obj <- Minimize(p_norm(R, 2) + alpha*p_norm(slack, 1))
  }else if (loss == "elasticnet"){
    obj <- Minimize(0.5*p_norm(R, 2) + 0.5*p_norm(R,1) + alpha*p_norm(slack, 1))
  }else{
    stop("loss must be one of lasso, ridge or elasticnet.")
  }

  constraints = list(diag(R) == 1,
                     diag(slack) == 0,
                     R <= pairwise_cor + constrained_overall_bound + slack,
                     R >= pairwise_cor - constrained_overall_bound - slack)
  prob <- Problem(obj, constraints)
  result <- solve(prob)
  R_hat = as.matrix(result$getValue(R))
  R_hat = cov2cor(R_hat)

  if(!is.null(colnames(data_with_missing))){
    rownames(R_hat) = colnames(data_with_missing)
    colnames(R_hat) = colnames(data_with_missing)
  }

  return(R_hat)
}
