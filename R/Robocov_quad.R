#' @title Robust correlation estimation using quadratic data fidelity function.
#' @description A robust estimation of correlation matrix for data with missing entries using
#' a local second order Taylor series approximation of the Fisher z-score values resulting in a
#' quadratic data fidelity function. This algorithm is mentioned in the Supplementary Note of the paper.
#'
#' @param data_with_missing Samples by features data matrix. May contain missing entries (NA) values.
#' @param alpha The tuning parameter for the L-1 scaling of the correlation values.
#' @param loss Specify if we minimize L-1 (`lasso`), L-2 (`ridge`) or elastic-net (`elasticnet`)
#'             loss functions.
#'
#' @examples
#' data("sample_by_feature_data")
#' out = Robocov_quad(sample_by_feature_data, alpha = 1)
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
#' @importFrom Matrix nearPD
#' @export


Robocov_quad <- function(data_with_missing,
                         alpha= 1,
                         loss = c("lasso", "ridge", "elasticnet")){

  #library(CVXR)
  ##################  Building matrix of common samples for pairwise comparisons  ####################

  ##  B: binary matrix B_{N x P} similar to X_{N x P}
  ##  B^{T}B (P x P ) matrix with each entry equal to n_{ij}
  ##  n_{ii} = 0 y construction

  binary_indicator = matrix(1, nrow(data_with_missing), ncol(data_with_missing))
  binary_indicator[is.na(data_with_missing)]= 0
  common_samples = t(binary_indicator)%*%binary_indicator
  diag(common_samples) = 0

  #################  Pairwise correlations computation  ###############################

  ## C = cor(X) pairwise sample cov. matrix:

  pairwise_cor = cor(data_with_missing, use = "pairwise.complete.obs")
  pairwise_cor[is.na(pairwise_cor)] = 0
  pairwise_cor[pairwise_cor > 0.95] = 0.95
  pairwise_cor[pairwise_cor < -0.95] = -0.95
  diag(pairwise_cor) = 1
  common_samples[common_samples <= 2] = 2

  ################# Computing sample Fisher Z scores   ###########################

  ## Compute Z = 0.5 log ((1+r)/(1-r))

  pairwise_zscores = apply(pairwise_cor, c(1,2), function(x) return (0.5*log((1+x)/(1-x))))
  diag(pairwise_zscores) = 0

  ################  Samples adjustment - bias and sd   ########################

  samples_adjustment = 1/(common_samples - 1)  + 2/(common_samples - 1)^2
  bias_cor = samples_adjustment*pairwise_cor*(1 + pairwise_cor)*(1 - pairwise_cor)
  sd_cor = sqrt(samples_adjustment)*(1+pairwise_cor)*(1-pairwise_cor)
  estimate_cor = pairwise_cor + bias_cor
  inverse_sd_cor = 1/sd_cor
  diag(inverse_sd_cor) = 0

  ################  Convec optimization problem  #################################

  R <- Semidef(dim(pairwise_cor)[1])

  if(loss == "lasso"){
    obj <- Minimize(alpha*p_norm(R, 1) +
                      p_norm(mul_elemwise(inverse_sd_cor, square(R - estimate_cor)), 1))
  }else if (loss == "ridge"){
    obj <- Minimize(alpha*p_norm(R, 2) +
                      p_norm(mul_elemwise(inverse_sd_cor, square(R - estimate_cor)), 1))
  }else if (loss == "elasticnet"){
    obj <- Minimize(alpha*0.5*p_norm(R, 1) + alpha*0.5*p_norm(R, 2) +
                      p_norm(mul_elemwise(inverse_sd_cor, square(R - estimate_cor)), 1))
  }else{
    stop("loss must be one of lasso, ridge or elasticnet.")
  }

  constraints = list(diag(R) == 1)
  prob <- Problem(obj, constraints)
  result <- solve(prob)
  R_hat = cov2cor(as.matrix(result$getValue(R)))

  if(!is.null(colnames(data_with_missing))){
    rownames(R_hat) = colnames(data_with_missing)
    colnames(R_hat) = colnames(data_with_missing)
  }

  return(R_hat)
}
