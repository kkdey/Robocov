
data("sample_by_feature_data")

NUM_SIM=30
ll <- vector(mode="list", length=NUM_SIM)

for(nsim in 1:NUM_SIM){
  rand = sample(1:nrow(sample_by_feature_data),
                floor(0.5*nrow(sample_by_feature_data)), replace = F)
  train_data = sample_by_feature_data[rand, ]
  test_data = sample_by_feature_data[-rand, ]

  data_missing = train_data
  corSigma = as.matrix(cor(test_data, use = "pairwise.complete.obs"))
  corSigma[which(is.na(corSigma))] = 0

  standard_cor = cor(data_missing, use = "pairwise.complete.obs")
  standard_cor[which(is.na(standard_cor))] = 0
  cov_sample_ML <-  CorShrinkData(data_missing, sd_boot = FALSE,
                                  ash.control = list())
  corshrink_cor = cov2cor(cov_sample_ML$cor)
  robocov_box_cor = Robocov_box(data_with_missing = data_missing)

  alpha_vec = c(1e-03, 1e-02, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100)
  nlik_list = rep(0, length(alpha_vec))
  for(m in 1:length(alpha_vec)){
    temp_cor = Robocov_box_slack(data_missing, alpha=alpha_vec[m])
    nlik_list[m] = nloglik(data_missing, temp_cor)
    cat("Passed alpha value: ", alpha_vec[m], "\n")
  }
  final_alpha = alpha_vec[which.min(nlik_list)]
  robocov_box_slack_cor =  Robocov_box_slack(data_missing, alpha=final_alpha)


  alpha_vec = c(1e-03, 1e-02, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100)
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

save(df, file = paste0("/Users/kushaldey/Documents/Robocov-pages/output/gtex_predictive_robocov.rda"))

arr = array(0, c(5, 5, 30))
for(m in 1:30){
  arr[,,m] = ll[[m]]
}

mean_metric = apply(arr, c(1,2), function(x) return(mean(x[!is.na(x)])))
sd_metric = apply(arr, c(1,2), function(x) return(sd(x[!is.na(x)])))
