

####################  Trials for Robocov-box  method   ###########################

data("sample_by_feature_data")
out = Robocov_box(sample_by_feature_data)
corrplot::corrplot(as.matrix(out), diag = FALSE,
        col = colorRampPalette(c("blue", "white", "red"))(200),
        tl.pos = "td", tl.cex = 0.4, tl.col = "black",
        rect.col = "white",na.label.col = "white",
        method = "color", type = "upper")

out2 = Robocov_box_slack(sample_by_feature_data, alpha=1)
corrplot::corrplot(as.matrix(out2), diag = FALSE,
                   col = colorRampPalette(c("blue", "white", "red"))(200),
                   tl.pos = "td", tl.cex = 0.4, tl.col = "black",
                   rect.col = "white",na.label.col = "white",
                   method = "color", type = "upper")


out = Robocov_local(sample_by_feature_data, alpha = 1)
corrplot::corrplot(as.matrix(out), diag = FALSE,
      col = colorRampPalette(c("blue", "white", "red"))(200),
       tl.pos = "td", tl.cex = 0.4, tl.col = "black",
       rect.col = "white",na.label.col = "white",
       method = "color", type = "upper")


out = Robocov_precision(sample_by_feature_data, alpha = 0.1)
corrplot::corrplot(as.matrix(out), diag = FALSE,
        col = colorRampPalette(c("blue", "white", "red"))(200),
        tl.pos = "td", tl.cex = 0.4, tl.col = "black",
        rect.col = "white",na.label.col = "white",
        method = "color", type = "upper")

