
data("sample_by_feature_data")
standard = cor(sample_by_feature_data, use = "pairwise.complete.obs")
corrplot::corrplot(standard, diag = TRUE,
                   col = colorRampPalette(c("lightblue4", "lightblue2", "white", "indianred1", "indianred3"))(200),
                   tl.pos = "ld", tl.cex = 0.2, tl.col = "black",
                   rect.col = "white",na.label.col = "white",
                   method = "color", type = "lower", tl.srt=45)

robocov = Robocov_cor(sample_by_feature_data)
colnames(robocov) = colnames(standard)
rownames(robocov) = rownames(standard)
corrplot::corrplot(robocov, diag = TRUE,
                   col = colorRampPalette(c("lightblue4", "lightblue2", "white", "indianred1", "indianred3"))(200),
                   tl.pos = "ld", tl.cex = 0.2, tl.col = "black",
                   rect.col = "white",na.label.col = "white",
                   method = "color", type = "lower", tl.srt=45)

probocov = Robocov_precision(sample_by_feature_data, alpha = 0.1)
colnames(probocov) = colnames(standard)
rownames(probocov) = rownames(standard)
corrplot::corrplot(probocov, diag = TRUE,
                   col = colorRampPalette(c("lightblue4", "lightblue2", "white", "indianred1", "indianred3"))(200),
                   tl.pos = "ld", tl.cex = 0.2, tl.col = "black",
                   rect.col = "white",na.label.col = "white",
                   method = "color", type = "lower", tl.srt=45)
