
data("sample_by_feature_data")
standard = cor(sample_by_feature_data, use = "pairwise.complete.obs")
png("/Users/kushaldey/Documents/Robocov/vignettes/standard.png", width = 4, height = 4, units = 'in', res = 2000)
corrplot::corrplot(standard, diag = TRUE,
                   col = colorRampPalette(c("lightblue4", "lightblue2", "white", "indianred1", "indianred3"))(200),
                   tl.pos = "ld", tl.cex = 0.2, tl.col = "black",
                   rect.col = "white",na.label.col = "white",
                   method = "color", type = "lower", tl.srt=45)
dev.off()

robocov = Robocov_cor(sample_by_feature_data)
colnames(robocov) = colnames(standard)
rownames(robocov) = rownames(standard)
png("/Users/kushaldey/Documents/Robocov/vignettes/robocov.png", width = 4, height = 4, units = 'in', res = 2000)
corrplot::corrplot(robocov, diag = TRUE,
                   col = colorRampPalette(c("lightblue4", "lightblue2", "white", "indianred1", "indianred3"))(200),
                   tl.pos = "ld", tl.cex = 0.2, tl.col = "black",
                   rect.col = "white",na.label.col = "white",
                   method = "color", type = "lower", tl.srt=45)
dev.off()

probocov = Robocov_precision(sample_by_feature_data, alpha = 0.1)
colnames(probocov) = colnames(standard)
rownames(probocov) = rownames(standard)
png("/Users/kushaldey/Documents/Robocov/vignettes/probocov.png", width = 4, height = 4, units = 'in', res = 2000)
corrplot::corrplot(probocov, diag = TRUE,
                   col = colorRampPalette(c("lightblue4", "lightblue2", "white", "indianred1", "indianred3"))(200),
                   tl.pos = "ld", tl.cex = 0.2, tl.col = "black",
                   rect.col = "white",na.label.col = "white",
                   method = "color", type = "lower", tl.srt=45)
dev.off()
