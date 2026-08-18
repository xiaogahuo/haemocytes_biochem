options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")) ##指定install.packages安装镜像，这个是清华镜像
install.packages("qqman") # location of installation can be changed but has to correspond with the library location
library("qqman")
# results_log <- read.table("logistic_results.assoc_2.logistic", head=TRUE)
# jpeg("QQ-Plot_logistic.jpeg")
# qq(results_log$P, main = "Q-Q plot of GWAS p-values : log")
# dev.off()

results_as <- read.table("logistic_results.assoc.logistic2", head=TRUE)
jpeg("QQ-Plot_assoc.jpeg")
qq(results_as$`P`, main = "Q-Q plot of GWAS with Logistic Regression P-values : log")
dev.off()

