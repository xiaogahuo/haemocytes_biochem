options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")) ##指定install.packages安装镜像，这个是清华镜像
install.packages("qqman") # location of installation can be changed but has to correspond with the library location
library("qqman")
# results_log <- read.table("logistic_results.assoc_2.logistic", head=TRUE)
# jpeg("Logistic_manhattan.jpeg")
# manhattan(results_log,chr="CHR",bp="BP",p="P",snp="SNP", main = "Manhattan plot: logistic")
# dev.off()

results_as <- read.table("logistic_results.assoc.logistic2", head=TRUE)
jpeg("assoc_manhattan.jpeg")
manhattan(results_as,chr="CHR",bp="BP",p="P",snp="SNP", col = c("blue", "orange"), main = "Manhattan plot: GWAS with Logistic Regression")
dev.off()  




