#!/usr/bin/env Rscript

# 加载包
suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(clusterProfiler)
})

# 定义 Ensembl ID
ensembl_id <- "ENSG00000127054"

# 使用 bitr 进行 ID 转换
gene_info <- bitr(
  ensembl_id,
  fromType = "ENSEMBL",
  toType = c("SYMBOL", "GENENAME"),
  OrgDb = org.Hs.eg.db
)

# 查看结果
print(gene_info)
