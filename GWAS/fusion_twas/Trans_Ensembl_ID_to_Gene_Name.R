# ============================================================
# 1. 加载依赖包
# ============================================================
if (!requireNamespace("readxl", quietly = TRUE))   install.packages("readxl")
if (!requireNamespace("biomaRt", quietly = TRUE))  install.packages("biomaRt")
if (!requireNamespace("writexl", quietly = TRUE))  install.packages("writexl")

library(readxl)
library(biomaRt)
library(writexl)

# ============================================================
# 2. 读取原始 xls 文件
#    原始文件第 1 行为 Supplementary Table 标题，第 2 行为列名，故 skip = 1
# ============================================================
input_file  <- "/home/user/pwas id.xls"          # 请根据实际路径修改
output_file <- "/home/user/pwas_id_with_gene_name.xlsx"

df <- read_xls(input_file, skip = 1)

# 查看列名确认结构
cat("Columns in file:\n")
print(colnames(df))

# ============================================================
# 3. 提取 Ensembl ID 并去除版本号
#    例如 ENSG00000155957.17 → ENSG00000155957
# ============================================================
ensembl_ids       <- df$ID
ensembl_ids_nover <- sub("\\..*", "", ensembl_ids)

cat("Total rows:", nrow(df), "\n")
cat("Unique Ensembl IDs (no version):", length(unique(ensembl_ids_nover)), "\n")

# ============================================================
# 4. 使用 biomaRt 查询 HGNC Gene Symbol
# ============================================================
cat("Querying Ensembl BioMart for gene symbols...\n")

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = unique(ensembl_ids_nover),
  mart       = mart
)

# 建立 Ensembl ID → Gene Symbol 的命名向量
id_map <- setNames(mapping$hgnc_symbol, mapping$ensembl_gene_id)

# 将空字符串转为 NA（表示无官方 Symbol）
id_map[id_map == ""] <- NA

# ============================================================
# 5. 将查询结果映射回数据框的 Gene Symbol 列
# ============================================================
df$`Gene Symbol` <- id_map[ensembl_ids_nover]

# 检查未匹配情况
unmatched <- unique(ensembl_ids_nover[is.na(df$`Gene Symbol`)])
if (length(unmatched) > 0) {
  cat("\nWarning:", length(unmatched), "unique ID(s) returned no HGNC symbol.\n")
  cat("Examples:\n")
  print(head(unmatched, 10))
} else {
  cat("\nAll IDs successfully mapped.\n")
}

# ============================================================
# 6. 保存结果
# ============================================================
# 方案 A：输出为 .xlsx（推荐，无需 Java，兼容 Excel 2007+）
write_xlsx(df, output_file)
cat("\nDone! File saved to:", normalizePath(output_file), "\n")

# ------------------------------------------------------------
# 方案 B：如需严格输出旧版 .xls（需安装 xlsx 包 + Java 环境）
# ------------------------------------------------------------
# if (!requireNamespace("xlsx", quietly = TRUE)) install.packages("xlsx")
# library(xlsx)
# write.xlsx(df, "pwas_id_with_gene_name.xls",
#            sheetName = "Sheet1", row.names = FALSE)