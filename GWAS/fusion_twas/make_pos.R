#!/usr/bin/env Rscript

# ============================================================
# 批量读取 PEC_TWAS/weights 下的 .RDat 权重文件，生成 FUSION .pos 文件
# 同时使用 org.Hs.eg.db + clusterProfiler 将 Ensembl ID 转换为
# 基因 SYMBOL 和 GENENAME。
#
# 注意：
# 1. 原始 ID 列必须保留，用于与权重文件名匹配。
# 2. SYMBOL 和 GENENAME 作为附加注释列写入输出文件。
# ============================================================

required_packages <- c("org.Hs.eg.db", "clusterProfiler")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "缺少以下 R/Bioconductor 包：",
    paste(missing_packages, collapse = ", "),
    "\n请先运行：\n",
    "if (!requireNamespace(\"BiocManager\", quietly = TRUE)) ",
    "install.packages(\"BiocManager\")\n",
    "BiocManager::install(c(\"org.Hs.eg.db\", ",
    "\"clusterProfiler\"))"
  )
}

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(clusterProfiler)
})


rdat_dir <- "/data/UKBioBank/TWAS_PWAS_Data/Weights/PEC_TWAS/weights"

out_pos <- "/data/UKBioBank/TWAS_PWAS_Data/Weights/PEC_TWAS/PEC_TWAS.pos"

# 基因 ID 转换结果单独保存，避免向 FUSION .pos 文件加入额外列
out_annotation <- "/data/UKBioBank/TWAS_PWAS_Data/Weights/PEC_TWAS/PEC_TWAS_gene_annotation.tsv"

if (!dir.exists(rdat_dir)) {
  stop("错误：权重文件夹不存在：", rdat_dir)
}

out_dir <- dirname(out_pos)

if (!dir.exists(out_dir)) {
  stop("错误：输出目录不存在：", out_dir)
}

rdat_files <- list.files(
  path = rdat_dir,
  pattern = "\\.RDat$",
  full.names = TRUE
)

if (length(rdat_files) == 0) {
  stop("错误：文件夹中没有找到 .RDat 文件：", rdat_dir)
}

cat("找到", length(rdat_files), "个 .RDat 文件\n")

process_one_rdat <- function(file, base_dir) {

  env <- new.env()

  ok <- tryCatch({
    load(file, envir = env)
    TRUE
  }, error = function(e) {
    cat("读取失败：", file, "\n")
    cat("原因：", conditionMessage(e), "\n")
    FALSE
  })

  if (!ok) {
    return(NULL)
  }

  obj_names <- ls(env)

  if (!("snps" %in% obj_names)) {
    cat("跳过：", file, "，没有找到 snps 对象\n")
    return(NULL)
  }

  snps <- get("snps", envir = env)

  if (!is.data.frame(snps)) {
    cat("跳过：", file, "，snps 不是 data.frame\n")
    return(NULL)
  }

  if (nrow(snps) == 0) {
    cat("跳过：", file, "，snps 为空\n")
    return(NULL)
  }

  if (!all(c("V1", "V4") %in% colnames(snps))) {
    cat("跳过：", file, "，snps 中没有 V1 或 V4\n")
    cat("snps列名：", paste(colnames(snps), collapse = ", "), "\n")
    return(NULL)
  }

  chr_values <- snps$V1
  pos_values <- suppressWarnings(as.numeric(snps$V4))

  chr_values <- chr_values[!is.na(chr_values)]
  pos_values <- pos_values[!is.na(pos_values)]

  if (length(chr_values) == 0 || length(pos_values) == 0) {
    cat("跳过：", file, "，CHR 或 POS 信息为空\n")
    return(NULL)
  }

  chr <- unique(chr_values)

  if (length(chr) > 1) {
    cat("警告：", file, " 中检测到多个染色体：", paste(chr, collapse = ","), "\n")
    cat("将使用出现次数最多的染色体\n")
    chr <- names(sort(table(chr_values), decreasing = TRUE))[1]
  } else {
    chr <- chr[1]
  }

  p0 <- min(pos_values, na.rm = TRUE)
  p1 <- max(pos_values, na.rm = TRUE)

  wgt_path <- file.path(
    basename(normalizePath(base_dir)),
    basename(file)
  )

  id <- basename(file)
  id <- sub("\\.wgt\\.RDat$", "", id)
  id <- sub("\\.RDat$", "", id)

  out <- data.frame(
    WGT = wgt_path,
    ID  = id,
    CHR = chr,
    P0  = p0,
    P1  = p1,
    stringsAsFactors = FALSE
  )

  return(out)
}

res_list <- list()

for (i in seq_along(rdat_files)) {
  f <- rdat_files[i]

  cat("[", i, "/", length(rdat_files), "] 处理：", basename(f), "\n", sep = "")

  tmp <- process_one_rdat(f, rdat_dir)

  if (!is.null(tmp)) {
    res_list[[length(res_list) + 1]] <- tmp
  }
}

if (length(res_list) == 0) {
  stop("没有成功提取任何 .RDat 文件的信息，未生成 .pos 文件")
}

pos_df <- do.call(rbind, res_list)

pos_df$CHR_sort <- suppressWarnings(
  as.numeric(gsub("^chr", "", pos_df$CHR, ignore.case = TRUE))
)

pos_df <- pos_df[order(pos_df$CHR_sort, pos_df$P0), ]

pos_df$CHR_sort <- NULL

# 暂不写出 .pos 文件。先将 Ensembl ID 转换成基因 SYMBOL，
# 再把第二列 ID 组装为“ENSEMBL.SYMBOL”。

# ============================================================
# Ensembl ID 转换为基因 SYMBOL 和 GENENAME
# ============================================================
annotation_df <- pos_df

# 处理 ENSG00000123456.7 这种带版本号的 Ensembl ID。
annotation_df$ENSEMBL <- sub("\\.[0-9]+$", "", annotation_df$ID)

# 当文件名还含有其他字符时，尝试从中提取 ENSG 编号。
extracted_ensembl <- regmatches(
  annotation_df$ENSEMBL,
  regexpr("ENSG[0-9]+", annotation_df$ENSEMBL, ignore.case = TRUE)
)

has_extracted_id <- nzchar(extracted_ensembl)
annotation_df$ENSEMBL[has_extracted_id] <- toupper(
  extracted_ensembl[has_extracted_id]
)

valid_ensembl <- unique(
  annotation_df$ENSEMBL[
    grepl("^ENSG[0-9]+$", annotation_df$ENSEMBL)
  ]
)

annotation_df$SYMBOL <- NA_character_
annotation_df$GENENAME <- NA_character_

if (length(valid_ensembl) > 0) {
  cat("开始转换", length(valid_ensembl), "个唯一 Ensembl ID...\n")

  gene_info <- suppressMessages(
    clusterProfiler::bitr(
      valid_ensembl,
      fromType = "ENSEMBL",
      toType = c("SYMBOL", "GENENAME"),
      OrgDb = org.Hs.eg.db
    )
  )

  if (!is.null(gene_info) && nrow(gene_info) > 0) {
    gene_info <- aggregate(
      cbind(SYMBOL, GENENAME) ~ ENSEMBL,
      data = gene_info,
      FUN = function(x) {
        paste(unique(x[!is.na(x) & nzchar(x)]), collapse = ";")
      }
    )

    symbol_map <- setNames(gene_info$SYMBOL, gene_info$ENSEMBL)
    genename_map <- setNames(gene_info$GENENAME, gene_info$ENSEMBL)

    annotation_df$SYMBOL <- unname(symbol_map[annotation_df$ENSEMBL])
    annotation_df$GENENAME <- unname(genename_map[annotation_df$ENSEMBL])
  }

  converted_num <- sum(
    !is.na(annotation_df$SYMBOL) | !is.na(annotation_df$GENENAME)
  )

  cat(
    "成功注释：", converted_num, "/", nrow(annotation_df), "行\n",
    sep = ""
  )
} else {
  warning(
    "未从权重文件名中识别出有效的 Ensembl ID，",
    "SYMBOL 和 GENENAME 将为空"
  )
}

# ============================================================
# 生成 .pos 第二列：ENSEMBL.SYMBOL
# 示例：ENSG00000066583.ISR2
# 如果某个 Ensembl ID 无法转换，则第二列保留原始 ID。
# ============================================================
annotation_df$POS_ID <- ifelse(
  !is.na(annotation_df$SYMBOL) & nzchar(annotation_df$SYMBOL),
  paste0(annotation_df$ENSEMBL, ".", annotation_df$SYMBOL),
  annotation_df$ID
)

# 写出 FUSION .pos 文件，仍然严格保持五列。
pos_output_df <- data.frame(
  WGT = annotation_df$WGT,
  ID  = annotation_df$POS_ID,
  CHR = annotation_df$CHR,
  P0  = annotation_df$P0,
  P1  = annotation_df$P1,
  stringsAsFactors = FALSE
)

write.table(
  pos_output_df,
  file = out_pos,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "NA"
)

annotation_df <- annotation_df[, c(
  "WGT", "ID", "POS_ID", "ENSEMBL", "SYMBOL", "GENENAME", "CHR", "P0", "P1"
)]

write.table(
  annotation_df,
  file = out_annotation,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "NA"
)

cat("\n完成！\n")
cat("成功生成带 ENSEMBL.SYMBOL ID 的 .pos：", out_pos, "\n")
cat("成功生成基因注释表：", out_annotation, "\n")
cat("共写入：", nrow(pos_df), "行\n")
