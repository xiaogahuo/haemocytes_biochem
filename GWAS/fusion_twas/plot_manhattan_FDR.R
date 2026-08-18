library(CMplot)
library(tidyverse)
library(biomaRt)

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) {
    stop(paste("Missing argument:", flag))
  }
  return(args[idx + 1])
}

input_file <- get_arg("--input")
out_file <- get_arg("--out")
fdr_threshold <- as.numeric(get_arg("--fdr"))

if (is.na(fdr_threshold)) {
  stop("FDR threshold must be numeric, for example: 0.1 or 0.2")
}

cat("Input file:", input_file, "\n")
cat("Output PDF:", out_file, "\n")
cat("FDR threshold:", fdr_threshold, "\n")

a_raw <- read.table(
  file = input_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

a_raw <- tibble::as_tibble(a_raw)
ensembl_ids_nover <- sub("\\..*$", "", a_raw$ID)

cat("Querying Ensembl for gene symbols...\n")
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = unique(ensembl_ids_nover),
  mart = mart
)
id_map <- setNames(mapping$hgnc_symbol, mapping$ensembl_gene_id)
id_map[id_map == ""] <- NA
a_raw$SYMBOL <- id_map[ensembl_ids_nover]

a <- data.frame(
  SNP = a_raw$SYMBOL,
  Chromosome = as.numeric(a_raw$CHR),
  Position = as.numeric(a_raw$P0),
  trait1 = as.numeric(a_raw$TWAS.P),
  FDR = as.numeric(a_raw$FDR)
)

a <- a %>%
  filter(
    !is.na(SNP),
    !is.na(Chromosome),
    !is.na(Position),
    !is.na(trait1),
    !is.na(FDR),
    trait1 > 0,
    Chromosome %in% 1:22
  )

fdr_hits <- a %>%
  filter(FDR < fdr_threshold) %>%
  arrange(FDR, trait1)

cat("Number of genes with FDR <", fdr_threshold, ":", nrow(fdr_hits), "\n")

# 自动设置纵坐标上限，避免 -log10(P) > 5 的点被截掉
max_y <- ceiling(max(-log10(a$trait1), na.rm = TRUE)) + 1
cat("Maximum -log10(P):", max(-log10(a$trait1), na.rm = TRUE), "\n")
cat("Y axis upper limit:", max_y, "\n")
out_dir <- dirname(out_file)
out_base <- tools::file_path_sans_ext(basename(out_file))
old_wd <- getwd()
setwd(out_dir)

if (nrow(fdr_hits) > 0) {
  fdr_p_cutoff <- max(fdr_hits$trait1)
  cat("Approximate TWAS.P cutoff for FDR <", fdr_threshold, ":", fdr_p_cutoff, "\n")
  fdr_col <- ifelse(fdr_hits$Chromosome %% 2 == 1, "#029658", "#fc6404")

  CMplot(
    a[, c("SNP", "Chromosome", "Position", "trait1")],
    file = "pdf",
    plot.type = "m",

    threshold = fdr_p_cutoff,
    threshold.col = "red",
    threshold.lty = 2,
    threshold.lwd = 1,

    col = c("#029658", "#fc6404"),

    cex = 0.6,
    signal.cex = 0.8,
    ylim = c(0, max_y),
    chr.labels = 1:22,
    # 所有 FDR < 阈值的点都高亮
    highlight = fdr_hits$SNP,
    highlight.col = fdr_col,
    highlight.cex = 1.0,

    # 所有 FDR < 阈值的点都标完整 ID
    highlight.text = fdr_hits$SNP,
    highlight.text.col = "black",
    highlight.text.cex = 0.7,

    file.name = out_base,
    dpi = 300,
    width = 14,
    height = 7
  )

} else {
  CMplot(
    a[, c("SNP", "Chromosome", "Position", "trait1")],
    file = "pdf",
    plot.type = "m",

    col = c("#029658", "#fc6404"),

    cex = 0.6,
    signal.cex = 0.8,
    ylim = c(0, max_y),
    chr.labels = 1:22,

    file.name = out_base,
    dpi = 300,
    width = 14,
    height = 7
  )
}

setwd(old_wd)
cat("Manhattan plot finished:", out_file, "\n")