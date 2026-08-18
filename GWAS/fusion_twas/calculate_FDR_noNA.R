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
sigout_file <- get_arg("--sigout")

if (is.na(fdr_threshold)) {
  stop("FDR threshold must be numeric, for example: 0.2")
}

cat("Input file:", input_file, "\n")
cat("Output file:", out_file, "\n")
cat("FDR threshold:", fdr_threshold, "\n")
cat("Significant output file:", sigout_file, "\n")

dat <- read.table(
  input_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  na.strings = c("NA", "NaN", "", " ")
)

dat$TWAS.P <- as.numeric(dat$TWAS.P)

dat2 <- dat[!is.na(dat$TWAS.P), ]

dat2$FDR <- p.adjust(dat2$TWAS.P, method = "BH")

write.table(
  dat2,
  file = out_file,
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

sig_dat <- dat2[dat2$FDR < fdr_threshold, ]

write.table(
  sig_dat,
  file = sigout_file,
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

cat("Total rows after removing NA TWAS.P:", nrow(dat2), "\n")
cat("Rows with FDR <", fdr_threshold, ":", nrow(sig_dat), "\n")
cat("FDR calculation finished.\n")