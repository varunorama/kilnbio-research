# One observation per donor per modality. Only the background molecular matrices
# are measured. Biology labels, controls, and handling effects below are planted.
args <- commandArgs(trailingOnly = TRUE)
source(args[1])
load(file.path(args[2], "CLL_data.RData"))
out <- args[3]
blocks <- CLL_data[c("mRNA", "Methylation")]
ids <- sort(Reduce(intersect, lapply(blocks, function(x) colnames(x)[colSums(is.finite(x)) > 0])))
set.seed(2718)
ids <- sample(ids, 80)
biology <- factor(rep(c("A", "B"), each = 40))
effect <- function(y) colMeans(y[biology == "B", , drop = FALSE]) - colMeans(y[biology == "A", , drop = FALSE])
rms <- function(x) sqrt(mean(x ^ 2))
rows <- list()
for (assay in names(blocks)) {
  block <- t(blocks[[assay]][, ids])
  ok <- colSums(is.finite(block)) == nrow(block)
  selected <- head(which(ok & apply(block, 2, sd, na.rm = TRUE) > 0), 60)
  stopifnot(length(selected) == 60)
  clean <- scale(block[, selected])
  # The first 10 positions within each block respond to the same biological
  # group. The next 10 positions respond only within the RNA block.
  clean[biology == "B", 1:10] <- clean[biology == "B", 1:10] + 1.4
  only_rna <- 11:20
  if (assay == "mRNA") clean[biology == "B", only_rna] <- clean[biology == "B", only_rna] + 1.2
  batch <- unlist(lapply(1:2, function(k) sample(rep(c("X", "Y"), each = 20))))
  drift <- as.numeric(scale(rnorm(80)))
  weights <- seq(-1.5, 1.5, length.out = ncol(clean))
  observed <- clean + outer(as.numeric(batch == "Y") - 0.5, rev(weights) * 2) + outer(drift, weights)
  controls <- sapply(seq_len(6), function(i) drift + rnorm(80, sd = 0.15))
  colnames(controls) <- paste0("spike_", seq_len(ncol(controls)))
  input <- cbind(observed, controls)
  fits <- list(unadjusted = observed,
               labeled_batch = adjust_assay(input, biology, batch)$adjusted,
               controls_only = adjust_assay(input, biology, controls = colnames(controls))$adjusted,
               batch_and_controls = adjust_assay(input, biology, batch, colnames(controls))$adjusted)
  for (method in names(fits)) {
    stopifnot(!is.null(fits[[method]]))
    estimate <- fits[[method]][, colnames(clean), drop = FALSE]
    rows[[length(rows) + 1L]] <- data.frame(assay, method, donors = length(ids),
      technical_RMSE = rms(estimate - clean),
      effect_RMSE = rms(effect(estimate) - effect(clean)),
      shared_block_contrast = mean(effect(estimate)[1:10]),
      oracle_shared_block_contrast = mean(effect(clean)[1:10]),
      RNA_specific_region_contrast = mean(effect(estimate)[only_rna]),
      oracle_RNA_specific_region_contrast = mean(effect(clean)[only_rna]))
  }
  stopifnot(is.null(adjust_assay(input, biology)$adjusted),
            is.null(adjust_assay(input, biology, as.character(biology))$adjusted))
}
write.table(do.call(rbind, rows), file.path(out, "noise_challenge.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
