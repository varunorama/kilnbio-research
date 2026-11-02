# Complete-feature implementation of the published protein and T2D projections.
# Invoke via clocks.py to verify the upstream resources before loading them.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% 4:5)
kind <- args[1]
cache <- args[2]
tab <- function(p) read.delim(p, check.names = FALSE, colClasses = "character")
beta <- tab(args[3])
ids <- names(beta)[-1]
stopifnot(names(beta)[1] == "probe", length(ids) > 0, all(nzchar(ids)), !anyDuplicated(ids),
          !anyNA(beta$probe), all(nzchar(beta$probe)), !anyDuplicated(beta$probe))
x <- as.matrix(beta[, -1, drop = FALSE])
rownames(x) <- beta$probe
storage.mode(x) <- "double"
if (any(!is.finite(x)) || any(x < 0 | x > 1)) stop("Expected finite methylation beta values in [0,1]")
out <- data.frame(sample = ids, check.names = FALSE)
if (kind == "proteins") {
  if (length(ids) < 2) stop("The authors' cohort-standardized projection requires multiple samples")
  k <- read.csv(file.path(cache, "KORA_weights_finalised_120121.csv"))
  l <- read.csv(file.path(cache, "LBC_weights_finalised_120121.csv"))
  weights <- rbind(data.frame(probe = k$CpG.site, w = k$Coefficient, key = paste0("KORA.", k$Predictor)),
                   data.frame(probe = l$CpG_Site, w = l$Coefficient, key = paste0("LBC.", l$Predictor)))
  stopifnot(length(unique(weights$key)) == 109, !anyDuplicated(weights[c("probe", "key")]))
  needed <- unique(weights$probe)
  if (!all(needed %in% rownames(x))) stop("The complete 109-score panel requires all published probes")
  z <- t(scale(t(x[needed, , drop = FALSE])))
  if (any(!is.finite(z))) stop("A required probe has zero between-sample variance")
  for (w in split(weights, weights$key)) {
    out[[w$key[1]]] <- drop(crossprod(w$w, z[w$probe, , drop = FALSE]))
  }
} else if (kind == "diabetes") {
  w <- readRDS(file.path(cache, "directEpiScoreCoefficients.rds"))
  stopifnot(length(w) == 145L, all(names(w) %in% rownames(x)))
  out$direct_T2D_EpiScore <- drop(crossprod(as.numeric(w), x[names(w), , drop = FALSE]))
  if (length(args) == 5) {
    meta <- tab(args[5])
    required <- c("age", "male", "bmi", "family_diabetes", "hypertension")
    stopifnot(all(c("sample", required) %in% names(meta)), !anyDuplicated(meta$sample), setequal(meta$sample, ids))
    m <- data.matrix(as.data.frame(lapply(meta[match(ids, meta$sample), required], as.numeric)))
    stopifnot(all(is.finite(m)), all(m[, "age"] > 0), all(m[, "bmi"] > 0),
              all(m[, c("male", "family_diabetes", "hypertension")] %in% 0:1))
    m <- cbind(m, direct = out$direct_T2D_EpiScore)
    center <- c(48.28064, 0.3796568, 26.83208, 0.2015488, 0.1391796, 12.43837)
    clinical <- readRDS(file.path(cache, "dModelCoefficients.rds"))
    keys <- c("Age", "sexM", "bmi", "family_diabetes", "high_BP", "dScore")
    stopifnot(all(keys %in% names(clinical)))
    linear <- drop(sweep(m, 2, center) %*% as.numeric(clinical[keys]))
    out$research_10y_T2D_probability <- -expm1(-0.02033703 * exp(linear))
  }
} else stop("Unknown projection")
write.table(out, args[4], sep = "\t", row.names = FALSE, quote = FALSE)
