# Per-assay exploratory adjustment on a continuous, already processed scale.
# The input rows are independent donors. Known batch levels or stable controls
# are the only technical anchors. Adjustments are centered to preserve the mean.
adjust_assay <- function(y, biology, batch = NULL, controls = character()) {
  stopifnot(is.matrix(y), is.numeric(y), all(is.finite(y)), !is.null(rownames(y)),
            !is.null(colnames(y)), !anyDuplicated(rownames(y)), !anyDuplicated(colnames(y)),
            length(biology) == nrow(y), !anyNA(biology), all(controls %in% colnames(y)),
            !anyDuplicated(controls))
  decline <- function(reason) list(status = reason, adjusted = NULL)
  if (length(unique(biology)) < 2) return(decline("biology_does_not_vary"))
  if (is.numeric(biology) && !all(is.finite(biology))) stop("Nonfinite biological covariate")
  protected <- model.matrix(~ biology)
  nuisance <- matrix(numeric(), nrow(y), 0)
  if (!is.null(batch)) {
    stopifnot(length(batch) == nrow(y), !anyNA(batch))
    batch <- factor(batch)
    if (nlevels(batch) > 1) {
      if (min(table(batch)) < 3L) return(decline("too_few_donors_per_batch"))
      if (is.numeric(biology)) {
        ranges <- vapply(split(biology, batch), range, numeric(2))
        if (max(ranges[1, ]) >= min(ranges[2, ])) return(decline("no_common_biological_range"))
      } else if (min(table(batch, biology)) < 2L) {
        return(decline("insufficient_biology_by_batch_overlap"))
      }
      nuisance <- model.matrix(~ batch)[, -1, drop = FALSE]
    }
  }
  if (length(controls)) {
    if (length(controls) < 3L) return(decline("fewer_than_three_control_features"))
    cmat <- y[, controls, drop = FALSE]
    if (any(apply(cmat, 2, sd) == 0)) return(decline("constant_control_feature"))
    control_fit <- prcomp(cmat, scale. = TRUE, rank. = 1)
    nuisance <- cbind(nuisance, control_PC1 = control_fit$x[, 1])
  }
  if (!ncol(nuisance)) return(decline("no_variable_technical_anchor"))
  nuisance <- scale(nuisance, center = TRUE, scale = FALSE)
  design <- cbind(protected, nuisance)
  if (nrow(y) - ncol(design) < 3L || qr(design)$rank < ncol(design)) {
    return(decline("not_enough_independent_design_information"))
  }
  # Unit-length columns make this screen invariant to the units of age or controls.
  standardized <- scale(design[, -1, drop = FALSE])
  if (kappa(standardized, exact = TRUE) > 10) return(decline("poorly_separated_design_terms"))
  fit <- lm.fit(design, y)
  idx <- ncol(protected) + seq_len(ncol(nuisance))
  shift <- nuisance %*% fit$coefficients[idx, , drop = FALSE]
  adjusted <- y - shift
  dimnames(adjusted) <- dimnames(y)
  projection <- NULL
  if (length(controls) && identical(colnames(nuisance), "control_PC1")) {
    projection <- list(controls = controls, center = control_fit$center, scale = control_fit$scale,
                       rotation = control_fit$rotation[, 1], score_mean = unname(attr(nuisance, "scaled:center")[1]),
                       coefficients = fit$coefficients[idx, , drop = TRUE])
  }
  list(status = "adjusted_exploratory", adjusted = adjusted, control_projection = projection,
       nuisance = nuisance, biology_coefficients = fit$coefficients[seq_len(ncol(protected)), , drop = FALSE])
}

# Apply a frozen control-only fit to unseen samples. No test labels or test-batch centering.
project_control_adjustment <- function(fit, y) {
  p <- fit$control_projection
  if (is.null(p)) stop("Requires a successful fit with only a control-feature nuisance term")
  stopifnot(is.matrix(y), is.numeric(y), all(is.finite(y)), !anyDuplicated(colnames(y)),
            setequal(colnames(y), names(p$coefficients)))
  z <- scale(y[, p$controls, drop=FALSE], center=p$center, scale=p$scale)
  score <- drop(z %*% p$rotation) - p$score_mean
  y - tcrossprod(as.numeric(score), p$coefficients[colnames(y)])
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!(length(args) %in% 4:5)) stop("Use: Rscript adjust.R matrix.tsv design.tsv categorical|continuous new_output_dir [controls.txt]")
  if (!args[3] %in% c("categorical", "continuous")) stop("Select categorical or continuous biology")
  read_tab <- function(path) read.delim(path, check.names = FALSE, colClasses = "character")
  data <- read_tab(args[1])
  meta <- read_tab(args[2])
  stopifnot(names(data)[1] == "donor", all(c("donor", "biology") %in% names(meta)),
            !anyNA(data$donor), all(nzchar(data$donor)), !anyDuplicated(data$donor),
            !anyDuplicated(meta$donor), setequal(meta$donor, data$donor))
  meta <- meta[match(data$donor, meta$donor), , drop = FALSE]
  y <- as.matrix(data[, -1, drop = FALSE])
  storage.mode(y) <- "double"
  rownames(y) <- data$donor
  stopifnot(!anyNA(meta$biology), all(nzchar(meta$biology)))
  if (!is.null(meta$batch)) stopifnot(!anyNA(meta$batch), all(nzchar(meta$batch)))
  bio <- if (args[3] == "continuous") as.numeric(meta$biology) else factor(meta$biology)
  control_ids <- if (length(args) == 5) readLines(args[5]) else character()
  fit <- adjust_assay(y, bio, meta$batch, control_ids)
  if (!dir.create(args[4], recursive = TRUE, showWarnings = FALSE)) stop("Use a new output directory")
  writeLines(fit$status, file.path(args[4], "status.txt"))
  if (!is.null(fit$adjusted)) {
    saveRDS(fit, file.path(args[4], "fit.rds"))
    write.table(data.frame(donor = rownames(y), fit$adjusted, check.names = FALSE),
                file.path(args[4], "exploratory_adjusted.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}
