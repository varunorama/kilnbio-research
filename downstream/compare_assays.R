# Four-arm comparison of two measured assays on a common development and holdout cohort.
fit_assays <- function(donors, response, covariates, blocks, recipes, train, test,
                       fold, endpoint, out, gender_coding=NULL) {
  write_tsv <- function(x, name) {
    write.table(x, file.path(out, name), sep="\t", quote=FALSE, row.names=FALSE, na="")
  }
  # Learn all transformations from the fit participants. Keep at most 24 features
  # per molecular block, ranked by association with the covariate-adjusted outcome.
  learn_block <- function(x, train_rows, residual_y) {
    train_x <- x[train_rows, , drop = FALSE]
    center <- colMeans(train_x, na.rm = TRUE)
    centered <- sweep(train_x, 2, center)
    observed <- colSums(is.finite(centered))
    centered[!is.finite(centered)] <- 0 # training-mean imputation
    scale <- sqrt(colSums(centered ^ 2) / pmax(1, observed - 1))
    keep <- which(observed >= ceiling(0.8 * length(train_rows)) & is.finite(scale) & scale > 0)
    stopifnot(length(keep) > 0)
    association <- abs(drop(crossprod(centered[, keep, drop = FALSE], residual_y)) /
                         sqrt(colSums(centered[, keep, drop = FALSE] ^ 2)))
    selected <- head(keep[order(-association, colnames(x)[keep])], 24)
    list(center = center[selected], scale = scale[selected])
  }

  design <- function(rows, transforms) {
    columns <- list(covariates[rows, , drop = FALSE])
    for (assay in names(transforms)) {
      tr <- transforms[[assay]]
      z <- sweep(blocks[[assay]][rows, names(tr$center), drop = FALSE], 2, tr$center)
      z[is.na(z)] <- 0
      columns[[assay]] <- sweep(z, 2, tr$scale, "/")
    }
    do.call(cbind, unname(columns))
  }

  learn_design <- function(rows, assays) {
    stopifnot(qr(covariates[rows, , drop = FALSE])$rank == ncol(covariates))
    residual_y <- lm.fit(covariates[rows, , drop = FALSE], response[rows])$residuals
    lapply(blocks[assays], learn_block, train_rows = rows, residual_y = residual_y)
  }

  ridge <- function(x, y, lambda) {
    # Minimize mean squared error + lambda * sum(molecular coefficients squared).
    penalty <- c(rep(0, ncol(covariates)), rep(nrow(x) * lambda, ncol(x) - ncol(covariates)))
    drop(solve(crossprod(x) + diag(penalty, nrow = ncol(x)), crossprod(x, y)))
  }

  pred <- data.frame(donor = donors[test], observed = unname(response[test]))
  models <- list()
  tuning <- list()
  for (method in names(recipes)) {
    assays <- recipes[[method]]
    grid <- if (length(assays)) c(0.01, 0.1, 1, 10, 100) else 0
    cv <- matrix(NA_real_, length(train), length(grid))
    for (k in unique(fold)) {
      fit_rows <- train[fold != k]
      transforms <- learn_design(fit_rows, assays)
      x <- design(fit_rows, transforms)
      v <- design(train[fold == k], transforms)
      for (j in seq_along(grid)) {
        cv[fold == k, j] <- v %*% ridge(x, response[fit_rows], grid[j])
      }
    }
    losses <- colMeans(abs(cv - response[train]))
    best <- which.min(losses)
    tuning[[method]] <- data.frame(model = method, lambda = grid, inner_MAE = losses,
                                   selected = seq_along(grid) == best)
    transforms <- learn_design(train, assays)
    x <- design(train, transforms)
    target <- design(test, transforms)
    coef <- ridge(x, response[train], grid[best])
    pred[[method]] <- drop(target %*% coef)
    # Independent least-squares solver for the fitted ridge coefficients.
    pseudo <- diag(sqrt(c(rep(0, ncol(covariates)),
                          rep(length(train) * grid[best], ncol(x) - ncol(covariates)))), nrow=ncol(x))
    check <- qr.solve(rbind(x, pseudo), c(response[train], rep(0, ncol(x))))
    stopifnot(max(abs(drop(target %*% check) - pred[[method]])) < 1e-8)
    models[[method]] <- list(endpoint = endpoint, transforms = transforms, coefficients = coef,
                             lambda = grid[best], covariate_names = colnames(covariates))
    if (!is.null(gender_coding)) models[[method]]$gender_coding <- gender_coding
  }
  write_tsv(do.call(rbind, tuning), "tuning.tsv")
  saveRDS(models, file.path(out, "development_models.rds"))
  write_tsv(pred, "predictions.tsv")
  errors <- abs(as.matrix(pred[names(recipes)]) - pred$observed)
  metrics <- data.frame(model = names(recipes), n = nrow(pred), MAE = colMeans(errors),
                        RMSE = sqrt(colMeans(errors ^ 2)))
  write_tsv(metrics, "performance.tsv")
  # Paired uncertainty conditional on the fitted models and this holdout.
  set.seed(17042)
  contrasts <- lapply(names(recipes)[1:3], function(baseline) {
    delta <- errors[, baseline] - errors[, names(recipes)[4]]
    samples <- replicate(2000, mean(sample(delta, replace = TRUE)))
    interval <- quantile(samples, c(0.025, 0.975), names = FALSE)
    data.frame(reference = baseline, improvement_MAE = mean(delta), lower = interval[1], upper = interval[2])
  })
  write_tsv(do.call(rbind, contrasts), "paired_differences.tsv")
  writeLines(capture.output(sessionInfo()), file.path(out, "R_session.txt"))
}

read_donors <- function(path) {
  d <- read.delim(path, check.names=FALSE, stringsAsFactors=FALSE, colClasses="character",
                  na.strings=c("", "NA"), comment.char="", quote="")
  if (names(d)[1] != "donor" || anyDuplicated(names(d)) || anyNA(names(d)) ||
      any(!nzchar(names(d))) || anyNA(d$donor) || any(d$donor != trimws(d$donor)) ||
      anyDuplicated(d$donor)) stop("Require a unique donor column and unique feature names in ", path)
  d
}

numeric_field <- function(x, label) {
  value <- suppressWarnings(as.numeric(x))
  if (any(!is.na(x) & !is.finite(value))) stop("Non-numeric or infinite value in ", label)
  value
}

run_comparison <- function(paths, dest) {
  if (file.exists(dest)) stop("Choose a new output directory")
  meta <- read_donors(paths[1])
  required <- c("donor", "observed", "partition", "inner_fold")
  if (!all(required %in% names(meta))) stop("Design needs donor, observed, partition, inner_fold")
  extra <- setdiff(names(meta), required)
  if (any(!startsWith(extra, "cov_"))) stop("Name each numeric baseline covariate with the prefix cov_")
  response <- numeric_field(meta$observed, "observed")
  if (any(!is.finite(response)) || !setequal(meta$partition, c("development", "holdout")))
    stop("Every donor needs a finite continuous outcome and a valid partition")
  train <- which(meta$partition=="development")
  test <- which(meta$partition=="holdout")
  fold <- meta$inner_fold[train]
  if (anyNA(fold) || length(unique(fold))<2 || any(!is.na(meta$inner_fold[test])) || length(test)<2)
    stop("Supply at least two development folds, no holdout fold labels, and at least two holdout donors")
  covariates <- matrix(1, nrow(meta), 1, dimnames=list(meta$donor, "intercept"))
  if (length(extra)) covariates <- cbind(covariates, as.matrix(data.frame(lapply(meta[extra], numeric_field, label="covariates"))))
  if (anyNA(covariates)) stop("Baseline covariates must be complete")
  for (k in unique(fold)) {
    tr <- train[fold!=k]
    if (length(tr)<=ncol(covariates) || qr(covariates[tr,,drop=FALSE])$rank<ncol(covariates) ||
        length(unique(response[tr]))<2) stop("A training fold lacks outcome variation or covariate support")
  }
  read_assay <- function(path) {
    data <- read_donors(path)
    if (ncol(data)<2 || !setequal(data$donor,meta$donor)) stop("Assay donor IDs must match the design: ",path)
    data <- data[match(meta$donor,data$donor),-1,drop=FALSE]
    values <- lapply(names(data), function(n) numeric_field(data[[n]], n))
    mat <- do.call(cbind, values)
    dimnames(mat) <- list(meta$donor, names(data))
    if (any(rowSums(is.finite(mat))==0)) stop("A donor has no observed features in ",path)
    mat
  }
  blocks <- list(A=read_assay(paths[2]), B=read_assay(paths[3]))
  recipes <- list(baseline=character(), assay_A="A", assay_B="B", combined=c("A","B"))
  dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
  staged <- tempfile(pattern=".kilnbio-fit-",tmpdir=dirname(dest))
  if (!dir.create(staged)) stop("Cannot stage output")
  on.exit(unlink(staged, recursive=TRUE), add=TRUE)
  fit_assays(meta$donor,response,covariates,blocks,recipes,train,test,fold,"observed",staged)
  write.table(meta[c("donor","partition","inner_fold")], file.path(staged,"partitions.tsv"),
              sep="\t",quote=FALSE,row.names=FALSE,na="")
  write.table(data.frame(role=c("design","assay_A","assay_B"),file=basename(paths),
                         md5=unname(tools::md5sum(paths))),file.path(staged,"inputs.tsv"),
              sep="\t",quote=FALSE,row.names=FALSE)
  if (file.exists(dest) || !file.rename(staged,dest)) stop("Could not publish new output directory")
}

if (sys.nframe()==0L) {
  args <- commandArgs(TRUE)
  if (identical(args,"--help")) {
    cat("kilnbio compare-assays design.tsv assay_A.tsv assay_B.tsv NEW_OUTPUT\n")
  } else {
    if (length(args)!=4) stop("Supply design.tsv assay_A.tsv assay_B.tsv NEW_OUTPUT")
    run_comparison(args[1:3],args[4])
  }
}
