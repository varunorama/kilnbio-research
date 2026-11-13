# Fit the same declared association or paired contrast in multiple assay matrices.
# Public limma / edgeR workflow; no adaptive test selection or trait permutation.
run_layers <- function(config_path, out) {
  extra <- Sys.getenv("KILNBIO_R_LIBRARY")
  if (nzchar(extra)) .libPaths(c(extra, .libPaths()))
  for (pkg in c("jsonlite", "limma")) if (!requireNamespace(pkg, quietly=TRUE)) stop("Missing R package ", pkg)
  if (file.exists(out)) stop("Choose a new output directory")
  cfg <- jsonlite::fromJSON(config_path, simplifyVector=FALSE)
  base <- dirname(normalizePath(config_path))
  path <- function(f) if (grepl("^(/|[A-Za-z]:[/\\\\])", f)) f else file.path(base, f)
  read <- function(f) read.delim(path(f), check.names=FALSE, colClasses="character", na.strings=c("", "NA"))
  numeric_matrix <- function(d) {
    x <- as.matrix(d)
    y <- suppressWarnings(matrix(as.numeric(x), nrow=nrow(x), dimnames=dimnames(x)))
    if (any(is.na(y) != is.na(x)) || any(is.infinite(y))) stop("Matrices require finite values or NA")
    y
  }
  donors <- read(cfg$donors)
  if (names(donors)[1] != "donor" || anyNA(donors$donor) || anyDuplicated(donors$donor)) stop("Need unique donor IDs")
  if (!cfg$estimand %in% c("association", "paired") || length(cfg$layers)<2) stop("Need association or paired and >=2 layers")
  if (!is.null(cfg$permutations) && cfg$permutations != 0) stop("Joint trait permutations with covariates were removed; omit permutations or set 0")
  covs <- unlist(cfg$covariates)
  if (anyDuplicated(covs) || (!is.null(cfg$trait) && cfg$trait %in% covs)) stop("Duplicate model term")
  cohort_policy <- if (is.null(cfg$cohort_policy)) "matched" else cfg$cohort_policy
  if (!cohort_policy %in% c("matched","available")) stop("cohort_policy must be matched or available")
  first_donors <- NULL
  if (cfg$estimand == "paired" && length(covs)) stop("Paired mode uses donor blocking; additional covariates need a study-specific model")
  names_layers <- vapply(cfg$layers, function(x) x$name, "")
  if (anyDuplicated(names_layers) || any(!grepl("^[A-Za-z][A-Za-z0-9_]*$", names_layers))) stop("Unique simple layer names required")
  dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
  stage <- tempfile(".layers-",tmpdir=dirname(out)); dir.create(stage)
  on.exit(unlink(stage,recursive=TRUE), add=TRUE)
  write <- function(d,n) write.table(d,file.path(stage,n),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
  manifest <- list(); audits <- list(); inputs <- c(config_path,path(cfg$donors))
  for (layer in cfg$layers) {
    tab <- read(layer$matrix); sp <- read(layer$specimens)
    inputs <- c(inputs,path(layer$matrix),path(layer$specimens))
    if (names(tab)[1]!="feature" || anyNA(tab$feature) || anyDuplicated(tab$feature) || anyDuplicated(names(tab))) stop("Unique feature matrix required")
    if (!all(c("specimen","donor") %in% names(sp)) || anyNA(sp[c("specimen","donor")]) || anyDuplicated(sp$specimen)) stop("Unique specimen mapping required")
    if (!setequal(names(tab)[-1],sp$specimen) || !all(sp$donor %in% donors$donor)) stop("Every matrix column needs an explicit specimen/donor map")
    if (is.null(first_donors)) first_donors <- unique(sp$donor)
    if (cohort_policy == "matched" && !setequal(sp$donor,first_donors)) stop("Matched cohort requires identical donor sets across layers")
    raw <- numeric_matrix(tab[,sp$specimen,drop=FALSE]); rownames(raw) <- tab$feature
    design_data <- donors[match(sp$donor, donors$donor),,drop=FALSE]
    if (cfg$estimand == "association") {
      if (anyDuplicated(sp$donor)) stop("Association mode requires one independent specimen per donor per layer")
      if (!all(c(cfg$trait,covs) %in% names(design_data))) stop("Missing trait or covariate")
      if (any(!grepl("^[A-Za-z][A-Za-z0-9_]*$",c(cfg$trait,covs)))) stop("Simple model column names required")
      for (name in c(cfg$trait,covs)) design_data[[name]] <- as.numeric(design_data[[name]])
      if (any(!is.finite(as.matrix(design_data[c(cfg$trait,covs)])))) stop("Supply complete finite numeric trait/covariates; encode categorical covariates explicitly")
      design <- model.matrix(reformulate(c(covs,cfg$trait)), design_data)
      coef <- match(cfg$trait,colnames(design))
    } else {
      cc <- cfg$condition
      if (is.null(cc) || !cc$column %in% names(sp) || identical(cc$test,cc$reference)) stop("Declare paired condition column, reference and test")
      if (!all(sp[[cc$column]] %in% c(cc$reference,cc$test)) || anyDuplicated(paste(sp$donor,sp[[cc$column]],sep="\037"))) stop("One specimen per donor and condition required")
      if (!all(table(sp$donor)==2)) stop("Paired mode requires complete pairs")
      design_data$donor <- factor(sp$donor)
      design_data$condition <- factor(sp[[cc$column]],levels=c(cc$reference,cc$test))
      design <- model.matrix(~ donor + condition,design_data); coef <- ncol(design)
    }
    if (qr(design)$rank != ncol(design) || nrow(design)-ncol(design)<2) stop("Confounded design or insufficient residual degrees of freedom")
    if (!layer$orientation %in% c(-1,1)) stop("Each layer needs orientation 1 or -1")
    scale <- layer$scale
    if (scale == "counts") {
      if (!requireNamespace("edgeR",quietly=TRUE)) stop("counts requires edgeR")
      if (anyNA(raw) || any(raw<0 | raw!=floor(raw))) stop("Counts must be complete nonnegative integers")
      y <- edgeR::DGEList(raw)
      keep <- edgeR::filterByExpr(y,design=design)
      if (sum(keep)<2) stop("Too few expressed count features")
      y <- edgeR::calcNormFactors(y[keep,,keep.lib.sizes=FALSE])
      v <- limma::voom(y,design,plot=FALSE)
      fit <- limma::eBayes(limma::lmFit(v,design))
    } else {
      if (!scale %in% c("beta","mvalue","log2")) stop("Supported scales: counts, beta, mvalue, log2")
      if (scale == "beta") {
        if (any(raw<0 | raw>1,na.rm=TRUE)) stop("Beta values outside [0,1]")
        raw <- log2(pmin(pmax(raw,1e-6),1-1e-6)/(1-pmin(pmax(raw,1e-6),1-1e-6)))
      }
      keep <- rowSums(!is.finite(raw))==0 & apply(raw,1,sd,na.rm=TRUE)>0
      if (sum(keep)<2) stop("Too few complete variable features")
      fit <- limma::eBayes(limma::lmFit(raw[keep,,drop=FALSE],design),trend=TRUE)
    }
    result <- data.frame(feature=rownames(raw),effect=NA_real_,p=NA_real_)
    result$effect[keep] <- fit$coefficients[,coef]
    result$p[keep] <- pmax(fit$p.value[,coef],1e-300)
    write(result,paste0(layer$name,".tsv"))
    fmap <- layer$feature_map
    if (!identical(fmap,"identity")) {
      map <- read(fmap)
      if (!all(c("feature","gene") %in% names(map))) stop("Map requires feature and gene")
      inputs <- c(inputs,path(fmap)); fmap <- paste0(layer$name,"_map.tsv")
      write(map[,c("feature","gene")],fmap)
    }
    manifest[[layer$name]] <- data.frame(layer=layer$name,results=paste0(layer$name,".tsv"),feature_map=fmap,orientation=layer$orientation,weight=1)
    audits[[layer$name]] <- data.frame(layer=layer$name,estimand=cfg$estimand,cohort_policy=cohort_policy,scale=scale,donors=length(unique(sp$donor)),features=nrow(raw),tested=sum(keep),residual_df=nrow(design)-ncol(design))
    write(data.frame(specimen=sp$specimen,donor=sp$donor,design,check.names=FALSE),paste0(layer$name,"_design.tsv"))
  }
  write(do.call(rbind,manifest),"manifest.tsv"); write(do.call(rbind,audits),"layers_design.tsv")
  write(data.frame(input=basename(inputs),md5=unname(tools::md5sum(inputs))),"inputs.tsv")
  writeLines(capture.output(sessionInfo()),file.path(stage,"session.txt"))
  if (file.exists(out) || !file.rename(stage,out)) stop("Cannot publish outputs")
}
if (sys.nframe()==0L) {
  a <- commandArgs(TRUE)
  if (identical(a,"--help")) cat("kilnbio multiomic-layers config.json NEW_OUTPUT\n") else if (length(a)!=2L) stop("Expected config.json NEW_OUTPUT") else run_layers(a[1],a[2])
}
