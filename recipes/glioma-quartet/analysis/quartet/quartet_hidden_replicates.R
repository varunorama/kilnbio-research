# Evaluate a frozen control-feature projection on previously unseen labs and replicates.
# Rscript quartet_hidden_replicates.R VERIFIED_QUARTET_CACHE NEW_OUTPUT
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE))))
source(file.path(here, "adjust.R"))
cache <- args[1]
out <- args[2]
files <- file.path(cache, paste0("quartet_", c("metadata","RNA","protein","metabolite"), ".csv"))
stopifnot(identical(unname(tools::md5sum(files)),c("dca5a516abbd81203adff88d0df5057f","ac5047aa080964bb4e08e08565475955","e5737b2bc28c5bb09993617371f7e67b","7ac78e34868f73d8197668f7d4974c45")))
stopifnot(dir.create(out,recursive=TRUE,showWarnings=FALSE))
meta <- read.csv(files[1],stringsAsFactors=FALSE)
materials <- c("D5","F7","M8")
rmse <- function(x) sqrt(mean(x^2))
write_tsv <- function(x, f) write.table(x,file.path(out,f),sep="\t",row.names=FALSE,quote=FALSE)
detail <- list()
folds <- list()
for (assay in c("RNA","protein","metabolite")) {
  m <- meta[tolower(meta$datatype)==tolower(assay),]
  stopifnot(!anyDuplicated(m$library), all(table(m$batch,m$sample,m$rep)==1L),
            all(tapply(m$lab,m$batch,function(a) length(unique(a)))==1L))
  tab <- read.csv(file.path(cache,paste0("quartet_",assay,".csv")),row.names=1,check.names=FALSE)
  x <- t(as.matrix(tab[,m$library,drop=FALSE]))
  stopifnot(is.numeric(x),all(is.finite(x)))
  bridge <- which(m$sample=="D6" & m$rep==1)
  ratio <- x - x[bridge[match(m$batch,m$batch[bridge])],,drop=FALSE]
  for (lab in sort(unique(m$lab))) {
    train <- m$lab!=lab & m$rep==1 & m$sample %in% materials
    b <- factor(m$sample[train])
    train_x <- x[train,,drop=FALSE]
    idx <- split(seq_len(nrow(train_x)),b)
    means <- t(vapply(idx,function(i) colMeans(train_x[i,,drop=FALSE]),numeric(ncol(x))))
    residual <- train_x - means[as.character(b),,drop=FALSE]
    within <- colMeans(residual^2)
    between <- colMeans(sweep(means,2,colMeans(means))^2)
    eligible <- which(within>0 & is.finite(within))
    n <- max(3L,min(20L,floor(ncol(x)/5)))
    stopifnot(length(eligible)>=n, n < ncol(x))
    controls <- colnames(x)[eligible[order(between[eligible]/within[eligible],colnames(x)[eligible])[seq_len(n)]]]
    scored <- setdiff(colnames(x),controls)
    fit <- adjust_assay(train_x,b,controls=controls)
    confounded <- adjust_assay(train_x,b,batch=b)
    stopifnot(is.null(confounded$adjusted))
    folds[[length(folds)+1L]] <- data.frame(assay,lab,control=controls,fit_status=fit$status,
                                           confounded_status=confounded$status,train_n=sum(train),score_features=length(scored))
    methods <- list(published=x,concurrent_D6=ratio)
    if (!is.null(fit$adjusted)) {
      stopifnot(isTRUE(all.equal(project_control_adjustment(fit,train_x),fit$adjusted,tolerance=1e-9)))
      methods$frozen_control_PC <- project_control_adjustment(fit,x)
    }
    hidden <- m$rep %in% c(2,3) & m$sample %in% materials
    other <- hidden & m$lab!=lab
    # Balance laboratory contributions irrespective of the number of batches per lab.
    templates <- function(mat) {
      vapply(materials,function(s) {
        ids <- which(other & m$sample==s)
        per_lab <- lapply(split(ids,m$lab[ids]),function(i) colMeans(mat[i,scored,drop=FALSE]))
        Reduce("+",per_lab)/length(per_lab)
      },numeric(length(scored)))
    }
    raw_template <- templates(x)
    for (method in names(methods)) {
      matrix <- methods[[method]]
      ref <- templates(matrix)
      for (batch in unique(m$batch[m$lab==lab])) for (replicate in c(2L,3L)) {
        ids <- which(m$batch==batch & m$rep==replicate & hidden)
        obs <- matrix[ids[match(materials,m$sample[ids])],scored,drop=FALSE]
        distances <- vapply(seq_along(materials),function(k) sqrt(colMeans((ref-obs[k,])^2)),numeric(3))
        pairs <- combn(1:3,2)
        contrasts <- vapply(seq_len(ncol(pairs)),function(k) {
          i <- pairs[1,k]; j <- pairs[2,k]
          a <- obs[i,]-obs[j,]; target <- raw_template[,i]-raw_template[,j]
          c(rmse=rmse(a-target),rho=cor(a,target,method="spearman"))
        },numeric(2))
        detail[[length(detail)+1L]] <- data.frame(assay,lab,batch,replicate,method,
            technical_RMSE=mean(diag(distances)),biological_contrast_RMSE=mean(contrasts[1,]),
            contrast_rho=mean(contrasts[2,]),identity_accuracy=mean(apply(distances,2,which.min)==1:3))
      }
    }
  }
}
d <- do.call(rbind,detail)
write_tsv(d,"hidden_replicates.tsv")
write_tsv(do.call(rbind,folds),"fold_audit.tsv")
bylab <- aggregate(d[6:9], d[c("assay","method","lab")], mean)
write_tsv(bylab,"held_out_labs.tsv")
summary <- aggregate(bylab[4:7], bylab[c("assay","method")], mean)
write_tsv(summary,"summary.tsv")
writeLines(capture.output(sessionInfo()),file.path(out,"R_session.txt"))
