# Compare methods on a fixed matrix, contrast, biological unit and feature family.
compare_script_path <- if (sys.nframe()>0L && !is.null(sys.frame(1)$ofile)) sys.frame(1)$ofile else sub("^--file=", "", grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
agreement_script <- file.path(dirname(normalizePath(compare_script_path)), "method_agreement.py")
compare_differential <- function(kind, matrix_path, design_path, reference, tested, out) {
  if (!kind %in% c("rna-counts","protein-log2") || reference==tested) stop("Choose a supported data type and two distinct groups")
  extra <- Sys.getenv("KILNBIO_R_LIBRARY")
  if (nzchar(extra)) .libPaths(c(extra,.libPaths()))
  pkgs <- if (kind=="rna-counts") c("limma","edgeR","DESeq2") else "limma"
  for (p in pkgs) if (!requireNamespace(p,quietly=TRUE)) stop("Missing R package: ",p)
  if (file.exists(out)) stop("Choose a new output directory")
  load_tsv <- function(p) read.delim(p,check.names=FALSE,stringsAsFactors=FALSE,colClasses="character",na.strings=c("","NA"),quote="",comment.char="")
  meta <- load_tsv(design_path); tab <- load_tsv(matrix_path)
  if (anyDuplicated(names(meta)) || !setequal(names(meta),c("sample","unit","group")) || anyNA(meta) || anyDuplicated(meta$sample) ||
      !setequal(meta$group,c(reference,tested)) || any(vapply(meta,function(v) any(trimws(v)!=v),logical(1)))) stop("Design requires exactly sample, unit, group and two specified levels")
  if (ncol(tab)<2 || names(tab)[1]!="feature" || anyNA(tab$feature) || anyDuplicated(tab$feature) ||
      anyDuplicated(names(tab)) || !setequal(names(tab)[-1],meta$sample)) stop("Require unique feature IDs and exact sample columns")
  raw <- as.matrix(tab[,meta$sample,drop=FALSE]); x <- matrix(suppressWarnings(as.numeric(raw)),nrow=nrow(raw),dimnames=list(tab$feature,meta$sample))
  if (any(is.na(x)!=is.na(raw)) || any(is.infinite(x))) stop("Malformed or nonfinite numeric measurement")
  meta$unit <- factor(meta$unit);meta$group <- factor(meta$group,levels=c(reference,tested))
  pairs <- table(meta$unit,meta$group)
  paired <- any(rowSums(pairs)>1)
  if (any(pairs>1L) || (paired && !all(pairs==1L)) || min(table(meta$group))<3L) stop("Need three independent units per arm, or at least three complete one-to-one pairs; aggregate technical repeats first")
  form <- if (paired) ~unit+group else ~group
  rownames(meta) <- meta$sample
  design <- model.matrix(form,meta)
  if (qr(design)$rank!=ncol(design) || nrow(design)<=ncol(design)) stop("Nonestimable design")
  coef <- ncol(design)
  if (kind=="rna-counts") {
    if (anyNA(x) || any(x<0 | x!=trunc(x) | x>.Machine$integer.max) || any(colSums(x)<=0)) stop("RNA requires nonnegative unnormalized integer fragment/read counts")
    keep <- rowSums(x>=10)>=min(table(meta$group))
    reason <- ifelse(keep,"retained","low_count")
  } else {
    finite <- rowSums(is.finite(x))==ncol(x)
    vary <- apply(x,1,function(z) length(unique(z[is.finite(z)]))>1)
    keep <- finite & vary
    reason <- ifelse(!finite,"incomplete",ifelse(!vary,"constant","retained"))
  }
  if (sum(keep)<2L) stop("Fewer than two eligible features")
  y <- x[keep,,drop=FALSE]
  methods <- list()
  arm <- function(effect,p) {
    p <- as.numeric(p); effect <- as.numeric(effect)
    if (length(p)!=nrow(y) || length(effect)!=nrow(y) || any(p<0|p>1,na.rm=TRUE)) stop("Unexpected engine output")
    q <- rep(NA_real_,length(p)); ix <- is.finite(p)
    if (any(ix)) q[ix] <- p.adjust(p[ix],"BH",n=nrow(y))
    data.frame(effect=effect,p=p,q=q)
  }
  if (kind=="rna-counts") {
    storage.mode(y) <- "integer"
    ds <- DESeq2::DESeqDataSetFromMatrix(y,meta,form)
    ds <- DESeq2::DESeq(ds,sfType="poscounts",minReplicatesForReplace=Inf,quiet=TRUE)
    dr <- as.data.frame(DESeq2::results(ds,contrast=c("group",tested,reference),independentFiltering=FALSE))
    stopifnot(identical(rownames(dr),rownames(y)))
    methods$DESeq2 <- arm(dr$log2FoldChange,dr$pvalue)
    dge <- edgeR::calcNormFactors(edgeR::DGEList(counts=y))
    fit <- edgeR::glmQLFit(edgeR::estimateDisp(dge,design),design)
    et <- edgeR::glmQLFTest(fit,coef=coef)$table
    stopifnot(identical(rownames(et),rownames(y)))
    methods$edgeR_QL <- arm(et$logFC,et$PValue)
    lm <- limma::eBayes(limma::lmFit(limma::voom(dge,design,plot=FALSE),design))
    methods$limma_voom <- arm(lm$coefficients[,coef],lm$p.value[,coef])
  } else {
    lm <- limma::eBayes(limma::lmFit(y,design),trend=TRUE)
    methods$limma_trend <- arm(lm$coefficients[,coef],lm$p.value[,coef])
    hi <- which(meta$group==tested);lo <- which(meta$group==reference)
    if (paired) lo <- lo[match(meta$unit[hi],meta$unit[lo])]
    A <- y[,hi,drop=FALSE]; B <- y[,lo,drop=FALSE]
    ptest <- function(fun,a,b) tryCatch(suppressWarnings(fun(a,b)$p.value),error=function(e) NA_real_)
    tt <- function(a,b)t.test(a,b,paired=paired)
    ww <- function(a,b)wilcox.test(a,b,paired=paired,exact=FALSE,correct=TRUE)
    vals <- t(vapply(seq_len(nrow(y)),function(i) {
      a <- A[i,];b <- B[i,]
      c(mean(a)-mean(b),ptest(tt,a,b),if(paired)median(a-b) else median(a)-median(b),ptest(ww,a,b))
    },numeric(4)))
    methods[[if(paired)"paired_t" else "Welch_t"]] <- arm(vals[,1],vals[,2])
    methods[[if(paired)"signed_rank" else "rank_sum"]] <- arm(vals[,3],vals[,4])
  }
  calls <- lapply(methods,function(z)!is.na(z$q)&z$q<=.05)
  summary <- do.call(rbind,lapply(combn(names(methods),2,simplify=FALSE),function(n){
    a <- methods[[n[1]]];b <- methods[[n[2]]];ka <- calls[[n[1]]];kb <- calls[[n[2]]]
    ok <- is.finite(a$effect)&is.finite(b$effect);u <- sum(ka|kb)
    data.frame(method_A=n[1],method_B=n[2],family=nrow(y),both_testable=sum(is.finite(a$p)&is.finite(b$p)),
      effect_pairs=sum(ok),effect_spearman=if(sum(ok)>=2)suppressWarnings(cor(a$effect[ok],b$effect[ok],method="spearman")) else NA_real_,
      A_calls=sum(ka),B_calls=sum(kb),shared_calls=sum(ka&kb),opposite_shared=sum(ka&kb&sign(a$effect)!=sign(b$effect),na.rm=TRUE),
      union_calls=u,jaccard=if(u)sum(ka&kb)/u else NA_real_)
  }))
  dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
  tmp <- tempfile(".differential-compare-",tmpdir=dirname(out));if (!dir.create(tmp)) stop("Cannot stage outputs")
  on.exit(unlink(tmp,recursive=TRUE),add=TRUE)
  write <- function(x,name)write.table(x,file.path(tmp,name),sep="\t",row.names=FALSE,quote=FALSE,na="NA")
  res <- data.frame(feature=rownames(y))
  for(n in names(methods)) for(k in names(methods[[n]])) res[[paste(n,k,sep="_")]] <- methods[[n]][[k]]
  write(res,"features.tsv");write(summary,"comparisons.tsv")
  write(data.frame(method=names(methods),study=unname(tools::md5sum(design_path)),contrast=paste(tested,"vs",reference),
                   feature_family=paste(unname(tools::md5sum(matrix_path)),sum(keep),sep=":"),feature_unit=if(kind=="rna-counts")"counted_feature" else "protein",
                   assembly="as_input",channel=kind),"method_manifest.tsv")
  write(data.frame(feature=rownames(x),eligibility=reason),"eligibility.tsv");write(meta,"design.tsv")
  write(data.frame(input=c("matrix","design"),md5=unname(tools::md5sum(c(matrix_path,design_path)))),"inputs.tsv")
  write(data.frame(type=kind,reference=reference,tested=tested,paired=paired,units=nrow(pairs),input_features=nrow(x),retained=nrow(y),alpha=.05),"analysis.tsv")
  writeLines(capture.output(sessionInfo()),file.path(tmp,"session.txt"))
  python <- Sys.getenv("KILNBIO_PYTHON",if(nzchar(Sys.which("python3"))) "python3" else "python")
  status <- system2(python,shQuote(c(agreement_script,file.path(tmp,"features.tsv"),file.path(tmp,"method_manifest.tsv"),"--out",file.path(tmp,"method_agreement.tsv"))))
  if(status!=0L || !file.exists(file.path(tmp,"method_agreement.tsv"))) stop("Per-feature agreement failed")
  if(file.exists(out)||!file.rename(tmp,out)) stop("Could not publish outputs")
}
if(sys.nframe()==0L){
  args <- commandArgs(TRUE)
  if(identical(args,"--help")) cat("kilnbio differential-compare {rna-counts|protein-log2} matrix.tsv design.tsv REFERENCE TEST NEW_OUTPUT\n")
  else if(length(args)!=6L) stop("Expected type, matrix, design, reference group, test group and new output path")
  else do.call(compare_differential,as.list(args))
}
