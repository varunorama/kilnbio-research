# Input matrices contain deposited continuous abundances (not raw RNA counts).
args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 1 && nzchar(args[2])) .libPaths(c(args[2], .libPaths()))
if (!requireNamespace("limma", quietly=TRUE)) stop("Install Bioconductor limma; see README.md")
if (args[1] == "--check") {
    cat(R.version.string, "| limma", as.character(packageVersion("limma")), "\n")
    quit(status=0L)
}
setwd(args[1])
readmat <- function(f) as.matrix(read.delim(f, row.names=1, check.names=FALSE, na.strings=c("NA","NaN","")))
write_tsv <- function(x,f) write.table(x,f,sep="\t",row.names=FALSE,quote=FALSE,na="NA")
roster <- read.delim("roster.tsv", check.names=FALSE)
label <- setNames(roster$fgfr3,roster$specimen)
rna <- readmat("rna_matrix.tsv"); protein <- readmat("protein_matrix.tsv")
stopifnot(is.numeric(rna),is.numeric(protein),!any(is.infinite(rna)),!any(is.infinite(protein)),
          all(rna>=0,na.rm=TRUE),identical(rownames(rna),rownames(protein)),
          identical(colnames(rna),roster$specimen),identical(colnames(protein),roster$specimen[roster$add_ok==1]))
rna <- log2(1+rna)
matrices <- list(A=rna,B=rna[,colnames(protein),drop=FALSE],P=protein)
pass <- function(x) {
    g <- label[colnames(x)]
    stopifnot(!anyNA(g),all(table(factor(g,levels=0:1))>=3))
    v <- apply(x,1,sd,na.rm=TRUE)
    rowMeans(is.finite(x[,g==0,drop=FALSE]))>=0.8 & rowMeans(is.finite(x[,g==1,drop=FALSE]))>=0.8 &
        is.finite(v) & v>0
}
mask <- lapply(matrices,pass)
eligible <- Reduce("&",mask)
write_tsv(data.frame(gene=rownames(rna),as.data.frame(mask),eligible=eligible),"eligibility.tsv")
if (sum(eligible)<15) stop("Fewer than 15 eligible genes")
fit_t <- function(x) {
    g <- unname(label[colnames(x)])
    fit <- limma::eBayes(limma::lmFit(x[eligible,,drop=FALSE],cbind(intercept=1,positive=g)))
    t <- fit$t[,"positive"]
    stopifnot(all(is.finite(t)))
    t
}
tstats <- lapply(matrices,fit_t)
signed <- function(t) sign(t)*rank(abs(t),ties.method="average")/(length(t)+1)
s <- lapply(tstats,signed)
s$C <- (s$B+s$P)/2
result <- do.call(rbind,lapply(names(s),function(arm){
    score <- s[[arm]]; ord <- order(-abs(score),names(score),method="radix")
    data.frame(arm=arm,rank=seq_along(ord),gene=names(score)[ord],direction=sign(score[ord]),
               score=unname(score[ord]),t=if (arm=="C") NA_real_ else unname(tstats[[arm]][ord]))
}))
write_tsv(result,"rankings.tsv")
write_tsv(result[result$rank<=15,],"shortlists.tsv")
capture.output(sessionInfo(),file="R-session.txt")
