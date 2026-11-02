# Rscript export_colon.R <new output>. Unsmoothed raw counts, three paired donors.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L,!dir.exists(args[1]))
suppressPackageStartupMessages(library(bsseq))
data(BS.cancer.ex,package='bsseqData')
b <- updateObject(BS.cancer.ex)
# Freeze chr21 as the bounded biological validation region; do not select on effect.
b <- b[as.character(seqnames(b))=='chr21',]
dir.create(args[1],recursive=TRUE)
meta <- as.data.frame(pData(b))
samples <- data.frame(sample=sampleNames(b),biological_id=meta$Pair,condition=meta$Type)
write.csv(samples,file.path(args[1],'samples.csv'),row.names=FALSE)
N <- as.matrix(getCoverage(b,type='Cov'));X <- as.matrix(getCoverage(b,type='M'))
stopifnot(all(X>=0),all(N>=X),all(X==floor(X)),all(N==floor(N)))
manifest <- data.frame(sample=sampleNames(b),counts=paste0(sampleNames(b),'.tsv'))
write.table(data.frame(chr=as.character(seqnames(b)),pos=start(b)),file.path(args[1],'positions.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
for(i in seq_len(ncol(b))) {
 tab <- data.frame(chr=as.character(seqnames(b)),pos=start(b),N=N[,i],X=X[,i])
 write.table(tab,file.path(args[1],manifest$counts[i]),sep='\t',quote=FALSE,row.names=FALSE)
}
write.csv(manifest,file.path(args[1],'methylation.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(args[1],'export_session.txt'))
cat(nrow(b),'chr21 CpGs;',ncol(b),'libraries; three paired patients\n')
