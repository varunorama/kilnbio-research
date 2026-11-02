args=commandArgs(trailingOnly=TRUE)
.libPaths(c(args[2],.libPaths()))
library(limma)
setwd(file.path(args[1],"analysis"))
readmat=function(f) as.matrix(read.delim(f,check.names=FALSE,row.names=1))
rna=readmat("RNA.tsv"); protein=readmat("protein.tsv")
roster=read.delim("roster.tsv",check.names=FALSE)
group=setNames(roster$egfr,roster$specimen)
obs=read.delim("source_manifest.tsv")
mixed=read.delim("planted_manifest.tsv")
withheld=read.delim("withheld_manifest.tsv")
tstat=function(x,ids){
    design=cbind(intercept=1,positive=as.numeric(group[ids]))
    result=eBayes(lmFit(x,design))$t[,2]
    stopifnot(all(is.finite(result)))
    result
}
signed_score=function(x) sign(x)*rank(abs(x),ties.method="average")/(length(x)+1)
single=function(x) signed_score(tstat(x,colnames(x)))
direct=intersect(colnames(rna),colnames(protein))
A=single(rna);B=single(rna[,direct]);P=single(protein[,direct]);C_standard=(B+P)/2
library(jsonlite)
get_pairs=function(file){
    d=fromJSON(file)$unique_pairs
    data.frame(rna=sub("^RNA_","",d$base_measurement),prot=sub("^protein_","",d$additional_measurement))
}
joint=function(p){
    x=p$rna; y=p$prot
    return((signed_score(tstat(rna[,x],x))+signed_score(tstat(protein[,y],x)))/2)
}
C=joint(get_pairs("clean_pairing.json"))
stopifnot(max(abs(C-C_standard))<1e-10)
pa=merge(mixed[mixed$assay=="RNA",c("donor","specimen")],mixed[mixed$assay=="protein",c("donor","specimen")],by="donor",suffixes=c("_rna","_protein"))
D=joint(setNames(pa[,c("specimen_rna","specimen_protein")],c("rna","prot")))
W=joint(get_pairs("withheld_pairing.json"))
S=list(A_RNA_all=A,B_RNA_paired=B,P_protein=P,C_joint=C,standard_exact=C_standard,planted_donor_join=D,planted_withhold=W)
rankings=do.call(rbind,lapply(names(S),function(method){
    x=S[[method]];o=order(-abs(x),names(x),method="radix")
    data.frame(arm=method,rank=seq_along(o),gene=names(x)[o],score=unname(x[o]),direction=as.integer(sign(x[o])))
}))
stopifnot(identical(subset(rankings,arm=="C_joint")$gene,subset(rankings,arm=="standard_exact")$gene),all(rankings$direction[rankings$rank<=15]!=0))
write.table(rankings,"rankings.tsv",sep="\t",row.names=FALSE,quote=FALSE)
write.table(subset(rankings,rank<=15),"frozen.tsv",sep="\t",row.names=FALSE,quote=FALSE)
capture.output(sessionInfo(),file="R-session.txt")

