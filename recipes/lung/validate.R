args=commandArgs(trailingOnly=TRUE)
.libPaths(c(args[2],.libPaths()));library(limma)
setwd(file.path(args[1],"analysis"))
x=as.matrix(read.delim("external_expression.tsv",check.names=FALSE,row.names=1))
patients=read.delim("external_patients.tsv",row.names=1)
g=patients[colnames(x),"egfr"]; nomin=read.delim("frozen.tsv")
arms=unique(nomin$arm)
fit=function(i){
    m=eBayes(lmFit(x[,i],cbind(1,g[i])))
    data.frame(gene=rownames(x),effect=m$coefficients[,2],p=m$p.value[,2],q=p.adjust(m$p.value[,2],"BH"))
}
yield=function(df){
    i=match(nomin$gene,df$gene)
    ok=(!is.na(i)) & df$q[i]<=0.05 & (df$effect[i]*nomin$direction>0)
    vapply(arms,function(a)sum(ok[nomin$arm==a],na.rm=TRUE),numeric(1))
}
obs=fit(seq_along(g))
write.table(obs,"external_tests.tsv",sep="\t",row.names=FALSE,quote=FALSE)
set.seed(621013)
boot=t(replicate(1000,{
    i=c(sample(which(g==0),sum(g==0),TRUE),sample(which(g==1),sum(g==1),TRUE));yield(fit(i))
}))
write.table(boot,"bootstrap.tsv",sep="\t",row.names=FALSE,quote=FALSE)
print(yield(obs))

