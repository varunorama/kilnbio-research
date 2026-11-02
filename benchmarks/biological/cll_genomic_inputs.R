# Source-specific preparation for the public MOFA CLL worked example.
args <- commandArgs(TRUE)
if(length(args)!=2L) stop("Rscript cll_genomic_inputs.R /path/CLL_data.RData NEW_DIR")
if(file.exists(args[2])) stop("Choose a new directory")
e <- new.env(parent=emptyenv());load(args[1],envir=e)
dat <- e$CLL_data
if(!all(c("Mutations","mRNA","Methylation")%in%names(dat)) || !"TP53"%in%rownames(dat$Mutations))stop("Unexpected source")
z<-dat$Mutations["TP53",];stopifnot(!anyDuplicated(names(z)),all(is.na(z)|z%in%c(0,1)))
dir.create(args[2],recursive=TRUE)
write <- function(x,n)write.table(x,file.path(args[2],n),sep="\t",quote=FALSE,row.names=FALSE,na="")
write(data.frame(event="TP53_gene_alteration",origin="somatic",donor=names(z),specimen="",state=ifelse(is.na(z),"unknown",ifelse(z==1,"altered","reference"))),"calls.tsv")
for (n in c("mRNA","Methylation")){
 x<-dat[[n]];stopifnot(!anyDuplicated(colnames(x)),!anyDuplicated(rownames(x)),setequal(colnames(x),names(z)))
 k<-colSums(is.finite(x))>0
 write(data.frame(sample=colnames(x)[k],donor=colnames(x)[k],specimen=""),paste0(n,"-samples.tsv"))
 write(data.frame(feature=rownames(x),x[,k,drop=FALSE],check.names=FALSE),paste0(n,"-traits.tsv"))
}
