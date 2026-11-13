# Direct R integration checks: real coefficients, paired blocking, counts and design rejection.
a <- commandArgs(TRUE)
if (length(a)!=2) stop("Expected source root and NEW output directory")
source(file.path(a[1],"downstream","multiomic_layers.R"))
if (file.exists(a[2])) stop("Choose a new directory")
dir.create(a[2],recursive=TRUE); root <- normalizePath(a[2])
write <- function(x,n) write.table(x,file.path(root,n),sep="\t",quote=FALSE,row.names=FALSE)
set.seed(391)
n <- 20; ids <- paste0("d",seq_len(n)); trait <- rep(c(0,1),each=n/2)
write(data.frame(donor=ids,trait=trait,age=rep(1:10,2)),"donors.tsv")
write(data.frame(specimen=ids,donor=ids),"specimens.tsv")
x <- matrix(rnorm(80*n),80,n)+seq_len(80)*.1
x[1,] <- x[1,]+5*trait
colnames(x)<-ids
write(data.frame(feature=paste0("G",1:80),x,check.names=FALSE),"protein.tsv")
counts <- matrix(rnbinom(80*n,mu=50,size=5),80,n); counts[1,trait==1] <- counts[1,trait==1]+150
colnames(counts)<-ids
write(data.frame(feature=paste0("G",1:80),counts,check.names=FALSE),"rna.tsv")
cfg <- list(donors="donors.tsv",estimand="association",trait="trait",covariates=list("age"),layers=list(
  list(name="protein",matrix="protein.tsv",specimens="specimens.tsv",scale="log2",feature_map="identity",orientation=1),
  list(name="rna",matrix="rna.tsv",specimens="specimens.tsv",scale="counts",feature_map="identity",orientation=1)))
config<-file.path(root,"config.json")
jsonlite::write_json(cfg,config,auto_unbox=TRUE,pretty=TRUE)
run_layers(config,file.path(root,"association"))
observed <- read.delim(file.path(root,"association","protein.tsv"))
expected <- coef(lm(x[1,]~rep(1:10,2)+trait))[3]
stopifnot(abs(observed$effect[1]-expected)<1e-10,observed$p[1]<1e-6)
bad <- cfg; bad$covariates<-list("trait")
jsonlite::write_json(bad,file.path(root,"bad.json"),auto_unbox=TRUE)
stopifnot(inherits(try(run_layers(file.path(root,"bad.json"),file.path(root,"bad")),silent=TRUE),"try-error"))
bad <- cfg; bad$permutations<-10
jsonlite::write_json(bad,file.path(root,"bad.json"),auto_unbox=TRUE)
stopifnot(inherits(try(run_layers(file.path(root,"bad.json"),file.path(root,"bad")),silent=TRUE),"try-error"))
# Paired data have large donor baselines; the fitted difference must remain 3.
paired_ids<-paste0("s",1:20); donor<-rep(paste0("p",1:10),2); condition<-rep(c("pre","post"),each=10)
write(data.frame(donor=unique(donor)),"paired_donors.tsv")
write(data.frame(specimen=paired_ids,donor=donor,visit=condition),"paired_specimens.tsv")
y<-matrix(rnorm(40*20,sd=.3),40,20)+rep(rep(1:10,2),each=40)
y[1,]<-rep(1:10,2)+rep(c(0,3),each=10)+rnorm(20,sd=.01);colnames(y)<-paired_ids
write(data.frame(feature=paste0("G",1:40),y,check.names=FALSE),"paired.tsv")
paired<-list(donors="paired_donors.tsv",estimand="paired",condition=list(column="visit",reference="pre",test="post"),layers=list(
  list(name="a",matrix="paired.tsv",specimens="paired_specimens.tsv",scale="log2",feature_map="identity",orientation=1),
  list(name="b",matrix="paired.tsv",specimens="paired_specimens.tsv",scale="log2",feature_map="identity",orientation=1)))
jsonlite::write_json(paired,file.path(root,"paired.json"),auto_unbox=TRUE)
run_layers(file.path(root,"paired.json"),file.path(root,"paired"))
observed<-read.delim(file.path(root,"paired","a.tsv"))
stopifnot(abs(observed$effect[1]-mean(y[1,11:20]-y[1,1:10]))<1e-10)
cat("PASS: limma coefficient, voom execution, invalid designs, and paired donor blocking\n")
