# Link one declared genomic event to one quantitative assay. Association, not causation.
args <- commandArgs(TRUE)
if(identical(args,"--help")) {cat("kilnbio genomic-associate calls.tsv event samples.tsv matrix.tsv {specimen|donor} NEW_OUTPUT\n");quit(status=0)}
if(length(args)!=6L)stop("See --help")
read <- function(f) read.delim(f,check.names=FALSE,colClasses="character",quote="",comment.char="",na.strings=c("","NA"))
calls <- read(args[1]); map <- read(args[3]); x <- read(args[4]); by <- args[5]; out <- args[6]
needed <- function(d,cols) all(cols%in%names(d)) && !anyDuplicated(names(d))
if(!needed(calls,c("event","origin","donor","specimen","state")) || !needed(map,c("sample","donor","specimen")) || names(x)[1]!="feature" || anyDuplicated(names(x)) || !by%in%c("specimen","donor") || file.exists(out))stop("Invalid schema/matching level or destination exists")
g <- calls[!is.na(calls$event) & calls$event==args[2],,drop=FALSE]
if(nrow(g)==0 || length(unique(g$origin))!=1 || !all(g$origin%in%c("germline","somatic")) || !all(g$state%in%c("altered","reference","unknown")))stop("Select one event and source class; unknown calls must be explicit")
key <- function(d) if(by=="donor")d$donor else paste(d$donor,d$specimen,sep="\037")
if(anyNA(g$donor)||anyNA(map$donor)||anyNA(map$sample)||anyDuplicated(map$sample)||anyDuplicated(map$donor)||anyDuplicated(key(g))||anyNA(x$feature)||anyDuplicated(x$feature)||!setequal(names(x)[-1],map$sample))stop("Ambiguous linkage; use one independent assay observation per donor and exact sample columns")
if(by=="specimen" && (anyNA(map$specimen)||anyNA(g$specimen)))stop("Explicit specimen IDs needed; do not assume matched tissue")
if(g$origin[1]=="germline" && anyDuplicated(g$donor))stop("Choose one germline genotype per donor")
raw <- as.matrix(x[,map$sample,drop=FALSE]); a <- matrix(suppressWarnings(as.numeric(raw)),nrow=nrow(raw),dimnames=list(x$feature,map$sample))
if(any(is.na(a)!=is.na(raw))||any(is.infinite(a)))stop("Need finite numerical traits or missing values")
j <- match(key(map),key(g)); state <- g$state[j];state[is.na(state)] <- "unmatched"
aa <- which(state=="altered");bb <- which(state=="reference")
if(min(length(aa),length(bb))<3L)stop("Fewer than three independent donors in a genomic group")
stat <- function(z) {
    u<-z[aa];v<-z[bb];u<-u[is.finite(u)];v<-v[is.finite(v)]
    if(min(length(u),length(v))<3)return(c(length(u),length(v),NA_real_,NA_real_))
    p<-if(length(unique(c(u,v)))<2)1 else suppressWarnings(wilcox.test(u,v,exact=FALSE,correct=TRUE)$p.value)
    c(length(u),length(v),median(u)-median(v),p)
}
res <- as.data.frame(t(apply(a,1,stat)));names(res)<-c("n_altered","n_reference","median_difference","p")
res$q<-NA_real_;good<-is.finite(res$p);if(any(good))res$q[good]<-p.adjust(res$p[good],"BH",n=nrow(res))
dir.create(dirname(out),showWarnings=FALSE,recursive=TRUE);stage<-tempfile(".genomics-",tmpdir=dirname(out));if(!dir.create(stage))stop("Cannot create staging directory")
write <- function(obj,n)write.table(obj,file.path(stage,n),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
write(data.frame(feature=x$feature,res),"association.tsv")
write(data.frame(map,genomic_state=state),"linkage.tsv")
write(data.frame(event=args[2],origin=g$origin[1],matching=by,source_calls=nrow(g),assay_donors=nrow(map),altered=length(aa),reference=length(bb),unknown=sum(state=="unknown"),unmatched=sum(state=="unmatched"),features=nrow(a),testable=sum(good),BH_family=nrow(a)),"audit.tsv")
write(data.frame(role=c("calls","design","traits"),md5=unname(tools::md5sum(args[c(1,3,4)]))),"inputs.tsv")
writeLines(capture.output(sessionInfo()),file.path(stage,"session.txt"));if(!file.rename(stage,out))stop("Could not publish results")
