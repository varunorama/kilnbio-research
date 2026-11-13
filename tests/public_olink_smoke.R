# Export the vendor's synthetic NPX fixture for an independent importer check.
# Source: Olink-Proteomics/OlinkRPackage commit 7385f477050b65e9470e709d30b1b23c6fcb5702
# OlinkAnalyze/data/npx_data1.rda; this is synthetic, not a biological validation.
a<-commandArgs(TRUE)
if(length(a)!=2)stop("Expected downloaded npx_data1.rda and NEW directory")
if(file.exists(a[2]))stop("Choose a new directory")
e<-new.env();load(a[1],envir=e);d<-as.data.frame(e$npx_data1)
d<-d[!grepl("^CONTROL",d$SampleID) & !is.na(d$Time) & d$Time=="Baseline",]
roster<-unique(d[,c("SampleID","Subject","Treatment")]);names(roster)<-c("sample","donor","treatment")
stopifnot(!anyDuplicated(roster$sample),!anyDuplicated(roster$donor))
dir.create(a[2],recursive=TRUE)
write.table(d,file.path(a[2],"export.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(roster,file.path(a[2],"samples.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
cat(nrow(d),"rows;",nrow(roster),"synthetic baseline samples\n")
