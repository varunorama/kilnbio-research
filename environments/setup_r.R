# Usage: Rscript setup_r.R /absolute/path/to/private-library
# Reviewed baseline: R 4.6.2 / Bioconductor 3.23. No global library changes.
args<-commandArgs(TRUE);if(length(args)!=1)stop('Supply a private library directory')
lib<-args[1];dir.create(lib,recursive=TRUE,showWarnings=FALSE);.libPaths(c(lib,.libPaths()))
options(repos=c(CRAN='https://cloud.r-project.org'),timeout=1200,Ncpus=2L)
if(!requireNamespace('BiocManager',quietly=TRUE))install.packages('BiocManager',lib=lib)
if(as.character(getRversion())!='4.6.2')stop('Baseline is R 4.6.2; evaluate other R releases explicitly')
BiocManager::install(version='3.23',ask=FALSE,update=FALSE,lib=lib)
packages<-c('DESeq2','edgeR','limma','DSS','bsseq','dmrseq','fgsea','tximport','jsonlite','anndataR','SingleCellExperiment')
BiocManager::install(packages,ask=FALSE,update=FALSE,lib=lib)
write.csv(data.frame(package=packages,version=vapply(packages,function(p)as.character(packageVersion(p)),character(1))),file.path(lib,'kilnbio_package_versions.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(lib,'kilnbio_sessionInfo.txt'))
# This bootstrap follows the stated Bioconductor release, not an immutable lock.
# Archive the installed packages/environment for exact reconstruction.
