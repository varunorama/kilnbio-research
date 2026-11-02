# Rscript export_airway.R <unpacked airway package or installed> <new output directory>
# Export original counts/metadata without installing or modifying the data package.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L, !dir.exists(args[2]))
suppressPackageStartupMessages(library(SummarizedExperiment))
data_env <- new.env()
if(args[1]=='installed') data('airway',package='airway',envir=data_env) else
  load(file.path(args[1], 'data', 'airway.rda'), envir=data_env)
se <- data_env$airway
stopifnot(ncol(se)==8L, all(assay(se)>=0), all(assay(se)==floor(assay(se))))
meta <- as.data.frame(colData(se))
stopifnot(all(c('cell','dex') %in% names(meta)))
samples <- data.frame(sample=colnames(se), biological_id=as.character(meta$cell),
                      condition=as.character(meta$dex))
stopifnot(all(table(samples$biological_id,samples$condition)==1L))
dir.create(args[2], recursive=TRUE)
write.table(data.frame(feature_id=rownames(se),assay(se),check.names=FALSE),
            file.path(args[2],'counts.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.csv(samples,file.path(args[2],'samples.csv'),row.names=FALSE)
write.csv(meta,file.path(args[2],'source_metadata.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(args[2],'export_session.txt'))
cat(nrow(se), 'genes;',ncol(se),'libraries;',length(unique(samples$biological_id)),'paired donors\n')
