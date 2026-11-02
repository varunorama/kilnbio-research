# Rscript export_kang.R <original EH2259 RDS> <new output>
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L,!dir.exists(args[2]))
suppressPackageStartupMessages(library(SingleCellExperiment))
s <- readRDS(args[1]); m <- as.data.frame(colData(s))
stopifnot(all(c('ind','stim','cell','multiplets') %in% names(m)))
dir.create(args[2],recursive=TRUE)
write.csv(as.data.frame(table(m$stim,m$multiplets)),file.path(args[2],'author_assignments.csv'),row.names=FALSE)
# Retain original genetic-demultiplexing singlets. No computational doublet
# algorithm is being evaluated by this example.
keep <- m$multiplets=='singlet' & !is.na(m$ind) & !is.na(m$cell)
s <- s[,keep];m <- m[keep,]
symbols <- rowData(s)$SYMBOL
ids <- trimws(rowData(s)$ENSEMBL)
stopifnot(length(symbols)==nrow(s),length(ids)==nrow(s),!anyDuplicated(ids),all(nzchar(ids)))
libs <- list(); cells <- list()
for(condition in sort(unique(m$stim))) {
  chosen <- which(m$stim==condition)
  folder <- file.path(args[2],condition);dir.create(folder)
  Matrix::writeMM(assay(s,'counts')[,chosen],file.path(folder,'matrix.mtx'))
  write.table(data.frame(ids,symbols,'Gene Expression'),file.path(folder,'features.tsv'),sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
  writeLines(colnames(s)[chosen],file.path(folder,'barcodes.tsv'))
  libs[[condition]] <- data.frame(sample=condition,matrix_dir=condition,condition=condition)
  cells[[condition]] <- data.frame(cell_id=paste0(condition,':',colnames(s)[chosen]),
                                  biological_id=as.character(m$ind[chosen]),cell_type=as.character(m$cell[chosen]))
}
write.csv(do.call(rbind,libs),file.path(args[2],'libraries.csv'),row.names=FALSE)
write.table(do.call(rbind,cells),file.path(args[2],'cell_metadata.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(args[2],'export_session.txt'))
cat(nrow(s),'genes;',ncol(s),'author singlets;',length(unique(m$ind)),'donors\n')
