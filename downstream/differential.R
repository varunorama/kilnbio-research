# KilnBio statistical engines. No normalized-input inference or automatic metadata omission.
args <- commandArgs(trailingOnly=TRUE)
extra <- Sys.getenv('KILNBIO_R_LIBRARY')
if(nzchar(extra)) .libPaths(c(extra,.libPaths()))
if(!requireNamespace('jsonlite',quietly=TRUE)) stop('Install jsonlite; see environments/setup_r.R')
cfg <- jsonlite::fromJSON(args[1],simplifyVector=FALSE); out <- args[2]
set.seed(20261031)
samples <- read.csv(cfg$samples,check.names=FALSE,colClasses='character',na.strings=c('NA',''))
rownames(samples) <- samples$sample
group <- cfg$group
samples[[group]] <- relevel(factor(samples[[group]]),ref=cfg$denominator)
covariates <- vapply(cfg$covariates,function(x)x$name,character(1))
for(term in cfg$covariates) samples[[term$name]] <- if(term$type=='factor')factor(samples[[term$name]]) else as.numeric(samples[[term$name]])
formula <- reformulate(c(covariates,group))
design <- model.matrix(formula,samples)
stopifnot(!anyNA(design),qr(design)$rank==ncol(design),nrow(design)>ncol(design))
# Identify the exact treatment-coded contrast by constructing design rows.
lo <- hi <- samples[1,,drop=FALSE]
lo[[group]] <- factor(cfg$denominator,levels=levels(samples[[group]]))
hi[[group]] <- factor(cfg$numerator,levels=levels(samples[[group]]))
contrast <- as.numeric(model.matrix(formula,hi,contrasts.arg=attr(design,'contrasts'))-
                       model.matrix(formula,lo,contrasts.arg=attr(design,'contrasts')))
stopifnot(sum(contrast!=0)==1, sum(contrast)==1)
coef_index <- which(contrast!=0)
write.table(design,file.path(out,'provenance/design.tsv'),sep='\t',quote=FALSE,col.names=NA)
writeLines(deparse(formula),file.path(out,'provenance/design_formula.txt'))
write_tsv <- function(x,name)write.table(x,file.path(out,'results',name),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
packages <- c('jsonlite')
if(cfg$mode=='de') {
  estimated <- cfg$count_source %in% c('salmon_tximport','rsem_tximport')
  if(estimated) {
    packages <- c(packages,'tximport','DESeq2')
    files <- unlist(cfg$quant_files)[samples$sample]
    if(cfg$count_source=='rsem_tximport') {
      read_rsem <- function(path) read.delim(path, check.names=FALSE, colClasses=c(gene_id='character'))
      txi <- tximport::tximport(files,type='rsem',txIn=FALSE,txOut=FALSE,importer=read_rsem)
    } else {
      tx2gene <- read.delim(cfg$tx2gene,check.names=FALSE,colClasses='character')
      txi <- tximport::tximport(files,type='salmon',tx2gene=tx2gene,countsFromAbundance='no')
    }
    counts <- txi$counts
  } else {
    tab <- read.delim(cfg$counts,check.names=FALSE,colClasses=c(feature_id='character'))
    counts <- as.matrix(tab[,samples$sample,drop=FALSE]);rownames(counts)<-tab$feature_id
    stopifnot(all(is.finite(counts)),all(counts>=0),all(counts==floor(counts)))
  }
  keep <- rowSums(counts>=cfg$min_count)>=cfg$min_samples
  filter_table <- data.frame(feature_id=rownames(counts),included=keep)
  if(estimated) {
    filter_table$valid_effective_length <- rowSums(!is.finite(txi$length) | txi$length<=0)==0
    if (cfg$count_source == 'rsem_tximport') {
      zero_length_ids <- unlist(cfg$rsem_zero_length_genes, use.names=FALSE)
      filter_table$valid_effective_length <- filter_table$valid_effective_length &
        !(filter_table$feature_id %in% zero_length_ids)
    }
    keep <- keep & filter_table$valid_effective_length
    filter_table$included <- keep
  }
  write_tsv(filter_table,'feature_filter.tsv')
  if(sum(keep)<10)stop('Fewer than ten retained features: insufficient for this dispersion workflow')
  if(any(colSums(counts[keep,,drop=FALSE])==0))stop('A retained library has zero counts')
  external_sf <- NULL
  if(!is.null(cfg$external_size_factors)) {
    sf <- read.csv(cfg$external_size_factors,colClasses=c(sample='character'))
    external_sf <- sf$size_factor[match(samples$sample,sf$sample)]
    stopifnot(!estimated,all(is.finite(external_sf)),all(external_sf>0))
    external_sf <- external_sf/exp(mean(log(external_sf)))
  }
  if(cfg$engine=='deseq2') {
    packages<-unique(c(packages,'DESeq2'))
    if(estimated) {
      for(component in c('abundance','counts','length')) txi[[component]] <- txi[[component]][keep,,drop=FALSE]
      dds<-DESeq2::DESeqDataSetFromTximport(txi,colData=samples,design=formula)
    } else {
      storage.mode(counts)<-'integer'
      dds<-DESeq2::DESeqDataSetFromMatrix(counts[keep,,drop=FALSE],colData=samples,design=formula)
    }
    if(!is.null(external_sf)) DESeq2::sizeFactors(dds) <- external_sf
    dds<-DESeq2::DESeq(dds,quiet=TRUE,sfType='poscounts',parallel=FALSE)
    result<-as.data.frame(DESeq2::results(dds,contrast=c(group,cfg$numerator,cfg$denominator),alpha=cfg$fdr))
    result<-data.frame(feature_id=rownames(result),result,check.names=FALSE)
    write_tsv(data.frame(feature_id=rownames(dds),DESeq2::counts(dds,normalized=TRUE),check.names=FALSE),'normalized_counts.tsv')
    saveRDS(dds,file.path(out,'results/model.rds'))
    write_tsv(result,'all_features.tsv')
    write_tsv(result[!is.na(result$padj)&result$padj<cfg$fdr,],'discoveries.tsv')
    diagnostic<-data.frame(sample=samples$sample)
    if(estimated) {
      diagnostic$median_gene_normalization_factor <- apply(DESeq2::normalizationFactors(dds),2,median)
    } else {
      diagnostic$size_factor <- DESeq2::sizeFactors(dds)
    }
    logcounts<-SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(dds,blind=FALSE))
  } else {
    packages<-c(packages,'edgeR')
    y<-edgeR::DGEList(counts=counts[keep,,drop=FALSE])
    if(is.null(external_sf)) {
      y<-edgeR::calcNormFactors(y,method='TMM')
    } else {
      # Effective library sizes are proportional to reference-derived exposure.
      nf <- external_sf/y$samples$lib.size
      y$samples$norm.factors <- nf/exp(mean(log(nf)))
    }
    # edgeR v4 adjusted-deviance QL path estimates dispersion when NULL.
    if(packageVersion('edgeR')<'4.2.0')stop('This preset requires edgeR >=4.2.0')
    fit<-edgeR::glmQLFit(y,design,robust=TRUE,legacy=FALSE)
    test<-edgeR::glmQLFTest(fit,contrast=contrast)
    result<-edgeR::topTags(test,n=Inf,sort.by='none')$table
    result<-data.frame(feature_id=rownames(result),result,check.names=FALSE)
    write_tsv(result,'all_features.tsv')
    write_tsv(result[!is.na(result$FDR)&result$FDR<cfg$fdr,],'discoveries.tsv')
    write_tsv(data.frame(feature_id=rownames(y),edgeR::cpm(y,normalized.lib.sizes=TRUE),check.names=FALSE),'normalized_cpm.tsv')
    saveRDS(list(dge=y,fit=fit,test=test),file.path(out,'results/model.rds'))
    diagnostic<-data.frame(sample=samples$sample,library_size=y$samples$lib.size,normalization_factor=y$samples$norm.factors)
    if(is.null(external_sf)) diagnostic$TMM_factor<-y$samples$norm.factors
    diagnostic$effective_library_size<-y$samples$lib.size*y$samples$norm.factors
    logcounts<-edgeR::cpm(y,log=TRUE,prior.count=2)
  }
  diagnostic$normalization_method <- if(!is.null(external_sf)) 'external_reference' else if(estimated) 'tximport_length_and_library' else if(cfg$engine=='edger') 'TMM' else 'poscounts'
  write_tsv(diagnostic,'sample_diagnostics.tsv')
  pdf(file.path(out,'results/QC.pdf'))
  heatmap(cor(logcounts),symm=TRUE,main='Library correlation (normalized transformed counts)')
  plot(density(logcounts[,1]),main='Normalized transformed distributions',xlab='VST (DESeq2) or log CPM (edgeR)')
  if(ncol(logcounts)>1)for(i in 2:ncol(logcounts))lines(density(logcounts[,i]),col=i)
  pcs<-prcomp(t(logcounts))$x
  plot(pcs[,1],pcs[,2],col=as.integer(samples[[group]]),pch=19,xlab='PC1',ylab='PC2',main='Exploratory sample PCA')
  text(pcs[,1],pcs[,2],labels=samples$sample,pos=3,cex=.7)
  lfc<-if(cfg$engine=='deseq2')result$log2FoldChange else result$logFC
  padj<-if(cfg$engine=='deseq2')result$padj else result$FDR
  plot(lfc,-log10(pmax(padj,.Machine$double.xmin)),pch=19,cex=.5,col=ifelse(!is.na(padj)&padj<cfg$fdr,'firebrick','gray50'),xlab='log2 fold change (numerator / denominator)',ylab='-log10 adjusted p',main='Differential features')
  dev.off()
} else {
  packages<-c(packages,'DSS','bsseq')
  input<-lapply(samples$sample,function(s){
    x<-read.delim(cfg$methylation_files[[s]],check.names=FALSE,colClasses=c('character','integer','numeric','numeric'))
    if(!identical(names(x),c('chr','pos','N','X')))stop('Methylation counts require chr,pos,N,X columns')
    if(anyNA(x)||any(x$pos<1)||any(x$N<0)||any(x$X<0)||any(x$X>x$N)||any(x$X!=floor(x$X))||any(x$N!=floor(x$N)))stop('Invalid methylation counts')
    if(anyDuplicated(x[,c('chr','pos')]))stop('Duplicate sites: combine strands or lanes before modeling')
    x[order(x$chr,x$pos),]
  })
  bs<-DSS::makeBSseqData(input,samples$sample)
  coverage<-bsseq::getCoverage(bs,type='Cov')
  pass<-rep(TRUE,nrow(coverage))
  for(level in levels(samples[[group]]))pass<-pass & rowSums(coverage[,samples[[group]]==level,drop=FALSE]>=cfg$min_coverage)>=2
  write_tsv(data.frame(chr=as.character(GenomicRanges::seqnames(bs)),pos=BiocGenerics::start(bs),included=pass),'locus_filter.tsv')
  if(sum(pass)<20)stop('Fewer than twenty sufficiently covered loci')
  bs<-bs[pass,]
  if(cfg$engine=='dmrseq') {
    packages<-c(packages,'dmrseq','BiocParallel')
    metadata<-samples
    metadata[[group]]<-as.integer(samples[[group]]==cfg$numerator)
    SummarizedExperiment::colData(bs)<-S4Vectors::DataFrame(metadata)
    regions<-dmrseq::dmrseq(bs,testCovariate=group,
      adjustCovariate=if(length(covariates))covariates else NULL,
      cutoff=cfg$region_cutoff,minNumRegion=max(3,cfg$dmr_min_cpg),
      maxPerms=cfg$max_permutations,smooth=cfg$smoothing,
      BPPARAM=BiocParallel::SerialParam())
    regions<-as.data.frame(regions)
    # Preserve nested region metadata in RDS; table columns must be scalar.
    saveRDS(regions,file.path(out,'results/dmrseq_regions.rds'))
    flat<-regions[,!vapply(regions,is.list,logical(1)),drop=FALSE]
    write_tsv(flat,'all_regions.tsv')
    discoveries<-if('qval' %in% names(flat))flat[!is.na(flat$qval)&flat$qval<cfg$fdr,,drop=FALSE] else flat[FALSE,,drop=FALSE]
    write_tsv(discoveries,'DMR.tsv')
    writeLines(c('dmrseq region-level q-values; review permutation exchangeability and attainable resolution.',
                 'Positive group coding is numerator=1 versus denominator=0.',
                 'No single-CpG DML calls are inferred from region results.'),file.path(out,'results/INTERPRETATION.txt'))
  } else {
  fit<-DSS::DMLfit.multiFactor(bs,design=samples,formula=formula,smoothing=cfg$smoothing)
  test<-DSS::DMLtest.multiFactor(fit,coef=coef_index)
  table<-as.data.frame(test)
  table$FDR_BH<-p.adjust(table$pvals,method='BH')
  # Group differences below are descriptive, unadjusted coverage-weighted fractions.
  X<-bsseq::getCoverage(bs,type='M');N<-bsseq::getCoverage(bs,type='Cov')
  means<-lapply(c(cfg$numerator,cfg$denominator),function(g){z<-samples[[group]]==g;rowSums(X[,z,drop=FALSE])/rowSums(N[,z,drop=FALSE])})
  table$unadjusted_delta_fraction<-means[[1]]-means[[2]]
  write_tsv(table,'all_loci.tsv')
  write_tsv(table[!is.na(table$FDR_BH)&table$FDR_BH<cfg$fdr,],'DML.tsv')
  regions<-DSS::callDMR(test,p.threshold=cfg$dmr_p_threshold,minlen=cfg$dmr_min_length,
                       minCG=cfg$dmr_min_cpg,dis.merge=cfg$dmr_merge_distance,pct.sig=cfg$dmr_fraction_significant)
  if(is.null(regions))regions<-data.frame(chr=character(),start=integer(),end=integer(),length=integer(),nCG=integer(),areaStat=numeric())
  write_tsv(regions,'DMR_candidates.tsv')
  writeLines(c('DSS multi-factor locus tests; BH adjustment applies to loci only.',
               'DMR_candidates.tsv contains DSS p-threshold/merging regions, NOT region-level FDR-controlled discoveries.',
               'Coordinates are 1-based closed; unadjusted_delta_fraction is descriptive, not the adjusted model coefficient.',
               'Modification-channel calling/conversion uncertainty is not modeled by DSS.'),file.path(out,'results/INTERPRETATION.txt'))
  saveRDS(list(fit=fit,test=test),file.path(out,'results/model.rds'))
  }
}
write.csv(data.frame(package=unique(packages),version=vapply(unique(packages),function(p)as.character(packageVersion(p)),character(1))),file.path(out,'provenance/package_versions.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'provenance/sessionInfo.txt'))
