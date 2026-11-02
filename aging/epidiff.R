# EpiDiff 0.1.0. Invoke via clocks.py to publish output atomically.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L, args[3] %in% c("complete", "upstream-mean"))
library(EpiDiff)
stopifnot(as.character(packageVersion("EpiDiff")) == "0.1.0")
model <- new.env()
data("rfTRdiff2", package = "EpiDiff", envir = model)
expected <- model$candDiffALL2.v
stopifnot(length(expected) == 1298L, !anyDuplicated(expected), model$rfTRdiff2.lo[[7]]$ntree == 400L)
d <- read.delim(args[1], check.names = FALSE, colClasses = "character", quote = "",
                comment.char = "", na.strings = c("", "NA"))
stopifnot(names(d)[1] == "probe", ncol(d) > 1, !anyDuplicated(names(d)), all(nzchar(names(d))),
          !anyNA(d$probe), !anyDuplicated(d$probe))
x <- as.matrix(d[-1])
storage.mode(x) <- "double"
rownames(x) <- d$probe
if (any(!is.finite(x)) || any(x < 0 | x > 1)) stop("Expected normalized finite beta values in [0,1]")
missing <- sum(!expected %in% rownames(x))
if (missing && args[3] == "complete") stop("Missing ", missing, " of 1298 scoring probes")
if (missing > 85) stop("Coverage below the bundled reference example (1213 of 1298)")
# The published implementation fills absent probe rows with mean(data.m)
# over this entire input. Retain all input rows for the explicit opt-in route.
scoring <- x
if (ncol(x) == 1L) {
    # DoDiffRF drops the sample dimension for one column. Identical copies
    # preserve that dimension and the upstream missing-row mean.
    scoring <- cbind(x, x)
    colnames(scoring) <- make.unique(rep(colnames(x), 2L))
}
res <- DoDiffRF(scoring)
pr <- as.matrix(res$prob)
stopifnot(identical(rownames(pr), colnames(scoring)), identical(colnames(pr), c("0","1","2")),
          all(is.finite(pr)), all(pr >= 0 & pr <= 1), all(abs(rowSums(pr)-1) < 1e-10),
          all(abs(res$di - (0.25 * pr[,2] + pr[,3])) < 1e-10))
if (ncol(x) == 1L) stopifnot(identical(unname(pr[1,]), unname(pr[2,])))
pr <- pr[colnames(x), , drop = FALSE]
res$di <- res$di[seq_len(ncol(x))]
res$class <- res$class[seq_len(ncol(x))]
out <- data.frame(sample = colnames(x), differentiation_index = unname(res$di),
                  p_pluripotent = pr[,1], p_multipotent = pr[,2], p_differentiated = pr[,3],
                  class = as.character(res$class), observed_probes = length(expected) - missing,
                  filled_probes = missing)
write.table(out, args[2], sep = "\t", row.names = FALSE, quote = FALSE)
