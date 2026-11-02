# Public MOFA CLL worked example. Run through run.py for input verification.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
input <- args[1]
out <- args[2]
load(file.path(input, "CLL_data.RData"))
load(file.path(input, "CLL_covariates.RData"))
write_tsv <- function(x, name) {
  write.table(x, file.path(out, name), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

# First study-coded drug, middle concentration step. No other response is a predictor.
endpoint <- "D_001_3"
donors <- sort(colnames(CLL_data$Drugs))
stopifnot(!anyDuplicated(donors), !anyNA(donors), endpoint %in% rownames(CLL_data$Drugs))
for (block in CLL_data[c("mRNA", "Methylation")]) {
  stopifnot(!anyDuplicated(colnames(block)), setequal(donors, colnames(block)),
            !anyDuplicated(rownames(block)), !any(is.infinite(block)))
}
stopifnot(!anyDuplicated(rownames(CLL_covariates)), setequal(donors, rownames(CLL_covariates)))
blocks <- lapply(CLL_data[c("mRNA", "Methylation")], function(x) t(x[, donors, drop = FALSE]))
response <- CLL_data$Drugs[endpoint, donors]
gender <- CLL_covariates[donors, "Gender"]
cohort <- data.frame(donor = donors, outcome = is.finite(response), gender = !is.na(gender),
                     RNA = rowSums(is.finite(blocks$mRNA)) > 0,
                     DNAm = rowSums(is.finite(blocks$Methylation)) > 0)
cohort$eligible <- with(cohort, outcome & gender & RNA & DNAm)
write_tsv(cohort, "cohort.tsv")
eligible <- which(cohort$eligible)
stopifnot(length(eligible) >= 40, all(levels(gender) == c("f", "m")))
covariates <- cbind(intercept = 1, source_gender_m = as.numeric(gender == "m"))
set.seed(17041)
test <- sort(sample(eligible, ceiling(length(eligible) / 4)))
train <- setdiff(eligible, test)
fold <- sample(rep(seq_len(5), length.out = length(train)))
write_tsv(data.frame(donor = donors[c(train, test)],
                     partition = rep(c("development", "holdout"), c(length(train), length(test))),
                     inner_fold = c(fold, rep(NA, length(test)))), "partitions.tsv")

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "..", "..", "downstream", "compare_assays.R"))
recipes <- list(gender_only=character(), RNA="mRNA", DNAm="Methylation", RNA_DNAm=c("mRNA","Methylation"))
fit_assays(donors, response, covariates, blocks, recipes, train, test, fold, endpoint, out, gender_coding="m=1, f=0")
