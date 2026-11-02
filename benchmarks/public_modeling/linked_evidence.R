# Repeated discovery screens with known truth. Base R, no participant data.
# All signs, equal weights, layers and loci are set before looking at outcomes.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || file.exists(args[1])) stop("Supply a new output directory")
set.seed(20261101)
k <- 4L
alpha <- 0.05
loci <- 1000L
screens <- 500L
rho <- 0.6
assay_p <- 0.08
signed_p <- function(z, rho) {
  pnorm(rowSums(z) / sqrt(k + k * (k - 1) * rho), lower.tail = FALSE)
}
correlations <- c(0, 0.5, 0.8, 1)
z_equal <- matrix(qnorm(assay_p, lower.tail = FALSE), 1, k)
illustration <- data.frame(rho = correlations, one_sided_input_p = assay_p,
    joint_global_p = vapply(correlations, function(r) signed_p(z_equal, r), numeric(1)),
    partial_at_least_two_bonferroni_p = (k - 1) * assay_p)
stopifnot(abs(tail(illustration$joint_global_p, 1) - assay_p) < 1e-12)

simulate_screen <- function(mixture) {
  z <- sqrt(1 - rho) * matrix(rnorm(loci * k), loci, k) +
    sqrt(rho) * rnorm(loci)
  type <- rep("null", loci)
  if (mixture) {
    # Equal numbers of singleton effects in each assay. All signals point in
    # the prespecified direction. Effect sizes are on the z-statistic scale.
    z[cbind(1:80, rep(seq_len(k), each = 20))] <- z[cbind(1:80, rep(seq_len(k), each = 20))] + 4
    z[81:160, ] <- z[81:160, ] + 1.9
    type[1:80] <- "one"
    type[81:160] <- "four"
  }
  p <- pnorm(z, lower.tail = FALSE)
  # Second order statistic of four p-values, using order statistics of pairs.
  a <- p[, 1]; b <- p[, 2]; c <- p[, 3]; d <- p[, 4]
  second <- pmin(pmax(a, b), pmax(c, d), pmax(pmin(a, b), pmin(c, d)))
  tests <- list(first_assay = p[, 1], any_assay_bonferroni = pmin(1, k * pmin(a, b, c, d)),
                signed_correct_covariance = signed_p(z, rho),
                signed_ignore_covariance = signed_p(z, 0),
                partial_two_bonferroni = pmin(1, (k - 1) * second))
  do.call(rbind, lapply(names(tests), function(method) {
    call <- p.adjust(tests[[method]], method = "BH") <= alpha
    null <- switch(method, partial_two_bonferroni = type != "four",
                   first_assay = type == "null" | (type == "one" & seq_len(loci) > 20),
                   type == "null")
    data.frame(method, fdp = sum(null & call) / max(1L, sum(call)), discoveries = sum(call),
               four_layer_power = if (mixture) mean(call[type == "four"]) else NA_real_,
               one_layer_call_rate = if (mixture) mean(call[type == "one"]) else NA_real_)
  }))
}
all_runs <- list()
for (scenario in c("all_null", "mixed")) {
  all_runs[[scenario]] <- do.call(rbind, lapply(seq_len(screens), function(i) {
    data.frame(scenario, screen = i, simulate_screen(scenario == "mixed"))
  }))
}
runs <- do.call(rbind, all_runs)
summary_rows <- do.call(rbind, lapply(split(runs, list(runs$scenario, runs$method)), function(x) {
  data.frame(scenario = x$scenario[1], method = x$method[1], screens = screens,
             mean_FDP = mean(x$fdp), MC_SE_FDR = sd(x$fdp) / sqrt(screens),
             mean_discoveries = mean(x$discoveries), four_layer_power = mean(x$four_layer_power),
             one_layer_call_rate = mean(x$one_layer_call_rate))
}))
dir.create(args[1], recursive = TRUE, showWarnings = FALSE)
save_table <- function(x, name) write.table(x, file.path(args[1], name), sep = "\t",
                                            quote = FALSE, row.names = FALSE, na = "NA")
save_table(illustration, "four_assays.tsv")
save_table(runs, "per_screen.tsv")
save_table(summary_rows, "calibration.tsv")
writeLines(capture.output(sessionInfo()), file.path(args[1], "R_session.txt"))
print(illustration, row.names = FALSE)
print(summary_rows, row.names = FALSE, digits = 4)
