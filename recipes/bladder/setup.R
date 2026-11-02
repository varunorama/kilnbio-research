args <- commandArgs(trailingOnly=TRUE)
if (getRversion()<"4.6.0" || getRversion()>="4.7.0") stop("For this pinned installation use R 4.6.x")
if (identical(args[1],"check")) quit(status=0L)
lib <- args[1]
.libPaths(c(lib,.Library))
for (p in args[3:4]) install.packages(p,repos=NULL,lib=lib,type=args[2],dependencies=FALSE)
for (pkg in c("statmod","limma")) {
    expected <- c(statmod="1.5.2",limma="3.68.5")
    stopifnot(as.character(packageVersion(pkg,lib.loc=lib))==unname(expected[pkg]))
}
cat("Installed verified statmod and limma in",lib,"\n")
