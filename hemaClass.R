################################################################################
#
# Master script for the hemaclass paper
#
#   Anders Ellern Bilgray & Steffen Falgreen
#
################################################################################

# NOTE: Running this script needs more approx 50 GBs of free disk space

# Initalization
rm(list = ls()) # Clear global enviroment
memory.limit(size = 60000)

# If any of the used packages are missing they will be installed.
pkgs <- c("devtools", "shiny", "matrixStats", "Rcpp", "RcppArmadillo",
          "RcppEigen", "testthat", "WriteXLS", "RLumShiny", "gdata",
          "Biobase", "affy", "affyio", "preprocessCore", "BiocInstaller",
          "AnnotationDbi", "GEOquery", "GEOquery", "oligo",
          "shinysky", "DLBCLdata", "hemaClass")
missing <- setdiff(pkgs, installed.packages()[, "Package"])

if (length(missing)) {
  # Packages from CRAN
  install.packages(pkgs[1:10])

  # Packages form Bioconductor
  source("http://bioconductor.org/biocLite.R")
  biocLite(pkgs[11:19])

  # Packages from github
  devtools::install_github("AnalytixWare/ShinySky")
  devtools::install_github("AEBilgrau/DLBCLdata")
  devtools::install_github("oncoclass/hemaClass", dependencies = TRUE)
}

# Load packages
library("hemaClass")  # Load the hemaclass package
library("DLBCLdata")  # Package for data handling and download
library("devtools")   # For source_url

# Global variables
recompute <- FALSE
verbose <- TRUE

# Load saved data
saved.file <- "saved.RData"
if (file.exists(saved.file)) { load(saved.file) }

################################################################################
# Auxiliary functions
################################################################################

# Logit function
logit <- function(p) log(p/(1-p))

# Load the resave function
source_url(
  "https://raw.githubusercontent.com/AEBilgrau/Bmisc/master/R/resave.R"
)

# Function for doing the one-by-one and study based reference normalization
normalizer <- function(study, gse, nsamples = 30, global = FALSE) {

  # Read files used in cohort
  files <- file.path("data", gse, colnames(rma$cohort[[study]]))
  affy.batch <- readCelfiles(files)

  if (global) {
    rma[["reference"]][["global"]] <- rmaPreprocessing(affy.batch)
  } else {
    rma[["onebyone"]][[study]] <-
      rmaReference(affy.batch, rma[["reference"]][["global"]])$exprs.sc
  }

  # Select the nsamples at random for reference
  ref.samples <- sample(colnames(affy.batch$exprs), nsamples)
  rma[["ref.samples"]][[study]] <- ref.samples

  affy.batch.ref       <- affy.batch
  affy.batch.ref$exprs <- affy.batch.ref$exprs[, rma$ref.samples[[study]]]

  affy.batch.refbased  <- affy.batch
  affy.batch.refbased$exprs <-
    affy.batch.refbased$exprs[, setdiff(colnames(affy.batch$exprs),ref.samples)]

  rma[["reference"]][[study]] <- rmaPreprocessing(affy.batch.ref)
  rma[["refbased"]][[study]]  <-
    rmaReference(affy.batch.refbased, rma[["reference"]][[study]])$exprs.sc

  rma <<- rma
}


################################################################################
# Download and prepare the data
#  + regular RMA normalization
################################################################################

# To get an overview of the data in DLBCLdata run:
data("DLBCL_overview")
# View(DLBCL_overview)
# ls("package:DLBCLdata")

studies <- DLBCL_overview[c(1, 5, 6, 7), 1:4]
rownames(studies) <- studies$GSE
dat <- list()
for (gse in as.character(studies$GSE)) {
  gse.file <- paste0("data/", gse, "/", gse,"_affy.Rds")

  if (!file.exists(gse.file)) { # Download and preprocess data
    downloadAndProcessGEO(geo_nbr = gse, destdir = "data", verbose = verbose)
    gc()  # Garbage collect (clear up some memory)
  }

  # Read data into list
  dat[[gse]] <- readRDS(gse.file)
}

# Organize data
if (!exists("rma", inherits = FALSE) || recompute) {
  rma <- list()
}
rma[["cohort"]][["LLMPPCHOP"]]  <- microarrayScale(exprs(dat$GSE10846$es$CHOP))
rma[["cohort"]][["LLMPPRCHOP"]] <- microarrayScale(exprs(dat$GSE10846$es$'R-CHOP'))
rma[["cohort"]][["MDFCI"]]      <- microarrayScale(exprs(dat$GSE34171$es$GPL570))
rma[["cohort"]][["IDRC"]]       <- microarrayScale(exprs(dat$GSE31312$es$Batch1))
rma[["cohort"]][["CHEPRETRO"]]  <- microarrayScale(exprs(dat$GSE56315$es$DLBCL))

################################################################################
# RMA normalization
# one-by-one (with and without study-based reference)
################################################################################

if (is.null(rma$reference) || recompute || TRUE) {

  # First the overall reference is made using the LLMPP CHOP data
  normalizer(study = "LLMPPCHOP", gse = "GSE10846", global = TRUE); gc()

  # Next the other datasets are nomalized according to LLMPPCHOP
  # and the 30 sample study based reference
  set.seed(1154538003)
  normalizer(study = "LLMPPRCHOP", gse = "GSE10846"); gc()
  normalizer(study = "MDFCI",      gse = "GSE34171"); gc()
  normalizer(study = "IDRC",       gse = "GSE31312"); gc()
  normalizer(study = "CHEPRETRO",  gse = "GSE56315"); gc()

  resave(rma, file = saved.file)
}

################################################################################
# Establish the results
################################################################################

results <- list()
for (study in c("CHEPRETRO", "MDFCI", "IDRC", "LLMPPRCHOP", "LLMPPCHOP")) {
  message("Creating results for", study)

  b <- (study != "LLMPPCHOP")

  # ABC/GCB
  results[["ABCGCB"]][[study]]$cohort   <- ABCGCB(rma$cohort[[study]])
  results[["ABCGCB"]][[study]]$refbased <- ABCGCB(rma$refbased[[study]])
  if (b) results[["ABCGCB"]][[study]]$onebyone <- ABCGCB(rma$onebyone[[study]])

  # BAGS
  results[["BAGS"]][[study]]$cohort   <- BAGS(rma$cohort[[study]])
  results[["BAGS"]][[study]]$refbased <- BAGS(rma$refbased[[study]])
  if (b) results[["BAGS"]][[study]]$onebyone <- BAGS(rma$onebyone[[study]])

  # REGS, the classifier for Cyclophosphamide, Doxorubicin, and Vincristine:
  RC <- ResistanceClassifier
  results[["REGS"]][[study]]$cohort   <- RC(rma$cohort[[study]])
  results[["REGS"]][[study]]$refbased <- RC(rma$refbased[[study]])
  if (b) results[["REGS"]][[study]]$onebyone <- RC(rma$onebyone[[study]])
}
rm(RC)


################################################################################
# Comparisons of the results
################################################################################

library("psych")
weight <- matrix(c(0, 1, 2, 1, 0, 1, 2, 1, 0), 3)/2

testfun <- function(x, y, dec = 2, weight = NULL){
  tab <- table(x, y)
  acc <- binom.test(c(sum(diag(tab)), sum(tab) - sum(diag(tab))))
  acc <- paste0(round(acc$estimate, dec),
                " (", paste(round(acc$conf.int, dec), collapse = ","), ")")
  kap <- cohen.kappa(data.frame(x,y), w = weight)
  conf.1 <- kap$confid[2,1]
  conf.2 <- min(kap$confid[2,3],1)
  kappa <- paste0(round(kap$weighted.kappa, dec),
                  " (", round(conf.1, dec), ",", round(conf.2, dec), ")")

  return(c("accuracy" = acc, "kappa" = kappa))
}

getLogit <- function(r, type.y, type.x = "cohort") {
  y <- logit(r[[type.y]]$prob)
  x <- logit(r[[type.x]]$prob)
  x <- x[rownames(y), , drop = FALSE]
  ans <- cbind(x, y)
  colnames(ans) <- c(type.x, type.y)
  return(as.data.frame(ans))
}

getClass <- function(r, type.y, type.x = "cohort") {
  y <- r[[type.y]]$class
  names(y) <- rownames(r[[type.y]]$prob)
  x <- r[[type.x]]$class
  names(x) <- rownames(r[[type.x]]$prob)
  x <- x[names(y)]
  ans <- data.frame(x, y)
  colnames(ans) <- c(type.x, type.y)
  return(ans)
}

plotline <- function(x, y, ...) {
  x[x == Inf | y == Inf] <- NaN
  y[y == Inf | x == Inf] <- NaN

  r <- prcomp(~x + y)
  slope <- r$rotation[2, 1]/r$rotation[1, 1]
  intercept <- r$center[2] - slope*r$center[1]
  abline(intercept, slope, ...)
}


addLegend <- function(r, type, cor.test = NULL, dec = 3) {
  tmp1 <- getClass(r, type)
  tmp2 <- testfun(tmp1[,1], tmp1[,2], weight = weight)
  if (!is.null(cor.test)) {
    tmp3 <- paste0(round(cor.test$estimate, dec), " (",
                   paste(round(cor.test$conf.int, dec), collapse = ", "), ")")
    tmp2 <- c("correlation" = tmp3, tmp2)
  }
  legend("topleft", legend = paste(names(tmp2), "=", tmp2))
}


#
# ABC/GCB
#

pdf("figures/results_overview_ABCGCB.pdf", width = 14)
for (study in c("CHEPRETRO", "MDFCI", "IDRC", "LLMPPRCHOP")) {
  res <- results[["ABCGCB"]][[study]]

  par(mfrow = c(1,2))

  # ONE BY ONE
  with(as.data.frame(getLogit(res, "onebyone")), {
    plot(cohort, onebyone)
    abline(0, 1, lty = 2, lwd = 2)
    abline(h = logit(c(0.1,0.9)), v = logit(c(0.1,0.9)), col = "darkgrey",
           lwd = 2, lty = 2)
    addLegend(res, "onebyone", cor.test(cohort, onebyone))
    plotline(cohort, onebyone, col = "red")
  })

  # REF BASED
  with(A <- as.data.frame(getLogit(res, "refbased")), {
    plot(cohort, refbased)
    abline(0, 1, lty = 2, lwd = 2)
    abline(h = logit(c(0.1,0.9)), v = logit(c(0.1,0.9)), col = "darkgrey",
           lwd = 2, lty = 2)
    addLegend(res, "refbased", cor.test(cohort, refbased))
    plotline(cohort, refbased, col = "red")
  })

  title(main = study, sub = "ABC/GCB classification", outer = TRUE, line = -1)
}
dev.off()










#
sink(file = "sessionInfo.txt")
print(sessionInfo())
sink()
