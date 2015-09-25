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
# memory.limit(size = 60000)  # If using a Windows machine

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
library("psych")
library("Hmisc")      # For latex

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
  affy.batch.ref$exprs <- affy.batch.ref$exprs[,  ref.samples]

  affy.batch.refbased       <- affy.batch
  affy.batch.refbased$exprs <-
    affy.batch.refbased$exprs[, setdiff(colnames(affy.batch$exprs),ref.samples)]

  # Build reference
  rma[["reference"]][[study]] <- rmaPreprocessing(affy.batch.ref)

  # RMA normalize using reference
  rma[["refbased"]][[study]]  <-
    rmaReference(affy.batch.refbased, rma[["reference"]][[study]])$exprs.sc
  # rmaReference(affy.batch, rma[["reference"]][[study]])$exprs.sc

  rma <<- rma
}


# Format confidence intervals
formatCI <- function(est, ci, dec = 2) {
  paste0(round(est, dec), " (", paste(round(ci, dec), collapse = ", "), ")")
}

# Functions for comparing results
testfun <- function(x, y, dec = 2, weight = NULL){
  tab <- table(x, y)
  acc <- binom.test(c(sum(diag(tab)), sum(tab) - sum(diag(tab))))
  acc <- formatCI(acc$estimate, acc$conf.int, dec)

  kap <- cohen.kappa(tab, w = weight)
  conf <- pmin(kap$confid[2, c("lower", "upper")], 1)
  kappa <- formatCI(kap$weighted.kappa, conf)

  return(c("accuracy" = acc, "kappa" = kappa))
}

summclasses <- function(xc, yc, x, y, weight = NULL, dec = 3) {
  ans <- testfun(xc, yc, weight = weight)
  if (!missing(x) & !missing(y)) {
    rho <- cor.test(x, y)
    ans <- c(ans, rho = formatCI(rho$estimate, rho$conf.int, dec))
  }
  return(ans)
}

# Add orthogonal regression (total least squares) line
plotline <- function(x, y, ...) {
  isna <- is.na(x) | is.na(y)
  x <- x[!isna]
  y <- y[!isna]
  x[x == Inf | y == Inf] <- NaN
  y[y == Inf | x == Inf] <- NaN

  r <- prcomp(~x + y)
  slope <- r$rotation[2, 1]/r$rotation[1, 1]
  intercept <- r$center[2] - slope*r$center[1]
  abline(intercept, slope, ...)
}

# Helper function for adding legends to plots
addLegend <- function(xc, yc, x, y, weight = NULL, dec = 3) {
  tmp <- summclasses(xc, yc, x, y, weight = weight, dec = dec)
  legend("topleft", legend = paste(names(tmp), "=", tmp))
}

# Function to create weightes for Cohen's weighted kappa
weightFun <- function(n) {
  weight <- matrix(1, n, n)
  weight[ncol(weight), ] <- 0.5
  weight[, nrow(weight)] <- 0.5
  diag(weight) <- 0
  return(weight)
}

# Rename and reorder ABC/GCB factors
reFactor <- function(x) {
  x <- as.character(x)
  x[x == "Unclassified"]  <- "NC"
  factor(x, levels = c("ABC", "NC", "GCB"))
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

if (is.null(rma$reference) || recompute) {

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

studies.vec <- c("CHEPRETRO", "MDFCI", "IDRC", "LLMPPRCHOP", "LLMPPCHOP")
names(studies.vec) <-  c(studies.vec[1:3], "LLMPP R-CHOP", "LLMPP CHOP")


pp <- function(x) {  # Prepend variable name to colnames
  colnames(x) <- paste0(deparse(substitute(x)), ".", colnames(x))
  return(x)
}

merge.by.rownames <- function(x, y) {
  ans <- merge(x, y, by = "row.names", all = TRUE)
  rownames(ans) <- ans$Row.names
  return(subset(ans, select = -Row.names))
}


results <- list()
for (study in studies.vec) {
  message("Creating results for ", study)

  b <- (study != "LLMPPCHOP")

  #
  # ABC/GCB
  #

  cohort          <- data.frame(ABCGCB(rma$cohort[[study]]))
  refbased        <- data.frame(ABCGCB(rma$refbased[[study]]))
  if (b) onebyone <- data.frame(ABCGCB(rma$onebyone[[study]]))

  abcgcb <- merge.by.rownames(pp(cohort), pp(refbased))
  if (b) abcgcb <- merge.by.rownames(abcgcb, pp(onebyone))
  results[["ABCGCB"]][[study]] <- abcgcb

  #
  # BAGS
  #

  results[["BAGS"]][[study]]$cohort   <- BAGS(rma$cohort[[study]])
  results[["BAGS"]][[study]]$refbased <- BAGS(rma$refbased[[study]])
  if (b) results[["BAGS"]][[study]]$onebyone <- BAGS(rma$onebyone[[study]])

  #
  # REGS, the classifier for Cyclophosphamide, Doxorubicin, and Vincristine:
  #

  RC <- ResistanceClassifier
  results[["REGS"]][[study]]$cohort   <- RC(rma$cohort[[study]])
  results[["REGS"]][[study]]$refbased <- RC(rma$refbased[[study]])
  if (b) results[["REGS"]][[study]]$onebyone <- RC(rma$onebyone[[study]])

}
rm(RC)


################################################################################
# Comparisons of the results
################################################################################

#
# Load classifications by Wright's method
#  and add to results-object
#

load("data/metadataCHEPRETRO.RData")
load("data/metadataMDFCI.RData")
load("data/metadataIDRC.RData")
load("data/metadataLLMPPRCHOP.RData")

metadataCHEPRETRO$WrightClass  <- reFactor(metadataCHEPRETRO$WrightClass)
metadataMDFCI$WrightClass      <- reFactor(metadataMDFCI$WrightClass)
metadataIDRC$WrightClass       <- reFactor(metadataIDRC$WrightClass)
metadataLLMPPRCHOP$WrightClass <- reFactor(metadataLLMPPRCHOP$WrightClass)

# Add these results to the results-object!
# Add wright to CHEPRETRO
tmp <- metadataCHEPRETRO[, c("WrightClass", "WrightProb")]
colnames(tmp) <- c("wright.class", "wright.prob")
# check order
tmp1 <- gsub("GSM[0-9]+|_|-|\\(|\\)| ", "", metadataCHEPRETRO$file)
tmp2 <- gsub("GSM[0-9]+|_|-|\\(|\\)| ", "", rownames(results$ABCGCB$CHEPRETRO))
stopifnot(all(tmp1 == tmp2))
rownames(tmp) <- rownames(results$ABCGCB$CHEPRETRO)
results$ABCGCB$CHEPRETRO <- merge.by.rownames(results$ABCGCB$CHEPRETRO, tmp)

# Add wright to MDFCI
tmp <- metadataMDFCI[, "WrightClass", drop = FALSE]
rownames(tmp)[!is.na(metadataMDFCI$HGU133Plus2)] <-
  paste0(na.omit(metadataMDFCI$HGU133Plus2), ".CEL")
colnames(tmp) <- "wright.class"
results$ABCGCB$MDFCI <- merge.by.rownames(results$ABCGCB$MDFCI, tmp)

# Add wright to IDRC
tmp <- metadataIDRC[, "WrightClass", drop = FALSE]
rownames(tmp) <- metadataIDRC$Array.Data.File
colnames(tmp) <- "wright.class"
results$ABCGCB$IDRC <- merge.by.rownames(results$ABCGCB$IDRC, tmp)

# Add wright to LLMPP R-CHOP
tmp <- metadataLLMPPRCHOP[, "WrightClass", drop = FALSE]
rownames(tmp) <- paste0(metadataLLMPPRCHOP$GEO.ID, ".cel")
colnames(tmp) <- "wright.class"
results$ABCGCB$LLMPPRCHOP <- merge.by.rownames(results$ABCGCB$LLMPPRCHOP, tmp)



################################################################################
# ABC/GCB
################################################################################

# TABLE 2 ######################################################################

ABCGCB.weights <- matrix(c(0,1,2,1,0,1,2,1,0), 3)/2

table2 <- matrix(nrow = 4, ncol = 2)
rownames(table2) <- names(studies.vec[-5])
colnames(table2) <- c("Rate of agreement", "Cohen's $\\kappa$")

for (i in seq_len(nrow(table2))) {
  table2[i, ] <-
    with(results[["ABCGCB"]][[studies.vec[i]]],
         testfun(cohort.class, wright.class, weight = ABCGCB.weights))
}

caption <- "Comparison of ABC/GCB classification performed using Wright's
method and the established elastic net classifier based on cohort normalisation
for both. The first column shows the rate of agreement between the classifiers
with $95\\%$ CI. The second column shows the Cohen's $\\kappa$ and $95\\%$ CI."

# Create latex table
w <- latex(table2, file = "tables/table2.tex",
           title = "",
           caption = caption,
           label = "tab:ABCGCBclassifier",
           size = "small")



# TABLE S1 #####################################################################

confuseABCGCB <- list();
for (study in studies.vec[-5]) {
  tmp <- results$ABCGCB[[study]]
  confuseABCGCB[[study]][["wright"]]   <- table(wright = tmp$wright.class,
                                                cohort = tmp$cohort.class)
  confuseABCGCB[[study]][["onebyone"]] <- table(wright = tmp$onebyone.class,
                                                cohort = tmp$cohort.class)
  confuseABCGCB[[study]][["refbased"]] <- table(wright = tmp$refbased.class,
                                                cohort = tmp$cohort.class)
}

tableS1 <- do.call(cbind, lapply(confuseABCGCB, function(l) do.call(rbind, l)))

caption <- "Confusion tables for the ABC/GCB classifiers.
The columns represent cohort based normalisation using the ABC/GCB classifier
based on elastic net.
The first part of the table compares Wright's method for ABC/GCB classification
with the elastic net based.
In the second and third part one-by-one and reference based normalisation is
compared to cohort based normalisation using the ABC/GCB classifier based on
elastic net."

w <- latex(tableS1,
           file = "tables/tableS1.tex",
           title = "",
           cgroup = names(studies.vec[-5]),
           rgroup = c("Wright's method", "One-by-one", "Reference based"),
           size = "footnotesize",
           label = "tab:confusionABCGCBHEMA",
           caption = caption)

# Figure 1 and overview ########################################################

# ABC/GCB Overview
pdf("figures/results_overview_ABCGCB.pdf", width = 7, height = 14)
par(mfrow = c(4,2))
for (study in studies.vec[-5]) {
  res <- results[["ABCGCB"]][[study]]
  with(res, {
    # ONE BY ONE
    plot(logit(cohort.prob), logit(onebyone.prob), main = study)
    abline(0, 1, lty = 2, lwd = 2)
    abline(h = logit(c(0.1,0.9)), v = logit(c(0.1,0.9)), col = "darkgrey",
           lwd = 2, lty = 2)
    addLegend(cohort.class, onebyone.class,
              logit(cohort.prob), logit(onebyone.prob), weight = ABCGCB.weights)
    plotline(logit(cohort.prob), logit(onebyone.prob), col = "red")

    # REF BASED
    plot(logit(cohort.prob), logit(refbased.prob), main = study)
    abline(0, 1, lty = 2, lwd = 2)
    abline(h = logit(c(0.1,0.9)), v = logit(c(0.1,0.9)), col = "darkgrey",
           lwd = 2, lty = 2)
    addLegend(cohort.class, refbased.class,
              logit(cohort.prob), logit(refbased.prob), weight = ABCGCB.weights)
    plotline(logit(cohort.prob), logit(refbased.prob), col = "red")
  })
}
dev.off()



# TABLE 3 ######################################################################

comp.matrix <-   # To be filled
  matrix("-", 4, 3, dimnames = list(studies.vec[-5], c("acc", "kappa", "rho")))
drugs <- c("Cyclophosphamide", "Doxorubicin", "Vincristine", "Combined")

# Add the weights for Cohens Kappa
weights <- list()
weights[["ABCGCB"]] <- ABCGCB.weights
weights[["BAGS"]]   <- weightFun(6)
weights[["REGS"]]   <- matrix(c(0,1,2,1,0,1,2,1,0), 3)/2

subtab <- list()
subtab[["ABCGCB"]][["onebyone"]] <- comp.matrix
subtab[["ABCGCB"]][["refbased"]] <- comp.matrix
subtab[["BAGS"]][["onebyone"]] <- comp.matrix
subtab[["BAGS"]][["refbased"]] <- comp.matrix
subtab[["REGS"]][["onebyone"]] <- comp.matrix
subtab[["REGS"]][["refbased"]] <- comp.matrix

for (study in studies.vec[-5]) {
  # ABC/GCB
  res <- results[["ABCGCB"]][[study]]

  subtab[["ABCGCB"]][["onebyone"]][study, ] <-
    with(res, summclasses(cohort.class, onebyone.class,
                          logit(cohort.prob), logit(onebyone.prob),
                          weight = weights$ABCGCB))
  subtab[["ABCGCB"]][["refbased"]][study, ] <-
    with(res, summclasses(cohort.class, refbased.class,
                          logit(cohort.prob), logit(refbased.prob),
                          weight = weights$ABCGCB))


  # BAGS
  res <- results[["BAGS"]][[study]]

  subtab[["BAGS"]][["onebyone"]][study, -3] <-
    testfun(res$cohort$class, res$onebyone$class, weight = weights$BAGS)
  subtab[["BAGS"]][["refbased"]][study, -3] <-
    testfun(res$cohort$class[names(res$refbased$class)], res$refbased$class,
            weight = weights$BAGS)

  # REGS
  res <- results[["REGS"]][[study]]

  subtab[["REGS"]][["onebyone"]][study, ] <-
    summclasses(unlist(res$cohort$class), unlist(res$onebyone$class),
                unlist(logit(res$cohort$prob)),unlist(logit(res$onebyone$prob)),
                weight = weights$REGS)
  coclass <- res$cohort$class[rownames(res$refbased$class), ]
  coprob  <- res$cohort$prob[rownames(res$refbased$prob), ]
  subtab[["REGS"]][["refbased"]][study, ] <-
    summclasses(unlist(coclass), unlist(res$refbased$class),
                unlist(logit(coprob)),unlist(logit(res$refbased$prob)),
                weight = weights$REGS)
}

table3 <- do.call(rbind, lapply(subtab, function(x) do.call(cbind, x)))
colnames(table3) <- rep(c(colnames(table2), "Pearson's $r$"), 2)
rownames(table3) <- rep(names(studies.vec[-5]), 3)

caption <- "Comparison of classifications obtained using cohort based
normalisation and \\hemaClass{}.
The classifications are compared in terms of accuracy, Cohen's weighted
$\\kappa$, and Pearson's correlation coefficient $r$ all supplied with $95\\%$
CIs. The comparisons in the first and last three columns are based on the
one-by-one normalisation method and the reference based normalisation method,
respectively."

w <- latex(table3,
           file = "tables/table3.tex",
           title = "",
           cgroup = gsub("ABCGCB", "ABC/GCB", names(results)),
           rgroup = c("One-by-one normalisation", "Reference based"),
           size = "scriptsize",
           label = "tab:classALL",
           caption = caption)


# TABLE S2 #####################################################################
# BAGS table

subtab <- list()
for (study in studies.vec[-5]) {
  # BAGS
  res <- results[["BAGS"]][[study]]

  subtab[[study]][["onebyone"]] <-
    table(cohort = res$cohort$class,
          onebyone = res$onebyone$class)
  subtab[[study]][["refbased"]] <-
    table(cohort = res$cohort$class[names(res$refbased$class)],
          refbased = res$refbased$class)
}


tableS2 <- do.call(rbind, lapply(subtab, function(x) do.call(cbind, x)))
abbrev <- c("Naive" = "N", "Centroblast" = "CB", "Centrocyte" = "CC",
            "Memory" = "M",  "Plasmablast" = "PB", "Unclassified" = "UC")
colnames(tableS2) <- abbrev[colnames(tableS2)]

flip <- function(x) {
  ans <- names(x)
  names(ans) <- x
  return(ans)
}

caption <- "Confusion tables for the BAGS classifier.} One-by-one and reference
based normalisation are shown in the columns and cohort normalisation in the
rows."

w <- latex(tableS2,
           file = "tables/tableS2.tex",
           title = "",
           rgroup = flip(studies.vec)[names(subtab)],
           cgroup = c("One-by-one normalisation", "Reference based"),
           size = "small",
           label = "tab:BAGShemaclass",
           caption = caption)




# TABLE S3, TABLE S4 ###########################################################
# REGS supp. tables

drugs <- c("Cyclophosphamide", "Doxorubicin", "Vincristine", "Combined")
subtab <- list()
for (study in studies.vec[-5]) {
  res <- results[["REGS"]][[study]]

  for (drug in drugs) {
    subtab[["onebyone"]][[drug]][[study]] <-
      table(cohort = res$cohort$class[, drug],
            onebyone = res$onebyone$class[, drug])
    subtab[["refbased"]][[drug]][[study]] <-
      table(cohort = res$cohort$class[rownames(res$refbased$class), drug],
            refbased = res$refbased$class[, drug])
  }
}

tableS3 <- do.call(rbind, lapply(subtab$onebyone, function(x) do.call(cbind,x)))
tableS4 <- do.call(rbind, lapply(subtab$refbased, function(x) do.call(cbind,x)))

abbrev2 <- c("Sensitive" = "Sen", "Intermediate" = "Int", "Resistant" = "Res")

colnames(tableS3) <- abbrev2[colnames(tableS3)]
colnames(tableS4) <- abbrev2[colnames(tableS4)]

captionS3 <- "Confusion tables for the REGS classifiers.
One-by-one normalisation are shown in the rows and cohort normalisation in the
columns."

captionS4 <- "Confusion tables for the REGS classifiers.
One-by-one normalisation are shown in the rows and cohort normalisation in the
columns."

w <- latex(tableS3,
           file = "tables/tableS3.tex",
           title = "",
           rgroup = names(subtab$onebyone),
           cgroup = flip(studies.vec)[names(subtab$onebyone[[1]])],
           size = "small",
           label = "tab:confusiondrugonebyone",
           caption = captionS3)

w <- latex(tableS4,
           file = "tables/tableS4.tex",
           title = "",
           rgroup = names(subtab$refbased),
           cgroup = flip(studies.vec)[names(subtab$refbased[[1]])],
           size = "small",
           label = "tab:confusiondrugreference",
           caption = captionS4)


################################################################################

sink(file = "sessionInfo.txt")
print(sessionInfo())
sink()
