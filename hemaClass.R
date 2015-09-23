################################################################################
#
# Master script for the hemaclass paper
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

# # Load the resave function
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

  rma[["ref.samples"]][[study]] <- sample(colnames(affy.batch$exprs), nsamples)
  affy.batch.ref <- affy.batch
  affy.batch.ref$exprs <- affy.batch.ref$exprs[, rma$ref.samples[[study]]]

  affy.batch.refbased <- affy.batch
  affy.batch.refbased$exprs <-
    affy.batch.refbased$exprs[, setdiff(colnames(affy.batch$exprs),
                                        rma$ref.samples[[study]])]

  rma[["reference"]][[study]] <- rmaPreprocessing(affy.batch.ref)
  rma[["refbased"]][[study]]  <-
    rmaReference(affy.batch.ref, rma[["reference"]][[study]])$exprs.sc

  rma <<- rma
}

#function for comparing results
testfun <- function(x, y, dec = 2, weight){
  a <- table(x, y)
  acc <- binom.test(c(sum(diag(a)),sum(a)-sum(diag(a))))
  acc <- paste(round(acc$estimate, dec),
               " (", paste(round(acc$conf.int, dec), collapse = ","), ")", sep = "")

  kap <- cohen.kappa(data.frame(x,y), w=weight)
  conf.1 <- kap$confid[2,1]
  conf.2 <- min(kap$confid[2,3],1)
  kappa <- paste(round(kap$weighted.kappa, dec),
                 " (", round(conf.1, dec), ",", round(conf.2, dec), ")", sep = "")

  c(acc, kappa)
}

weightFun <- function(n){
  weight <- matrix(1, n, n)
  weight[ncol(weight), ] <- 0.5
  weight[, nrow(weight)] <- 0.5
  diag(weight) <- 0
  return(weight)
}

testfunCor <- function(x, y){

  logit <- function(p) log(p/(1-p))

  cort <- cor.test(logit(x), logit(y))

  leg <- paste0(round(cort$estimate, 2),
                " (", paste(round(cort$conf.int, 2), collapse = ","), ")")
  return(leg)
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
  # and a 30 sample study based reference
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
names(studies.vec) <- c("CHEPRETRO", "MDFCI", "IDRC", "LLMPP R-CHOP", "LLMPP CHOP")

results <- list()
for (study in studies.vec) {

  if (study != "LLMPPCHOP") {

    results[["ABCGCB"]][[study]]$onebyone <- ABCGCB(rma$onebyone[[study]])
    results[["BAGS"]][[study]]$onebyone   <- BAGS(rma$onebyone[[study]])
    results[["REGS"]][[study]]$onebyone   <- ResistanceClassifier(rma$onebyone[[study]])

  }

  # ABC/GCB
  results[["ABCGCB"]][[study]]$refbased <- ABCGCB(rma$refbased[[study]])
  results[["ABCGCB"]][[study]]$cohort   <- ABCGCB(rma$cohort[[study]])

  # BAGS
  results[["BAGS"]][[study]]$refbased <- BAGS(rma$refbased[[study]])
  results[["BAGS"]][[study]]$cohort   <- BAGS(rma$cohort[[study]])

  # REGS
  # The classifier for Cyclophosphamide, Doxorubicin, and Vincristine:
  results[["REGS"]][[study]]$refbased <- ResistanceClassifier(rma$refbased[[study]])
  results[["REGS"]][[study]]$cohort   <- ResistanceClassifier(rma$cohort[[study]])

}



################################################################################
# Comparisons of the results
################################################################################
comparison <- list()

comp.matrix <- matrix("-", nrow = length(studies.vec), ncol = 3,
                      dimnames = list(studies.vec,
                                      c("Accuracy", "Kappa", "rho")))

# Add the weights for Cohens Kappa
comparison[["ABCGCB"]]$weight <- matrix(c(0,1,2,1,0,1,2,1,0), 3) / 2
comparison[["BAGS"]]$weights  <- weightFun(6)
comparison[["REGS"]]$weights  <- weightFun(3)

comparison[["ABCGCB"]][["onebyone"]]$result <- comp.matrix
comparison[["ABCGCB"]][["refbased"]]$result <- comp.matrix

comparison[["BAGS"]][["onebyone"]]$result <- comp.matrix
comparison[["BAGS"]][["refbased"]]$result <- comp.matrix

comparison[["REGS"]][["onebyone"]]$result <- comp.matrix
comparison[["REGS"]][["refbased"]]$result <- comp.matrix


for(study in studies.vec){
  comparison[["ABCGCB"]][["onebyone"]][[study]]<- table(c("ABC", "NC", "GCB"), c("ABC", "NC", "GCB"))
  comparison[["ABCGCB"]][["refbased"]][[study]] <- table(c("ABC", "NC", "GCB"), c("ABC", "NC", "GCB"))

  comparison[["BAGS"]][["onebyone"]][[study]] <- table(c("N", "CB", "CC", "M", "PB"), c("N", "CB", "CC", "M", "PB"))
  comparison[["BAGS"]][["refbased"]][[study]] <- table(c("N", "CB", "CC", "M", "PB"), c("N", "CB", "CC", "M", "PB"))

  for(drug in c("Cyclophosphamide", "Doxorubicin", "Vincristine", "Combined")){
    comparison[["REGS"]][["onebyone"]][[study]][[drug]] <- table(c("Sensitive", "Intermediate", "Resistant"), c("Sensitive", "Intermediate", "Resistant"))
    comparison[["REGS"]][["refbased"]][[study]][[drug]] <- table(c("Sensitive", "Intermediate", "Resistant"), c("Sensitive", "Intermediate", "Resistant"))
  }
}










sink(file = "sessionInfo.txt")
print(sessionInfo())
sink()
