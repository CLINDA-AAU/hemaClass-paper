################################################################################
#
# Master script for the hemaclass paper
#
#   Anders Ellern Bilgrau & Steffen Falgreen
#
################################################################################

# NOTE: Running this script needs approx 50 GBs of free disk space

# Initalization
set.seed(220744703) # Set random seed
rm(list = ls())     # Clear global enviroment
# memory.limit(size = 60000)  # If using a Windows machine

# If any of the used packages are missing they will be installed.
pkgs <- c("glmnet", "rgdal", "psych", "devtools", "Hmisc", "shiny",
          "matrixStats", "Rcpp",
          "RcppArmadillo", "RcppEigen", "testthat", "WriteXLS", "RLumShiny",
          "gdata", "Biobase", "affy", "affyio", "preprocessCore",
          "BiocInstaller", "AnnotationDbi", "GEOquery", "GEOquery", "oligo",
          "annotate", "hgu133plus2.db",
          "shinysky", "DLBCLdata", "hemaClass")
missing <- setdiff(pkgs, installed.packages()[, "Package"])

if (length(missing)) {
  # Packages from CRAN
  install.packages(pkgs[1:13])

  # Packages form Bioconductor
  source("http://bioconductor.org/biocLite.R")
  biocLite(pkgs[14:23])

  # Packages from github
  devtools::install_github("AnalytixWare/ShinySky")
  devtools::install_github("AEBilgrau/DLBCLdata")
  devtools::install_github("oncoclass/hemaClass", dependencies = TRUE)
}

# Load needed packages
library("hemaClass")      # Load the hemaclass package
library("DLBCLdata")      # Package for data handling and download
library("devtools")       # For source_url()
library("psych")          # For cohen.kappa()
library("Hmisc")          # For latex()
library("glmnet")         # For fitting ABC/GCB classifier
library("annotate")       # For annotation
library("hgu133plus2.db") # For annotation

# Global variables
recompute <- FALSE
verbose   <- TRUE

twocol    <- 16 / 2.54 # Width of two column figures
onecol    <-  8 / 2.54 # Width of one column figures
maxheight <- 23 / 2.54 # Maximum figure height
pointsize <-  8        # Pointsize in figures


# Load saved data
saved.file <- "saved.RData"
if (file.exists(saved.file)) { load(saved.file) }

################################################################################
# Auxiliary functions / objects
################################################################################

# Colour information
col.data <-
  data.frame(col  = c("#3257A4", "#99991E", "#7F170F", "#71C7D9", "#7EBB46",
                      "#A44392", "#F59E20", "#F3E628", "#E9521C",
                      "red", "blue", "grey"),
             cell = c("Immature", "Pre-BI", "Pre-BII", "Naive","Centroblast",
                      "Centrocyte", "Memory",  "Plasmablast", "Plasmacell",
                      "ABC", "GCB", "Unclassified"),
             abr  = c("I", "PreB1", "PreBII", "N", "CB", "CC", "M",  "PB", "PC",
                      "ABC", "GCB", "UC"),
             order = 1:12, stringsAsFactors = FALSE)
rownames(col.data) <- col.data$abr


# Logit function
logit <- function(p) {
  return(log(p/(1 - p)))
}

# Load the resave function
source_url(
  "https://raw.githubusercontent.com/AEBilgrau/Bmisc/master/R/resave.R"
)

# Function for doing the one-by-one and study based reference normalization
normalizer <- function(study, gse, nsamples = 30, global = FALSE) {
  if (verbose) message("RMA normalizing ", gse, " (", study, ")")

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
    affy.batch.refbased$exprs[, setdiff(colnames(affy.batch$exprs), ref.samples)]

  # Build reference
  rma[["reference"]][[study]] <- rmaPreprocessing(affy.batch.ref)

  # RMA normalize using reference
  rma[["refbased"]][[study]]  <-
    rmaReference(affy.batch.refbased, rma[["reference"]][[study]])$exprs.sc

  rma <<- rma
}

# Format confidence intervals
formatCI <- function(est, ci, dec = 2) {
  est <- sprintf(paste0("%.", dec, "f"), round(est, dec))
  ci  <- sprintf(paste0("%.", dec, "f"), round(ci, dec))
  return(paste0("$", est, "~(", paste(ci, collapse = ", "), ")$"))
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
  tmp <- gsub("\\$|~", "", tmp)
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
  x[x == "Unclassified" | x == "UC"]  <- "NC"
  return(factor(x, levels = c("ABC", "NC", "GCB")))
}

plotCoef <- function(cv.fit,
                     col.chosen = 1, col.discarded = "grey",
                     ylab = NULL, xlab = NULL,
                     max.names = 20, label.pct = 33, labels = TRUE,
                     main = NULL,
                     ...){
  fit <- cv.fit$glmnet.fit
  index <- log(cv.fit$lambda)
  beta.which <- NULL
  COEF <- coef(fit)
  if(class(fit)[1] == "multnet"){
    beta.which <- names(COEF)[1]
    COEF <- COEF[[beta.which]]
  }
  COEF <- COEF[-1,]
  COEF <- COEF[-nrow(COEF), ]
  beta <- as.matrix(COEF[rowSums(COEF) != 0,,drop = FALSE])

  wh <- length(index) - which.min(rev(cv.fit$cvm)) + 1

  COEF.chosen <- rownames(beta[beta[,wh] != 0, , drop = FALSE])
  not.chosen  <- setdiff(names(beta[,wh]), COEF.chosen)
  ylim <- range(beta)
  xlim <- (rev(range(index)))
  xlim2 <- xlim
  if(labels){
    xlim2 <- xlim
    xlim[2] <- (xlim[2]-xlim[1]*label.pct/100) / (1-label.pct/100)
    #xlim[2] <- xlim[2] - abs(xlim[1] - xlim[2]) / l.ratio

    # xlim <- xlim2
  }
  if(is.null(xlab))
    xlab <- expression(paste("Regularization Parameter - log(", lambda, ")"))
  if(is.null(ylab) & !is.null(beta.which))
    ylab <- as.expression(substitute(paste("Coefficients: Response ", a),
                                     list(a= beta.which)))
  if(is.null(ylab))
    ylab <- as.expression(substitute(paste("Coefficients: Response")))


  plot(0, 0,
       xlim = xlim, ylim = ylim, type = "n", axes = FALSE,
       xlab = xlab,
       ylab = ylab)

  get_axp <- function(x) {
    c1 <- ifelse(x[1] < 0, floor(x[1]), ceiling(x[1]))
    c2 <- ifelse(x[2] > 0, floor(x[2]), ceiling(x[2]))
    c(c1, c2)
  }
  get_axp <- function(x) {
    c1 <- ifelse(x[1] < 0, ceiling(x[1]), floor(x[1]))
    c2 <- ifelse(x[2] > 0, ceiling(x[2]), floor(x[2]))
    c(c1, c2)
  }


  #unique((at.x <- round(axTicks(side=1, usr=xlim2,
  #                       axp=c(get_axp(xlim2), n=6), log=FALSE, nintLog=5))))
  at.x <- axTicks(side = 1, usr = xlim2,
                  axp = c(get_axp(xlim2), n = 6), log = FALSE, nintLog = 5)
  at.x <- unique( round(at.x ,1))
  (at.y <- axTicks(side=2))
  axis(side = 1, at = at.x)
  axis(side = 2, las = 2)
  q <- cv.fit$nzero
  q.round <- vector()
  for (a in 1:length(at.x)) {
    q.round[a] <- round(q)[which.min(abs(index - at.x[a]))]
  }
  #axis(side = 3, at = index[seq(1, length(index), length.out = 5)],
  #     labels = round(q)[seq(1, length(index), length.out = 5)])
  axis(side = 3, at = at.x,
       labels = q.round)

  graphics:::box()

  grid2 <- function(xlim, ylim, at.y, at.x, lwd = par("lwd"),
                    col = "grey", lty = "dotted"){

    for(y in at.y)
      segments(xlim[1], y, xlim[2], y, lty = lty,
               col = col, lwd = lwd)
    for(x in at.x)
      segments(x, ylim[1], x, ylim[2], at.y, lty = lty,
               col = col, lwd = lwd)
  }

  grid2(xlim = c(xlim2[1] + 1, xlim2[2]),
        ylim = c(ylim[1]-1, ylim[2]+1),
        at.y = at.y, at.x = at.x )


  if(length(not.chosen) > 0)
    matplot(index, t(beta[not.chosen,  , drop = FALSE]), lty = 1,
            type = "l",col = col.discarded, add = TRUE)# ...)
  if(length(COEF.chosen) > 0)
    matplot(index, t(beta[COEF.chosen, , drop = FALSE]), lty = 1,
            type = "l", col = col.chosen, add = TRUE)#, ...)

  graphics:::box()

  abline(v = index[wh], lty = 2, col = 2)
  if(!is.null(COEF.chosen)){
    if(length(COEF.chosen) == 1){
      which2 <- COEF.chosen
    }else{
      which2 <- names(sort(abs(beta[COEF.chosen, wh]),decreasing=TRUE))
    }
    if(length(which2) > max.names)
      which2 <- names(sort(abs(beta[COEF.chosen, wh]),
                           decreasing=TRUE)[1:max.names])

    genes <- unlist(lookUp(which2, "hgu133plus2", "SYMBOL"))

    if(any(is.na(genes)))
      genes[is.na(genes)] <- names(genes[is.na(genes)])

    if(labels){
      y1.cords <- beta[which2, ncol(beta)]
      names(y1.cords) <- which2
      y2.cords <- seq(min(ylim), max(ylim), length.out = length(which2))
      if(length(which2) == 1)
        y2.cords <- y1.cords
      names(y2.cords) <- names(y1.cords)[order(y1.cords)]
      x1.cords <- min(index)
      x2.cords <- min(index) - (xlim[1] - xlim[2])/10

      for(i in which2){
        segments(x1.cords, y1.cords[i], x2.cords, y2.cords[i], lty = 1,
                 col = "grey", lwd = par("lwd"))
        #if(labels & !is.null(array))
        #  text(x2.cords, y2.cords[i], genes[i], pos = 4)
        #if(labels & is.null(array))
        text(x2.cords, y2.cords[i], genes[i], pos = 4, cex = 0.8)
      }
    }
  }
  if (is.null(main)) {
    main <- paste("Regularization Curves for", beta.which)
  }
  title(main, font.main = 1, line = 3)
}

plotCV <- function(x, y, lo, up,
                   xlab = NULL, ylab = NULL,
                   col = "#333333", log = "y",
                   col.alpha = 20, ...){

  if(grepl("y", log)){
    y <- log(y)
    lo <- log(lo)
    up <- log(up)
  }
  if(grepl("x", log)){
    x <- log(x)
  }

  matplot(x, cbind(y, lo, up), col = 0, ylab = ylab, xlab = xlab,...)


  polygon(x= c(x, rev(x)),
          y = c(lo, rev(up)),
          col = paste(col, col.alpha, sep = ""), border = 0)
  lines(x, y, lwd= 2, col = col)
}

plotCValpha <- function(x, y, lo, up, nzero, parameter = "alpha",
                        xlab = NULL,
                        ylab = NULL,
                        main = NULL,
                        main1 = NULL, ...){
  if(is.null(main))
    if(is.null(main1)){
      main <- as.expression(
        substitute(paste("Cross validation: ", alpha, " = ", b) ,
                   list(b= round(a[which.min(mat.list$AUC[,c(1)])],2))))
    }else{
      main <- as.expression(
        substitute(paste(m," (CV: ", alpha, " = ", b, ")") ,
                   list(m = main1, p = parameter, b= round(x[which.min(y)],2))))
    }
  if(is.null(xlab))
    xlab <- expression(paste("Model parameter: ", alpha, sep = ""))

  if(is.null(ylab))
    ylab <- "Mean-Squared Error"

  matplot(x, cbind(y, lo, up), col = 0,
          xlab = xlab, ylab = ylab )


  polygon(x= c(x, rev(x)),
          y = c(lo, rev(up)),
          col = "#33333320", border = 0)
  lines(x, y, lwd= 2, col = "#333333")

  abline(v=x[which.min(y)], lty = 2)


  axis(side = 3, at = x[seq(1, length(x), length.out = 5)],
       labels = round(nzero)[seq(1, length(x), length.out = 5)])
  title(main, line = 3)
}

plotCVlambda <- function(x, y, lo, up, nzero, parameter = "alpha",
                         xlab = NULL,
                         ylab = NULL,
                         main = NULL,
                         main1 = NULL, ...){
  if(is.null(main))
    if(is.null(main1)){
      main <- as.expression(substitute(
        paste("Cross validation: ", lambda, " = ",
              b) , list(b= round(a[which.min(mat.list$AUC[,c(1)])],2))))
    }else{
      main <- as.expression(substitute(paste(m," (CV: log(", lambda, ") = ",
                                             b, ")") ,
                                       list(m = main1, p = parameter,
                                            b= round((x[which.min(y)]),2))))
    }
  if(is.null(xlab))
    xlab <- expression(paste("Penalty: log(", lambda, ")", sep = ""))

  if(is.null(ylab))
    ylab <- "Mean-Squared Error"

  matplot(x, cbind(y, lo, up), col = 0,
          xlab = xlab, ylab = ylab )


  polygon(x= c(x, rev(x)),
          y = c(lo, rev(up)),
          col = "#33333320", border = 0)
  lines(x, y, lwd= 2, col = "#333333")
  min <- length(y) -which.min(rev(y))
  abline(v=x[min], lty = 2)


  axis(side = 3, at = x[seq(1, length(x), length.out = 5)],
       labels = round(nzero)[seq(1, length(x), length.out = 5)])

  title(main, line = 3)
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
  set.seed(1154538003)

  # First the overall reference is made using the LLMPP CHOP data
  normalizer(study = "LLMPPCHOP", gse = "GSE10846", global = TRUE); gc()

  # Next the other datasets are nomalized according to LLMPPCHOP
  # and the 30 sample study based reference
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
# Load classifications by Wright's naive Bayes method
#  and add to results-object
#

# Add these results to the results-object!
# Add wright to CHEPRETRO
metadataCHEPRETRO <-
  dat$GSE56315$metadata[, "characteristics_ch1.5", drop = FALSE]
names(metadataCHEPRETRO) <- "wright.class"
levels(metadataCHEPRETRO$wright.class) <- c("", "ABC", "GCB", "Unclassified")
metadataCHEPRETRO <- subset(metadataCHEPRETRO, wright.class != "")
metadataCHEPRETRO$wright.class <- reFactor(metadataCHEPRETRO$wright.class)

# check order
tmp <- metadataCHEPRETRO
tmp1 <- rownames(metadataCHEPRETRO)
tmp2 <- gsub("(GSM[0-9]+)_.+", "\\1", rownames(results$ABCGCB$CHEPRETRO))
stopifnot(all(tmp1 == tmp2))  # Make sure the sorting is the same
rownames(tmp) <- rownames(results$ABCGCB$CHEPRETRO)

results$ABCGCB$CHEPRETRO <- merge.by.rownames(results$ABCGCB$CHEPRETRO, tmp)

# Add wright to MDFCI
metadataMDFCI <- read.csv(file = "data/metadataMDFCI.csv", row.names = 1)
metadataMDFCI$wright.class <- reFactor(metadataMDFCI$wright.class)
tmp <- metadataMDFCI[, "wright.class", drop = FALSE]
notna <- !is.na(metadataMDFCI$HGU133Plus2)
tmp <- tmp[notna, , drop = FALSE]
rownames(tmp) <- paste0(metadataMDFCI$HGU133Plus2[notna], ".CEL")

results$ABCGCB$MDFCI <- merge.by.rownames(results$ABCGCB$MDFCI, tmp)

# Add wright to IDRC
metadataIDRC <- dat$GSE31312$metadata[, "characteristics_ch1", drop = FALSE]

names(metadataIDRC) <- "wright.class"
metadataIDRC$wright.class <- gsub("gene expression profiling subgroup: ", "",
                                 metadataIDRC$wright.class)
metadataIDRC$wright.class <- reFactor(metadataIDRC$wright.class)

tmp <- metadataIDRC[, "wright.class", drop = FALSE]
rownames(tmp) <- paste0(rownames(tmp), ".CEL")
results$ABCGCB$IDRC <- merge.by.rownames(results$ABCGCB$IDRC, tmp)

# Add wright to LLMPP R-CHOP
metadataLLMPPRCHOP <- dat$GSE10846$metadata[, "WrightClass", drop = FALSE]
metadataLLMPPRCHOP$wright.class <- reFactor(metadataLLMPPRCHOP$WrightClass)

tmp <- metadataLLMPPRCHOP[, "wright.class", drop = FALSE]
rownames(tmp) <- gsub(".CEL", ".cel", rownames(metadataLLMPPRCHOP))
results$ABCGCB$LLMPPRCHOP <- merge.by.rownames(results$ABCGCB$LLMPPRCHOP, tmp)

# Remove
na <- rowSums(!is.na(subset(results$ABCGCB$LLMPPRCHOP, select = -wright.class)))
results$ABCGCB$LLMPPRCHOP <- results$ABCGCB$LLMPPRCHOP[na != 0, ]



################################################################################
# Create output
################################################################################

# TABLE 2 ######################################################################

ABCGCB.weights <- matrix(c(0,1,2,1,0,1,2,1,0), 3)/2

table2 <- matrix(nrow = 4, ncol = 2)
rownames(table2) <- names(studies.vec[-5])
colnames(table2) <- c("Agreement", "Cohen's $\\kappa$")

for (i in seq_len(nrow(table2))) {
  table2[i, ] <-
    with(results[["ABCGCB"]][[studies.vec[i]]],
         testfun(cohort.class, wright.class, weight = ABCGCB.weights))
}

caption <- "Comparison of ABC/GCB classification performed using Wright's naive
Bayes classifier \\citep{Wright2003} and the established elastic net classifier
based on cohort normalization
for both. The second column shows the rate of agreement (accuracy) between the
classifiers with $95\\%$ CI. The third column shows the Cohen's weighted
$\\kappa$ with $95\\%$ CI."

# Create latex table
w <- latex(table2, file = "tables/table2.tex",
           title = "Dataset",
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
The columns represent cohort based normalisztion using the ABC/GCB classifier
based on elastic net.
The first part of the table compares Wright's method for ABC/GCB classification
with the elastic net based.
In the second and third part one-by-one and reference based normalization is
compared to cohort based normalization using the ABC/GCB classifier based on
elastic net."

w <- latex(tableS1,
           file = "tables/tableS1.tex",
           title = "",
           cgroup = names(studies.vec[-5]),
           rgroup = c("Wright's method", "One-by-one", "Reference based"),
           size = "footnotesize",
           label = "tab:confusionABCGCBHEMA",
           caption = caption)

# ABC/GCB Overview ############################################################

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
normalization and one-by-one normalization as well as reference based
normalization.
The classifications are compared in terms of rate of argreement (accuracy),
Cohen's weighted $\\kappa$, and Pearson's correlation coefficient $r$ all
supplied with $95\\%$ CIs. The comparisons in the first and last three columns
are based on the one-by-one normalization method and the reference based
normalization method, respectively."

# Remove leading zeros:
table3 <- gsub("0\\.", ".", table3)
table3 <- gsub("1\\.0", "1.", table3)

w <- latex(table3,
           file = "tables/table3.tex",
           title = "",
           rgroup = gsub("ABCGCB", "ABC/GCB", names(results)),
           cgroup = c("One-by-one pre-processing", "Reference based pre-processing"),
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

caption <- "Confusion tables for the BAGS classifier. One-by-one and reference
based normalization are shown in the columns and cohort normalization in the
rows."

w <- latex(tableS2,
           file = "tables/tableS2.tex",
           title = "",
           rgroup = flip(studies.vec)[names(subtab)],
           cgroup = c("One-by-one normalization", "Reference based normalization"),
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
One-by-one normalization are shown in the rows and cohort normalization in the
columns."

captionS4 <- "Confusion tables for the REGS classifiers.
Reference based normalization are shown in the rows and cohort normalization in
the columns. Note, 30 samples were used as reference data and hence not present
in this table."

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


# FIGURE 2 FIGURE 3 ############################################################

subplot <- function(x, y,
                    cut.x = NULL,
                    cut.y = cut.x,
                    col1 = "#90EE9060",
                    col2 = "#FFA0A060",
                    ...) {
  tmp <- c(x, y)
  tmp <- tmp[is.finite(tmp)]
  rng <- range(tmp, na.rm = TRUE)
  plot(x, y, type = "n", xlab = "", ylab = "", axes = FALSE,
       xlim = rng, ylim = rng,
       ...)

  grid()

  rect(rng[1] - abs(rng[1]*10), rng[1] - abs(rng[1]*10),
       rng[2] + abs(rng[2]*10), rng[2] + abs(rng[2]*10), col = "grey95")

  if(!is.null((cut.x))){
    rect(max(cut.x), max(cut.y),  100,  100, col = col1, border = NA)
    rect(min(cut.x), min(cut.y), -100, -100, col = col1, border = NA)
    rect(min(cut.x), max(cut.y), -100,  100, col = col2, border = NA)
    rect(max(cut.x), min(cut.y),  100, -100, col = col2, border = NA)
  }

  abline(0, 1, lty = 2, lwd = 2, col = "darkgrey")

  axis(1)
  axis(2)
  box()

  points(x, y, pch = 16, col = "#000000A0")

  cc <- cor.test(x, y)
  rval <- formatCI(cc$estimate, cc$conf.int, dec = 4)
  rval <- gsub("~|\\$", " ", rval)
  leg <- bquote(italic(r) == .(rval))
  legend("bottomright", legend = leg, bty = "n")
  plotline(x, y, col = "black", lwd = 1)
}

myplot <- function(x1, y1, x2, y2, panel = c("A", "B"),
                   cut.x = NULL, cut.y1 = cut.x, cut.y2 = cut.y1, ...) {
  # ONE BY ONE
  subplot(x1, y1, cut.x = cut.x, cut.y = cut.y1, ...)
  mtext(panel[1], font = 2, adj = -0.1, line = 0.5, cex = 1.2)

  # REF BASED
  subplot(x2, y2, cut.x = cut.x, cut.y = cut.y2, ...)
  mtext(panel[2], font = 2, adj = -0.1, line = 0.5, cex = 1.2)
}



#
# Plot FIGURE 2
# ABC/GCB + REGS
#

f <- 0.5
pdf("figures/figure2.pdf", height = 5*7*f, width = 2*7*f)
{

  par(mfrow = c(5, 2), mar = c(2,4,2.2,0.2) + 0.1, oma = c(2,0,0,0))

  # ABC/GCB:
  res <- results$ABCGCB$CHEPRETRO
  x <-  logit(res$cohort.prob)
  y1 <- logit(res$onebyone.prob)
  y2 <- logit(res$refbased.prob)
  myplot(x1 = x, y1, x2 = x, y2, panel = c("A", "B"), main = "ABC/GCB",
         cut.x = logit(c(0.1, 0.9)), asp = 1)

  # REGS
  res <- results$REGS$CHEPRETRO
  i <- 3
  for (drug in drugs) { # skal vendes om
    x1 <- logit(res$cohort$prob[, drug])
    y1 <- logit(res$onebyone$prob[, drug])
    y2 <- logit(res$refbased$prob[, drug])
    x2 <- x1[names(y2)]
    cut.x <- res$cohort$cut[[drug]]
    cut.y <- res$refbased$cut[[drug]]
    stopifnot(all.equal(cut.y, res$onebyone$cut[[drug]]))
    myplot(x1, y1, x2, y2, panel = LETTERS[i + 0:1],
           main = paste0("REGS (", drug, ")"),
           cut.x = cut.x, cut.y1 = cut.y, asp = 1)
    i <- i + 2
  }

  mtext("Cohort based classification", side = 1, outer = TRUE)
  mtext("One-by-one based classification", side = 2, outer = TRUE, line = -1.5)
  mtext("Reference based classification", side = 2, outer = TRUE, line = -28)
}
dev.off()


#
# Plot FIGURE 3
# BAGS
#

pdf("figures/figure3.pdf", height = 5*7*f, width = 2*7*f)
{

  par(mfrow = c(5, 2), mar = c(2,4,2.2,0.2) + 0.1, oma = c(2,0,0,0))

  # BAGS
  res <- results$BAGS$CHEPRETRO
  i <- 1
  for (subtype in names(abbrev)[-6]) {
    x1 <- logit(res$cohort$prob.mat[, subtype])
    y1 <- logit(res$onebyone$prob.mat[, subtype])
    y2 <- logit(res$refbased$prob.mat[, subtype])
    x2 <- x1[names(y2)]

    cut.x <- c(-100, res$cohort$cut)
    cut.y1 <- c(-100, res$refbased$cut)
    cut.y2 <- c(-100, res$onebyone$cut)

    myplot(x1, y1, x2, y2, panel = LETTERS[i + 0:1],
           main = subtype,
           col1 = paste0(col.data$col, "60")[col.data$cell == subtype],
           col2 = "White", asp = 1,
           cut.x = cut.x, cut.y1 = cut.y1, cut.y2 = cut.y2)
    i <- i + 1
  }

  mtext("Cohort based classification", side = 1, outer = TRUE)
  mtext("One-by-one based classification", side = 2, outer = TRUE, line = -1.5)
  mtext("Reference based classification", side = 2, outer = TRUE, line = -28)
}
dev.off()


################################################################################

################################################################################
# Comparisons of the pre-processing results
################################################################################

for (type in c("ABCGCB", "BAGS", "REGS")) {

  pdf(paste0("figures/figure4", type,".pdf"), height = 4*7*f, width = 2*7*f)


  par(mfrow = c(4, 2), mar = c(2,4,2.2,0.2) + 0.1, oma = c(2,0,0,0))

  if (type == "ABCGCB") {
    probe.order <- rownames(hemaClass::readABCGCBCoef())[-1]
  }
  if (type == "BAGS") {
    probe.order <- rownames(hemaClass::readBAGSCoef())[-1]
  }
  if (type == "REGS") {
    probe.order <- rownames(hemaClass::readClasCoef())[-1]
  }
  sample <- colnames(rma[["refbased"]][["CHEPRETRO"]])[1]

  x <-  microarrayScale(exprs(dat$GSE56315$es$DLBCL))[probe.order ,sample]
  y1 <- rma[["onebyone"]][["CHEPRETRO"]][probe.order, sample]
  y2 <- rma[["refbased"]][["CHEPRETRO"]][probe.order, sample]

  myplot(x, y1, x, y2, panel = c("A", "B"),
         main = "CHEPRETRO", asp = 1)

  sample <- colnames(rma[["refbased"]][["MDFCI"]])[1]

  x1 <-  microarrayScale(exprs(dat$GSE34171$es$GPL570))[probe.order ,sample]
  y1 <- rma[["onebyone"]][["MDFCI"]][probe.order, sample]
  y2 <- rma[["refbased"]][["MDFCI"]][probe.order, sample]
  x2 <- x1[names(y2)]

  myplot(x1, y1, x2, y2, panel = c("C", "D"),
         main = "MDFCI", asp = 1)

  sample <- colnames(rma[["refbased"]][["IDRC"]])[1]

  x1 <-  microarrayScale(exprs(dat$GSE31312$es$Batch1))[probe.order ,sample]
  y1 <- rma[["onebyone"]][["IDRC"]][probe.order, sample]
  y2 <- rma[["refbased"]][["IDRC"]][probe.order, sample]
  x2 <- x1[names(y2)]

  myplot(x1, y1, x2, y2, panel = c("E", "F"),
         main = "IDRC", asp = 1)


  # Compare onebyone and refbased pre-processing to cohort for LLMPP R-CHOP
  sample <- colnames(rma[["refbased"]][["LLMPPRCHOP"]])[1]

  x1 <-  microarrayScale(exprs(dat$GSE10846$es$`R-CHOP`))[probe.order ,sample]
  y1 <- rma[["onebyone"]][["LLMPPRCHOP"]][probe.order, sample]
  y2 <- rma[["refbased"]][["LLMPPRCHOP"]][probe.order, sample]
  x2 <- x1[names(y2)]

  myplot(x1, y1, x2, y2, panel = c("G", "H"),
         main = "LLMPP R-CHOP", asp = 1)

  mtext("Cohort based Pre-processing", side = 1, outer = TRUE)
  mtext("One-by-one based Pre-processing", side = 2, outer = TRUE, line = -1.5)
  mtext("Reference based Pre-processing", side = 2, outer = TRUE, line = -28)


  dev.off()
}

################################################################################
# Establish the ABC/GCB classifier
################################################################################

# In order to make the LLMPP CHOP data comparable to other datasets
# it is median centered.
GEPLLMPPCHOP.med <- microarrayScale(exprs(dat$GSE10846$es$CHOP), scale = FALSE)
file.ABCGCB.class.cv <- file.path("GeneratedData", "ABCGCB.class.cv.RData")

alphas <- seq(0.1, 1, 0.025)

if (file.exists(file.ABCGCB.class.cv)) {
  load(file.ABCGCB.class.cv)
} else {

  dir.create(dirname(file.ABCGCB.class.cv),
             showWarnings = FALSE, recursive = TRUE)
  mat <- matrix(NaN, ncol = 4, nrow = length(alphas))
  colnames(mat) <- c("ABCGCB", "ABCGCB.lo", "ABCGCB.up", "ABCGCB.n")
  row.names(mat) <- alphas

  ABCGCB.fit.cv <- list()

  min.error <- Inf

  metadata <- dat$GSE10846$metadata[
    dat$GSE10846$metadata$microarray.diagnosis %in% c("ABC", "GCB") &
      dat$GSE10846$metadata$chemo == "CHOP" , ]

  metadata$microarray.diagnosis <-
    factor(metadata$microarray.diagnosis, levels = c("GCB", "ABC"))

  set.seed(1000)

  sample <- sample(rep(seq(10), length = nrow(metadata)))

  keep <- rownames(GEPLLMPPCHOP.med)[!grepl("AFFX", rownames(GEPLLMPPCHOP.med))]

  x <- t((GEPLLMPPCHOP.med[keep, paste0(metadata$GEO.ID, ".cel")]))
  y <- metadata$microarray.diagnosis
  for (alpha in alphas) {
    cat(alpha, "\n")

    fit <- cv.glmnet(x, y,
                     standardize = FALSE, family = "binomial",
                     foldid = sample,
                     lambda = exp(seq(-10, 2, length.out = 200)),
                     alpha  = alpha, type.measure = "deviance",
                     keep   = TRUE)

    wh.min <- length(fit$cvm) - which.min(rev(fit$cvm))
    mat[paste(alpha), "ABCGCB"] <- min(fit$cvm)
    mat[paste(alpha), paste("ABCGCB", ".lo", sep = "")] <-
      fit$cvlo[wh.min]
    mat[paste(alpha), paste("ABCGCB", ".up", sep = "")] <-
      fit$cvup[wh.min]
    mat[paste(alpha), paste("ABCGCB", ".n", sep = "")] <-
      fit$nzero[wh.min]

    if (min(fit$cvm) < min.error) {
      cat("new.min =", min(fit$cvm), "\n")
      min.error <- min(fit$cvm)
      ABCGCB.fit <- fit
    }

  }
  LLMPPCHOP.sd <- matrixStats::rowSds(GEPLLMPPCHOP.med[keep, metadata$GEO.ID])
  names(LLMPPCHOP.sd ) <- keep
  save(ABCGCB.fit, mat, sample, LLMPPCHOP.sd, file = file.ABCGCB.class.cv)
}


# Plot the result

pdf(file.path("figures/figureS1.pdf"),
    pointsize = pointsize ,
    width = twocol, height = twocol/3)
{
  par(mfrow = c(1, 3))
  plotCValpha(x  = as.numeric(row.names(mat)),
              y  = mat[,"ABCGCB"],
              lo = mat[, paste("ABCGCB", ".lo", sep = "")],
              up = mat[, paste("ABCGCB", ".up", sep = "")],
              nzero = mat[,paste("ABCGCB", ".n", sep = "")],
              main1 = "ABC/GCB",
              ylim  = c(0,0.08),
              ylab  = "Logistic Deviance")
  mtext(LETTERS[1], 3, line = 1.5, adj = -0.145, cex = 10/pointsize)


  plotCVlambda(x  = log(ABCGCB.fit$lambda),
               y  = ABCGCB.fit$cvm,
               lo = ABCGCB.fit$cvlo,
               up = ABCGCB.fit$cvup,
               nzero = ABCGCB.fit$nzero,
               main1 = "ABC/GCB",
               ylab  = "Logistic Deviance")
  mtext(LETTERS[2], 3, line = 1.5, adj = -0.145, cex = 10/pointsize)

  plotCoef(ABCGCB.fit,
           label.pct = 35,
           main = paste("Regularization Curves"),
           ylab = "Coeffecients: ABC")
  mtext(LETTERS[3], 3, line = 1.5, adj = -0.145, cex = 10/pointsize)
}
dev.off()

# The coefficients are created for the hemaclass package
coef <- coef(ABCGCB.fit, s = "lambda.min")
ABCGCB.coef <- as.matrix(coef)[as.matrix(coef) != 0, , drop = FALSE]
ABCGCB.coef <- ABCGCB.coef * c(1, LLMPPCHOP.sd[rownames(ABCGCB.coef)[-1]])
nrow(ABCGCB.coef) - 1

# The established coefficients are saved
dir.create("hemaclass.org/ABCGCB", showWarnings = FALSE, recursive = TRUE)
saveRDS(ABCGCB.coef, file = "hemaclass.org/ABCGCB/ABCGCB.coef.rds")

################################################################################
# Write session info
################################################################################

sink(file = "sessionInfo.txt")
  print(sessionInfo())
sink()
