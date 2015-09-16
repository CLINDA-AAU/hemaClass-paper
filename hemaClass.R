################################################################################
#
# Master script for the hemaclass paper
#
################################################################################

rm(list = ls())

# Initialization
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

# Load the resave function
source_url(
  "https://raw.githubusercontent.com/AEBilgrau/Bmisc/master/R/resave.R"
)

# Other stuff to come...

################################################################################
# Download and prepare the data
#  + regular RMA normalization
################################################################################

# To get an overview of the data in DLBCLdata run:
data("DLBCL_overview")
# View(DLBCL_overview)
# ls("package:DLBCLdata")

studies <- DLBCL_overview[c(1, 5, 6, 7), 1:4]
for (gse in as.character(studies$GSE)) {
  downloadAndProcessGEO(gse, destdir = "data", cdf = "affy", verbose = verbose)
}


################################################################################
# RMA normalization
#  one-by-one (with and without reference)
################################################################################





################################################################################
# ABC/GCB comparison
################################################################################




################################################################################
# BAGS
################################################################################



################################################################################
# REGS
################################################################################

