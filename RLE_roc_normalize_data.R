##################################
#### Generate data for RLE ROC ###
##################################
## hemaclass.R must be run first
## to download data for the script

library(hemaClass)
setwd("h:/hema_temp/hemaClass-paper")

### Read references for RMA normalization
CHEP=readCHEPRETROreference()
CHOP=readLLMPPCHOPreference()
RCHOP=readLLMPPRCHOPreference()
IDRC=readIDRCreference()
MDFCI=readMDFCIreference()

### CHEPRETRO
if(file.exists("data/chep_rma.rds")){
  chep_rma=readRDS("data/chep_rma.rds")
} else{
  ### Read processed file to get sample names
  chep=readRDS("data/GSE56315/GSE56315_affy.Rds")
  chepFiles=colnames(exprs(chep$es$DLBCL))
  rm(chep)

  ### Choose a random InLab reference
  set.seed(43)
  chep_userRefFiles=sample(chepFiles, 30)
  chep_sampleFiles=chepFiles[!chepFiles%in%chep_userRefFiles]

  ### Read data
  chep_userRef.cel=readCelfiles(file.path("data/GSE56315/",chep_userRefFiles))
  chep_sample.cel=readCelfiles(file.path("data/GSE56315/",chep_sampleFiles))

  ### RMA normalize CHEPRETRO dataset
  chep_userRef.rma=rmaPreprocessing(chep_userRef.cel)
  chep_rma=list()
  
  chep_rma[["cohort"]]=rmaPreprocessing(chep_sample.cel)
  chep_rma[["InLab"]]=rmaReference(chep_sample.cel,chep_userRef.rma)
  chep_rma[["LLMPP CHOP"]]=rmaReference(chep_sample.cel,CHOP)
  chep_rma[["LLMPP R-CHOP"]]=rmaReference(chep_sample.cel,RCHOP)
  chep_rma[["IDRC"]]=rmaReference(chep_sample.cel,IDRC)
  chep_rma[["MDFCI"]]=rmaReference(chep_sample.cel, MDFCI)
  
  saveRDS(chep_rma, "data/chep_rma.rds")
  rm(chep_sample.cel, chep_userRef.cel, chep_userRef.rma, chep_rma)
}


### R-CHOP
if(file.exists("data/RCHOP_rma.rds")){
  RCHOP_rma=readRDS("data/RCHOP_rma.rds")
} else{
  ### Read processed file to get sample names
  LLMPPRCHOP=readRDS("data/GSE10846/GSE10846_affy.Rds")
  id=colnames(exprs(LLMPPRCHOP$es$`R-CHOP`))
  rm(LLMPPRCHOP)
  
  set.seed(42)
  ### Choose random 30 as userRef
  userRef=sample(1:length(id),30)
  userRefID=id[userRef]
  sampleID=id[-userRef]
  
  RCHOP_userRef.cel=readCelfiles(file.path("data/GSE10846/", userRefID))
  RCHOP_SampleFiles=readCelfiles(file.path("data/GSE10846/", sampleID))
  
  ### Normalize with different references
  RCHOP_userRef.rma=rmaPreprocessing(RCHOP_userRef.cel)
  RCHOP_rma=list()
  
  RCHOP_rma[["cohort"]]=rmaPreprocessing(RCHOP_SampleFiles)
  RCHOP_rma[["InLab"]]=rmaReference(RCHOP_SampleFiles, RCHOP_userRef.rma)
  RCHOP_rma[["LLMPPP CHOP"]]=rmaReference(RCHOP_SampleFiles, CHOP)
  RCHOP_rma[["CHEPRETRO"]]=rmaReference(RCHOP_SampleFiles,CHEP)
  RCHOP_rma[["IDRC"]]=rmaReference(RCHOP_SampleFiles,IDRC)
  RCHOP_rma[["MDFCI"]]=rmaReference(RCHOP_SampleFiles,MDFCI)
  
  saveRDS(RCHOP_rma,"data/RCHOP_rma.rds")
  rm(RCHOP_SampleFiles, RCHOP_userRef.cel, RCHOP_userRef.rma, RCHOP_rma)
}

### CHOP
if(file.exists("data/CHOP_rma.rds")){
  RCHOP_rma=readRDS("data/CHOP_rma.rds")
} else{
  ### Read processed file to get sample names
  LLMPPCHOP=readRDS("data/GSE10846/GSE10846_affy.Rds")
  id=colnames(exprs(LLMPPCHOP$es$CHOP))
  rm(LLMPPCHOP)
  
  set.seed(42)
  ### Choose random 30 as userRef
  userRef=sample(1:length(id),30)
  userRefID=id[userRef]
  sampleID=id[-userRef]
  
  CHOP_userRef.cel=readCelfiles(file.path("data/GSE10846/", userRefID))
  CHOP_SampleFiles=readCelfiles(file.path("data/GSE10846/", sampleID))
  
  ### Normalize with different references
  CHOP_userRef.rma=rmaPreprocessing(CHOP_userRef.cel)
  CHOP_rma=list()
  
  CHOP_rma[["cohort"]]=rmaPreprocessing(CHOP_SampleFiles)
  CHOP_rma[["InLab"]]=rmaReference(CHOP_SampleFiles, CHOP_userRef.rma)
  CHOP_rma[["LLMPPP R-CHOP"]]=rmaReference(CHOP_SampleFiles, RCHOP)
  CHOP_rma[["CHEPRETRO"]]=rmaReference(CHOP_SampleFiles,CHEP)
  CHOP_rma[["IDRC"]]=rmaReference(CHOP_SampleFiles,IDRC)
  CHOP_rma[["MDFCI"]]=rmaReference(CHOP_SampleFiles,MDFCI)
  
  saveRDS(CHOP_rma,"data/CHOP_rma.rds")
  rm(CHOP_SampleFiles, CHOP_userRef.cel, CHOP_userRef.rma, CHOP_rma)
}


### IDRC
if(file.exists("data/IDRC_rma.rds")){
  IDRC_rma=readRDS("data/IDRC_rma.rds")
} else{
  IDRC=readRDS("data/GSE31312/GSE31312_affy.Rds")
  IDRCFiles=colnames(exprs((IDRC$es$Batch1)))
  rm(IDRC)
  
  set.seed(42)
  ### Choose random 30 as userRef
  userRef=sample(1:length(IDRCFiles),30)
  userRefID=IDRCFiles[userRef]
  sampleID=IDRCFiles[-userRef]
  
  IDRC_userRef.cel=readCelfiles(file.path("data/GSE31312/",userRefID))
  IDRC_SampleFiles=readCelfiles(file.path("data/GSE31312/",sampleID))
  
  IDRC_userRef.rma=rmaPreprocessing(IDRC_userRef.cel)
  IDRC_rma=list()
  
  IDRC_rma[["cohort"]]=rmaPreprocessing(IDRC_SampleFiles)
  IDRC_rma[["InLab"]]=rmaReference(IDRC_SampleFiles, IDRC_userRef.rma)
  IDRC_rma[["LLMPP CHOP"]]=rmaReference(IDRC_SampleFiles,CHOP)
  IDRC_rma[["LLMPP R-CHOP"]]=rmaReference(IDRC_SampleFiles, RCHOP)
  IDRC_rma[["CHEPRETRO"]]=rmaReference(IDRC_SampleFiles,CHEP)
  IDRC_rma[["MDFCI"]]=rmaReference(IDRC_SampleFiles,MDFCI)
  
  saveRDS(IDRC_rma,"data/IDRC_rma.rds")
  rm(IDRC_SampleFiles, IDRC_userRef.cel, IDRC_userRef.rma, IDRC_rma)
}

### MDFCI
if(file.exists("data/MDFCI_rma.rds")){
  MDFCI_rma=readRDS("data/MDFCI_rma.rds")
} else{
  MDFCI=readRDS("data/GSE34171/GSE34171_affy.Rds")
  MDFCIFiles=colnames(exprs(MDFCI$es$GPL570))
  rm(MDFCI)
  
  set.seed(42)
  ### Choose random 30 as userRef
  userRef=sample(1:length(MDFCIFiles),30)
  userRefID=MDFCIFiles[userRef]
  sampleID=MDFCIFiles[-userRef]
  
  MDFCI_userRef.cel=readCelfiles(file.path("data/GSE34171/",userRefID))
  MDFCI_SampleFiles=readCelfiles(file.path("data/GSE34171/",sampleID))
  
  MDFCI_userRef.rma=rmaPreprocessing(MDFCI_userRef.cel)
  MDFCI_rma=list()
  
  MDFCI_rma[["cohort"]]=rmaPreprocessing(MDFCI_SampleFiles)
  MDFCI_rma[["InLab"]]=rmaReference(MDFCI_SampleFiles, MDFCI_userRef.rma)
  MDFCI_rma[["LLMPP CHOP"]]=rmaReference(MDFCI_SampleFiles,CHOP)
  MDFCI_rma[["LLMPP R-CHOP"]]=rmaReference(MDFCI_SampleFiles, RCHOP)
  MDFCI_rma[["CHEPRETRO"]]=rmaReference(MDFCI_SampleFiles,CHEP)
  MDFCI_rma[["IDRC"]]=rmaReference(MDFCI_SampleFiles,IDRC)
  
  saveRDS(MDFCI_rma,"data/MDFCI_rma.rds")
  rm(MDFCI_rma, MDFCI_userRef.cel, MDFCI_SampleFiles, MDFCI_userRef.rma)
}