require(knitr)
library(hemaClass)
library(ggbiplot)
library(matrixStats)
library(pROC)
library(xtable)

### Plot output
tiff_out=FALSE
outDir="H:/hema_temp/hemaClass-paper/figures"
tableOut="H:/hema_temp/hemaClass-paper/tables"

### Function for plotting RLE values
plot_RLE=function(rma_data, statistic){
  if(!statistic %in% c("mean", "median", "iqr","total")){
    cat('error statistic must be "mean", "median" or "iqr"')
  } else{
      ## List references
      refs=names(rma_data)
      nRef=length(refs)
      
      ## Set column number in RLE.stats matrix
      if( statistic=="mean" ){ i=1 }
      if( statistic=="median" ){ i=2 }
      if( statistic=="iqr" ){ i=3 }
      if( statistic=="total"){i=5
        for(ref in refs){
          rma_data[[ref]]$RLE.stats[,i]=abs(rma_data[[ref]]$RLE.stats[,2])+rma_data[[ref]]$RLE.stats[,3]
        }
      }
      
      ## Find max and min
      temp=c()
      for(ref in refs){
        temp=c(temp,abs(rma_data[[ref]]$RLE.stats[,i]))
      }
      minY=min(temp)
      maxY=max(temp)+0.2*(max(temp)-min(temp))
      
      ## Plot RLE values
      plot(abs(rma_data[[refs[1]]]$RLE.stats[,i]), ylim=c(minY,maxY), pch=16, ylab=paste("RLE",statistic), xlab="Sample no.")
      #abline(a=quantile(abs(rma_data[[refs[1]]]$RLE.stats[,i]), 0.025), b=0)
      #abline(a=quantile(abs(rma_data[[refs[1]]]$RLE.stats[,i]), 0.975), b=0)

      if(nRef>1){
        for(j in 2:nRef){
          points(abs(rma_data[[refs[j]]]$RLE.stats[,i]), col=j, pch=16)
        }
      }
      legend("top", refs, col=1:nRef, pch=rep(20,nRef), cex=0.8, pt.cex=1, bty="n", ncol=3)
    }
}

calc_roc=function(rma_data,value,statistic){
  ## List references
  refs=names(rma_data)[-1]
  nRef=length(refs)
  
  ## Calculate statistic used for ROC
  temp=list()
  for(ref in refs){
    temp[[ref]]=apply(rma_data[[ref]][[value]],2,statistic)
  }
  
  ## Calculate ROC objects for InLab vs ExLab
  nSamples=length(temp[[1]])
  response=c(rep("User",nSamples),rep("notUser",nSamples))
  temp2=list()
  for(ref in refs[-1]){
    temp2[[ref]]=roc(response,abs(c(temp[[1]], temp[[ref]])))
  }
  return(temp2)  
}

plot_roc=function(roc_list,statistic){
  plot.roc(roc_list[[1]], col=3, main=paste("InLab vs ExLab with", statistic))
  plot.roc(roc_list[[2]], col=4, add=T)
  plot.roc(roc_list[[3]], col=5, add=T)
  plot.roc(roc_list[[4]], col=6, add=T)
  auc=round(c(roc_list[[1]]$auc, roc_list[[2]]$auc,roc_list[[3]]$auc,roc_list[[4]]$auc),2)
  legend("bottomright", paste(names(roc_list)," - AUC",auc), fill=3:6, cex=0.8, bty="n")

}

optimize_roc=function(roc_object){
  ## Calculate optimal threshold
  ## Using Youden index
  youden=roc_object$sensitivities+roc_object$specificities
  y=which.max(youden)
  
  optimal=list("Threshold"=roc_object$thresholds[y],
               "Sensitiviy"=roc_object$sensitivities[y],
               "Specificity"=roc_object$specificities[y])
  return(optimal)  
}

REGS=function(rma_data){
  ResistanceClassifier(rma_data$exprs.sc)
}

BAGS2=function(rma_data){
  BAGS(rma_data$exprs.sc)
}

ABCGCB2=function(rma_data){
  ABCGCB(rma_data$exprs.sc)
}

###################################
#### Plots for CHEPRETRO ##########
###################################
chep_rma=readRDS("h:/hema_temp/RLE/data/chep_rma.rds")
names(chep_rma)[2]="InLab"

### Calculate ROC curves
chep_roc_rle_median=calc_roc(chep_rma,"RLE",median)
chep_roc_rle_iqr=calc_roc(chep_rma,"RLE",iqr)

if(tiff_out){
  tiff(file.path(outDir,"figureS2.tiff"), height=8, width=8, units="in", res=300, compression="lzw")
} else{
pdf(file.path(outDir,"chep_rle.pdf"), height=8, width=8)
}
  
par(mfcol=c(2,2))
### Scatter plots
plot_RLE(chep_rma, "median")
mtext("A", font = 2, adj = -0.1, line = 2.5, cex = 1.2)
plot_RLE(chep_rma, "iqr")
mtext("B", font = 2, adj = -0.1, line = 2.5, cex = 1.2)

### Plot ROC curves
plot_roc(chep_roc_rle_median, "RLE Median")
mtext("C", font = 2, adj = -0.1, line = 2.5, cex = 1.2)
plot_roc(chep_roc_rle_iqr, "RLE IQR")
mtext("D", font = 2, adj = -0.1, line = 2.5, cex = 1.2)
dev.off()

###############################
#### Plots for RCHOP ##########
###############################

RCHOP_rma=readRDS("h:/hema_temp/RLE/data/RCHOP_rma.rds")
names(RCHOP_rma)[2]="InLab"
names(RCHOP_rma)[4]="CHEPRETRO"

### Calculate ROC curves
RCHOP_roc_rle_median=calc_roc(RCHOP_rma,"RLE",median)
RCHOP_roc_rle_iqr=calc_roc(RCHOP_rma,"RLE",iqr)

if(tiff_out){
  tiff(file.path(outDir,"figureS3.tiff"), height=8, width=8, units="in", res=300, compression="lzw")
} else{
pdf(file.path(outDir,"RCHOP_rle.pdf"), height=8, width=8)
}
par(mfcol=c(2,2))
### Scatter plots
plot_RLE(RCHOP_rma, "median")
mtext("A", font = 2, adj = -0.1, line = 2.5, cex = 1.2)
plot_RLE(RCHOP_rma, "iqr")
mtext("B", font = 2, adj = -0.1, line = 2.5, cex = 1.2)

### Plot ROC curves
plot_roc(RCHOP_roc_rle_median, "RLE Median")
mtext("C", font = 2, adj = -0.1, line = 2.5, cex = 1.2)
plot_roc(RCHOP_roc_rle_iqr, "RLE IQR")
mtext("D", font = 2, adj = -0.1, line = 2.5, cex = 1.2)
dev.off()

##############################
#### Plots for CHOP ##########
##############################

CHOP_rma=readRDS("h:/hema_temp/RLE/data/CHOP_rma.rds")
names(CHOP_rma)[2]="InLab"
names(CHOP_rma)[4]="CHEPRETRO"

### Calculate ROC curves
CHOP_roc_rle_median=calc_roc(CHOP_rma,"RLE",median)
CHOP_roc_rle_iqr=calc_roc(CHOP_rma,"RLE",iqr)

if(tiff_out){
  tiff(file.path(outDir,"figureS4.tiff"), height=8, width=8, units="in", res=300, compression="lzw")
} else{
pdf(file.path(outDir,"CHOP_rle.pdf"), height=8, width=8)
}

par(mfcol=c(2,2))
### Scatter plots
plot_RLE(CHOP_rma, "median")
mtext("A", font = 2, adj = -0.1, line = 2.5, cex = 1.2)
plot_RLE(CHOP_rma, "iqr")
mtext("B", font = 2, adj = -0.1, line = 2.5, cex = 1.2)

### Plot ROC curves
plot_roc(CHOP_roc_rle_median, "RLE Median")
mtext("C", font = 2, adj = -0.1, line = 2.5, cex = 1.2)
plot_roc(CHOP_roc_rle_iqr, "RLE IQR")
mtext("D", font = 2, adj = -0.1, line = 2.5, cex = 1.2)
dev.off()

##############################
#### Plots for IDRC ##########
##############################

IDRC_rma=readRDS("h:/hema_temp/RLE/data/IDRC_rma.rds")
names(IDRC_rma)[2]="InLab"
names(IDRC_rma)[5]="CHEPRETRO"

### Calculate ROC curves
IDRC_roc_rle_median=calc_roc(IDRC_rma,"RLE",median)
IDRC_roc_rle_iqr=calc_roc(IDRC_rma,"RLE",iqr)

if(tiff_out){
  tiff(file.path(outDir,"figureS5.tiff"), height=8, width=8, units="in", res=300, compression="lzw")
} else{
pdf(file.path(outDir,"IDRC_rle.pdf"), height=8, width=8)
}
par(mfcol=c(2,2))
### Scatter plots
plot_RLE(IDRC_rma, "median")
mtext("A", font = 2, adj = -0.1, line = 2.5, cex = 1.2)
plot_RLE(IDRC_rma, "iqr")
mtext("B", font = 2, adj = -0.1, line = 2.5, cex = 1.2)

### Plot ROC curves
plot_roc(IDRC_roc_rle_median, "RLE Median")
mtext("C", font = 2, adj = -0.1, line = 2.5, cex = 1.2)
plot_roc(IDRC_roc_rle_iqr, "RLE IQR")
mtext("D", font = 2, adj = -0.1, line = 2.5, cex = 1.2)
dev.off()

rm(IDRC_rma)

##############################
#### Plots for MDFCI##########
##############################

MDFCI_rma=readRDS("h:/hema_temp/RLE/data/MDFCI_rma.rds")
names(MDFCI_rma)[2]="InLab"
names(MDFCI_rma)[5]="CHEPRETRO"

### Calculate ROC curves
MDFCI_roc_rle_median=calc_roc(MDFCI_rma,"RLE",median)
MDFCI_roc_rle_iqr=calc_roc(MDFCI_rma,"RLE",iqr)

if(tiff_out){
  tiff(file.path(outDir,"figureS6.tiff"), height=8, width=8, units="in", res=300, compression="lzw")
} else{
pdf(file.path(outDir,"MDFCI_rle.pdf"), height=8, width=8)
}
  
par(mfcol=c(2,2))
### Scatter plots
plot_RLE(MDFCI_rma, "median")
mtext("A", font = 2, adj = -0.1, line = 2.5, cex = 1.2)
plot_RLE(MDFCI_rma, "iqr")
mtext("B", font = 2, adj = -0.1, line = 2.5, cex = 1.2)

### Plot ROC curves
plot_roc(MDFCI_roc_rle_median, "RLE Median")
mtext("C", font = 2, adj = -0.1, line = 2.5, cex = 1.2)
plot_roc(MDFCI_roc_rle_iqr, "RLE IQR")
mtext("D", font = 2, adj = -0.1, line = 2.5, cex = 1.2)
dev.off()

################################################
### Calculate optimal thresholds for RLE IQR ###
################################################

chep_roc_rle_iqr_optimal=lapply(chep_roc_rle_iqr, optimize_roc)
RCHOP_roc_rle_iqr_optimal=lapply(RCHOP_roc_rle_iqr, optimize_roc)
CHOP_roc_rle_iqr_optimal=lapply(CHOP_roc_rle_iqr, optimize_roc)
IDRC_roc_rle_iqr_optimal=lapply(IDRC_roc_rle_iqr, optimize_roc)
MDFCI_roc_rle_iqr_optimal=lapply(MDFCI_roc_rle_iqr, optimize_roc)

### Format as table
dataset=c(rep("CHEPRETRO", 4), rep("RCHOP",4), rep("CHOP", 4), rep("IDRC", 4), rep("MDFCI", 4), "Median")
ref=c(names(chep_roc_rle_iqr), names(RCHOP_roc_rle_iqr), names(CHOP_roc_rle_iqr), names(IDRC_roc_rle_iqr), names(MDFCI_roc_rle_iqr), "-")

x1=matrix(unlist(chep_roc_rle_iqr_optimal), ncol=3, byrow=T)
x2=matrix(unlist(RCHOP_roc_rle_iqr_optimal), ncol=3, byrow=T)
x3=matrix(unlist(CHOP_roc_rle_iqr_optimal), ncol=3, byrow=T)
x4=matrix(unlist(IDRC_roc_rle_iqr_optimal), ncol=3, byrow=T)
x5=matrix(unlist(MDFCI_roc_rle_iqr_optimal), ncol=3, byrow=T)
dat=rbind(x1,x2,x3,x4,x5)
dat=round(rbind(dat,colMedians(dat)),2)

optimal_threshold_rle=data.frame(dataset,ref,dat)
names(optimal_threshold_rle)=c("Dataset", "RMA reference", "Threshold", "Sensitivity", "Specificity")
optimal_threshold_rle[,2]=gsub("CHEP", "CHEPRETRO",optimal_threshold_rle[,2])
optimal_threshold_rle

### Save table with results
sink(file.path(tableOut,"RLE_table.tex"))
print(xtable(optimal_threshold_rle,
             caption = "Optimal thresholds for RLE IQR",
             label="rleTable"),
      include.rownames=FALSE)
sink()

############################################
### RLE IQR vs Classification: CHEPRETRO ###
############################################

chep_REGS=lapply(chep_rma, REGS)
chep_BAGS=lapply(chep_rma, BAGS2)
chep_ABCGCB=lapply(chep_rma, ABCGCB2)

### RLE IQR vs BAGS accuracy
result_bags=list()
threshold=seq(0.3,1, by=0.01)

references=names(chep_rma)[-1]
nSamples=ncol(chep_rma[["cohort"]]$exprs)

for(ref in references){
  i=1
  for(t in threshold){
    confusionMatrix=table(chep_BAGS[["cohort"]]$class[abs(chep_rma[[ref]]$RLE.stats[,3])<t]
                          , chep_BAGS[[ref]]$class[abs(chep_rma[[ref]]$RLE.stats[,3])<t])
    
    result_bags[[ref]]$threshold_total_prop[i]=sum(confusionMatrix)/nSamples
    result_bags[[ref]]$threshold_diag_prop[i]=sum(diag(confusionMatrix))/sum(confusionMatrix)
    i=i+1
  }
}
par(mfcol=c(2,5))
for(ref in references){
  plot(threshold, result_bags[[ref]]$threshold_total_prop, xlab="RLE IQR", ylab="% Samples", main=ref, ylim=c(0,1))
  abline(v=0.6, col="red")
  plot(threshold, result_bags[[ref]]$threshold_diag_prop, xlab="RLE IQR", ylab="Accuracy", main=ref, ylim=c(0,1))
  abline(v=0.6, col="red")
}

if(tiff_out){
  tiff(file.path(outDir,"figureS7.tiff"), height=8.75, width=7.5, units="in", res=300, compression="lzw")
} else{
pdf(file.path(outDir,"chep_rle_classification_bags.pdf"), height=8.75, width=7.5)
}
par(mfrow=c(3,2))
panel=1
for(ref in c("InLab","RCHOP","CHOP")){
  plot(threshold, result_bags[[ref]]$threshold_total_prop, xlab="RLE IQR", ylab="% Samples", main=ref, ylim=c(0,1))
  abline(v=0.6, col="red")
  mtext(LETTERS[panel], font = 2, adj = -0.1, line = 2.5, cex = 1.2)
  plot(threshold, result_bags[[ref]]$threshold_diag_prop, xlab="RLE IQR", ylab="Accuracy", main=ref, ylim=c(0,1))
  abline(v=0.6, col="red")
  mtext(LETTERS[panel+1], font = 2, adj = -0.1, line = 2.5, cex = 1.2)
  panel=panel+2
}
dev.off()


### RLE IQR vs ABC/GCB accuracy
result_abcgcb=list()
threshold=seq(0.3,1, by=0.01)

references=names(chep_rma)[-1]
nSamples=ncol(chep_rma[["cohort"]]$exprs)

for(ref in references){
  i=1
  for(t in threshold){
    confusionMatrix=table(chep_ABCGCB[["cohort"]]$class[abs(chep_rma[[ref]]$RLE.stats[,3])<t]
                          , chep_ABCGCB[[ref]]$class[abs(chep_rma[[ref]]$RLE.stats[,3])<t])
    
    result_abcgcb[[ref]]$threshold_total_prop[i]=sum(confusionMatrix)/nSamples
    result_abcgcb[[ref]]$threshold_diag_prop[i]=sum(diag(confusionMatrix))/sum(confusionMatrix)
    i=i+1
  }
}
par(mfcol=c(2,5))
for(ref in references){
  plot(threshold, result_abcgcb[[ref]]$threshold_total_prop, xlab="RLE IQR", ylab="% Samples", main=ref, ylim=c(0,1))
  abline(v=0.6, col="red")
  plot(threshold, result_abcgcb[[ref]]$threshold_diag_prop, xlab="RLE IQR", ylab="Accuracy", main=ref, ylim=c(0,1))
  abline(v=0.6, col="red")
}

if(tiff_out){
  tiff(file.path(outDir,"figureS8.tiff"), height=8.75, width=7.5, units="in", res=300, compression="lzw")
} else{
pdf(file.path(outDir,"chep_rle_classification_abcgcb.pdf"), height=8.75, width=7.5)
}
par(mfrow=c(3,2))
panel=1
for(ref in c("InLab","RCHOP","CHOP")){
  plot(threshold, result_abcgcb[[ref]]$threshold_total_prop, xlab="RLE IQR", ylab="% Samples", main=ref, ylim=c(0,1))
  abline(v=0.6, col="red")
  mtext(LETTERS[panel], font = 2, adj = -0.1, line = 2.5, cex = 1.2)
  plot(threshold, result_abcgcb[[ref]]$threshold_diag_prop, xlab="RLE IQR", ylab="Accuracy", main=ref, ylim=c(0,1))
  abline(v=0.6, col="red")
  mtext(LETTERS[panel+1], font = 2, adj = -0.1, line = 2.5, cex = 1.2)
  panel=panel+2
}
dev.off()


### RLE IQR vs REGS accuracy
result_regs=list()
threshold=seq(0.3,1, by=0.01)

references=names(chep_rma)[-1]
nSamples=ncol(chep_rma[["cohort"]]$exprs)

for(ref in references){
  i=1
  for(t in threshold){
    confusionMatrix=table(chep_REGS[["cohort"]]$class[,4][abs(chep_rma[[ref]]$RLE.stats[,3])<t]
                          , chep_REGS[[ref]]$class[,4][abs(chep_rma[[ref]]$RLE.stats[,3])<t])
    
    result_regs[[ref]]$threshold_total_prop[i]=sum(confusionMatrix)/nSamples
    result_regs[[ref]]$threshold_diag_prop[i]=sum(diag(confusionMatrix))/sum(confusionMatrix)
    i=i+1
  }
}
par(mfcol=c(2,5))
for(ref in references){
  plot(threshold, result_regs[[ref]]$threshold_total_prop, xlab="RLE IQR", ylab="% Samples", main=ref, ylim=c(0,1))
  abline(v=0.6, col="red")
  plot(threshold, result_regs[[ref]]$threshold_diag_prop, xlab="RLE IQR", ylab="Accuracy", main=ref, ylim=c(0,1))
  abline(v=0.6, col="red")
}

if(tiff_out){
  tiff(file.path(outDir,"figureS9.tiff"), height=8.75, width=7.5, units="in", res=300, compression="lzw")
} else{
pdf(file.path(outDir,"chep_rle_classification_regs.pdf"), height=8.75, width=7.5)
}
par(mfrow=c(3,2))
panel=1
for(ref in c("InLab","RCHOP","CHOP")){
  plot(threshold, result_regs[[ref]]$threshold_total_prop, xlab="RLE IQR", ylab="% Samples", main=ref, ylim=c(0,1))
  abline(v=0.6, col="red")
  mtext(LETTERS[panel], font = 2, adj = -0.1, line = 2.5, cex = 1.2)
  plot(threshold, result_regs[[ref]]$threshold_diag_prop, xlab="RLE IQR", ylab="Accuracy", main=ref, ylim=c(0,1))
  abline(v=0.6, col="red")
  mtext(LETTERS[panel+1], font = 2, adj = -0.1, line = 2.5, cex = 1.2)
  panel=panel+2
}
dev.off()

########################################
### RLE IQR vs Classification: RCHOP ###
########################################

RCHOP_REGS=lapply(RCHOP_rma, REGS)
RCHOP_BAGS=lapply(RCHOP_rma, BAGS2)
RCHOP_ABCGCB=lapply(RCHOP_rma, ABCGCB2)

### RLE IQR vs BAGS accuracy
result_bags=list()
threshold=seq(0.3,1, by=0.01)

references=names(RCHOP_rma)[-1]
nSamples=ncol(RCHOP_rma[["cohort"]]$exprs)

for(ref in references){
  i=1
  for(t in threshold){
    confusionMatrix=table(RCHOP_BAGS[["cohort"]]$class[abs(RCHOP_rma[[ref]]$RLE.stats[,3])<t]
                          , RCHOP_BAGS[[ref]]$class[abs(RCHOP_rma[[ref]]$RLE.stats[,3])<t])
    
    result_bags[[ref]]$threshold_total_prop[i]=sum(confusionMatrix)/nSamples
    result_bags[[ref]]$threshold_diag_prop[i]=sum(diag(confusionMatrix))/sum(confusionMatrix)
    i=i+1
  }
}
par(mfcol=c(2,5))
for(ref in references){
  plot(threshold, result_bags[[ref]]$threshold_total_prop, xlab="RLE IQR", ylab="% Samples", main=ref, ylim=c(0,1))
  abline(v=0.6, col="red")
  plot(threshold, result_bags[[ref]]$threshold_diag_prop, xlab="RLE IQR", ylab="Accuracy", main=ref, ylim=c(0,1))
  abline(v=0.6, col="red")
}

if(tiff_out){
  tiff(file.path(outDir,"figureS10.tiff"), height=8.75, width=7.5, units="in", res=300, compression="lzw")
} else{
pdf(file.path(outDir,"RCHOP_rle_classification_bags.pdf"), height=8.75, width=7.5)
}
par(mfrow=c(3,2))
panel=1
for(ref in c("InLab","CHEPRETRO","CHOP")){
  plot(threshold, result_bags[[ref]]$threshold_total_prop, xlab="RLE IQR", ylab="% Samples", main=ref, ylim=c(0,1))
  abline(v=0.6, col="red")
  mtext(LETTERS[panel], font = 2, adj = -0.1, line = 2.5, cex = 1.2)
  plot(threshold, result_bags[[ref]]$threshold_diag_prop, xlab="RLE IQR", ylab="Accuracy", main=ref, ylim=c(0,1))
  abline(v=0.6, col="red")
  mtext(LETTERS[panel+1], font = 2, adj = -0.1, line = 2.5, cex = 1.2)
  panel=panel+2
}
dev.off()


### RLE IQR vs ABC/GCB accuracy
result_abcgcb=list()
threshold=seq(0.3,1, by=0.01)

references=names(RCHOP_rma)[-1]
nSamples=ncol(RCHOP_rma[["cohort"]]$exprs)

for(ref in references){
  i=1
  for(t in threshold){
    confusionMatrix=table(RCHOP_ABCGCB[["cohort"]]$class[abs(RCHOP_rma[[ref]]$RLE.stats[,3])<t]
                          , RCHOP_ABCGCB[[ref]]$class[abs(RCHOP_rma[[ref]]$RLE.stats[,3])<t])
    
    result_abcgcb[[ref]]$threshold_total_prop[i]=sum(confusionMatrix)/nSamples
    result_abcgcb[[ref]]$threshold_diag_prop[i]=sum(diag(confusionMatrix))/sum(confusionMatrix)
    i=i+1
  }
}
par(mfcol=c(2,5))
for(ref in references){
  plot(threshold, result_abcgcb[[ref]]$threshold_total_prop, xlab="RLE IQR", ylab="% Samples", main=ref, ylim=c(0,1))
  abline(v=0.6, col="red")
  plot(threshold, result_abcgcb[[ref]]$threshold_diag_prop, xlab="RLE IQR", ylab="Accuracy", main=ref, ylim=c(0,1))
  abline(v=0.6, col="red")
}

if(tiff_out){
  tiff(file.path(outDir,"figureS11.tiff"), height=8.75, width=7.5, units="in", res=300, compression="lzw")
} else{
pdf(file.path(outDir,"RCHOP_rle_classification_abcgcb.pdf"), height=8.75, width=7.5)
}
par(mfrow=c(3,2))
panel=1
for(ref in c("InLab","CHEPRETRO","CHOP")){
  plot(threshold, result_abcgcb[[ref]]$threshold_total_prop, xlab="RLE IQR", ylab="% Samples", main=ref, ylim=c(0,1))
  abline(v=0.6, col="red")
  mtext(LETTERS[panel], font = 2, adj = -0.1, line = 2.5, cex = 1.2)
  plot(threshold, result_abcgcb[[ref]]$threshold_diag_prop, xlab="RLE IQR", ylab="Accuracy", main=ref, ylim=c(0,1))
  abline(v=0.6, col="red")
  mtext(LETTERS[panel+1], font = 2, adj = -0.1, line = 2.5, cex = 1.2)
  panel=panel+2
}
dev.off()


### RLE IQR vs REGS accuracy
result_regs=list()
threshold=seq(0.3,1, by=0.01)

references=names(RCHOP_rma)[-1]
nSamples=ncol(RCHOP_rma[["cohort"]]$exprs)

for(ref in references){
  i=1
  for(t in threshold){
    confusionMatrix=table(RCHOP_REGS[["cohort"]]$class[,4][abs(RCHOP_rma[[ref]]$RLE.stats[,3])<t]
                          , RCHOP_REGS[[ref]]$class[,4][abs(RCHOP_rma[[ref]]$RLE.stats[,3])<t])
    
    result_regs[[ref]]$threshold_total_prop[i]=sum(confusionMatrix)/nSamples
    result_regs[[ref]]$threshold_diag_prop[i]=sum(diag(confusionMatrix))/sum(confusionMatrix)
    i=i+1
  }
}
par(mfcol=c(2,5))
for(ref in references){
  plot(threshold, result_regs[[ref]]$threshold_total_prop, xlab="RLE IQR", ylab="% Samples", main=ref, ylim=c(0,1))
  abline(v=0.6, col="red")
  plot(threshold, result_regs[[ref]]$threshold_diag_prop, xlab="RLE IQR", ylab="Accuracy", main=ref, ylim=c(0,1))
  abline(v=0.6, col="red")
}

if(tiff_out){
  tiff(file.path(outDir,"figureS12.tiff"), height=8.75, width=7.5, units="in", res=300, compression="lzw")
} else{
pdf(file.path(outDir,"RCHOP_rle_classification_regs.pdf"), height=8.75, width=7.5)
}  
par(mfrow=c(3,2))
panel=1
for(ref in c("InLab","CHEPRETRO","CHOP")){
  plot(threshold, result_regs[[ref]]$threshold_total_prop, xlab="RLE IQR", ylab="% Samples", main=ref, ylim=c(0,1))
  abline(v=0.6, col="red")
  mtext(LETTERS[panel], font = 2, adj = -0.1, line = 2.5, cex = 1.2)
  plot(threshold, result_regs[[ref]]$threshold_diag_prop, xlab="RLE IQR", ylab="Accuracy", main=ref, ylim=c(0,1))
  abline(v=0.6, col="red")
  mtext(LETTERS[panel+1], font = 2, adj = -0.1, line = 2.5, cex = 1.2)
  panel=panel+2
}
dev.off()


##################################
### Annotate classifier probes ###
##################################

### Read Affy Annotation
affyAnno=read.csv("//ithfil08/HAELAB125/ForskningsLAB/Statlab/Database2/Annotation/hgu133plus2/Affymetrix/HG-U133_Plus_2.na33.annot.csv", comment.char = "#", na.string="---")

### Function for annotating probe IDs
list_anno=function(x){
  n=length(x)
  nHugo=length(unique(affyAnno$Gene.Symbol[affyAnno$Probe.Set.ID%in%x]))
  nEnSembl=length(unique(affyAnno$Ensembl[affyAnno$Probe.Set.ID%in%x]))
  return(list("nProbes"=n, "nHGNC"=nHugo, "nEnsembl"=nEnSembl))
}

classGenes=list()

### ABCGCB probes
classGenes[["ABC/GCB"]]=rownames(readABCGCBCoef())[-1]

### BAGS probes
classGenes[["BAGS"]]=rownames(readBAGSCoef())[-1]

## Regs Probes
classGenes[["Vincristine Classifier"]]=names(readVincristineClasCoef())
classGenes[["Vincristine Predictor"]]=names(readVincristinePredCoef())

classGenes[["Cyclophosphamide Classifier"]]=names(readCyclophosphamideClasCoef())
classGenes[["Cyclophosphamide Predictor"]]=names(readCyclophosphamidePredCoef())


classGenes[["Doxorubicine Classifier"]]=names(readDoxorubicinClasCoef())
classGenes[["Doxorubicine Predictor"]]=names(readDoxorubicinPredCoef())

classGenes[["Combined Classifier"]]=rownames(readClasCoef())
classGenes[["Combined Predictor"]]=rownames(readPredCoef())

### Annotate
probeTable=t(sapply(classGenes,list_anno))
probeTable2=data.frame("Classifier"=row.names(probeTable),probeTable)

sink(file.path(tableOut,"probeTable.tex"))
print(xtable(probeTable2,
             caption = "Number of probes used in the classifiers and the number of corresponding HGNC and Ensembl gene IDs",
             label="probeTable",
             digits=0),
      include.rownames=FALSE)
sink()