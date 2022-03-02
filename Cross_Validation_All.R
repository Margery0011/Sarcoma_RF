#-----------------------------------------------------------------------------------
# nested cross-validation 
#                                                                   
rm(list=ls())

library(randomForest)
library(parallel)

library(minfi)
library(limma)

ntrees <- 500
cores <- 4
seed <- 180314
p <- 10000
folds <- 3

message("loading filtered Mset ...",Sys.time())
load(file.path("results_all/","Mset_filtered.RData"))

anno_45 <- read.csv("Merge_45.csv")
anno_45['platform'] <- '450k'

anno_epic <- read.csv("Merge_Epic.csv")
anno_epic['platform']<-'850k'
anno <- rbind(anno_45,anno_epic)


y <-as.factor(anno$Methylation.Class.Name)
#y <- as.factor(anno$`methylation class:ch1`)
batch <- as.factor(anno$platform)

source(file.path("R","makefolds.R"))
source(file.path("R","train.R"))
source(file.path("R","CombineCalculateCVfold.R"))
source(file.path("R","CombineBatchadjust.R"))

if(!file.exists(file.path("CV_alls","nfolds.RData"))){
  dir.create("CV_alls",showWarnings = FALSE)
  nfolds <- makenestedfolds(y,folds)
  save(nfolds,file=file.path("CV_alls","nfolds.RData"))
}
load(file.path("CV_alls","nfolds.RData"))

message("performing nested CV ...", Sys.time())
message("check minimal class sizes for inner training loops")

# check minimal class sizes for inner training loops
minclasssize <- matrix(0,ncol=length(nfolds),nrow=length(nfolds))
for(i in 1:length(nfolds)){
  for(j in 1:length(nfolds))
    minclasssize[i,j]  <- min(table(y[nfolds[[i]][[2]][[j]]$train]))
}
colnames(minclasssize) <- paste0("innfold",1:folds)
rownames(minclasssize) <- paste0("fold",1:folds)
print(minclasssize)

for(K in 1:folds){
  
  for(k in 0:folds){
    
    if(k>0){  message("calculateing fold ",K,".",k,"  ...",Sys.time())
      fold <- nfolds[[K]][[2]][[k]]
    }else{
      message("calculateing outer fold ",K,"  ...",Sys.time())
      fold <- nfolds[[K]][[1]][[1]]
    }
    
    rf.scores <- CombinecalcultateCVfold(Mset,y,batch,fold,p,cores,ntrees)
    
    fname <- paste("CVfold",K,k,"RData",sep=".")
    save(rf.scores,file=file.path("CV_alls",fname))
    
    rm(rf.scores)
    gc()
  }
}
message("finished ...",Sys.time())
