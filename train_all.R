library(randomForest)
library(parallel)

ntrees <- 500  # 10000 in the paper, here 500 to speed up the example
cores <- 4
seed <- 180314
p <- 10000 

source(file.path("R","train.R"))

y <- as.factor(anno$Methylation.Class.Name)

# sd pre filtering to 20k probes, to speed up the example
load("/Volumes/lyt/Data_Practice_Sarcoma/results_all/betas.ba.RData")
betas <- betas[,order(-apply(betas,2,sd))[1:20000]]

set.seed(seed,kind ="L'Ecuyer-CMRG") 
message("seed: ",seed)
message("cores: ",cores)
message("ntrees: ",ntrees)  
message("n: ",nrow(betas))

message("p: ",ncol(betas))  


rf.varsel <- rfp(betas,
                 y,
                 mc=cores,
                 ntree=ntrees,
                 sampsize=rep(min(table(y)),length(table(y))),
                 importance=TRUE)
# get permutation variable importance
imp.meandecrease <- rf.varsel$importance[,dim(rf.varsel$importance)[2]-1]
# save selection forest
save(rf.varsel,file=file.path("results_all/","varsel.RData"))
save(rf.varsel,file=file.path("results_all/","varsel.try.RData"))
save(rf.varsel,file=file.path("results_all/","varsel.try2.RData"))

rm(rf.varsel)

# reduce data matrix
or <- order(imp.meandecrease,decreasing=T)
p <- 10000
betasy <- betas[,or[1:p]]
gc()

message("single core")
message("ntrees: ",ntrees)  
message("n: ",nrow(betasy))
message("p: ",ncol(betasy))

rf.pred <- randomForest(betasy,
                        y,
                        #mc=cores,
                        ntree=ntrees,
                        #strata=y,
                        #mtry=sqrt(ncol(betas)),
                        sampsize=rep(min(table(y)),length(table(y))),
                        proximity=TRUE,
                        oob.prox=TRUE,
                        importance=TRUE,
                        keep.inbag=TRUE,
                        do.trace=FALSE,
                        seed=seed
)
save(rf.pred,file=file.path("results_all/","rf.pred.RData"))
pred <- as.data.frame(rf.pred$confusion)
