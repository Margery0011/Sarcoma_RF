gc()
rm(list=ls())
library(minfi)
library(GEOquery)
library(limma)

source(file.path("R","MNPprocessIDAT_functions.R"))
dir.create("results_all")

anno_45 <- read.csv("Merge_45.csv")

anno_45['platform'] <- '450k'


# read raw data downloaded from GEO and extracted in GSE140686_RAW
filepath <- file.path("GSE140686_RAW/",gsub("_Grn.*","",gsub(".*suppl/","",anno_45$supplementary_file)))
RGset_45 <- read.metharray(filepath,verbose=TRUE,force=TRUE)
anno_epic <- read.csv("Merge_Epic.csv")
anno_epic['platform']<-'850k'
anno <- rbind(anno_45,anno_epic)

filepath <- file.path("GSE140686_RAW/",gsub("_Grn.*","",gsub(".*suppl/","",anno_epic$supplementary_file)))
RGset_45 <- read.metharray(filepath,verbose=TRUE,force=TRUE)
RGset_epic <- read.metharray(filepath,verbose=TRUE,force=TRUE)

RGset <- combineArrays(RGset_45, RGset_epic,
                       outType = c("IlluminaHumanMethylation450k",
                                   "IlluminaHumanMethylationEPIC"),
                       verbose = TRUE)
gc()

save(RGset,file=file.path("results_all/","RGset.RData"))
save(RGset_45,file=file.path("results_all/","RGset_45.RData"))
rm(RGset_45)
rm(RGset_epic)
save(RGset_epic,file=file.path("results_all/","RGset_epic.RData"))
Mset <- MNPpreprocessIllumina(RGset)

message("probe filtering ...",Sys.time())
amb.filter <- read.table(file.path("filter","amb_3965probes.vh20151030.txt"),header=F)
#epic.filter <- read.table(file.path("filter","epicV1B2_32260probes.vh20160325.txt"),header=F)
snp.filter <- read.table(file.path("filter","snp_7998probes.vh20151030.txt"),header=F)
xy.filter <- read.table(file.path("filter","xy_11551probes.vh20151030.txt"),header=F)
rs.filter <- grep("rs",rownames(Mset))
ch.filter <- grep("ch",rownames(Mset))

remove <- unique(c(match(amb.filter[,1], rownames(Mset)),
                   match(snp.filter[,1], rownames(Mset)),
                   match(xy.filter[,1], rownames(Mset)),
                   rs.filter,
                   ch.filter))


Mset_filtered <- Mset[-remove[!is.na(remove)],]
save(Mset,Mset_filtered,file=file.path("results_all/","Mset_filtered.RData"))
rm(RGset)
rm(amb.filter,snp.filter,xy.filter)

methy <- getMeth(Mset_filtered)
unmethy <- getUnmeth(Mset_filtered)

four5k <- anno$platform
batch <- ifelse(four5k == '450k',2,1)

# get FFPE/Frozen type
#ffpe <- anno$material.preparation.ch1
#batch2 <- ifelse(ffpe == "FFPE", 3, 4)

# remove batch effects by linear model
library(limma)

methy.ba <- 2^removeBatchEffect(log2(methy +1), batch)

gc()
rm(Mset)
gc()
unmethy.ba <- 2^removeBatchEffect(log2(unmethy +1), batch)

# extract effects to adjust diagnostic samples
s.850k <- min(which(batch == 1))
s.450k <- min(which(batch == 2))
methy.coef <- unmethy.coef <- list()
methy.coef[["850k"]] <- log2(methy.ba[, s.850k]) - log2(methy[, s.850k] +1)
methy.coef[["450k"]] <- log2(methy.ba[, s.450k]) - log2(methy[, s.450k] +1)
unmethy.coef[["850k"]] <- log2(unmethy.ba[, s.850k]) - log2(unmethy[, s.850k] +1)
unmethy.coef[["450k"]] <- log2(unmethy.ba[, s.450k]) - log2(unmethy[, s.450k] +1)

save(methy.coef,unmethy.coef,file=file.path("results_all/","ba.coef.RData"))
# recalculate betas, illumina like
betas <- methy.ba / (methy.ba +unmethy.ba +100)
betas <- as.data.frame(t(betas))
save(betas,anno,file=file.path("results_all/","betas.ba.RData"))  
save(methy,methy.ba,file=file.path("results_all","methy.RData"))
save(unmethy,unmethy.ba,file=file.path("results_all/","unmethy.RData"))
gc()
rm(methy,methy.ba,unmethy,unmethy.ba)
gc()


