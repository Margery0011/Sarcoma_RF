#-----------------------------------------------------------------------------------
# script to include the batchadjustment into the cross validation
#                                                                     
# Martin Sill
# m.sill@dkfz.de                                                                  
# 
# 2018-03-14 UTC
#------------------------------------------------------------------------------------                   
#options(max.print = 1000)
#options(stringsAsFactors = FALSE)
#options(scipen = 999)


Combinebatchadjust <- function(Mset,batchall,fold){
  methy <- getMeth(Mset[,fold$train])
  unmethy <- getUnmeth(Mset[,fold$train])
  batch <- batchall[fold$train]
  methy.ba <- 2^removeBatchEffect(log2(methy +1), batch)
  unmethy.ba <- 2^removeBatchEffect(log2(unmethy +1), batch)
  # recalculate betas, illumina like
  betas <- methy.ba / (methy.ba +unmethy.ba +100)
  betas <- as.data.frame(t(betas))
  # extract coefficients that can be applied to diagnostic samples
  s.850k <- min(which(batch == '850k'))
  s.450k <- min(which(batch == '450k'))
  methy.coef <- unmethy.coef <- list()
  methy.coef[["850k"]] <- log2(methy.ba[, s.850k]) - log2(methy[, s.850k] +1)
  methy.coef[["450k"]] <- log2(methy.ba[, s.450k]) - log2(methy[, s.450k] +1)
  unmethy.coef[["850k"]] <- log2(unmethy.ba[, s.850k]) - log2(unmethy[, s.850k] +1)
  unmethy.coef[["450k"]] <- log2(unmethy.ba[, s.450k]) - log2(unmethy[, s.450k] +1)
  methy <- getMeth(Mset[,fold$test])
  unmethy <- getUnmeth(Mset[,fold$test])
  batch <- batchall[fold$test]
  # perform batch adjustment
  methy.b <- log2(methy +1) + matrix(unlist(methy.coef[match(batch,names(methy.coef))]),ncol=length(batch))
  unmethy.b <- log2(unmethy +1) + matrix(unlist(unmethy.coef[match(batch,names(unmethy.coef))]),ncol=length(batch))
  methy.b[methy.b < 0] <- 0
  unmethy.b[unmethy.b < 0] <- 0
  methy.ba <- 2^methy.b
  unmethy.ba <- 2^unmethy.b
  # illumina-like beta values
  betas.test <- methy.ba / (methy.ba +unmethy.ba +100)
  betas.test <- as.data.frame(t(betas.test))
  return(list(betas.train=betas,betas.test=betas.test))
}