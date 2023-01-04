# clear ws
rm(list = ls(all = TRUE))

# library
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(data.table)
library(xtable)

# source file
source("varname.R")
source("fun.R")

# Load data
load("info_reco.RData")
load("./data/fore_csreco.RData")
load("./data/fore_ctreco.RData")
load("./data/fore_treco.RData")
load("./data/fore_seq.RData")
load("./data/fore_ka.RData")
load("./data/fore_ite.RData")
load("./data/fore_ite_tol.RData")
load("./data/fore_seqL2.RData")
load("./data/fore_base_test.RData")

## nRMSE ----
RMSE <- function(actual, predict){ sqrt(mean((actual-predict)^2))/mean(actual)*100 }

kid <- rep(thf_info$kset, h*thf_info$m/thf_info$kset)
opday <- unlist(sapply(thf_info$m/thf_info$kset, function(x) rep(1:h, each =x)))

nRMSEtibble <- tibble(k = as.numeric(),
                      id = as.character(),
                      comb =  as.character(),
                      method =  as.character(),
                      type = as.character(),
                      value = as.numeric(), 
                      h = as.numeric())

## cross-sectional ----
cont <- c("cols", "cshr", "cstruc", "cwls")
rownames(hts_info$S) <- varname
for(i in cont){
  # free
  tilde <- simplify2array(get(i))
  
  # sntz
  tilde_bts <- tilde[,-c(1:6),]
  tilde_bts[tilde_bts<0] <- 0
  tmp <- apply(tilde_bts, 3, function(x) t(as.matrix(hts_info$S %*% t(x))), simplify = FALSE)
  tilde0 <- simplify2array(tmp)
  
  # osqp
  tilde_nn <- simplify2array(get(paste0(i, "_nn")))
  
  for(k in thf_info$kset){
    for(j in 1:NCOL(tilde)){
      nRMSEtibble <- add_row(nRMSEtibble, 
                             k = k,
                             id = colnames(tilde)[j],
                             comb =  str_sub(i, 2, -1),
                             method =  "cs",
                             type = c("tilde", "tilde0", "osqp"),
                             value = c(RMSE(test[kid==k,j,], tilde[kid==k,j,]),
                                       RMSE(test[kid==k,j,], tilde0[kid==k,j,]),
                                       RMSE(test[kid==k,j,], tilde_nn[kid==k,j,])),
                             h = 0
      )
      
      nRMSEtibble <- add_row(nRMSEtibble, 
                             k = k,
                             id = colnames(tilde)[j],
                             comb =  str_sub(i, 2, -1),
                             method =  "cs",
                             type = c("tilde", "tilde0", "osqp"),
                             value = c(RMSE(test[kid==k & opday == 1,j,], 
                                            tilde[kid==k & opday == 1,j,]),
                                       RMSE(test[kid==k & opday == 1,j,], 
                                            tilde0[kid==k & opday == 1,j,]),
                                       RMSE(test[kid==k & opday == 1,j,], 
                                            tilde_nn[kid==k & opday == 1,j,])),
                             h = 1
      )
      
      nRMSEtibble <- add_row(nRMSEtibble, 
                             k = k,
                             id = colnames(tilde)[j],
                             comb =  str_sub(i, 2, -1),
                             method =  "cs",
                             type = c("tilde", "tilde0", "osqp"),
                             value = c(RMSE(test[kid==k & opday == 2,j,], 
                                            tilde[kid==k & opday == 2,j,]),
                                       RMSE(test[kid==k & opday == 2,j,], 
                                            tilde0[kid==k & opday == 2,j,]),
                                       RMSE(test[kid==k & opday == 2,j,], 
                                            tilde_nn[kid==k & opday == 2,j,])),
                             h = 2
      )
    }
    cat(k, " ")
  }
  cat(i, "\n")
}

## cross-temporal ----
tempocont <- c("ctols", "ctbdshr", "ctstruc", "ctwlsv")
for(i in tempocont){
  # free
  tilde <- simplify2array(get(i))
  
  # sntz
  tilde0_bts <- tilde[rep(thf_info$kset, h*thf_info$m/thf_info$kset) == 1,-c(1:6),]
  tilde0_bts[tilde0_bts<0] <- 0
  tmp <- apply(tilde0_bts, 3, function(x) t(as.matrix(hts_info$S %*% t(x))), simplify = FALSE)
  tmp <- apply(simplify2array(tmp), 3, function(x) as.matrix(thf_info$R %*% x), simplify = FALSE)
  tilde0 <- simplify2array(tmp)
  
  # osqp
  tilde_nn <- simplify2array(get(paste0(i, "_nn")))
  
  for(k in thf_info$kset){
    for(j in 1:NCOL(tilde)){
      nRMSEtibble <- add_row(nRMSEtibble, 
                             k = k,
                             id = colnames(tilde)[j],
                             comb =  str_sub(i, 3, -1),
                             method =  "ct",
                             type = c("tilde", "tilde0", "osqp"),
                             value = c(RMSE(test[kid==k,j,], tilde[kid==k,j,]),
                                       RMSE(test[kid==k,j,], tilde0[kid==k,j,]),
                                       RMSE(test[kid==k,j,], tilde_nn[kid==k,j,])),
                             h = 0
      )
      
      nRMSEtibble <- add_row(nRMSEtibble, 
                             k = k,
                             id = colnames(tilde)[j],
                             comb =  str_sub(i, 3, -1),
                             method =  "ct",
                             type = c("tilde", "tilde0", "osqp"),
                             value = c(RMSE(test[kid==k & opday == 1,j,], 
                                            tilde[kid==k & opday == 1,j,]),
                                       RMSE(test[kid==k & opday == 1,j,], 
                                            tilde0[kid==k & opday == 1,j,]),
                                       RMSE(test[kid==k & opday == 1,j,], 
                                            tilde_nn[kid==k & opday == 1,j,])),
                             h = 1
      )
      
      nRMSEtibble <- add_row(nRMSEtibble, 
                             k = k,
                             id = colnames(tilde)[j],
                             comb =  str_sub(i, 3, -1),
                             method =  "ct",
                             type = c("tilde", "tilde0", "osqp"),
                             value = c(RMSE(test[kid==k & opday == 2,j,], 
                                            tilde[kid==k & opday == 2,j,]),
                                       RMSE(test[kid==k & opday == 2,j,], 
                                            tilde0[kid==k & opday == 2,j,]),
                                       RMSE(test[kid==k & opday == 2,j,], 
                                            tilde_nn[kid==k & opday == 2,j,])),
                             h = 2
      )
    }
    cat(k, " ")
  }
  cat(i, "\n")
}

## base ----
for(k in thf_info$kset){
  for(j in 1:NCOL(base)){
    nRMSEtibble <- add_row(nRMSEtibble, 
                           k = k,
                           id = colnames(base)[j],
                           comb =  "base",
                           method =  "none",
                           type = "none",
                           value = c(RMSE(test[kid==k,j,], base[kid==k,j,])),
                           h = 0
    )
    
    nRMSEtibble <- add_row(nRMSEtibble, 
                           k = k,
                           id = colnames(base)[j],
                           comb =  "base",
                           method =  "none",
                           type = "none",
                           value = c(RMSE(test[kid==k & opday == 1,j,], 
                                          base[kid==k & opday == 1,j,])),
                           h = 1
    )
    
    nRMSEtibble <- add_row(nRMSEtibble, 
                           k = k,
                           id = colnames(base)[j],
                           comb =  "base",
                           method =  "none",
                           type = "none",
                           value = c(RMSE(test[kid==k & opday == 2,j,], 
                                          base[kid==k & opday == 2,j,])),
                           h = 2
    )
  }
  cat(k, " ")
}

## temporal ----
tempo <- c("tols", "tsar1", "tstruc", "twlsv")
for(i in tempo){
  # free
  tilde <- simplify2array(get(i))
  
  # z
  ddot <- tilde
  ddot[ddot<0] <- 0
  
  # yang
  yang <- ddot
  yang[test==0] <- 0
  
  # sntz
  tilde0_1 <- tilde[rep(thf_info$kset, h*thf_info$m/thf_info$kset) == 1,,]
  tilde0_1[tilde0_1<0] <- 0
  tmp <- apply(tilde0_1, 3, function(x) as.matrix(thf_info$R %*% x), simplify = FALSE)
  tilde0 <- simplify2array(tmp)
  
  # osqp
  tilde_nn <- simplify2array(get(paste0(i, "_nn")))
  
  for(k in thf_info$kset){
    for(j in 1:NCOL(tilde)){
      nRMSEtibble <- add_row(nRMSEtibble, 
                             k = k,
                             id = colnames(tilde)[j],
                             comb =  str_sub(i, 2, -1),
                             method =  "t",
                             type = c("tilde", "ddot", "yang","tilde0", "osqp"),
                             value = c(RMSE(test[kid==k,j,], tilde[kid==k,j,]),
                                       RMSE(test[kid==k,j,], ddot[kid==k,j,]),
                                       RMSE(test[kid==k,j,], yang[kid==k,j,]),
                                       RMSE(test[kid==k,j,], tilde0[kid==k,j,]),
                                       RMSE(test[kid==k,j,], tilde_nn[kid==k,j,])),
                             h = 0
      )
      
      nRMSEtibble <- add_row(nRMSEtibble, 
                             k = k,
                             id = colnames(tilde)[j],
                             comb =  str_sub(i, 2, -1),
                             method =  "t",
                             type = c("tilde", "ddot", "yang", "tilde0", "osqp"),
                             value = c(RMSE(test[kid==k & opday == 1,j,], 
                                            tilde[kid==k & opday == 1,j,]),
                                       RMSE(test[kid==k & opday == 1,j,], 
                                            ddot[kid==k & opday == 1,j,]),
                                       RMSE(test[kid==k & opday == 1,j,], 
                                            yang[kid==k & opday == 1,j,]),
                                       RMSE(test[kid==k & opday == 1,j,], 
                                            tilde0[kid==k & opday == 1,j,]),
                                       RMSE(test[kid==k & opday == 1,j,], 
                                            tilde_nn[kid==k & opday == 1,j,])),
                             h = 1
      )
      
      nRMSEtibble <- add_row(nRMSEtibble, 
                             k = k,
                             id = colnames(tilde)[j],
                             comb =  str_sub(i, 2, -1),
                             method =  "t",
                             type = c("tilde", "ddot", "yang", "tilde0", "osqp"),
                             value = c(RMSE(test[kid==k & opday == 2,j,], 
                                            tilde[kid==k & opday == 2,j,]),
                                       RMSE(test[kid==k & opday == 2,j,], 
                                            ddot[kid==k & opday == 2,j,]),
                                       RMSE(test[kid==k & opday == 2,j,], 
                                            yang[kid==k & opday == 2,j,]),
                                       RMSE(test[kid==k & opday == 2,j,], 
                                            tilde0[kid==k & opday == 2,j,]),
                                       RMSE(test[kid==k & opday == 2,j,], 
                                            tilde_nn[kid==k & opday == 2,j,])),
                             h = 2
      )
    }
    cat(k, " ")
  }
  cat(i, "\n")
}


## BU ----
BU_1 <- base[rep(thf_info$kset, h*thf_info$m/thf_info$kset) == 1,,]
tmp <- apply(BU_1, 3, function(x) as.matrix(thf_info$R %*% x), simplify = FALSE)
BU <- simplify2array(tmp)
for(k in thf_info$kset){
  for(j in 1:NCOL(BU)){
    nRMSEtibble <- add_row(nRMSEtibble, 
                           k = k,
                           id = colnames(BU)[j],
                           comb =  "bu",
                           method =  "t",
                           type = "none",
                           value = c(RMSE(test[kid==k,j,], BU[kid==k,j,])),
                           h = 0
    )
    
    nRMSEtibble <- add_row(nRMSEtibble, 
                           k = k,
                           id = colnames(BU)[j],
                           comb =  "bu",
                           method =  "t",
                           type = "none",
                           value = c(RMSE(test[kid==k & opday == 1,j,], 
                                          BU[kid==k & opday == 1,j,])),
                           h = 1
    )
    
    nRMSEtibble <- add_row(nRMSEtibble, 
                           k = k,
                           id = colnames(BU)[j],
                           comb =  "bu",
                           method =  "t",
                           type = "none",
                           value = c(RMSE(test[kid==k & opday == 2,j,], 
                                          BU[kid==k & opday == 2,j,])),
                           h = 2
    )
  }
  cat(k, " ")
}

## pers ----
meas <- as.data.frame(fread("./meas.csv", header = TRUE, sep = ","))
split_row <- rep(1:365, each = 24)

fmeas <- lapply(1:365, function(x) as.matrix(meas[split_row%in%c(x, x+1),]))
fmeas <- fmeas[c(13:362)]
fmeas <- simplify2array(fmeas)

tmp <- apply(fmeas, 3, function(x){
  out <- t(as.matrix(hts_info$S %*% t(x)))
  colnames(out) <- varname
  rownames(out) <- NULL
  out
}, simplify = FALSE)
tmp <- apply(simplify2array(tmp), 3, function(x) as.matrix(thf_info$R %*% x), simplify = FALSE)
pers <- simplify2array(tmp)
for(k in thf_info$kset){
  for(j in 1:NCOL(pers)){
    nRMSEtibble <- add_row(nRMSEtibble, 
                           k = k,
                           id = colnames(pers)[j],
                           comb =  "pers",
                           method =  "bf",
                           type = "none",
                           value = c(RMSE(test[kid==k,j,], pers[kid==k,j,])),
                           h = 0
    )
    
    nRMSEtibble <- add_row(nRMSEtibble, 
                           k = k,
                           id = colnames(pers)[j],
                           comb =  "pers",
                           method =  "bf",
                           type = "none",
                           value = c(RMSE(test[kid==k & opday == 1,j,], 
                                          pers[kid==k & opday == 1,j,])),
                           h = 1
    )
    
    nRMSEtibble <- add_row(nRMSEtibble, 
                           k = k,
                           id = colnames(pers)[j],
                           comb =  "pers",
                           method =  "bf",
                           type = "none",
                           value = c(RMSE(test[kid==k & opday == 2,j,], 
                                          pers[kid==k & opday == 2,j,])),
                           h = 2
    )
  }
  cat(k, " ")
}

## cross-sectional + bu ----
cont <- c("cols", "cshr", "cstruc", "cwls")
rownames(hts_info$S) <- varname
for(i in cont){
  # free
  tilde <- simplify2array(get(i))
  varname <- colnames(tilde)
  tilde_bts <- tilde[rep(thf_info$kset, h*thf_info$m/thf_info$kset) == 1,-c(1:6),]
  tmp <- apply(tilde_bts, 3, function(x) t(as.matrix(hts_info$S %*% t(x))), simplify = FALSE)
  tmp <- apply(simplify2array(tmp), 3, function(x) as.matrix(thf_info$R %*% x), simplify = FALSE)
  tilde <- simplify2array(tmp)
  colnames(tilde) <- varname
  
  # sntz
  tilde0_bts <- tilde[rep(thf_info$kset, h*thf_info$m/thf_info$kset) == 1,-c(1:6),]
  tilde0_bts[tilde0_bts<0] <- 0
  tmp <- apply(tilde0_bts, 3, function(x) t(as.matrix(hts_info$S %*% t(x))), simplify = FALSE)
  tmp <- apply(simplify2array(tmp), 3, function(x) as.matrix(thf_info$R %*% x), simplify = FALSE)
  tilde0 <- simplify2array(tmp)
  colnames(tilde0) <- varname
  
  for(k in thf_info$kset){
    for(j in 1:NCOL(tilde)){
      nRMSEtibble <- add_row(nRMSEtibble, 
                             k = k,
                             id = colnames(tilde)[j],
                             comb =  paste0("cs-", str_sub(i, 2, -1)),
                             method =  "ctbu",
                             type = c("tilde", "tilde0"),
                             value = c(RMSE(test[kid==k,j,], tilde[kid==k,j,]),
                                       RMSE(test[kid==k,j,], tilde0[kid==k,j,])),
                             h = 0
      )
      
      nRMSEtibble <- add_row(nRMSEtibble, 
                             k = k,
                             id = colnames(tilde)[j],
                             comb =  paste0("cs-", str_sub(i, 2, -1)),
                             method =  "ctbu",
                             type = c("tilde", "tilde0"),
                             value = c(RMSE(test[kid==k & opday == 1,j,], 
                                            tilde[kid==k & opday == 1,j,]),
                                       RMSE(test[kid==k & opday == 1,j,], 
                                            tilde0[kid==k & opday == 1,j,])),
                             h = 1
      )
      
      nRMSEtibble <- add_row(nRMSEtibble, 
                             k = k,
                             id = colnames(tilde)[j],
                             comb =  paste0("cs-", str_sub(i, 2, -1)),
                             method =  "ctbu",
                             type = c("tilde", "tilde0"),
                             value = c(RMSE(test[kid==k & opday == 2,j,], 
                                            tilde[kid==k & opday == 2,j,]),
                                       RMSE(test[kid==k & opday == 2,j,], 
                                            tilde0[kid==k & opday == 2,j,])),
                             h = 2
      )
    }
    cat(k, " ")
  }
  cat(i, "\n")
}

## temporal + bu ----
tempo <- c("tols", "tsar1", "tstruc", "twlsv")
for(i in tempo){
  # free
  tilde <- simplify2array(get(i))
  varname <- colnames(tilde)
  tilde_bts <- tilde[rep(thf_info$kset, h*thf_info$m/thf_info$kset) == 1,-c(1:6),]
  tmp <- apply(tilde_bts, 3, function(x) t(as.matrix(hts_info$S %*% t(x))), simplify = FALSE)
  tmp <- apply(simplify2array(tmp), 3, function(x) as.matrix(thf_info$R %*% x), simplify = FALSE)
  tilde <- simplify2array(tmp)
  colnames(tilde) <- varname
  
  # sntz
  tilde0_bts <- tilde[rep(thf_info$kset, h*thf_info$m/thf_info$kset) == 1,-c(1:6),]
  tilde0_bts[tilde0_bts<0] <- 0
  tmp <- apply(tilde0_bts, 3, function(x) t(as.matrix(hts_info$S %*% t(x))), simplify = FALSE)
  tmp <- apply(simplify2array(tmp), 3, function(x) as.matrix(thf_info$R %*% x), simplify = FALSE)
  tilde0 <- simplify2array(tmp)
  colnames(tilde0) <- varname
  
  for(k in thf_info$kset){
    for(j in 1:NCOL(tilde)){
      nRMSEtibble <- add_row(nRMSEtibble, 
                             k = k,
                             id = colnames(tilde)[j],
                             comb =  paste0("t-", str_sub(i, 2, -1)),
                             method =  "ctbu",
                             type = c("tilde", "tilde0"),
                             value = c(RMSE(test[kid==k,j,], tilde[kid==k,j,]),
                                       RMSE(test[kid==k,j,], tilde0[kid==k,j,])),
                             h = 0
      )
      
      nRMSEtibble <- add_row(nRMSEtibble, 
                             k = k,
                             id = colnames(tilde)[j],
                             comb =  paste0("t-", str_sub(i, 2, -1)),
                             method =  "ctbu",
                             type = c("tilde", "tilde0"),
                             value = c(RMSE(test[kid==k & opday == 1,j,], 
                                            tilde[kid==k & opday == 1,j,]),
                                       RMSE(test[kid==k & opday == 1,j,], 
                                            tilde0[kid==k & opday == 1,j,])),
                             h = 1
      )
      
      nRMSEtibble <- add_row(nRMSEtibble, 
                             k = k,
                             id = colnames(tilde)[j],
                             comb =  paste0("t-", str_sub(i, 2, -1)),
                             method =  "ctbu",
                             type = c("tilde", "tilde0"),
                             value = c(RMSE(test[kid==k & opday == 2,j,], 
                                            tilde[kid==k & opday == 2,j,]),
                                       RMSE(test[kid==k & opday == 2,j,], 
                                            tilde0[kid==k & opday == 2,j,])),
                             h = 2
      )
    }
    cat(k, " ")
  }
  cat(i, "\n")
}

## 3TIER ----
TIER_1 <- base[rep(thf_info$kset, h*thf_info$m/thf_info$kset) == 1,-c(1:6),]
tmp <- apply(TIER_1, 3, function(x) t(as.matrix(hts_info$S %*% t(x))), simplify = FALSE)
tmp <- apply(simplify2array(tmp), 3, function(x) as.matrix(thf_info$R %*% x), simplify = FALSE)
TIER <- simplify2array(tmp)
varname <- colnames(base)
colnames(TIER) <- varname

for(k in thf_info$kset){
  for(j in 1:NCOL(TIER)){
    nRMSEtibble <- add_row(nRMSEtibble, 
                           k = k,
                           id = colnames(TIER)[j],
                           comb =  "3TIER",
                           method =  "ctbu",
                           type = "none",
                           value = c(RMSE(test[kid==k,j,], TIER[kid==k,j,])),
                           h = 0
    )
    
    nRMSEtibble <- add_row(nRMSEtibble, 
                           k = k,
                           id = colnames(TIER)[j],
                           comb =  "3TIER",
                           method =  "ctbu",
                           type = "none",
                           value = c(RMSE(test[kid==k & opday == 1,j,], 
                                          TIER[kid==k & opday == 1,j,])),
                           h = 1
    )
    
    nRMSEtibble <- add_row(nRMSEtibble, 
                           k = k,
                           id = colnames(TIER)[j],
                           comb =  "3TIER",
                           method =  "ctbu",
                           type = "none",
                           value = c(RMSE(test[kid==k & opday == 2,j,], 
                                          TIER[kid==k & opday == 2,j,])),
                           h = 2
    )
  }
  cat(k, " ")
}

## ka ----
ka <- c("twlsv_wls", "twlsv_shr", "tstruc_wls", "tstruc_shr",
         "cwlsv_wls", "cwlsv_shr", "cstruc_wls", "cstruc_shr")
for(i in ka){
  # free
  tilde <- simplify2array(get(i))
  
  # sntz
  tilde0_1 <- tilde[rep(thf_info$kset, h*thf_info$m/thf_info$kset) == 1,,]
  tilde0_1[tilde0_1<0] <- 0
  tmp <- apply(tilde0_1, 3, function(x) as.matrix(thf_info$R %*% x), simplify = FALSE)
  tilde0 <- simplify2array(tmp)
  
  for(k in thf_info$kset){
    for(j in 1:NCOL(tilde)){
      nRMSEtibble <- add_row(nRMSEtibble, 
                             k = k,
                             id = colnames(tilde)[j],
                             comb =  str_sub(i, 2, -1),
                             method =  ifelse(str_sub(i, 1, 1) == "t", "ka (tcs)", "ka (cst)"),
                             type = c("tilde", "tilde0"),
                             value = c(RMSE(test[kid==k,j,], tilde[kid==k,j,]),
                                       RMSE(test[kid==k,j,], tilde0[kid==k,j,])),
                             h = 0
      )
      
      nRMSEtibble <- add_row(nRMSEtibble, 
                             k = k,
                             id = colnames(tilde)[j],
                             comb =  str_sub(i, 2, -1),
                             method =  ifelse(str_sub(i, 1, 1) == "t", "ka (tcs)", "ka (cst)"),
                             type = c("tilde", "tilde0"),
                             value = c(RMSE(test[kid==k & opday == 1,j,], 
                                            tilde[kid==k & opday == 1,j,]),
                                       RMSE(test[kid==k & opday == 1,j,], 
                                            tilde0[kid==k & opday == 1,j,])),
                             h = 1
      )
      
      nRMSEtibble <- add_row(nRMSEtibble, 
                             k = k,
                             id = colnames(tilde)[j],
                             comb =  str_sub(i, 2, -1),
                             method =  ifelse(str_sub(i, 1, 1) == "t", "ka (tcs)", "ka (cst)"),
                             type = c("tilde", "tilde0"),
                             value = c(RMSE(test[kid==k & opday == 2,j,], 
                                            tilde[kid==k & opday == 2,j,]),
                                       RMSE(test[kid==k & opday == 2,j,], 
                                            tilde0[kid==k & opday == 2,j,])),
                             h = 2
      )
    }
    cat(k, " ")
  }
  cat(i, "\n")
}


## ite ----
ite <- c("tols_struc", "tstruc_ols", "twlsv_wls", "twlsv_shr", "tstruc_wls", "twlsv_struc", "tstruc_shr",
           "cols_struc", "cstruc_ols", "cwlsv_wls", "cwlsv_shr", "cstruc_wls", "cwlsv_struc", "cstruc_shr")
for(i in ite){
  # free
  tilde <- simplify2array(get(i))
  
  # sntz
  tilde0_bts <- tilde[kid == 1,-c(1:6),]
  tilde0_bts[tilde0_bts<0] <- 0
  tmp <- apply(tilde0_bts, 3, function(x) t(as.matrix(hts_info$S %*% t(x))), simplify = FALSE)
  tmp <- apply(simplify2array(tmp), 3, function(x) as.matrix(thf_info$R %*% x), simplify = FALSE)
  tilde0 <- simplify2array(tmp)
  
  for(k in thf_info$kset){
    for(j in 1:NCOL(tilde)){
      nRMSEtibble <- add_row(nRMSEtibble, 
                             k = k,
                             id = colnames(tilde)[j],
                             comb =  str_sub(i, 2, -1),
                             method =  ifelse(str_sub(i, 1, 1) == "t", "ite (tcs)", "ite (cst)"),
                             type = c("tilde", "tilde0"),
                             value = c(RMSE(test[kid==k,j,], tilde[kid==k,j,]),
                                       RMSE(test[kid==k,j,], tilde0[kid==k,j,])),
                             h = 0
      )
      
      nRMSEtibble <- add_row(nRMSEtibble, 
                             k = k,
                             id = colnames(tilde)[j],
                             comb =  str_sub(i, 2, -1),
                             method =  ifelse(str_sub(i, 1, 1) == "t", "ite (tcs)", "ite (cst)"),
                             type = c("tilde", "tilde0"),
                             value = c(RMSE(test[kid==k & opday == 1,j,], 
                                            tilde[kid==k & opday == 1,j,]),
                                       RMSE(test[kid==k & opday == 1,j,], 
                                            tilde0[kid==k & opday == 1,j,])),
                             h = 1
      )
      
      nRMSEtibble <- add_row(nRMSEtibble, 
                             k = k,
                             id = colnames(tilde)[j],
                             comb =  str_sub(i, 2, -1),
                             method =  ifelse(str_sub(i, 1, 1) == "t", "ite (tcs)", "ite (cst)"),
                             type = c("tilde", "tilde0"),
                             value = c(RMSE(test[kid==k & opday == 2,j,], 
                                            tilde[kid==k & opday == 2,j,]),
                                       RMSE(test[kid==k & opday == 2,j,], 
                                            tilde0[kid==k & opday == 2,j,])),
                             h = 2
      )
    }
    cat(k, " ")
  }
  cat(i, "\n")
}

## ite tol ----
ite <- c("twlsv_wls_tol1", "twlsv_wls_tol2", "twlsv_wls_tol3", 
         "cwlsv_wls_tol1", "cwlsv_wls_tol2", "cwlsv_wls_tol3")
for(i in ite){
  # free
  tilde <- simplify2array(get(i))
  
  # sntz
  tilde0_bts <- tilde[kid == 1,-c(1:6),]
  tilde0_bts[tilde0_bts<0] <- 0
  tmp <- apply(tilde0_bts, 3, function(x) t(as.matrix(hts_info$S %*% t(x))), simplify = FALSE)
  tmp <- apply(simplify2array(tmp), 3, function(x) as.matrix(thf_info$R %*% x), simplify = FALSE)
  tilde0 <- simplify2array(tmp)
  
  for(k in thf_info$kset){
    for(j in 1:NCOL(tilde)){
      nRMSEtibble <- add_row(nRMSEtibble, 
                             k = k,
                             id = colnames(tilde)[j],
                             comb =  str_sub(i, 2, -1),
                             method =  ifelse(str_sub(i, 1, 1) == "t", "ite (tcs)", "ite (cst)"),
                             type = c("tilde", "tilde0"),
                             value = c(RMSE(test[kid==k,j,], tilde[kid==k,j,]),
                                       RMSE(test[kid==k,j,], tilde0[kid==k,j,])),
                             h = 0
      )
      
      nRMSEtibble <- add_row(nRMSEtibble, 
                             k = k,
                             id = colnames(tilde)[j],
                             comb =  str_sub(i, 2, -1),
                             method =  ifelse(str_sub(i, 1, 1) == "t", "ite (tcs)", "ite (cst)"),
                             type = c("tilde", "tilde0"),
                             value = c(RMSE(test[kid==k & opday == 1,j,], 
                                            tilde[kid==k & opday == 1,j,]),
                                       RMSE(test[kid==k & opday == 1,j,], 
                                            tilde0[kid==k & opday == 1,j,])),
                             h = 1
      )
      
      nRMSEtibble <- add_row(nRMSEtibble, 
                             k = k,
                             id = colnames(tilde)[j],
                             comb =  str_sub(i, 2, -1),
                             method =  ifelse(str_sub(i, 1, 1) == "t", "ite (tcs)", "ite (cst)"),
                             type = c("tilde", "tilde0"),
                             value = c(RMSE(test[kid==k & opday == 2,j,], 
                                            tilde[kid==k & opday == 2,j,]),
                                       RMSE(test[kid==k & opday == 2,j,], 
                                            tilde0[kid==k & opday == 2,j,])),
                             h = 2
      )
    }
    cat(k, " ")
  }
  cat(i, "\n")
}

## seqL2 ----
ite <- c("tols2_ols", "tstruc2_struc", "tols2_struc", "tstruc2_ols")
rownames(hts_info$S) <- varname
for(i in ite){
  # free
  tilde <- simplify2array(get(i))
  
  # sntz
  tilde_bts <- tilde[,-c(1:6),]
  tilde_bts[tilde_bts<0] <- 0
  tmp <- apply(tilde_bts, 3, function(x) t(as.matrix(hts_info$S %*% t(x))), simplify = FALSE)
  tilde0 <- simplify2array(tmp)
  
  for(k in thf_info$kset){
    for(j in 1:NCOL(tilde)){
      nRMSEtibble <- add_row(nRMSEtibble, 
                             k = k,
                             id = colnames(tilde)[j],
                             comb =  str_sub(i, 2, -1),
                             method =  ifelse(str_sub(i, 1, 1) == "t", "ite (tcs)", "ite (cst)"),
                             type = c("tilde", "tilde0"),
                             value = c(RMSE(test[kid==k,j,], tilde[kid==k,j,]),
                                       RMSE(test[kid==k,j,], tilde0[kid==k,j,])),
                             h = 0
      )
      
      nRMSEtibble <- add_row(nRMSEtibble, 
                             k = k,
                             id = colnames(tilde)[j],
                             comb =  str_sub(i, 2, -1),
                             method =  ifelse(str_sub(i, 1, 1) == "t", "ite (tcs)", "ite (cst)"),
                             type = c("tilde", "tilde0"),
                             value = c(RMSE(test[kid==k & opday == 1,j,], 
                                            tilde[kid==k & opday == 1,j,]),
                                       RMSE(test[kid==k & opday == 1,j,], 
                                            tilde0[kid==k & opday == 1,j,])),
                             h = 1
      )
      
      nRMSEtibble <- add_row(nRMSEtibble, 
                             k = k,
                             id = colnames(tilde)[j],
                             comb =  str_sub(i, 2, -1),
                             method =  ifelse(str_sub(i, 1, 1) == "t", "ite (tcs)", "ite (cst)"),
                             type = c("tilde", "tilde0"),
                             value = c(RMSE(test[kid==k & opday == 2,j,], 
                                            tilde[kid==k & opday == 2,j,]),
                                       RMSE(test[kid==k & opday == 2,j,], 
                                            tilde0[kid==k & opday == 2,j,])),
                             h = 2
      )
    }
    cat(k, " ")
  }
  cat(i, "\n")
}

nRMSEdata <- nRMSEtibble %>%
  mutate(group = ifelse(id %in% varname[c(1)], "L1",
                        ifelse(id %in% varname[c(2:6)], "L2", "L3")))
save(nRMSEdata, file = "./data/nRMSEdata.RData")

