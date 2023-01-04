# clear ws
rm(list = ls(all = TRUE))

# Load data
load("info_reco.RData")

# source file
source("varname.R")
source("fun.R")

# base and test ----
fls <- list.files("./Results_Time", full.names = TRUE)
base <- array(NA, dim = c(thf_info$kt*h, length(varname), 350), 
              dimnames = list(NULL, varname, NULL))
test <- array(NA, dim = c(thf_info$kt*h, length(varname), 350), 
              dimnames = list(NULL, varname, NULL))
for(j in fls){
  all <- mget(load(j, envir=(NE. <- new.env())), envir=NE.)
  rm(NE.)
  base[,strsplit(basename(j), "--")[[1]][1],] <- apply(all$Y.hat, 2, 
                                                       function(x) unname(unlist(rev(split(matrix(x, ncol = 2), 
                                                                                           rep(thf_info$kset, 24/thf_info$kset))))))
  test[,strsplit(basename(j), "--")[[1]][1],] <- as.matrix(thf_info$R %*% all$Y.obs)
  cat(".")
}
save(base, test, file = "./data/fore_base_test.RData")

# temporal ----
fls <- list.files('./results_thf', full.names = TRUE)
tols <- vector(mode = "list", length = 350)
tols_nn <- vector(mode = "list", length = 350)
tstruc <- vector(mode = "list", length = 350)
tstruc_nn <- vector(mode = "list", length = 350)
twlsv <- vector(mode = "list", length = 350)
twlsv_nn <- vector(mode = "list", length = 350)
tsar1 <- vector(mode = "list", length = 350)
tsar1_nn <- vector(mode = "list", length = 350)
for(j in fls){
  all <- mget(load(j, envir=(NE. <- new.env())), envir=NE.)
  rm(NE.)
  tmp <- all$free[['ols']]
  tols[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- tmp
  tmp[,!is.na(all$index_nn[['ols']])] <- all$nn[['ols']][,!is.na(all$index_nn[['ols']])]
  tols_nn[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- tmp
  
  tmp <- all$free[['struc']]
  tstruc[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- tmp
  tmp[,!is.na(all$index_nn[['struc']])] <- all$nn[['struc']][,!is.na(all$index_nn[['struc']])]
  tstruc_nn[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- tmp
  
  tmp <- all$free[['wlsv']]
  twlsv[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- tmp
  tmp[,!is.na(all$index_nn[['wlsv']])] <- all$nn[['wlsv']][,!is.na(all$index_nn[['wlsv']])]
  twlsv_nn[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- tmp
  
  tmp <- all$free[['sar1']]
  tsar1[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- tmp
  tmp[,!is.na(all$index_nn[['sar1']])] <- all$nn[['sar1']][,!is.na(all$index_nn[['sar1']])]
  tsar1_nn[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- tmp
  cat(".")
}
save(tols, tols_nn, tstruc, tstruc_nn, twlsv, twlsv_nn, tsar1, tsar1_nn, 
     file = "./data/fore_treco.RData")

# cross-sectional ----
fls <- list.files('./results_hts', full.names = TRUE)
cols <- vector(mode = "list", length = 350)
cols_nn <- vector(mode = "list", length = 350)
cstruc <- vector(mode = "list", length = 350)
cstruc_nn <- vector(mode = "list", length = 350)
cwls <- vector(mode = "list", length = 350)
cwls_nn <- vector(mode = "list", length = 350)
cshr <- vector(mode = "list", length = 350)
cshr_nn <- vector(mode = "list", length = 350)
for(j in fls){
  all <- mget(load(j, envir=(NE. <- new.env())), envir=NE.)
  rm(NE.)
  tmp <- all$free[['ols']]
  cols[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- tmp
  tmp[!is.na(all$nn[['ols']])] <- all$nn[['ols']][!is.na(all$nn[['ols']])]
  cols_nn[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- tmp
  
  tmp <- all$free[['struc']]
  cstruc[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- tmp
  tmp[!is.na(all$nn[['struc']])] <- all$nn[['struc']][!is.na(all$nn[['struc']])]
  cstruc_nn[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- tmp
  
  tmp <- all$free[['wls']]
  cwls[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- tmp
  tmp[!is.na(all$nn[['wls']])] <- all$nn[['wls']][!is.na(all$nn[['wls']])]
  cwls_nn[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- tmp
  
  tmp <- all$free[['shr']]
  cshr[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- tmp
  tmp[!is.na(all$nn[['shr']])] <- all$nn[['shr']][!is.na(all$nn[['shr']])]
  cshr_nn[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- tmp
  cat(".")
}
save(cols, cols_nn, cstruc, cstruc_nn, cwls, cwls_nn, cshr, cshr_nn, 
     file = "./data/fore_csreco.RData")

# cross-temporal ----
fls <- list.files('./results_oct', full.names = TRUE)
ctols <- vector(mode = "list", length = 350)
ctols_nn <- vector(mode = "list", length = 350)
ctstruc <- vector(mode = "list", length = 350)
ctstruc_nn <- vector(mode = "list", length = 350)
ctwlsv <- vector(mode = "list", length = 350)
ctwlsv_nn <- vector(mode = "list", length = 350)
ctbdshr <- vector(mode = "list", length = 350)
ctbdshr_nn <- vector(mode = "list", length = 350)
for(j in fls){
  all <- mget(load(j, envir=(NE. <- new.env())), envir=NE.)
  rm(NE.)
  tmp <- all$free[['ols']]
  ctols[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- tmp
  if(!all(is.na(all$nn[['ols']]))){
    tmp <- all$nn[['ols']]
  }
  ctols_nn[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- tmp
  
  tmp <- all$free[['struc']]
  ctstruc[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- tmp
  if(!all(is.na(all$nn[['struc']]))){
    tmp <- all$nn[['struc']]
  }
  ctstruc_nn[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- tmp
  
  tmp <- all$free[['wlsv']]
  ctwlsv[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- tmp
  if(!all(is.na(all$nn[['wlsv']]))){
    tmp <- all$nn[['wlsv']]
  }
  ctwlsv_nn[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- tmp
  
  tmp <- all$free[['bdshr']]
  ctbdshr[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- tmp
  if(!all(is.na(all$nn[['bdshr']]))){
    tmp <- all$nn[['bdshr']]
  }
  ctbdshr_nn[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- tmp
  cat(".")
}
save(ctols, ctols_nn, ctstruc, ctstruc_nn, ctwlsv, ctwlsv_nn, ctbdshr, ctbdshr_nn, 
     file = "./data/fore_ctreco.RData")

# ite ----
fls <- list.files('./results_ite', full.names = TRUE)
tols_struc <- tstruc_ols <- vector(mode = "list", length = 350)
twlsv_wls <- twlsv_shr <- vector(mode = "list", length = 350)
tstruc_wls <- twlsv_struc <- tstruc_shr <- vector(mode = "list", length = 350)
cols_struc <- cstruc_ols <- vector(mode = "list", length = 350)
cwlsv_wls <- cwlsv_shr <- vector(mode = "list", length = 350)
cstruc_wls <- cwlsv_struc <- cstruc_shr <- vector(mode = "list", length = 350)
ite <- vector(mode = "list", length = 350)
diff <- vector(mode = "list", length = 350)
for(j in fls){
  all <- mget(load(j, envir=(NE. <- new.env())), envir=NE.)
  rm(NE.)
  tols_struc[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$tcs[['ols-struc']]
  tstruc_ols[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$tcs[['struc-ols']]
  twlsv_wls[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$tcs[['wlsv-wls']]
  twlsv_shr[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$tcs[['wlsv-shr']]
  tstruc_wls[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$tcs[['struc-wls']]
  twlsv_struc[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$tcs[['wlsv-struc']]
  tstruc_shr[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$tcs[['struc-shr']]
  
  cols_struc[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$cst[['ols-struc']]
  cstruc_ols[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$cst[['struc-ols']]
  cwlsv_wls[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$cst[['wlsv-wls']]
  cwlsv_shr[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$cst[['wlsv-shr']]
  cstruc_wls[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$cst[['struc-wls']]
  cwlsv_struc[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$cst[['wlsv-struc']]
  cstruc_shr[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$cst[['struc-shr']]
  ite[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- cbind(reshape::melt(all$ite), id = as.numeric(strsplit(basename(j), "--")[[1]][2]))
  diff[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- simplify2array(as.matrix(all$diff))
  cat(".")
}
save(tols_struc, tstruc_ols, twlsv_wls, twlsv_shr, tstruc_wls, twlsv_struc, tstruc_shr,
     cols_struc, cstruc_ols, cwlsv_wls, cwlsv_shr, cstruc_wls, cwlsv_struc, cstruc_shr, file = "./data/fore_ite.RData")

# KA ----
fls <- list.files('./results_ka', full.names = TRUE)
tstruc_wls <- tstruc_shr <- vector(mode = "list", length = 350)
twlsv_wls <- twlsv_shr <- vector(mode = "list", length = 350)
cstruc_wls <- cstruc_shr <- vector(mode = "list", length = 350)
cwlsv_wls <- cwlsv_shr <- vector(mode = "list", length = 350)
ite <- vector(mode = "list", length = 350)
diff <- vector(mode = "list", length = 350)
for(j in fls){
  all <- mget(load(j, envir=(NE. <- new.env())), envir=NE.)
  rm(NE.)
  twlsv_wls[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$tcs[['wlsv-wls']]
  twlsv_shr[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$tcs[['wlsv-shr']]
  tstruc_wls[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$tcs[['struc-wls']]
  tstruc_shr[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$tcs[['struc-shr']]
  
  cwlsv_wls[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$cst[['wlsv-wls']]
  cwlsv_shr[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$cst[['wlsv-shr']]
  cstruc_wls[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$cst[['struc-wls']]
  cstruc_shr[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$cst[['struc-shr']]
  diff[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- simplify2array(as.matrix(all$diff))
  cat(".")
}
save(twlsv_wls, twlsv_shr, tstruc_wls, tstruc_shr,
     cwlsv_wls, cwlsv_shr, cstruc_wls, cstruc_shr, file = "./data/fore_ka.RData")

# seqL2 ----
fls <- list.files('./results_seqL2', full.names = TRUE)
tols2_ols <- tstruc2_struc <- vector(mode = "list", length = 350)
tols2_struc <- tstruc2_ols <- vector(mode = "list", length = 350)
for(j in fls){
  all <- mget(load(j, envir=(NE. <- new.env())), envir=NE.)
  rm(NE.)
  tols2_ols[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$tsr[['ols-ols']]
  tstruc2_struc[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$tsr[['struc-struc']]
  tols2_struc[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$tsr[['ols-struc']]
  tstruc2_ols[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$tsr[['struc-ols']]
  cat(".")
}
save(tols2_ols, tstruc2_struc, tols2_struc, tstruc2_ols, file = "./data/fore_seqL2.RData")


# ite tol ----
fls <- list.files('./results_ite_tol', full.names = TRUE)
twlsv_wls_tol1 <- twlsv_wls_tol2 <- twlsv_wls_tol3 <- vector(mode = "list", length = 350)
cwlsv_wls_tol1 <- cwlsv_wls_tol2 <- cwlsv_wls_tol3 <- vector(mode = "list", length = 350)
ite <- vector(mode = "list", length = 350)
diff <- vector(mode = "list", length = 350)
for(j in fls){
  all <- mget(load(j, envir=(NE. <- new.env())), envir=NE.)
  rm(NE.)
  twlsv_wls_tol1[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$tcs$`1`
  twlsv_wls_tol2[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$tcs$`2`
  twlsv_wls_tol3[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$tcs$`3`
  
  cwlsv_wls_tol1[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$cst$`1`
  cwlsv_wls_tol2[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$cst$`2`
  cwlsv_wls_tol3[[as.numeric(strsplit(basename(j), "--")[[1]][2])]] <- all$cst$`3`
  
  cat(".")
}
save(twlsv_wls_tol1, twlsv_wls_tol2, twlsv_wls_tol3, 
     cwlsv_wls_tol1, cwlsv_wls_tol2, cwlsv_wls_tol3, file = "./data/fore_ite_tol.RData")

