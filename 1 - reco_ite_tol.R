#Clear all workspace
rm(list = ls(all = TRUE))
libs <- c("FoReco", "doSNOW", "progress", "doParallel")
invisible(lapply(libs, library, character.only = TRUE))

source("fun.R")
load("info_reco.RData")

ncores <- detectCores()-1
cat("\n\n\n#############################################\n\n",
    "Cores that we will use: ", ncores,"\n",
    "\n#############################################\n\n\n", sep = "")

cl <- makeCluster(ncores)
registerDoSNOW(cl)

rep = 1:350
#rep = 1:2
comb <- as.character(1:3)
tol <- c(1e-6, 1e-8, 1e-10)
pb <- txtProgressBar(max = max(rep), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
results <- foreach(rp = rep, .options.snow = opts) %dopar% {
  all <- load_replication(i = rp)
  tcs <- list()
  tcs[comb] <- rep(list(matrix(NA, ncol = NCOL(all$Yhat), nrow = NROW(all$Yhat), 
                               dimnames = dimnames(all$Yhat))), length(comb))
  cst <- tcs
  
  time <- list()
  time[["tcs"]][comb] <- rep(list(NA), length(comb))
  time[["cst"]][comb] <- rep(list(NA), length(comb))
  
  ite <- list()
  ite[["tcs"]][comb] <- rep(list(NA), length(comb))
  ite[["cst"]][comb] <- rep(list(NA), length(comb))
  
  diff <- list()
  for(j in comb){
    Start <- Sys.time()
    tmp <- FoReco::iterec(basef = t(all$Yhat), start_rec = "thf",
                          thf_comb = "wlsv",
                          hts_comb = "wls", 
                          tol = tol[as.numeric(j)],
                          m = thf_info$m, note = FALSE,
                          C = hts_info$C, res = t(all$E), keep = "recf")
    tcs[[j]] <- drop_zeros(t(tmp$recf))
    time[["tcs"]][[j]] <- difftime(Sys.time(), Start, units = "secs")
    ite[["tcs"]][[j]] <- tmp$iter
    
    Start <- Sys.time()
    tmp1 <- FoReco::iterec(basef = t(all$Yhat), start_rec = "hts",
                           thf_comb = "wlsv",
                           hts_comb = "wls", 
                           tol = tol[as.numeric(j)],
                           m = thf_info$m, note = FALSE, 
                           C = hts_info$C, res = t(all$E), keep = "recf")
    cst[[j]] <- drop_zeros(t(tmp1$recf))
    time[["cst"]][[j]] <- difftime(Sys.time(), Start, units = "secs")
    ite[["cst"]][[j]] <- tmp1$iter
    
    diff[[j]] <- tmp1$recf-tmp$recf
  }
  file.name <- paste("./results_ite_tol/rep--", rp, "--ite_tol--reco.RData", sep = "")
  save(cst, tcs, time, ite, diff, file = file.name)
}