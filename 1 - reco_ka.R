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
comb <- c("struc-wls", "struc-shr", "wlsv-wls", "wlsv-shr")
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
  time[c("cst", "tcs")] <- NULL
  time[["tcs"]][comb] <- rep(list(NA), length(comb))
  time[["cst"]][comb] <- rep(list(NA), length(comb))
  
  diff <- list()
  for(j in comb){
    splitcomb <- strsplit(j, "-")[[1]]
    Start <- Sys.time()
    tmp <- FoReco::tcsrec(basef = t(all$Yhat),
                          thf_comb = splitcomb[1],
                          hts_comb = splitcomb[2], 
                          m = thf_info$m, 
                          C = hts_info$C, res = t(all$E))
    tcs[[j]] <- drop_zeros(t(tmp$recf))
    time[["tcs"]][[j]] <- difftime(Sys.time(), Start, units = "secs")
    
    Start <- Sys.time()
    tmp1 <- FoReco::cstrec(basef = t(all$Yhat),
                           thf_comb = splitcomb[1],
                           hts_comb = splitcomb[2], 
                           m = thf_info$m, 
                           C = hts_info$C, res = t(all$E))
    cst[[j]] <- drop_zeros(t(tmp1$recf))
    time[["cst"]][[j]] <- difftime(Sys.time(), Start, units = "secs")
    
    diff[[j]] <- tmp1$recf-tmp$recf
  }
  file.name <- paste("./results_ka/rep--", rp, "--ka--reco.RData", sep = "")
  save(cst, tcs, time, diff, file = file.name)
}