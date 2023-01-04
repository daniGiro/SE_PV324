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
comb <- c("ols", "struc", "wlsv", "bdshr")
pb <- txtProgressBar(max = max(rep), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
results <- foreach(rp = rep, .options.snow = opts) %dopar% {
  all <- load_replication(i = rp)
  free <- list()
  free[comb] <- rep(list(matrix(NA, ncol = NCOL(all$Yhat), nrow = NROW(all$Yhat), 
                                dimnames = dimnames(all$Yhat))), length(comb))
  nn <- free
  
  index_nn <- list()
  index_nn[comb] <- rep(list(NA), length(comb))
  
  pri_res <- list()
  pri_res[comb] <- rep(list(NA), length(comb))
  
  time <- list()
  time[c("nn", "free")] <- NULL
  time[["nn"]][comb] <- rep(list(NA), length(comb))
  time[["free"]][comb] <- rep(list(NA), length(comb))
  for(j in comb){
    Start <- Sys.time()
    free[[j]] <- drop_zeros(t(FoReco::octrec(basef = t(all$Yhat), comb = j, m = thf_info$m,
                        C = hts_info$C, res = t(all$E), keep = "recf", 
                        type = ifelse(j %in% c("ols", "struc", "wlsv"), "M", "S"))))
    time[["free"]][[j]] <- difftime(Sys.time(), Start, units = "secs")
    
    if(any(free[[j]]<0)){
      Start <- Sys.time()
      tmp <- FoReco::octrec(basef = t(all$Yhat), comb = j, m = thf_info$m,
                    C = hts_info$C, res = t(all$E), keep = "recf", sol = "osqp", 
                    type = ifelse(j %in% c("ols", "struc", "wlsv"), "M", "S"),
                    nn = TRUE)
      time[["nn"]][[j]] <- difftime(Sys.time(), Start, units = "secs")
      
      if(is.list(tmp)){
        nn[[j]] <- drop_zeros(t(tmp$recf), tol = 10*max(tmp$info[,"pri_res"]))
        index_nn[[j]] <- ifelse(any(tmp$info[,"status"]!=1), 
                                tmp$info[tmp$info[,"status"]!=1,"status"][1], 1)
        pri_res[[j]] <- max(tmp$info[,"pri_res"])
      }else{
        nn[[j]][yid==k,] <- drop_zeros(t(tmp))
        index_nn[[j]] <- 99
      }
    }
  }
  file.name <- paste("./results_oct/rep--", rp, "--oct--reco.RData", sep = "")
  save(free, nn, index_nn, time, pri_res, file = file.name)
}