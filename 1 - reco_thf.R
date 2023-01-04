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
comb <- c("ols", "struc", "wlsv", "sar1")
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
  index_nn[comb] <- rep(list(rep(NA, NCOL(all$Yhat))), length(comb))
  
  pri_res <- list()
  pri_res[comb] <- rep(list(rep(NA, NCOL(all$Yhat))), length(comb))
  
  time <- list()
  time[c("nn", "free")] <- NULL
  time[["nn"]][comb] <- rep(list(rep(NA, NCOL(all$Yhat))), length(comb))
  time[["free"]][comb] <- rep(list(rep(NA, NCOL(all$Yhat))), length(comb))
  for(j in 1:NCOL(all$Yhat)){
    for(k in comb){
      Start <- Sys.time()
      free[[k]][,j] <- drop_zeros(FoReco::thfrec(basef = all$Yhat[,j], m = thf_info$m, 
                                                 comb = k, res = all$E[,j], keep = "recf"))
      time[["free"]][[k]][j] <- difftime(Sys.time(), Start, units = "secs")
      
      if(any(free[[k]][,j]<0)){
        Start <- Sys.time()
        tmp <- FoReco::thfrec(basef = all$Yhat[,j], m = thf_info$m, comb = k, 
                              res = all$E[,j], keep = "recf", nn = TRUE)
        time[["nn"]][[k]][j] <- difftime(Sys.time(), Start, units = "secs")
        if(is.list(tmp)){
          nn[[k]][,j] <- drop_zeros(tmp$recf, tol = 10*max(tmp$info[,"pri_res"]))
          index_nn[[k]][j] <- ifelse(any(tmp$info[,"status"]!=1), tmp$info[tmp$info[,"status"]!=1,"status"][1], 1)
          pri_res[[k]][j] <- max(tmp$info[,"pri_res"])
        }else{
          nn[[k]][,j] <- drop_zeros(tmp)
          index_nn[[k]][j] <- 99
        }
      }
    }
  }
  file.name <- paste("./results_thf/rep--", rp, "--thf--reco.RData", sep = "")
  save(free, nn, index_nn, time, pri_res, file = file.name)
}
