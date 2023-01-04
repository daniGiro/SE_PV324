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
comb <- c("ols", "struc", "wls", "shr")
pb <- txtProgressBar(max = max(rep), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
results <- foreach(rp = rep, .options.snow = opts) %dopar% {
  all <- load_replication(i = rp)
  free <- list()
  free[comb] <- rep(list(matrix(NA, ncol = NCOL(all$Yhat), nrow = NROW(all$Yhat), 
                                dimnames = dimnames(all$Yhat))), length(comb))
  nn <- free
  
  eid <- rep(thf_info$kset, (NROW(all$E)/thf_info$kt)*thf_info$m/thf_info$kset)
  yid <- rep(thf_info$kset, h*thf_info$m/thf_info$kset)
  
  index_nn <- list()
  index_nn[comb] <- rep(list(setNames(rep(NA, length(thf_info$kset)), thf_info$kset)), length(comb))
  
  pri_res <- list()
  pri_res[comb] <- rep(list(setNames(rep(NA, length(thf_info$kset)), thf_info$kset)), length(comb))
  
  time <- list()
  time[c("nn", "free")] <- NULL
  time[["nn"]][comb] <- rep(list(setNames(rep(NA, length(thf_info$kset)), thf_info$kset)), length(comb))
  time[["free"]][comb] <- rep(list(setNames(rep(NA, length(thf_info$kset)), thf_info$kset)), length(comb))
  for(k in thf_info$kset){
    for(j in comb){
      Start <- Sys.time()
      free[[j]][yid==k,] <- drop_zeros(FoReco::htsrec(basef = all$Yhat[yid==k,,drop = FALSE], 
                                                      comb = j, C = hts_info$C, 
                                                      res = all$E[eid==k,,drop = FALSE], 
                                                      keep = "recf"))
      time[["free"]][[j]][thf_info$kset==k] <- difftime(Sys.time(), Start, units = "secs")
      
      if(any(free[[j]][yid==k,]<0)){
        Start <- Sys.time()
        tmp <- FoReco::htsrec(basef = all$Yhat[yid==k,,drop = FALSE], comb = j,
                              C = hts_info$C, res = all$E[eid==k,,drop = FALSE], 
                              keep = "recf", nn = TRUE)
        time[["nn"]][[j]][thf_info$kset==k] <- difftime(Sys.time(), Start, units = "secs")
        
        if(is.list(tmp)){
          nn[[j]][yid==k,] <- drop_zeros(tmp$recf, tol = 10*max(tmp$info[,"pri_res"]))
          index_nn[[j]][thf_info$kset==k] <- ifelse(any(tmp$info[,"status"]!=1), 
                                                    tmp$info[tmp$info[,"status"]!=1,"status"][1], 1)
          pri_res[[j]][thf_info$kset==k] <- max(tmp$info[,"pri_res"])
        }else{
          nn[[j]][yid==k,] <- drop_zeros(tmp)
          index_nn[[j]][yid==k] <- 99
        }
      }
    }
  }
  file.name <- paste("./results_hts/rep--", rp, "--hts--reco.RData", sep = "")
  save(free, nn, index_nn, pri_res, time, file = file.name)
}