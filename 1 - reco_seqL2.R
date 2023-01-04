#Clear all workspace
rm(list = ls(all = TRUE))
libs <- c("FoReco", "doSNOW", "progress", "doParallel")
invisible(lapply(libs, library, character.only = TRUE))

source("fun.R")
source("varname.R")
load("info_reco.RData")

rownames(hts_info$S) <- varname
ncores <- detectCores()-1
cat("\n\n\n#############################################\n\n",
    "Cores that we will use: ", ncores,"\n",
    "\n#############################################\n\n\n", sep = "")

cl <- makeCluster(ncores)
registerDoSNOW(cl)
rep = 1:350
#rep = 1:2
comb <- c("ols-ols", "ols-struc", "struc-struc", "struc-ols")
pb <- txtProgressBar(max = max(rep), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
results <- foreach(rp = rep, .options.snow = opts) %dopar% {
  all <- load_replication(i = rp)
  tsr <- list()
  tsr[comb] <- rep(list(matrix(NA, ncol = NCOL(all$Yhat), nrow = NROW(all$Yhat), 
                               dimnames = dimnames(all$Yhat))), length(comb))
  
  time <- list()
  time[["tsr"]][comb] <- rep(list(NA), length(comb))
  
  for(j in comb){
    splitcomb <- strsplit(j, "-")[[1]]
    
    Start <- Sys.time()
    step1 <- apply(all$Yhat[,-c(1:6)], 2, FoReco::thfrec, comb = splitcomb[1], m = thf_info$m, keep = "recf")
    step1 <- cbind(all$Yhat[,c(1:6)], step1)
    step2 <- FoReco::htsrec(step1, comb = splitcomb[2], C = hts_info$C, keep = "recf")
    
    tsr[[j]] <- drop_zeros(step2)
    time[["tsr"]][[j]] <- difftime(Sys.time(), Start, units = "secs")
  }
  file.name <- paste("./results_seqL2/rep--", rp, "--seqL2--reco.RData", sep = "")
  save(tsr, time, file = file.name)
}