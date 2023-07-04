Add_noise <- function(mtx, noise){
  mtx_noise <- matrix(rnorm(length(mtx), 0, noise), nrow = nrow(mtx))
  mtx <- mtx + mtx_noise
  for(i in 1:nrow(mtx)){
    for(j in 1:ncol(mtx)){
      if(mtx[i,j]<0) mtx[i,j] <- 0
    }
  }
  return(mtx)
}


Sens.Deconv <- function(bulk, ref, B, cell_type, type_order, real_prop,
                        lambda.set = seq(10,200,10), Pi.start=NULL, max.iter=1e4, tol=1e-6,
                        test = "wilcox", nmrk = 50, seed = 123, noises = seq(0,2.4,0.1), parallel = T, os = "win", ncore = 60 ){
  source("/storage/windows-backup/D drive/Penghui/Sensitivity Analysis/HiDecon all fun.R")
  source("/storage/windows-backup/D drive/Penghui/Sensitivity Analysis/Other methods.R")
  library(parallel)
  
  if(parallel){
    
    pb <- progress_bar$new(
      format = "Current Sample : :current [:bar] :elapsed | percent: :percent",
      total = length(noises),
      clear = FALSE,
      force = TRUE,
      width = 60)
    
    progress_letter <- rep(1:10, 10)  # token reported in progress bar
    
    # allowing progress bar to be used in foreach -----------------------------
    progress <- function(n){
      pb$tick(tokens = list(letter = progress_letter[n]))
    }
    
    
    opts <- list(progress = progress)
    
    if(os == "win"){
      cl = makeCluster(ncore, outfile="")
    }else{
      cl = makeCluster(ncore, setup_strategy = "sequential")
    }
    
    
    registerDoSNOW(cl)
    clusterCall(cl, function(x) .libPaths(x), .libPaths())
    result_list = foreach(i = 1:length(noises),.options.snow = opts, .errorhandling='pass') %dopar% {
      source("/storage/windows-backup/D drive/Penghui/Sensitivity Analysis/HiDecon all fun.R")
      source("/storage/windows-backup/D drive/Penghui/Sensitivity Analysis/Other methods.R")
      
      Other_res <- Other_methods(bulk = bulk, ref = ref, cell_type = cell_type,
                                 real_prop = real_prop, noise = noises[i], seed = seed)
      HiDecon_res <- cv_HiDecon(bulk = bulk, ref = ref, B = B, cell_type = cell_type, type_order, noise = noises[i], seed = seed)$cv_res
      rownames(HiDecon_res) <- rownames(Other_res[[1]])
      colnames(HiDecon_res) <- colnames(Other_res[[1]])
      All_res <- list("HiDecon" = HiDecon_res, "CIBERSORT" = Other_res$CIBERSORT, "dtangle" = Other_res$dtangle,
                      "MuSiC" = Other_res$MuSiC)
      return(All_res)
    }
    stopCluster(cl)
  }else{
    result_list <- list()
    for(i in 1:length(noises)){
      cat("\n Noise sd:" , noises[i])
      Other_res <- Other_methods(bulk = bulk, ref = ref, cell_type = cell_type,
                                 real_prop = real_prop, noise = noises[i], seed = seed)
      HiDecon_res <- cv_HiDecon(bulk = bulk, ref = ref, B = B, cell_type = cell_type, type_order, noise = noises[i], seed = seed)$cv_res
      rownames(HiDecon_res) <- rownames(Other_res[[1]])
      colnames(HiDecon_res) <- colnames(Other_res[[1]])
      All_res <- list("HiDecon" = HiDecon_res, "CIBERSORT" = Other_res$CIBERSORT, "dtangle" = Other_res$dtangle,
                      "MuSiC" = Other_res$MuSiC)
      result_list <- c(result_list, list(All_res))
    }
  }
  return(result_list)
}










