#' Estimate drop-out probability of scRNA-seq expression data
#' 
#' @description Estimate drop-out probability of scRNA-seq expression 
#'              data (raw count or TPM values) via mixture model. 
#'              Poisson-NB model for raw count; Gamma-Normal model for 
#'              TPM (or other normalized expression values). 
#'              
#' @param input matrix containing scRNA-seq expression values,
#'              raw count or TPM values. Row is gene, col is sample.
#' @param is.count logical; if TRUE the input is raw count,
#'                 otherwise, it is TPM.
#' @param offsets offsets required for Poisson-NB model, the defaul is 
#'                0 for all samples. To account for the influence of 
#'                library size on the estimates of drop-out probability, 
#'                user could use log transformed library size (total counts) 
#'                as offsets.
#' @param mcore number of cores used in calculations. The parallel 
#'              computation depends on package foreach and doParallel.
#' @param ... further arguments passed to param.est and drop.prob functions.
#' 
#' @return a matrix as the same dimension of the input and 
#'         containing drop-out probability for each gene at each sample
#'         
#' @examples
#' # raw count example
#' data(zebrafish)
#' 
#' offsets <- as.numeric(log(colSums(zebrafish)))
#' count <- zebrafish[rowSums(zebrafish > 5) > 4, ]
#' 
#' dp.mat <- prob.dropout(input = count, is.count = T, offsets = offsets, mcore = 3)
#' 
#' 
#' # TPM example
#' data(lung)
#' 
#' tpm <- lung[rowSums(lung>5) >4, ]
#' 
#' dp.mat <- prob.dropout(input = tpm, is.count = F, mcore = 3)
#' 
#' @export
#' 


prob.dropout <- function(input, is.count = T, offsets = rep(0, ncol(input)), mcore = 1, ...) {
  
  tic <- proc.time()
  
  if (is.count) {
    
    cat(paste0("Input is raw count which has ", nrow(input), " genes and ", ncol(input), " samples.\n"))
    
    cat("   Start calculation of drop-out probability via Poisson-NB mixture model.\n")
    
    if (any(offsets > 0)) {
      cat("   log transformed library size is provided.\n")
    } else {
      cat("   log transformed library size is not provided and default value 0 will be used for all samples.\n")
    }
    
    if (mcore > 1) {
      
      require(foreach)
      require(doParallel)
      
      # setup parallel backend to use many processors
      cores <- detectCores()
      
      mcores <- min(mcore, (cores[1]-1))
      
      # not to overload your computer
      cl <- makeCluster(mcores) 
      
      cat(paste0("   Multiple cores were requested and ", mcores, " cores are using.\n"))
      
      # register the parallel backend
      registerDoParallel(cl)
      
      tempMatrix <- foreach(i=1:nrow(input), .combine=rbind, 
                            .packages = c("scDoc", "MASS")) %dopar% {
        
        xtmp <- as.numeric(input[i, ])
        
        rs <- try(param.est.pnb(y=xtmp, offsets = offsets, ...))
        
        if(!inherits(rs, "try-error")) {
          finalMatrix <- drop.prob.pnb(pars=c(rs$gamma, rs$alpha, rs$delta), 
                                  s=rs$s, y = xtmp, offset = offsets, ...)
        }
        
        finalMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
      }
      
      #stop cluster
      stopCluster(cl)
      
    } else {
      
      tempMatrix <- matrix(NA, nrow=nrow(input), ncol=ncol(input))
      
      for (i in 1:nrow(tempMatrix)) {
        xtmp <- as.numeric(input[i, ])
        
        rs <- try(param.est.pnb(y=xtmp, offsets = offsets, ...))
        
        if(!inherits(rs, "try-error")) {
          tempMatrix[i,] <- drop.prob.pnb(pars=c(rs$gamma, rs$alpha, rs$delta), 
                                      s=rs$s, y = xtmp, offset = offsets, ...)
        } else {
          print("Error!!!!!!")
        }
        
      }
      
    }
    
  } else {
    
    cat(paste0("Input is normalized expression value which has ", nrow(input), " genes and ", ncol(input), " samples.\n"))
    
    cat("   Start calculation of drop-out probability via Gamma-Normal mixture model.\n")
    
    xdata <- log10(input + 1.01)
    point <- log10(1.01)
    
    if (mcore > 1) {
      
      require(foreach)
      require(doParallel)
      
      # setup parallel backend to use many processors
      cores <- detectCores()
      
      mcores <- min(mcore, (cores[1]-1))
      
      # not to overload your computer
      cl <- makeCluster(mcores) 
      
      cat(paste0("   Multiple cores were requested and ", mcores, " cores are using.\n"))
      
      # register the parallel backend
      registerDoParallel(cl)
      
      tempMatrix <- foreach(i=1:nrow(xdata), .combine=rbind, 
                            .packages = c("scDoc", "MASS")) %dopar% {
        
        xtmp <- as.numeric(xdata[i, ])
        
        rs <- try(get_mix(xdata=xtmp, point=point))
        
        if(!inherits(rs, "try-error")) {
          finalMatrix <- dmix(x=xtmp, pars=rs)
        }
        
        finalMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
      }
      
      #stop cluster
      stopCluster(cl)
      
    } else {
      
      tempMatrix <- matrix(NA, nrow=nrow(xdata), ncol=ncol(xdata))
      
      for (i in 1:nrow(tempMatrix)) {
        
        xtmp <- as.numeric(xdata[i, ])
        
        rs <- try(get_mix(xdata=xtmp, point=point))
        
        if(!inherits(rs, "try-error")) {
          tempMatrix[i,] <- dmix(x=xtmp, pars=rs)
        } else {
          print("Error!!!!!!")
        }
        
      }
      
    }
    
  }
  
  
  rownames(tempMatrix) <- rownames(input)
  colnames(tempMatrix) <- colnames(input)
  
  toc <- proc.time() - tic
  
  cat("   End calculation and", toc[3], "seconds elapsed.\n")
  
  return(tempMatrix)
  
}
