#' KNN imputation
#' 
#' @description Impute dropout events via KNN double weighted imputation
#' 
#' @param input raw count matrix, row is gene and col is cell
#' @param dmat matrix containing dropout probability
#' @param sim.mat cell-to-cell similarity matrix
#' @param k number of nearest neighbours used for imputation
#' @param dprob threshold to determine whether an observation 
#'              is a dropout event. If observation >= dprob, 
#'              it will be considered as a dropout. The default 
#'              is 0.5.
#' @param sim.cut cut-point used to determine what is the cutoff 
#'                for a trustable similarity. The cutoff is the 
#'                "sim.cut"*100th quantile of all cell-to-cell 
#'                similarities. Any cell having similarity greater 
#'                than cutoff to a target cell will be considered 
#'                to use in imputation.
#' @param seed seed to generate reproducible results
#' 
#' @return a list containing imputation results:
#' \itemize{
#'   \item output.imp matrix containing imputed values
#'   \item num.imp number of observations imputed
#' }
#' 
#' @examples
#' data(zebrafish)
#' 
#' offsets <- as.numeric(log(colSums(zebrafish)))
#' count <- zebrafish[rowSums(zebrafish > 5) > 4, ]
#' 
#' # matrix of dropout probability
#' dp.mat <- prob.dropout(input = count, is.count = T, offsets = offsets, mcore = 3)
#' 
#' # similarity matrix
#' sim.mat <- sim.calc(log2(count+1), dp.mat)
#' 
#' # KNN imputation
#' impute.mat <- impute.knn(input = count, dmat = dp.mat, sim.mat = sim.mat, 
#'                          k = 6, sim.cut = 1e-4)
#' 
#' @export
#' 


impute.knn <- function(input, dmat, sim.mat, k=5, dprob=0.5, sim.cut=1e-4, seed=2018) {
  
  if (!all.equal(sort(rownames(input)), sort(rownames(dmat)))) {
    stop("Row names are not idential between input and dmat")
  }
  
  if (!all.equal(sort(colnames(input)), sort(colnames(dmat)))) {
    stop("Column names are not idential between input and dmat")
  }
  
  dmat <- dmat[rownames(input), colnames(input)]
  
  if (!all.equal(sort(colnames(input)), sort(colnames(dmat)))) {
    stop("Samples names are not idential between input and dmat")
  }
  
  sim.mat <- sim.mat[colnames(input), colnames(input)]
  
  cat("Imputing drop-outs using KNN...\n")
  
  # replace zero with double weighted average
  nc <- ncol(input)
  
  # to avoid impute value to rare cell type
  # cell-to-cell similarity must be greater than certain citeria
  # for instance, greater than median cell-to-cell similarity 
  # greater than 75% quantile (sim.cut=0.75)
  cutoff <- quantile(sim.mat[sim.mat!=0][lower.tri(sim.mat, diag = FALSE)], sim.cut) 
  cutoff.E <- quantile(sim.mat[sim.mat!=0][lower.tri(sim.mat, diag = FALSE)], 0.9) 
  
  # find k nearest neighbour for each cell
  nn <- matrix(NA, nrow=nc, ncol=k)
  outlier <- NULL
  for (i in 1: nc) {
    nn[i, ] <- order(sim.mat[i, -i], decreasing = T)[1:k]
    if (!(any(sim.mat[i,nn[i, ]] > cutoff))) {
      outlier <- c(outlier, i)
    }
  }
  
  dat.imp <- cbind(input, dmat)
  
  # imputation
  set.seed(seed)
  dat.imp.o <- t(apply(dat.imp, 1, gimpute, 
                       sim.mat=sim.mat,
                       nc=nc,
                       nn=nn,
                       dprob=dprob,
                       cutoff=cutoff,
                       cutoff.E=cutoff.E))

  count.imp <- dat.imp.o[,1:nc]
  dmat.imp <-  dat.imp.o[,-c(1:nc)]
  
  colnames(count.imp) <- colnames(dmat.imp) <- colnames(input)
  rownames(count.imp) <- rownames(dmat.imp) <- rownames(input)
  
  num.imp <- sum(dmat != dmat.imp) / prod(dim(dmat))
  
  return(list(output.imp = count.imp,
              num.imp = num.imp,
              outlier=outlier))
  
}



# impute function
gimpute <- function(x, sim.mat, nc, nn, dprob, cutoff, cutoff.E) {
  
  x <- as.numeric(x)
  count.tmp <- x[1:nc]
  dmat.tmp  <- x[-(1:nc)]
  
  rnc <- sample(1:nc)
  
  if (!anyNA(x)) {
    for (j in rnc) {
      
      if (dmat.tmp[j] > dprob) {
        
        prob  <- dmat.tmp[nn[j,]]
        value <- count.tmp[nn[j,]]
        dd <- sim.mat[j, nn[j,]]
        
        if (any(dd>cutoff)) {
          
          value <- value[dd>cutoff]
          prob <- prob[dd>cutoff]
          dd <- dd[dd>cutoff]
          
          if (sum(prob < dprob) >= 2) {
            
            # count.out[j] <- sum(value*prob*dd) / sum(prob*dd)
            count.tmp[j] <- sum(value*dd) / sum(dd)
            
            dmat.tmp[j] <- 1e-4
            
          } else {
            
            if (sum(prob < dprob) == 1) {
              
              if (dd[prob < dprob] > cutoff.E) {
                
                count.tmp[j] <- dd[prob < dprob]
                
                dmat.tmp[j] <- 1e-4
                
              } else {
                
                count.tmp[j] <- max(1, median(value*dd) / median(dd))
                
                dmat.tmp[j] <- dprob
                
              }
            }
          }
        }
      }
    }
  }
  
  return(c(ceiling(count.tmp), dmat.tmp))
}
