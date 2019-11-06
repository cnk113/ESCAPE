#' Similarity calculation via Weigheted Cosine similarity
#' 
#' @description Calculate weighted cosine similarity between cells 
#'              based on both input data and dropout probability
#'
#' @param input log2 transformation of raw count matrix, row is gene and col is sample
#' @param dmat dropout probability matrix, each entry is the probability of 
#'             observation for gene i in cell j is a dropout event
#' @param dprob threshold to determine whether an observation is a dropout 
#'              event. If its dropout probability > dprob, it's a statistical 
#'              dropout event and will be given weight. The default is 0.5.
#' 
#' @return a matrix containing cell-to-cell weighted consine similarity
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
#' @export
#' 


sim.calc <- function(input, dmat, dprob=0.5) {
  
  if (!all.equal(sort(rownames(input)), sort(rownames(dmat)))) {
    stop("Row names are not idential between input and dmat")
  }
  
  if (!all.equal(sort(colnames(input)), sort(colnames(dmat)))) {
    stop("Column names are not idential between input and dmat")
  }
  
  dmat <- dmat[rownames(input), colnames(input)]
  
  cat("Calculating of weighted Cosine similarity...\n")
  
  nc <- dim(input)[2]
  ng <- dim(input)[1]
  
  w.prior <- rowSums(dmat > dprob) / nc

  sim.o <- matrix(0, nrow=nc, ncol=nc)
  rownames(sim.o) <- colnames(sim.o) <- colnames(input)
  
  for (i in 1:nc) {
    for (j in 1:nc) {
      
      if (i < j) {
        
        prob0 <- rep(1, ng)
        prob0[(dmat[,i]>dprob) | (dmat[,j]>dprob)] <- -1 
        prob0[prob0 == -1] <- w.prior[prob0 == -1]
        
        A <- input[,i] * sqrt(prob0)
        B <- input[,j] * sqrt(prob0)
        
        sim.o[i, j] <- sim.o[j, i] <- sum(A*B) / sqrt(sum(A^2)*sum(B^2))

      }
    }
  }
  
  out <- sim.o / max(sim.o)
  
  return(out)
}
