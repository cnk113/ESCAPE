# Run DeepImpute
DeepImpute.EnImpute = function(count){
  source('deepimp.py')
  impute <- deepimp(count)
}


# Run I-Impute
iImpute.EnImpute = function(count){
  source('iimpute.py')
  impute <- iimpute_norm(count)
}


# Run netNMF-sc
netNMF.EnImpute = function(count){

}


# Run scNPF
scNPF.EnImpute = function(count){

}


# BISCUIT
BISCUIT.EnImpute = function(count){

}


# Rescue
rescue.EnImpute = function(count){
  impute <- rescue::bootstrapImputation(expression_matrix = expression_dropout@assays$RNA@data,
             number_pcs = 30, cores = 2, use_mclapply = FALSE, snn_resolution = 1, proportion_genes = 0.6,
             bootstrap_samples = 100)
}


# Run bayNorm
bayNorm.EnImpute = function(count){

  impute <- bayNorm::bayNorm(count)
  impute
}


# Run AutoImpute



# Run scDoc - Considers Dropout
scDoc.EnImpute = function(count,){


}


# Run scRecover - Considers Dropout 
scRecover.EnImpute = function(count, kc){

  scRecover::scRecover(count, Kcluster=kc, outputdir='./outDir_scRecover/', parallel=T)
  impute <- read.table()
}


# Run SCRIBE
SCRIBE.EnImpute = function(count, batch.ind, bio.ind){

  impute <- fitSCRIBE(count, batch.ind, bio.ind)
  impute
} 


# Run Scrabble?

# Run ALRA
ALRA.EnImpute = function(count, k = 0, q = 10, mkl = F){

  count.t = t(count)
  # library and log normalization
  count_norm = normalize_data(count.t)
  # Impute using the function alra
  count.ALRA = alra(count_norm, k=k, q=q, use.mkl=mkl)[[2]]
  count.ALRA = t(count.ALRA)

  count.ALRA = exp(count.ALRA)- 1

  row.names(count.ALRA) = row.names(count)
  colnames(count.ALRA) = colnames(count)
  count.ALRA = as.matrix(count.ALRA)
  count.ALRA
}


# Run DCA
DCA.EnImpute= function(count, normtype = "zheng", type = "zinb-conddisp",
                       l2 = 0, l1 = 0, l2enc = 0, l1enc = 0, ridge = 0,
                       gradclip = 5, activation = "relu", hiddensize = "64,32,64",
                       hyper = FALSE, hypern = 1000){

  count = round(count)
  dir.create("./DCA_result")

  command = paste("dca ./tmp/count.csv ./DCA_result/result",
                 "--normtype", eval(normtype),
                 "--type", eval(type),
                 "--l2", eval(l2),
                 "--l1", eval(l1),
                 "--l2enc", eval(l2enc),
                 "--l1enc", eval(l1enc),
                 "--ridge", eval(ridge),
                 "--gradclip", eval(gradclip),
                 "--activation", eval(activation),
                 "--hiddensize", eval(hiddensize))
  if (hyper==TRUE){
    command =  paste(command, "--hyper", "--hypern", eval(hypern))
  }

  # Impute using DCA
  system(eval(command))

  count.DCA = utils::read.csv("./DCA_result/result/mean.tsv", sep="\t",header = TRUE, row.names = NULL)
  unlink("./DCA_result", recursive=TRUE)

  count.DCA = as.matrix(count.DCA)
  count.DCA.final <- count.DCA[,-1]
  rownames(count.DCA.final) <- count.DCA[,1]
  rm(count.DCA)
  gc()
  count.DCA.final
}


# Run DrImpute
DrImpute.EnImpute = function(count, ks = 10:15, dists = c("spearman", "pearson"), method = "mean",
                             cls = NULL){
  start.time <- Sys.time()
  # Preprocess gene expression matrix
  count = DrImpute::preprocessSC(count, min.expressed.gene = 0, min.expressed.cell = 0, max.expressed.ratio = 1,
                                 normalize.by.size.effect = TRUE)
  # log transformation
  logcount = log(count+1)
  # Impute using the function DrImpute
  count_DrImpute = DrImpute::DrImpute(logcount, ks = ks, dists = dists, method = method,
                                      cls = cls)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(paste0("DrImpute time taken: ", time.taken))
  count_DrImpute = exp(count_DrImpute)-1

  row.names(count_DrImpute) = row.names(count)
  colnames(count_DrImpute) = colnames(count)
  count_DrImpute = as.matrix(count_DrImpute)
  count_DrImpute
}


# Run MAGIC
MAGIC.EnImpute = function(count, k = 10, alpha = 15, t = "auto", npca = 20,
                          t.max = 20, knn.dist.method = "euclidean", n.jobs = 1){
  count.t = t(count)
  # Library size normalization
  count.normalized = Rmagic::library.size.normalize(count.t)
  # Log normalization
  count.log = log(count.normalized + 1)
  # Impute using the function magic
  count_MAGIC = Rmagic::magic(count.log, k = k, alpha = alpha, t = t, npca = npca,
                              t.max = t.max, knn.dist.method = knn.dist.method, n.jobs = n.jobs)

  count_MAGIC = t(as.matrix(count_MAGIC))
  count_MAGIC[count_MAGIC<0]=0
  count_MAGIC = exp(count_MAGIC)-1

  row.names(count_MAGIC) = row.names(count)
  colnames(count_MAGIC) = colnames(count)
  count_MAGIC = as.matrix(count_MAGIC)
  count_MAGIC
}


# Run SAVER
SAVER.EnImpute = function(count, do.fast = TRUE, ncores = 1, size.factor = NULL,
                          npred = NULL, null.model = FALSE, mu = NULL){
  start.time <- Sys.time()
  # Impute using the function saver
  count_SAVER= SAVER::saver(count, do.fast = do.fast, ncores = ncores, size.factor = size.factor,
                            npred = npred, null.model = null.model, mu=mu)$estimate
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(paste0("SAVER time taken: ", time.taken))
  row.names(count_SAVER) = row.names(count)
  colnames(count_SAVER) = colnames(count)
  count_SAVER = as.matrix(count_SAVER)
  count_SAVER
}


# Run scImpute
scImpute.EnImpute = function(count, drop_thre = 0.5, Kcluster = 10, labeled = FALSE,
                             labels = NULL, genelen = NULL, ncores = 1){

  dir.create("./scImpute_result")
  saveRDS(count, file ="./scImpute_result/count.rds")
  # Run scImpute using the function scimpute
  start.time <- Sys.time()
  out = scImpute::scimpute(count_path = "./scImpute_result/count.rds",
                           infile = "rds",
                           outfile = "rds",
                           type = "count",
                           out_dir = "./scImpute_result/",
                           labeled = labeled,
                           drop_thre = drop_thre,
                           Kcluster = Kcluster,
                           genelen = genelen,
                           ncores = ncores)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(paste0("scImpute time taken: ", time.taken))
  count_scImpute = readRDS("./scImpute_result/scimpute_count.rds")
  unlink("./scImpute_result", recursive=TRUE)

  row.names(count_scImpute) = row.names(count)
  colnames(count_scImpute) = colnames(count)

  count_scImpute = as.matrix(count_scImpute)
  count_scImpute

}


# Run scRMD
scRMD.EnImpute = function(count, tau = NULL, lambda = NULL, candidate = 0.05){

  # library and log normalization
  totalUMIPerCell = colSums(count)
  count_norm = log10(sweep(count, 2, totalUMIPerCell/1000000, '/')+1);

  count.t = t(count_norm)
  cutoff = quantile(count.t[count.t>0], candidate)
  # Impute using the function rmd
  start.time <- Sys.time()
  count.scRMD  = scRMD::rmd(count.t, candidate = cutoff)$exprs
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(paste0("scRMD time taken: ", time.taken))
  count.scRMD = t(count.scRMD)

  count.scRMD = 10^count.scRMD-1

  row.names(count.scRMD) = row.names(count.scRMD)
  colnames(count.scRMD) = colnames(count.scRMD)
  count.scRMD = as.matrix(count.scRMD)
  count.scRMD
}



#' Run Empute on a raw read count matrix
#' @param count raw read count matrix. The rows correspond to genes and the columns correspond to cells.
#' @param scale.factor scale factor used to re-scale the imputed results generated by
#' different individual methods. Default is 10000.
#' @param trim  specifies the fraction (between 0 and 0.5)  of observations to be trimmed
#' from each end before the mean is computed. Default is 0.3.
#' @param ALRA  a boolean variable that defines whether to impute the raw data using the ALRA method.
#' Default is TRUE.
#' @param DCA a boolean variable that defines whether to impute the raw data using the DCA method.
#' Default is TRUE.
#' @param DrImpute a boolean variable that defines whether to impute the raw data using the DrImpute method.
#' Default is "TRUE".
#' @param MAGIC a boolean variable that defines whether to impute the raw data using the MAGIC method.
#' Default is TRUE.
#' @param SAVER  a boolean variable that defines whether to impute the raw data using the SAVER method.
#' Default is TRUE.
#' @param scImpute a boolean variable that defines whether to impute the raw data using the scImpute method.
#' Default is TRUE.
#' @param scRMD a boolean variable that defines whether to impute the raw data using the scRMD method.
#' Default is TRUE.
#' @param ALRA.k the rank of the rank-k approximation in ALRA. Set to 0 for automated choice of k.
#' Default is 0.
#' @param ALRA.q the number of power iterations in randomized SVD used by ALRA. Default is 10.
#' @param DCA.normtype a string variable specifying the type of size factor estimation in DCA.
#' Possible values: "deseq", "zheng". Default is "zheng".
#' @param DCA.type a string variable specifying type of autoencoder in DCA. Possible values:
#' "normal", "poisson", "nb", "nb-shared", "nb-conddisp", "nb-fork", "zinb", "zinb-shared", "zinb-conddisp",
#' "zinb-fork". Default is "zinb-conddisp".
#' @param DCA.l2 a real number specifying the L2 regularization coefficient in DCA.  Default is 0.
#' @param DCA.l1 a real number specifying the L1 regularization coefficient in DCA.  Default is 0.
#' @param DCA.l2enc a real number specifying the encoder-specific L2 regularization coefficient in DCA.
#' Default is 0.
#' @param DCA.l1enc a real number specifying the encoder-specific L1 regularization coefficient in DCA.
#' Default is 0.
#' @param DCA.ridge a real number specifying the L2 regularization coefficient for dropout probabilities
#' in DCA. Default is 0.
#' @param DCA.gradclip a real number specifying the Clip grad values in DCA. Default is 5.
#' @param DCA.activation a string value specifying the activation function of hidden unit in DCA. Default is "relu".
#' @param DCA.hiddensize a string vector specifying the size of hidden layers in DCA. Default is "64,32,64".
#' @param DCA.hyper a logical value specifying whether hyperparameter search is performed in DCA.
#' @param DCA.hypern an integer specifying the number of samples drawn from hyperparameter distributions
#' during optimization in DCA. Default is 1000.
#' @param DrImpute.ks an integer vector specifying the number of cell clustering groups in DrImpute.
#' Default is 10:15.
#' @param DrImpute.dists a string vector specifying the distance metrics in DrImpute. Default is
#' c("spearman", "pearson").
#' @param DrImpute.method a string specifying the method used for imputation in DrImpute. Use "mean"
#' for mean imputation, "med" for median imputation.
#' @param DrImpute.cls a matrix specifying the clustering information manually provided by users in DrImpute.
#' The rows represent different clusterings, and the columns represent cells. Default is NULL,
#' which means the user do not provide the clustering information.
#' @param MAGIC.k an integer specifying the number of nearest neighbors on which to build kernel in MAGIC.
#' Default is 10.
#' @param MAGIC.alpha an integer specifying the decay rate of kernel tails in MAGIC. Default is 15.
#' @param MAGIC.t an integer specifying the diffusion time for the Markov Affinity Matrix in MAGIC.
#' Default is "auto". For detail about the approach to set paramter t automatically,
#' please refer to the reference.
#' @param MAGIC.npca  an integer specifying the number of PCA components in MAGIC.
#' Default is 20.
#' @param MAGIC.t.max an integer specifying the maximum value of t to test for automatic t selection in MAGIC.
#' Default is 20.
#' @param MAGIC.knn.dist.method  a string value specifying the metric for building kNN graph in MAGIC.
#' Recommended values: "euclidean", "cosine". Default is "euclidean".
#' @param MAGIC.n.jobs an integer specifying the number of jobs used for computation in MAGIC. If -1 all CPUs are used.
#' If 1 is given, no parallel computing code is used at all. For n.jobs below -1, (n.cpus + 1 + n.jobs)
#' are used. Thus for n.jobs = -2, all CPUs but one are used.
#' @param SAVER.do.fast a boolean variable specifying whether the prediction step is
#' approximated in SAVER. Default is TRUE.
#' @param SAVER.ncores number of cores to use in SAVER. Default is 1.
#' @param SAVER.size.factor a vector of cell size specifying the normalization factors in SAVER.
#' If the data is already normalized or normalization is not desired, set size.factor = 1.
#' Default uses mean library size normalization.
#' @param SAVER.npred number of genes for regression prediction in SAVER. Selects the top npred genes in
#' terms of mean expression for regression prediction. Default is all genes.
#' @param SAVER.null.model a boolean variable specifying whether to use mean gene expression as prediction
#' in SAVER. Default is FALSE
#' @param SAVER.mu matrix of prior means in SAVER.
#' @param scImpute.drop_thre  a number (between 0 and 1) specifying the threshold on dropout probability in scImpute.
#' Default is 0.5.
#' @param scImpute.Kcluster an integer specifying the number of cell subpopulations in scImpute. Default is 10.
#' @param scImpute.labeled  a boolean variable indicating whether cell type information is given in scImpute. Default is FALSE.
#' @param scImpute.labels  a character vector specifying the cell type in scImpute. Only needed when \code{labeled = TRUE}.
#' Default is NULL
#' @param scImpute.genelen an integer vector giving the length of each gene in scImpute.  Default is NULL.

#' @param scImpute.ncores an integer specifying the number of cores used for parallel computation in scImpute. Default is 1.
#' @param scRMD.tau a non-negative real number specifying the tuning parameter to penalize the sparse term. Default is NULL.
#' @param scRMD.lambda a non-negative real number specifying the tuning parameter to penalize the row rank term. Default is NULL.
#' @param scRMD.candidate a real number (0 to 1) specifying the cutoff for candidate drop out. Default is 0.05.

#'
#' @return a list with the following components
#' \item{\code{count.EnImpute.log}}{Imputed count matrix generated by EnImpute (log scale).}
#' \item{\code{count.EnImpute.exp}}{Imputed count matrix generated by EnImpute (exp scale).}
#' \item{\code{count.imputed.individual.exp}}{Imputed count matrices generated by different individual imputation methods (exp scale).}
#' \item{\code{Methods.used}}{The individual methods used by EnImpute.}
#'
#' @export
#' @import DrImpute Rmagic SAVER scImpute scRMD Seurat rsvd
#'
#' @author Xiao-Fei Zhang  <zhangxf@mail.ccnu.edu.cn>
#'
EnImpute = function(count, scale.factor = 10000, trim = 0.3,
                  ALRA = TRUE, DCA = TRUE, DrImpute = TRUE, MAGIC = TRUE, SAVER = TRUE,
                  scImpute = TRUE, scRMD = TRUE, ALRA.k = 0, ALRA.q = 10, ALRA.mkl = TRUE,
                  DCA.normtype = "zheng", DCA.type = "zinb-conddisp",
                  DCA.l2 = 0, DCA.l1 =0, DCA.l2enc = 0, DCA.l1enc = 0, DCA.ridge = 0,
                  DCA.gradclip = 5, DCA.activation = "relu", DCA.hiddensize = "64,32,64",
                  DCA.hyper = FALSE, DCA.hypern = 1000,
                  DrImpute.ks = 10:15, DrImpute.dists = c("spearman", "pearson"),
                  DrImpute.method = "mean", DrImpute.cls = NULL,
                  MAGIC.k = 10, MAGIC.alpha = 15, MAGIC.t = "auto", MAGIC.npca = 30,
                  MAGIC.t.max = 20, MAGIC.knn.dist.method = "euclidean", MAGIC.n.jobs = 1,
                  SAVER.do.fast = TRUE, SAVER.ncores = 2, SAVER.size.factor = NULL,
                  SAVER.npred = NULL, SAVER.null.model = FALSE, SAVER.mu = NULL,
                  scImpute.drop_thre = 0.5, scImpute.Kcluster = 5, scImpute.labeled = FALSE,
                  scImpute.labels = NULL, scImpute.genelen = NULL, scImpute.ncores = 1,
                  scRMD.tau = NULL, scRMD.lambda = NULL, scRMD.candidate = 0.05){
  
  Methods = c("ALRA", "DCA", "DrImpute", "MAGIC", "SAVER", "scImpute","scRMD")
  Methods.idx = c(ALRA, DCA,  DrImpute, MAGIC, SAVER, scImpute, scRMD)
  
  if(sum(Methods.idx)==0)
    stop("You need choose at least one individual imputation method.")
  Methods.used = Methods[Methods.idx]

  K = length(Methods.used)

  p = dim(count)[1]
  n = dim(count)[2]

  count.imputed.individual = array(0, dim=c(p,n,K))
  dimnames(count.imputed.individual)[[1]] = rownames(count)
  dimnames(count.imputed.individual)[[2]] = colnames(count)
  dimnames(count.imputed.individual)[[3]]= Methods.used
  
  dir.create("./tmp")
  utils::write.csv(count,"./tmp/count.csv")
  options(future.globals.maxSize = 100000 * 1024^2)
  future::plan("multicore")
  imputed.list <- future.apply::future_lapply(Methods.used, function(i){
    if(i == "ALRA"){
      ALRA.EnImpute(count, k = ALRA.k, q = ALRA.q, mkl = ALRA.mkl)
    }
    else if(i == "DCA"){
      DCA.EnImpute(count, normtype = DCA.normtype, type = DCA.type,
                                                 l2 = DCA.l2, l1 = DCA.l1, l2enc = DCA.l2enc, l1enc = DCA.l1enc, ridge = DCA.ridge,
                                                 gradclip = DCA.gradclip, activation = DCA.activation, hiddensize = DCA.hiddensize,
                                                 hyper = DCA.hyper, hypern = DCA.hypern)
    }
    else if(i == "DrImpute"){
      DrImpute.EnImpute(count, ks = DrImpute.ks, dists = DrImpute.dists, method = DrImpute.method, cls = DrImpute.cls)
    }
    else if(i == "MAGIC"){
      MAGIC.EnImpute(count, k = MAGIC.k, alpha = MAGIC.alpha, t = MAGIC.t, npca = MAGIC.npca, t.max = MAGIC.t.max,
                                                   knn.dist.method = MAGIC.knn.dist.method, n.jobs = MAGIC.n.jobs)
    }
    else if(i == "SAVER"){
      SAVER.EnImpute(count, do.fast = SAVER.do.fast, ncores = SAVER.ncores, size.factor = SAVER.size.factor, npred = SAVER.npred,
                                                   null.model = SAVER.null.model, mu = SAVER.mu)
    }
    else if(i == "scImpute"){
      scImpute.EnImpute(count, drop_thre = scImpute.drop_thre, Kcluster = scImpute.Kcluster,
                                                      labeled = scImpute.labeled, labels = scImpute.labels,
                                                      genelen = scImpute.genelen, ncores = scImpute.ncores)
    }
    else if(i == "scRMD"){
      scRMD.EnImpute(count, tau = scRMD.tau, lambda = scRMD.lambda, candidate = scRMD.candidate)
    }
    else if(i == "iImpute"){
      iimpute.EnImpute("./tmp/count.csv")
    }
    else if(i == "DeepImpute"){
    }
    else if(i == "IImpute"){
    }
    else if(i == ""){
    }
  })
  print("===========================Starting Consensus=========================")
  count.imputed.individual <- abind::abind(imputed.list, along = 3)
  mode(count.imputed.individual) <- 'numeric'
#  count.imputed.individual =  floor(count.imputed.individual)
  count.imputed.individual[count.imputed.individual<=0] = 0
  count.imputed.individual[is.na(count.imputed.individual)] = 0

  # Rescale the imputed count matrices
  count.imputed.individual.rescaled = count.imputed.individual
  for (k in 1:K){
    totalUMIPerCell = colSums(count.imputed.individual[,,k])
    count.imputed.individual.rescaled[,,k] = sweep(count.imputed.individual[,,k] , 2, totalUMIPerCell/scale.factor, '/');
  }
  count.imputed.individual.rescaled = log(count.imputed.individual.rescaled+1)

  count.EnImpute.log = apply(count.imputed.individual.rescaled, 1:2, mean,  trim = trim)
  rownames(count.EnImpute.log) = rownames(count)
  colnames(count.EnImpute.log) = colnames(count)
  count.EnImpute.exp = exp(count.EnImpute.log) - 1


  result = list(count.EnImpute.log = count.EnImpute.log, count.EnImpute.exp = count.EnImpute.exp,
                count.imputed.individual = count.imputed.individual, count.imputed.individual.rescaled = count.imputed.individual.rescaled,
                Methods.used = Methods.used)
  result
}
