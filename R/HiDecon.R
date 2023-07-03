
#' @title reference.
#' @description Calculate signature matrix (reference) for some layer.
#'
#' @param ref matrix: Single cell data matrix (cells * genes) after count per million normalization and log2 fold change.
#' @param ref.type vector: cell type label of reference.
#' @param B list: contains cell type mapping matrices (same length as tree depth).
#' @param type_order vector: specify cell type order in the bottom layer of tree.
#' @param layer integer: the layer needed to be calculated.
#'
#' @return matrix: signature matrix for the specified layer
#' @export
#'
reference <- function(ref,ref.type,B,type_order,layer){
  n.level <- length(B)
  if(layer == n.level){
    type.relationship <- diag(1,nrow=length(type_order))
  }
  if(layer < n.level){
    type.relationship <- diag(1,nrow=nrow(B[[layer+1]]))
    for(i in (layer+1):n.level){
      type.relationship <- type.relationship %*% B[[i]]
    }
  }
  n.type <- nrow(type.relationship)
  reference <- matrix(nrow = ncol(ref),ncol = n.type)
  for(i in 1:n.type){
    types <- type.relationship[i,]
    types <- type_order[which(types == 1)]
    in.cluster <- ref.type %in% types
    reference[,i] <- apply(ref[in.cluster,, drop=F] ,MARGIN = 2, mean)
  }
  return(reference)
}


#' @title GeneratesigReference.
#' @description General signature matrix calculator.
#' @param ref matrix (genes * cells): single cell reference.
#' @param cell.types vector: cell type label of reference.
#'
#' @return matrix: signature matrix
#' @export
#' @importFrom matrixStats rowMeans2
GeneratesigReference <- function(ref,cell.types) {
  cell.labels <- factor(cell.types)
  all.cell.types <- levels(cell.labels)
  aggr.fn <- function(cell.type) {
    rowMeans2(ref[,cell.labels == cell.type, drop=F])
  }
  template <- numeric(nrow(ref))
  sig <- vapply(all.cell.types, aggr.fn, template)
  rownames(sig) = rownames(ref)
  return(sig)
}



#' @title Create.A.
#' @description Create the A.tilde in the HiDecon paper.
#'
#' @param Signature.list a length L list of #markers_l * #cell_types signature matrices for the L layers.
#' @param Cell.size.list a length L list of cell size vectors. Each vector is of length K.l=#cell_types in layer l.
#' @param standardize if the loss should be standardized by the number of marker genes. Default is TRUE.
#'
#' @return A list to be used in function Est.AllPi.
#' @export
#'
Create.A <- function(Signature.list, Cell.size.list=NULL, standardize=T) {
  P <- sapply(Signature.list,nrow)
  K <- sapply(Signature.list,ncol)
  L <- length(Signature.list)
  A <- matrix(0,nrow=sum(P),ncol=sum(K))
  for (l in 1:L) {
    ind.cols <- ifelse(l==1,1,1+sum(K[1:(l-1)])):sum(K[1:l])
    ind.rows <- ifelse(l==1,1,1+sum(P[1:(l-1)])):sum(P[1:l])
    A[ind.rows,ind.cols] <- Signature.list[[l]]
    if (!is.null(Cell.size.list)) {
      A[ind.rows,ind.cols] <- A[ind.rows,ind.cols]%*%diag(Cell.size.list[[l]],nrow=length(Cell.size.list[[l]]))
    }
    if (standardize) {
      A[ind.rows,ind.cols] <- 1/sqrt(P[l])*A[ind.rows,ind.cols]
    }
  }
  return(list(A=A,K.L=K[L],include.cell.size=!is.null(Cell.size.list),standardized=standardize))
}




#' @title Create.B.
#' @description Create the B.tilde in the HiDecon paper.
#'
#' @param B.list list: contains cell type mapping matrices (from B1,2).
#' @param standardize if the loss should be standardized by the number of cell types. Default is true.
#'
#' @return A matrix to be used in function Est.AllPi.
#' @export
#'
Create.B <- function(B.list, standardize=T) {
  K <- c(nrow(B.list[[1]]),sapply(B.list,ncol))
  L <- length(K)
  if (standardize) {
    B.ident <- diag(rep(1/sqrt(K[1:(L-1)]),times=K[1:(L-1)]))
  } else {
    B.ident <- diag(1,sum(K[1:(L-1)]))
  }
  B.ident <- cbind(B.ident, matrix(0,nrow=sum(K[1:(L-1)]),ncol=K[L]))
  B.pen <- matrix(0,nrow=nrow(B.ident),ncol=ncol(B.ident))
  for (l in 1:(L-1)) {
    ind.rows <- (1 + ifelse(l==1,0,sum(K[1:(l-1)]))):sum(K[1:l])
    ind.cols <- (1 + sum(K[1:l])):sum(K[1:(l+1)])
    B.pen[ind.rows,ind.cols] <- B.list[[l]]
    if (standardize) {
      B.pen[ind.rows,ind.cols] <- 1/sqrt(K[l])*B.pen[ind.rows,ind.cols]
    }
  }
  return(B.ident - B.pen)
}



#' @title Est.AllPi
#' @description Estimate cellular fractions for all samples.
#'
#' @param Y.list a length L list of #individuals x #markers bulk data matrices. The dimension of the lth element of this list is #individuals x #markers_l.
#' @param A the output from Create.A.
#' @param B the output from Create.B.
#' @param lambda numerical, the penalty parameter lambda in the HiDecon paper.
#' @param Pi.start vector, non-negative starting point for the estimation of sample i. Default is NULL, then function will assign the starting point.
#' @param max.iter the maximum number of iterations, default is 1e4.
#' @param tol tolerance for the KKT conditions, default is 1e-6.
#'
#' @return a list with elements:
#'    * Pi: a #individuals x K.L matrix of cellular fractions
#'    * n.iter: number of iterations for each individual
#'    * out: whether the optimization converged. 0:yes, 1:no
#' @export
#' @importFrom rlist list.cbind
Est.AllPi <- function(Y.list, A, B, lambda, Pi.start=NULL, max.iter=1e4, tol=1e-6) {
  K.L <- A$K.L
  if (A$standardized) {Y.list <- lapply(1:length(Y.list),function(l){1/sqrt(nrow(Y.list[[l]]))*Y.list[[l]]})}
  A <- A$A
  Y <- list.cbind(Y.list); rm(Y.list)
  H <- t(A)%*%A + lambda*t(B)%*%B
  H.inv <- solve(H)
  n <- nrow(Y)
  out <- list()
  out$Pi <- matrix(NA,nrow=n,ncol=K.L)
  out$n.iter <- rep(NA,n)
  out$out <- rep(NA,n)
  YA <- Y%*%A; rm(Y,A)
  for (i in 1:n) {
    pi.i <- EstPi(H = H, H.inv = H.inv, Aty = YA[i,], K.L = K.L, pi.start = NULL, max.iter = max.iter, tol = tol)
    out$n.iter[i] <- pi.i$n.iter
    out$out[i] <- pi.i$out
    out$Pi[i,] <- pi.i$pi[(length(pi.i$pi)-K.L+1):length(pi.i$pi)]
    out$Pi[i,] <- out$Pi[i,]/sum(out$Pi[i,])
  }
  return(out)
}



#' @title EstPi.
#' @description Estimate cellular fractions for one sample.
#' @param H matrix, Hessian matrix in the HiDecon paper.
#' @param H.inv inverse of H.
#' @param Aty the ith row of matrix Y*A.
#' @param K.L number of cell types in layer L.
#' @param pi.start vector, non-negative starting point for the estimation of sample i. Default is NULL, then function will assign the starting point.
#' @param max.iter the maximum number of iterations, default is 1e4.
#' @param tol tolerance for the KKT conditions, default is 1e-6.
#'
#' @return a list with elements:
#'    * Pi: cellular fractions for one sample
#'    * n.iter: number of iterations to convergence
#'    * out: whether the optimization converged. 0:yes, 1:no
#' @export
EstPi <- function(H, H.inv, Aty, K.L, pi.start=NULL, max.iter=1e4, tol=1e-6) {
  test <- c(H.inv%*%Aty)  #The un-constrained solution
  if (all(test>=0)) {return(list(pi=test,n.iter=0,out=0))}
  if (is.null(pi.start)) {
    pi.start <- pmax(test,0.01)
    pi.start[(length(pi.start)-K.L+1):length(pi.start)] <- pi.start[(length(pi.start)-K.L+1):length(pi.start)]/sum(pi.start[(length(pi.start)-K.L+1):length(pi.start)])
  }
  rm(test)
  pi <- pi.start; rm(pi.start)
  K <- length(pi)
  for (i in 1:max.iter) {
    for (k in 1:K) {
      pi.minusk <- pi; pi.minusk[k] <- 0
      pi[k] <- max(0, (Aty[k] - sum(H[k,]*pi.minusk))/H[k,k])
    }
    if (Check.KKT(grad = c(H%*%pi - Aty), pi = pi, tol = tol)) {return(list(pi=pi,n.iter=i,out=0))}
  }
  return(list(pi=pi,n.iter=i,out=1))
}



#' @title Check.KKT.
#' @description Check the KKT conditions as described in the HiDecon paper.
#'
#' @param grad gradient of the objective function.
#' @param pi the update of cellular fraction estimates.
#' @param tol tolerance for the KKT conditions, default is 1e-6.
#'
#' @return logic.If KKT conditions satisfied, TRUE. Otherwise, FALSE.
#' @export
Check.KKT <- function(grad, pi, tol=1e-6) {
  out <- (abs(grad)<=tol) | (pi==0 & grad>=-tol)
  return(all(out))
}




#' @title HiDecon_marker
#' @description Select marker genes for each layer of the hierarchical tree.
#'
#' @param ref matrix: Single cell data matrix (cells * genes) after count per million normalization.
#' @param cell_type vector: cell type label of single cell reference.
#' @param B list: contains cell type mapping matrices (same length as tree depth).
#' @param type_order vector: specify cell type order in the bottom layer of tree.
#' @param test the "test.type" in scran::findMarkers. The test used to find markers. test  = c("t", "wilcox", "binom"). Default is "wilcox".
#' @param nmrk numerical, number of marker genes chosen for each cell type, default is 50.
#'
#' @return a list in which elements are marker gene names for each layer of tree.
#' @export
#' @importFrom scran findMarkers
HiDecon_marker <- function(ref, cell_type, B, type_order,  test="wilcox", nmrk = 50){
  ref <- log2(ref+1)
  CT.mapping <- diag(1, nrow = ncol(B[[length(B)]]), ncol = ncol(B[[length(B)]]))
  CTs <- type_order
  hierarchical_markers <- list()
  for(i in length(B):1){
    cell_type_2 <- cell_type
    if(i < length(B)) CT.mapping <- B[[i+1]] %*% CT.mapping
    for(j in 1:nrow(CT.mapping)) {
      CT.in.cluster <- CTs[as.logical(CT.mapping[j,])]
      cell_type_2[cell_type_2 %in% CT.in.cluster] <- paste0("cluster",j)
    }
    out <-  findMarkers(ref, groups=cell_type_2,direction="up",
                        test=test,pval.type ="all")
    markers <- unique(unlist(lapply(out, function(x) x@rownames[1:nmrk])))
    hierarchical_markers <- c(list(markers), hierarchical_markers)
  }
  return(hierarchical_markers)
}




#' @title HiDecon_input.
#' @description Calculate inputs for HiDecon and
#'
#' @param bulk matrix, bulk data of genes by samples. It should be count per million normalized.
#' @param ref matrix: Single cell data matrix (genes by cells) after count per million normalization.
#' @param cell_type vector: cell type label of single cell reference.
#' @param B list: contains cell type mapping matrices (same length as tree depth).
#' @param type_order vector: specify cell type order in the bottom layer of tree.
#' @param test the "test.type" in scran::findMarkers. The test used to find markers. test  = c("t", "wilcox", "binom"). Default is "wilcox".
#' @param nmrk numerical, number of marker genes chosen for each cell type, default is 50.
#'
#' @return a list with input elements for function HiDecon and select_HiDecon.
#' @export
#'
HiDecon_input <- function(bulk, ref,  cell_type, B, type_order, test = "wilcox", nmrk = 50){
  markers <- HiDecon_marker(ref = ref, cell_type = cell_type, B = B,
                            type_order = type_order, test = test, nmrk = 50)
  Signature.list <- list()
  for(i in 1:length(B)){
    Signature.list <- c(Signature.list ,
                        list(reference(ref = t(log2(ref[markers[[i]],] + 1)),
                                       ref.type = cell_type, B = B,
                                       type_order = type_order, layer = i)))
  }
  Cell.size.list <- lapply(Signature.list, colMeans)
  A.tilde <- Create.A(Signature.list, Cell.size.list, standardize=T)
  B.tilde <- Create.B(B.list = B[-1], standardize=T)
  Y.list <- list()
  for(i in 1:length(B)){
    Y.list <- c(Y.list, list(t(log2(bulk[markers[[i]],]+1))))
  }
  return(list("hierarchical_marker" = markers, "sig_list" = Signature.list,
              "cellsize_list" = Cell.size.list, "A.tilde" = A.tilde,
              "B.tilde" = B.tilde, "Y.list" = Y.list))
}


#' @title HiDecon.
#' @description Estimate cellular fractions for bulk tissue data using single cell reference with hierarchical tree.
#'
#' @param bulk matrix, bulk data of genes by samples. It should be count per million normalized.
#' @param ref matrix: Single cell data matrix (genes by cells) after count per million normalization.
#' @param B list: contains cell type mapping matrices (same length as tree depth).
#' @param cell_type vector: cell type label of single cell reference.
#' @param type_order vector: specify cell type order in the bottom layer of tree.
#' @param lambda numerical, the penalty parameter lambda in the HiDecon paper.
#' @param Pi.start vector, non-negative starting point for the estimation of sample i. Default is NULL, then function will assign the starting point.
#' @param max.iter the maximum number of iterations, default is 1e4.
#' @param tol tolerance for the KKT conditions, default is 1e-6.
#' @param test the "test.type" in scran::findMarkers. The test used to find markers. test  = c("t", "wilcox", "binom"). Default is "wilcox".
#' @param nmrk numerical, number of marker genes chosen for each cell type, default is 50.
#'
#' @return matrix, HiDecon cellular fraction estimates for samples.
#' @export
#'
HiDecon <- function(bulk, ref, B, cell_type, type_order,
                    lambda = 40, Pi.start=NULL, max.iter=1e4, tol=1e-6,
                    test = "wilcox", nmrk = 50){
  ind = intersect(rownames(ref),rownames(bulk))
  bulk <- bulk[ind,]
  ref <- ref[ind,]
  input_dat <- HiDecon_input(bulk = bulk, ref = ref,  cell_type = cell_type,
                             B = B, type_order = type_order, test = test,
                             nmrk = nmrk)
  HiDecon_est <- Est.AllPi(Y.list = input_dat$Y.list, A = input_dat$A.tilde,
                           B = input_dat$B.tilde, lambda = lambda,
                           Pi.start=NULL, max.iter=max.iter, tol=tol)
  colnames(HiDecon_est$Pi) <- type_order
  rownames(HiDecon_est$Pi) <- colnames(bulk)
  return(HiDecon_est$Pi)
}



#' @title GeneratePseudoBulk.
#' @description Generate bulk data surrogate for tuning parameter selection by resampling from single cell reference based on NNLS fraction estimates.
#'
#' @param bulk matrix, bulk data of genes by samples. It should be count per million normalized.
#' @param ref matrix: Single cell data matrix (genes by cells) after count per million normalization.
#' @param cell_type vector: cell type label of single cell reference.
#' @param type_order vector: specify cell type order in the bottom layer of tree.
#' @param test the "test.type" in scran::findMarkers. The test used to find markers. test  = c("t", "wilcox", "binom"). Default is "wilcox".
#' @param nmrk numerical, number of marker genes chosen for each cell type, default is 50.
#' @param seed random seed for resampling. Default is 123.
#'
#' @return a list with elements:
#'    * bulk.sim: bulk data simulated using fractions estimated by NNLS
#'    * frac.sim: fractions estimated by NNLS
#' @export
#' @importFrom scran findMarkers
#' @importFrom glmnet glmnet
#' @importFrom stats coef
GeneratePseudoBulk <- function(bulk, ref, cell_type, type_order, test = "wilcox", nmrk = 50, seed = 123){
  set.seed(seed)
  out <-  findMarkers(log2(ref + 1), groups = cell_type, direction = "up", test=test, pval.type ="all")
  markers <-  unlist(lapply(out, function(x) x@rownames[1:nmrk]))
  sig <- GeneratesigReference (log2(ref[markers,]+1),cell_type)
  sig <- sig[,type_order]
  sim_prop <-  apply(t(log2(bulk[markers,]+1)),1,function(x)coef(glmnet(sig[markers,],x,lambda = 0, lower.limits = 0,intercept = FALSE,standardize = FALSE))[-1,])
  sim_prop <-  t(sim_prop)/colSums(sim_prop)
  sim <- matrix(0, nrow = nrow(bulk), ncol = ncol(bulk))
  sim.bulk.frac <- sim_prop[1:ncol(sim),]
  sim.bulk.frac <- round(1000*sim.bulk.frac) + 2
  colnames(sim.bulk.frac) <- colnames(sim_prop)
  for(i in 1:nrow(sim.bulk.frac)){
    for(j in 1:length(unique(cell_type))){
      sim[,i] <- sim[,i] + rowSums(ref[,sample(which(cell_type==colnames(sim.bulk.frac)[j]),size = sim.bulk.frac[i,j], replace = TRUE)])
    }
    sim[,i] <- sim[,i]/sum(sim.bulk.frac[i,])
  }
  sim.bulk.frac <- t(apply(sim.bulk.frac, MARGIN = 1, function(x){x/sum(x)}))
  bulk.sim <- list("bulk.sim" = sim, "frac.sim" = sim.bulk.frac)
}


#' @title select_HiDecon
#' @description Estimate cellular fractions for bulk tissue data using single cell reference with hierarchical tree. Use parameter selected by the resampling parameter selection method.
#'
#' @param bulk bulk matrix, bulk data of genes by samples. It should be count per million normalized.
#' @param ref matrix: Single cell data matrix (genes by cells) after count per million normalization.
#' @param B list: contains cell type mapping matrices (same length as tree depth).
#' @param cell_type vector: cell type label of single cell reference.
#' @param type_order vector: specify cell type order in the bottom layer of tree.
#' @param lambda.set vector: parameter set used for parameter selection. Default is seq(10,200,10).
#' @param Pi.start vector, non-negative starting point for the estimation of sample i. Default is NULL, then function will assign the starting point.
#' @param max.iter the maximum number of iterations, default is 1e4.
#' @param tol tolerance for the KKT conditions, default is 1e-6.
#' @param test the "test.type" in scran::findMarkers. The test used to find markers. test  = c("t", "wilcox", "binom"). Default is "wilcox".
#' @param nmrk numerical, number of marker genes chosen for each cell type, default is 50.
#' @param seed random seed for resampling. Default is 123.
#'
#' @return a list with elements:
#'    * res: cellular fraction estimates by HiDecon using the selected lambda
#'    * lambda: the selected lambda
#'    * mCCC: mean Lin's concordance correlation coefficients under different lambda
#' @export
#'
select_HiDecon <- function(bulk, ref, B, cell_type, type_order,
                       lambda.set = seq(10,200,10), Pi.start=NULL, max.iter=1e4, tol=1e-6,
                       test = "wilcox", nmrk = 50, seed = 123){
  ind = intersect(rownames(ref),rownames(bulk))
  bulk <- bulk[ind,]
  ref <- ref[ind,]
  bulk.sim <- GeneratePseudoBulk(bulk = bulk, ref = ref, cell_type = cell_type,
                                 type_order = type_order, test = test,
                                 nmrk = nmrk, seed = seed)
  sim <- bulk.sim$bulk.sim
  rownames(sim) <- ind
  frac.sim <- bulk.sim$frac.sim
  sim.input <- HiDecon_input(bulk = sim, ref = ref,  cell_type = cell_type,
                             B = B, type_order = type_order, test = test,
                             nmrk = nmrk)
  mCCC <- rep(0,length(lambda.set))
  names(mCCC) <- lambda.set
  for(i in 1:length(lambda.set)){
    lambda <- lambda.set[i]
    frac <- Est.AllPi(Y.list = sim.input$Y.list, A = sim.input$A.tilde,
                                     B = sim.input$B.tilde, lambda = lambda,
                                     Pi.start=NULL, max.iter=max.iter, tol=tol)$Pi
    mCCC[i] <- my.CCC(frac, frac.sim)[length(type_order)+1]
  }
  lambda <- lambda.set[which.max(mCCC)]
  input_dat <- HiDecon_input(bulk = bulk, ref = ref,  cell_type = cell_type,
                             B = B, type_order = type_order, test = test,
                             nmrk = nmrk)
  HiDecon_est <- Est.AllPi(Y.list = input_dat$Y.list, A = input_dat$A.tilde,
                           B = input_dat$B.tilde, lambda = lambda,
                           Pi.start=NULL, max.iter=max.iter, tol=tol)
  colnames(HiDecon_est$Pi) <- type_order
  rownames(HiDecon_est$Pi) <- colnames(bulk)
  return(list("res" = HiDecon_est$Pi, "lambda" = lambda, "mCCC" = mCCC))
}


#' @title my.CCC.
#' @description Calculate Lin's concordance correlation coefficients between estimates and ground truth for cell types and the mean CCC.
#'
#' @param est matrix, estimated cellular fractions (samples * cell types)
#' @param tr matrix, ground truth (sample * cell types)
#'
#' @return CCC and mean CCC (the last element of the vector)
#' @export
#' @importFrom DescTools CCC
my.CCC <- function(est, tr){
  CT.CCC <- c()
  for(i in 1:ncol(est)){
    CT.CCC <- c(CT.CCC, unlist(CCC(est[,i], tr[,i])$rho.c[1]))
  }
  tmp <- CT.CCC
  tmp <- replace(tmp, is.na(tmp), 0)
  CT.CCC <- c(CT.CCC, mean(tmp))

  return(CT.CCC)
}



























