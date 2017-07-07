ebmme.FitNewTissue <- function(data, cat, block_id, nclust, levels,
                    fs_binned, inner.iters){

  bins <- matrix(NA, nrow=nrow(data), ncol=ncol(data))
  for (i in 1:ncol(data)){
    levs <- unlist(strsplit(levels[, i], ","))
    levs <- gsub("\\(", "", levs)
    levs <- gsub("]", "", levs)
    levs <- sort(as.numeric(unique(levs)))
    levs[1] <- min(levs[1], min(data[, i]))
    levs[length(levs)] <- max(levs[length(levs)], max(data[, i]))
    bins[, i] <- cut(data[, i], breaks=levs, include.lowest=T)
    bins[, i] <- bins[, i] - 1
  }
  fs_binned_num <- as.numeric(fs_binned)
  fs_binned_num[fs_binned_num==0] = min(fs_binned_num[fs_binned_num>0])
  dims <- dim(fs_binned)
  fs_binned <- array(fs_binned_num, dim=dims)
  alpha<-vector('numeric',nclust)
  for (k in 1:nclust){
    alpha[k] <- 1;
  }
  newtissue <- newtissue(cat, block_id, nclust, bins, fs_binned, alpha, inner.iters)
  p <- newtissue$p
  a <- newtissue$a 
  list(p=p, a=a)
}

ebmme.Predict <- function(data, cat, block_id, nclust, levels, 
                    fs_binned, a){
  
  bins <- matrix(NA, nrow=nrow(data), ncol=ncol(data))
  for (i in 1:ncol(data)){
    levs <- unlist(strsplit(levels[, i], ","))
    levs <- gsub("\\(", "", levs)
    levs <- gsub("]", "", levs)
    levs <- sort(as.numeric(unique(levs)))
    levs[1] <- min(levs[1], min(data[, i]))
    levs[length(levs)] <- max(levs[length(levs)], max(data[, i]))
    bins[, i] <- cut(data[, i], breaks=levs, include.lowest=T)
    bins[, i] <- bins[, i] - 1
  }
  fs_binned_num <- as.numeric(fs_binned)
  fs_binned_num[fs_binned_num==0] = min(fs_binned_num[fs_binned_num>0])
  dims <- dim(fs_binned)
  fs_binned <- array(fs_binned_num, dim=dims)
  prediction <- predictlogsum(cat, block_id, nclust, bins, fs_binned, a)
  p <- prediction$p
  a <- prediction$a
  f <- prediction$f
  list(p=p, a=a, f=f)
}
 
# cat must start with 0 and consist of consecutive integers for the different categories
ebmme.lda <- function(data, cat, block.id, iters = 100, 
                            inner.iters = 100,
                            nclust=10, kde.nbins=NULL)
{
  set.seed(1)
  data <- as.matrix(data)	
  m <- nrow(data)
  k <- ncol(data)
  B <- length(unique(block.id))
  
  H.inv <- matrix(0, k, k)
  for (l in 1:k){
    data.l <- as.matrix( data[, l] )
    h <- 0.9*m^-0.2*min(stats::sd(data.l), stats::IQR(data.l)/1.34)
    H.inv[l, l] <- 1/h
  }
  
  H1.inv <- H.inv
  rm(H.inv)
  # initialize p
  p <- matrix(0, nclust, m)
  for (k in 1:nclust)
    for (i in 1:m)
      p[k,i]=stats::runif(1);
  for (k in 1:nclust)
    for (i in 1:m)
      p[k,i]=p[k,i]/sum(p[,i]);
  
  
  # initialize alpha
  alpha<-vector('numeric',nclust)
  for (k in 1:nclust){
    alpha[k] <- 1;
  }
  
  bins <- matrix(NA, nrow=nrow(data), ncol=ncol(data))
  bin.data <- matrix(NA, nrow=kde.nbins, ncol=ncol(data))
  levels.data <- matrix(NA, nrow=kde.nbins, ncol=ncol(data))
  colnames(levels.data) <- colnames(data)
  for (i in 1:ncol(data)){
    # factor information is lost when put into matrix
    cc <- cut(data[, i], breaks=kde.nbins)
    bins[, i] <- cc
    bottoms <- unlist(strsplit(levels(cc), ","))[c(T,F)]
    tops <- unlist(strsplit(levels(cc), ","))[c(F,T)]
    bottoms <- as.numeric(gsub("\\(", "", bottoms))
    tops <- as.numeric(gsub("]", "", tops))
    bin.data[, i] <- (bottoms + tops) / 2
    levels.data[, i] <- levels(cc)
  }
  # need to change array indexing for C++, bin number and category number should start from 0
  bins <- bins - 1
  fit <- ebmme_cpp_binned(data, cat, block.id, H1.inv, p, alpha, iters, 
                          inner.iters, nclust, 
                          bins, bin.data, kde.nbins)
  p <- as.matrix( fit$p )
  alpha <- as.vector( fit$alpha )
  a <- as.matrix( fit$a )
  f <- as.matrix( fit$f )
  return(list(p=p, alpha=alpha,
              a=a, f=f, fs_binned=fit$fs_binned, 
              bins=bins, bin.data=bin.data, kde.nbins=kde.nbins, 
              data=data, cat=cat, block.id=block.id, H1.inv=H1.inv, 
              iters=iters, inner.iters=inner.iters, 
              nclust=nclust,  
              levels.data=levels.data, all_as=fit$all_as))
}
