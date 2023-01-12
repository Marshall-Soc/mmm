
# NOTE: Adapted from https://rdrr.io/cran/cluster/src/R/clusGap.R

clusGap2 <- function (x, FUNcluster, K.max, B = 100, verbose = interactive(),
          do_parallel = FALSE, ...)
{
  stopifnot(is.function(FUNcluster), length(dim(x)) == 2, K.max >= 1,
            (n <- nrow(x)) >= 1, (p <- ncol(x)) >= 1)
  if(B != (B. <- as.integer(B)) || (B <- B.) <= 0)
    stop("'B' has to be a positive integer")
  
  apply_fun <- lapply
  
  if (do_parallel)
  {
    apply_fun <- parallel::mclapply
  }
  
  if(is.data.frame(x))
    x <- as.matrix(x)
  ii <- seq_len(n)
  W.k <- function(X, kk) {
    clus <- if(kk >= 1) FUNcluster(X, kk, ...)$cluster else rep.int(1L, nrow(X))
    ##                 ---------- =  =       -------- kmeans() has 'cluster'; pam() 'clustering'
    0.5* sum(vapply(split(ii, clus),
                    function(I) { xs <- X[I,, drop=FALSE]
                    sum((1-(cca(xs, filter.significance=T, filter.value=0.05, zero.action=c("drop"), 
                                verbose=F)$cormat))/nrow(xs)) }, 0.))
  }
  
  E.logW <- SE.sim <- numeric(K.max)
  
  if(verbose) cat("Clustering with modularity steps ar k = 0,1,..., K.max (= ",K.max,"): .. ", sep='')
  
  logW <- unlist(apply_fun(0:K.max, function(k) log(W.k(x, k))))
  
  if(verbose) cat("done\n")
  
  ## Scale 'x' into "hypercube" -- we later fill with H0-generated data
  xs <- scale(x, center=TRUE, scale=FALSE)
  m.x <- rep(attr(xs,"scaled:center"), each = n)# for back transforming
  V.sx <- svd(xs)$v
  rng.x1 <- apply(xs %*% V.sx, # = transformed(x)
                  2, range)
  
  if(verbose) cat("Bootstrapping, b = 1,2,..., B (= ", B,
                  ")  [one \".\" per sample]:\n", sep="")
  
  logWksList <- apply_fun(1:B,
                          function(b)
                          {
                            ## Generate "H0"-data as "parametric bootstrap sample" :
                            z1 <- apply(rng.x1, 2,
                                        function(M, nn) runif(nn, min=M[1], max=M[2]),
                                        nn=n)
                            z <- tcrossprod(z1, V.sx) + m.x # back transformed
                            curLogWks <- unlist(lapply(0:K.max, function(k) log(W.k(z, k))))
                            if(verbose && !do_parallel) cat(".", if(b %% 50 == 0) paste(b,"\n"))
                            
                            curLogWks
                          })
  
  logWks <- matrix(unlist(logWksList),
                   byrow = TRUE,
                   nrow = B)
  
  if(verbose && (B %% 50 != 0)) cat("",B,"\n")
  E.logW <- colMeans(logWks)
  SE.sim <- sqrt((1 + 1/B) * apply(logWks, 2, var))
  structure(class = "clusGap",
            list(Tab = cbind(logW, E.logW, gap = E.logW - logW, SE.sim),
                 ## K.max == nrow(T)
                 n = n, B = B, FUNcluster=FUNcluster))
}