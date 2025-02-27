
kernel <- function(x1, x2, sigma = 1, lengthscale = c(0.1,0.1)) {
  sigma**2 * exp(-0.5*sum(((x1-x2)/lengthscale)**2))
}

which.maxn <- function(x,n=1){
  if (n==1)
    which.max(x)
  else
  {
    if (n>1){
      ii <- order(x,decreasing=TRUE)[1:min(n,length(x))]
      ii[!is.na(x[ii])]
    }
    else {
      stop("n must be >=1")
    }
  }
}


getmode <- function(v, dist) {
  tab <- table(v)
  #if tie, break to shortest distance
  if(sum(tab == max(tab)) > 1){
    tied <- names(tab)[tab == max(tab)]
    sub  <- dist[v %in% tied]
    names(sub) <- v[v %in% tied]
    return(names(sub)[which.min(sub)])
  } else {
    return(names(tab)[which.max(tab)])
  }
}

# Loop using a squared exponential kernel
# S <- matrix(as.numeric(NA), nrow=N, ncol=N)
# for (i in 1:nrow(S)) {
#   print(i)
#   for (j in 1:ncol(S)) {
#     S[i,j] <- kernel(foo[i,],foo[j,])
#   }
# }


# linear kernel. IT DOESNT WORK, WHY? UNCENTERED DATA???
# S <- spatial.coordinates%*%t(spatial.coordinates)
