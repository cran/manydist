mahalanobis_like <- function(x, validate_x=NULL,  sig = NULL) {
  x=as.matrix(x)
  
  if (is.null(sig)) {
    sig <- solve(cov(x))
  }
  if (is.null(validate_x)) {
    validate_x <- x
  } else {
    validate_x=as.matrix(validate_x)
  }
  

  if (det(sig) == 0) {
    stop("The covariance matrix is singular and cannot be inverted.")
  }
  
  dist <- matrix(0, nrow = nrow(x), ncol = nrow(validate_x))
  for (i in 1:nrow(x)) {
    for (j in 1:nrow(validate_x)) {
      diff <- x[i, ] - validate_x[j, ]
      dist[i, j] <- sqrt(t(diff) %*% sig %*% diff)
    }
  }
  
  dist=t(dist)
  return(dist)
}
