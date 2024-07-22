BGLlinearization_orthogonal <- function(Qarray, inv_error_variances, projected_data) {
  p <- length(inv_error_variances)
  l <- dim(Qarray)[3]
  m <- ncol(projected_data)
  
  idiag <- diag(inv_error_variances)
  S_star_array <- Qarray
  
  for (i in 1:l) {
    blockindices <- ((i-1) * p + 1):(i * p)
    
    ll <- solve(Qarray[, , i] + idiag)
    tt <- ll %*% projected_data[blockindices, ]
    S_star_array[, , i] <- ll + (tt %*% t(tt)) / m
  }
  
  S_star_array
}
