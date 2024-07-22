BGLblockdiag_nopenalty <- function(inv_error_variances, projected_data) {
  p <- length(inv_error_variances)
  l <- nrow(projected_data) / p
  
  Qarray <- array(0, dim = c(p, p, l))
  for (i in 1:l) {
    Qarray[, , i] <- diag(p)
  }
  
  tol <- 0.05
  my_norm <- 1
  counter <- 1
  max_iterations <- 1000
  
  while (my_norm >= tol && counter <= max_iterations) {
    
    cat(sprintf("Starting Iter: %d", counter))
    cat("...linearization...")
    
    S_star_array <- BGLlinearization_orthogonal(Qarray, inv_error_variances, projected_data)
    
    sum_sq <- numeric(l)
    sum_sq_diff <- numeric(l)
    
    cat("...BGL...\n")
    
    for (i in 1:l) {
      Qguess <- solve(S_star_array[, , i])
      sum_sq_diff[i] <- norm(Qguess - Qarray[, , i], type = "F")^2
      sum_sq[i] <- norm(Qarray[, , i], type = "F")^2
      Qarray[, , i] <- Qguess
    }
    
    my_norm <- sqrt(sum(sum_sq_diff)) / sqrt(sum(sum_sq))
    
    cat("\n")
    cat(sprintf("Ending Iter: %d, norm: %.12e\n", counter, my_norm))
    cat("\n")
    
    counter <- counter + 1
    
  }
  
  list(Qarray = Qarray, counter = counter)
}
