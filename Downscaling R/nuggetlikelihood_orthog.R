nuggetlikelihood_orthog <- function(x, Phi_S_Phi, trS, n) {
  l <- nrow(Phi_S_Phi)
  
  epsilon <- .Machine$double.eps 
  Phi_S_Phi_over_nugsq <- Phi_S_Phi / (x[1]^2 + epsilon)
  
  #original code
  #Phi_S_Phi_over_nugsq <- Phi_S_Phi / (x[1]^2)
  Phi_S_Phi_over_nugsq2 <- Phi_S_Phi_over_nugsq
  
  for (j in 1:l) {
    
    p2 <- x[2] + 1/(x[1] + epsilon)
    Phi_S_Phi_over_nugsq2[, j] <- Phi_S_Phi_over_nugsq[, j]/p2
  }
  
  #for (j in 1:l) {
  #  Phi_S_Phi_over_nugsq3[, j] <- Phi_S_Phi_over_nugsq[, j] * 1 / (x[2] + 1 / x[1])
  #}
  
  logdetpart <- l * log(x[2] + 1 / x[1]) - l * log(x[2])
  tracepart <- sum(diag(Phi_S_Phi_over_nugsq2))
  
  constantpart <- n * log(x[1]) + trS / x[1]
  
  val <- logdetpart - tracepart + constantpart
  
  return(val)
}
