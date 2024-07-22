library(pracma) 
library(Matrix)  
library(nloptr)
source(file.path(path,"nuggetlikelihood_orthog.R"))
source(file.path(path,"BGLblockdiag_nopenalty.R"))
source(file.path(path,"BGLlinearization_orthogonal.R"))

# Read Obs_data.csv and convert to matrix
Obs <- as.matrix(read.csv(file.path(path, 'Obs_data.csv')))
Obs_pixels <- Obs[, 1:2]
Obs <- Obs[, -c(1, 2)]

# Load stand downscaled training data from Stand_Downscaled_training.csv
Stand_SD <- as.matrix(read.csv(file.path(path, 'Stand_Downscaled_training.csv')))
Stand_SD <- Stand_SD[, -c(1, 2)]

# Calculate dimensions
N <- nrow(Obs)
T_o <- ncol(Obs)

# Load Mu1 and Mu2 data
Mu1 <- as.matrix(read.csv(file.path(path, 'Mu1.csv')))
Mu2 <- as.matrix(read.csv(file.path(path, 'Mu2.csv')))

# Subtract the mean from the processes
Z1 <- Stand_SD[, 1:T_o] - Mu1[, 1:T_o]
Z2 <- Obs - Mu2[, 1:T_o]

# Remove unnecessary columns
Stand_SD <- Stand_SD[, -(1:T_o)]
Mu1 <- Mu1[, -(1:T_o)]
Mu2 <- Mu2[, -(1:T_o)]

# Calculate Z1 downscale and Z2 mean downscale
Z1_downscale <- Stand_SD - Mu1
Z2_mean_downscale <- Mu2

# Define some constants and dimensions
t_downscale <- ncol(Z1_downscale)
N <- nrow(Obs)

# Setup Data

# Create a data_array with two variables, N locations, and T_o realizations
data_array <- array(NA, dim = c(2, N, T_o))
data_array[1, , ] <- Z1
data_array[2, , ] <- Z2
p <- dim(data_array)[1]
n <- dim(data_array)[2]
m <- dim(data_array)[3]

# Print the results (optional)
cat("Number of variables (p):", p, "\n")
cat("Number of locations (n):", n, "\n")
cat("Number of realizations (m):", m, "\n")

####################################################################################
############################BASIS FUNCTIONS############################
####################################################################################


# Calculate orthogonal basis matrix using EOFs
reshaped_data_array <- cbind(Z1, Z2)

# Perform Singular Value Decomposition (SVD)
svd_result <- svd(reshaped_data_array)
#U <- svd_result$u
U <- -svd_result$u
S <- diag(svd_result$d)

# Calculate cumulative variance
S_squared <- diag(S)^2
cumulative_variance <- cumsum(S_squared) / sum(S_squared)

# Find number of components explaining 99.99% variance
Var99_explained <- which(cumulative_variance >= 0.9999)[1]
U <- U[, 1:Var99_explained]

# Retain basis functions for the stochastic term
smallPhi <- U
l <- ncol(smallPhi)






# Nugget Estimate

# Initialize variables
x0 <- c(1, 1)
lb <- c(0.0001,0.0001)
ub <- c(50, 50)
nuggetestimates <- numeric(p)
smoothnessestimates <- numeric(p)

################################################################################################################################
for (i in 1:p) {
  tt <- t(smallPhi) %*% data_array[i, , ]
  #tt <- t(smallPhi) %*% array(data_array[i, , ], dim = c(ncol(data_array), ncol(data_array)))
  #smallPhi_S_Phi <- tt %*% t(tt) / m
  smallPhi_S_Phi <- tt %*% t(tt) / m
  trS <- norm(data_array[i, , ], "F")^2 / m
  #result <- fmincon(x0 = x0, fn = nuggetlikelihood_orthog, lb = lb, ub = ub,Phi_S_Phi = smallPhi_S_Phi, trS = trS, n = n)
  #result <- optim(par = x0, fn = nuggetlikelihood_orthog, Phi_S_Phi = smallPhi_S_Phi, trS = trS, n = n,method="L-BFGS-B", lower = lb, upper = ub)
  
  
  opts <- list("algorithm"="NLOPT_LN_COBYLA",             "xtol_rel"=1.0e-10, maxeval=1000, "xtol_abs"=1.0e-10)
  eval_g0 <- function( x,Phi_S_Phi , trS , n ) {    return( x )}
  eval_jac_g0 <- function( x,Phi_S_Phi , trS , n ) {    return( rbind( c( 1, 0 ),                    c( 0, 1.0 ) ) )}
  result <- nloptr( x0=c(0.0001,0.001),                eval_f=nuggetlikelihood_orthog,                lb=c(0.000001, 0.00001),                ub = c(1,1),                eval_g_ineq = eval_g0,               eval_jac_g_ineq =  eval_jac_g0,               opts=opts,                 Phi_S_Phi = smallPhi_S_Phi,                trS = trS,                n = n)
  result$solution
  
  
  
  
  nuggetestimates[i] <- result$solution[1]
  smoothnessestimates[i] <- result$solution[2]
}



################################################################################################################################
# matlab fmincon result
###############################################################################################################
#nuggetestimates = c(1.7190e-05, 1.8638e-04)
#smoothnessestimates = c(2.8473e-04,1.2555e-04)

################################################################################################################################
# R nloptr result
###############################################################################################################

#nuggetestimates = c(9.999909e-06, 1.000002e-05)
#smoothnessestimates = 1.000002e-05 1.000121e-05

inv_error_variances <- 1 / nuggetestimates

# Projecting Data and Creating Matrix
projected_data <- matrix(0, nrow = p * l, ncol = m)

for (i in 1:m) {
  projected_data[, i] <- as.vector((inv_error_variances * data_array[, , i]) %*% smallPhi)
}

# Main Algorithms



mleguess <- BGLblockdiag_nopenalty(inv_error_variances, projected_data)

# Calculate linear transformation matrices A1 and Q
A1_start <- diag(l)
A1 <- matrix(NA, nrow = l, ncol = p * l)
for (r in 1:l) {
  col1 <- 2 * r - 1
  col2 <- 2 * r
  A1[1:l, col1:col2] <- cbind(rep(0, l), A1_start[1:l, r])
}

Q <- bdiag(lapply(1:l, function(r) mleguess$Qarray[, , r]))

# Calculate Q_inv
mleguess_inv <- array(NA, dim = c(p, p, l))
for (r in 1:l) {
  mleguess_inv[, , r] <- solve(mleguess$Qarray[, , r])
}
Q_inv <- bdiag(lapply(1:l, function(r) mleguess_inv[, , r]))

Predicted_temp_all <- matrix(NA, nrow = N, ncol = t_downscale + 2)
Predicted_sd_all <- matrix(NA, nrow = N, ncol = t_downscale + 2)

Predicted_temp_all[, 1:2] <- Obs_pixels
Predicted_sd_all[, 1:2] <- Obs_pixels

# Calculate C and Sigma for future prediction
C <- diag(l)
Sigma <- diag(l)
for (r in 1:l) {
  C[r, r] <- mleguess_inv[1, 2, r] / mleguess_inv[1, 1, r]
  Sigma[r, r] <- mleguess_inv[2, 2, r] - (mleguess_inv[1, 2, r])^2 / mleguess_inv[1, 1, r]
}

# Begin downscaling
Start <- proc.time()
for (month in 1:t_downscale) {
  
  # New Z1 and the corresponding mean for future prediction
  Z1_new <- Z1_downscale[, month]
  Mean_new <- Z2_mean_downscale[, month]
  
  SC <- smallPhi %*% C
  Inv <- solve(solve(A1 %*% Q_inv %*% t(A1)) + diag(l) / nuggetestimates[1])
  
  SY <- t(smallPhi) %*% Z1_new
  Predicted_temp <- Mean_new + (SC %*% Inv %*% SY) / nuggetestimates[1]
  
  Matrix <- Sigma + C %*% Inv %*% t(C)
  Predicted_sd <- apply(smallPhi, 1, function(row) sqrt(t(row) %*% Matrix %*% row + nuggetestimates[2]))
  Predicted_sd_l  <- sapply(Predicted_sd , function(x) as.numeric(x))
  
  Predicted_temp_all[, month + 2] <- Predicted_temp@x
  Predicted_sd_all[, month + 2] <- Predicted_sd_l
  
  cat(paste("t", month, " done!\n", sep = ""))
}
cat("\n")

# Save the results to CSV files
write.csv(Predicted_temp_all, file.path(path, 'BGL_Downscaled.csv'), row.names = FALSE)
write.csv(Predicted_sd_all, file.path(path, 'BGL_Downscaled_UQ.csv'), row.names = FALSE)

# Delete unnecessary CSV files
file.remove(file.path(path, 'Stand_Downscaled_training.csv'))
file.remove(file.path(path, 'Mu1.csv'))
file.remove(file.path(path, 'Mu2.csv'))






