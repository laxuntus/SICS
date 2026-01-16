library(elasticnet)

# Function for fixing signs columnwise
fix_signs <- function(B){
  K <- ncol(B)
  p <- nrow(B)
  
  for(i in 1:K){
    my_col <- B[, i]
    nonzero_ind <- which(abs(my_col) >= 1e-4)[1]
    B[, i] <- my_col*sign(my_col[nonzero_ind])
  }
  
  B
}

# The main algorithm for sparse invariant coordinate selection
# input:
# S1, S2 are (scatter) matrices
# K is the number of extracted invariant coordinates
# varnum_vec is a vector of non-zero components for each coordinate
#or a number corresponding to same number for every coordinate
# output:
# B is the ICS loading matrix

SICS = function(S1, S2, K = ncol(S1), varnum_vec = NULL, maximiter = 500){
  p = ncol(S1)
  
  if(length(varnum_vec) == 1){
    varnum_vec <- rep(varnum_vec, K)
  }
  if(is.null(varnum_vec)){
    varnum_vec <- rep(p, K)
  }
  
  S1_inv = solve(S1)
  eig_S1 = eigen(S1, symmetric = TRUE)
  S1_sqrt = eig_S1$vectors %*% diag(sqrt(eig_S1$values)) %*% t(eig_S1$vectors)
  S1_invsqrt <- solve(S1_sqrt)
  
  eig_S2 = eigen(S2, symmetric = TRUE)
  R = eig_S2$vectors %*% diag(sqrt(eig_S2$values)) %*% t(eig_S2$vectors)
  
  A <- eigen(S1_invsqrt %*% S2 %*% S1_invsqrt)$vectors[, 1:K, drop = FALSE]
  B <- S1_invsqrt %*% A
  B <- fix_signs(B)
  
  B_temp <- matrix(0, p, K)
  
  iterat = 0
  while (sum((B_temp-B)^2) > 1e-12) {
    iterat = iterat+1
    
    B_temp = B
    
    for (i in 1:K) {
      B[, i] <- solvebeta(R, R%*%S1_invsqrt%*%A[, i], paras = c(1e-8, varnum_vec[i]), sparse = "varnum")
    }
    B <- fix_signs(B)
    
    svd_A = svd(S1_invsqrt%*%S2%*%B)
    A_0 = svd_A$u %*% t(svd_A$v)
    
    O = eigen(t(A_0)%*%S1_invsqrt%*%S2%*%S1_invsqrt%*%A_0)$vectors
    A = A_0%*%O
    
    if(iterat == maximiter){
      break
    }
  }
  
  return(list(B=B, A=A, iters=iterat))
}

###
# Example
###

library(ICS)

# Unmixing matrix to be estimated
A_inv = matrix(c(1, 0, -2, 0,
                 2, 0, 0, -2,
                 0, -1, 0, -4,
                 0, 0, 1, 1), 4, 4, byrow = TRUE)

# Generating the data
A <- solve(A_inv)
n <- 100000
X = cbind(rpois(n, 1)-1, (rpois(n, 2)-2)/sqrt(2), rnorm(n), (runif(n)-0.5)*sqrt(12)) %*% t(A)

# Scatter matrices: Covariance matrix and the scatter matrix based on the 4th moments
S1 = cov(X)
S2 = cov4(X)

# Estimating four independent components with three non-zero coefficients
model_SICS = SICS(S1, S2, 4, 3)

t(model_SICS$B) # Estimate
A_inv # Real value
