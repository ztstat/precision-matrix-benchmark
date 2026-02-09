######################
rm(list = ls())
######################

################## Packages ##################
library(QUIC)
library(dpglasso)
#library(sglasso)
library(glasso)
library(JGL)
#library(tailoredGlasso)
library(hglasso)
library(CVglasso)
library(MASS)

library(foreach)
library(doParallel)



##################  dpglasso.new  ############
dpglasso.new <- function (Sigma, X = NULL, invX = NULL, rho, outer.Maxiter = 100, 
                          obj.seq = FALSE, outer.tol = 10^-5) 
{
  if (is.null(Sigma)) {
    print("Sigma is required as input")
    return()
  }
  else {
    p = nrow(Sigma)
  }
  check1 <- (nrow(Sigma) == ncol(Sigma)) && (mean(Sigma == t(Sigma))==1)
  if (check1 == 0) {
    print("Sigma sould be square and symmetric")
    return()
  }
  if (is.null(X)) 
    X = diag(1/(rep(rho, p) + diag(Sigma)))
  if (is.null(invX)) {
    U.mat = diag(rep(rho, p))
    invX <- Sigma + U.mat
  }
  else {
    U.mat <- invX - Sigma
  }
  p1 <- nrow(X)
  p11 <- ncol(X)
  p2 <- nrow(invX)
  p22 <- ncol(invX)
  if (p1 != p11) {
    print("check dimensions of X")
    return()
  }
  if (p2 != p22) {
    print("check dimensions of invX")
    return()
  }
  check.dim <- (p1 == p2) && (p1 == p)
  if (check.dim == 0) {
    print("check dimensions of X,invX,Sigma")
    return()
  }
  ids <- rep(1:p, outer.Maxiter)
  time.counter.QP <- array(0, dim = c(length(ids), 3))
  rel.err <- rep(0, outer.Maxiter)
  obj.vals <- rep(0, outer.Maxiter)
  sparse.nos <- rep(0, outer.Maxiter)
  ii <- 0
  kk <- 0
  tol = 10
  vec.diagsold <- X
  for (cd.iter in ids) {
    ii <- ii + 1
    t <- proc.time()
    obj <- box_qp_f(X[-cd.iter, -cd.iter], u = U.mat[-cd.iter, 
                                                     cd.iter], b = Sigma[-cd.iter, cd.iter], rho, Maxiter = 1000, 
                    tol = 10^-7)
    t <- (proc.time() - t)
    time.counter.QP[ii, ] <- as.numeric(c(t[1], t[2], t[3]))
    theta.hat <- -obj$grad_vec * 0.5/(Sigma[cd.iter, cd.iter] + 
                                        rho)
    theta.hat[abs(obj$u) < rho * 0.99999] = 0
    theta.hat[abs(theta.hat) < 10^-5] = 0
    X[cd.iter, -cd.iter] <- theta.hat
    X[-cd.iter, cd.iter] <- X[cd.iter, -cd.iter]
    X[cd.iter, cd.iter] = sum((obj$u + Sigma[cd.iter, -cd.iter]) * 
                                X[cd.iter, -cd.iter])
    X[cd.iter, cd.iter] <- (1 - X[cd.iter, cd.iter])/(Sigma[cd.iter, 
                                                            cd.iter] + rho)
    U.mat[-cd.iter, cd.iter] <- as.numeric(obj$u, ncol = 1)
    U.mat[cd.iter, -cd.iter] <- U.mat[-cd.iter, cd.iter]
    U.mat[cd.iter, cd.iter] = rho
    if (cd.iter == p) {
      kk <- kk + 1
      vec.diagsnew <- X
      tol = max(abs(vec.diagsnew - vec.diagsold))/max(abs(vec.diagsold))
      rel.err[kk] <- tol
      vec.diagsold <- vec.diagsnew
      sparse.nos[kk] <- sum(abs(X) <= 10^-9)
      if (obj.seq == TRUE) 
        obj.vals[kk] <- -sum(log(abs(diag(qr(X)$qr)))) + 
        sum(Sigma * X) + rho * sum(abs(X))
      if ((tol < outer.tol) & (kk > 1)) 
        break
    }
  }
  time.counter.QP = time.counter.QP[1:ii, ]
  time.counter.QP = colSums(time.counter.QP)
  if (obj.seq == TRUE) {
    return(list(time.counter.QP = time.counter.QP, rel.err = rel.err[1:kk], 
                obj.vals = obj.vals[1:kk], X = X, invX = Sigma + 
                  U.mat, sparse.nos = sparse.nos[1:kk]))
  }
  else {
    return(list(time.counter.QP = time.counter.QP, rel.err = rel.err[1:kk], 
                X = X, invX = Sigma + U.mat, sparse.nos = sparse.nos[1:kk]))
  }
}




####################################  Settings 

##################  Setting 1: Banded ################## ##################
Banded <- function(p, diag, offdiag) {
  Sigma_inv <- matrix(0, nrow = p, ncol = p)
  diag(Sigma_inv) <- diag
  for (i in 2:p) {
    Sigma_inv[i, i - 1] <- offdiag
    Sigma_inv[i - 1, i] <- offdiag
  }
  return(Sigma_inv)
}

##################  Setting 2: Dense ################## ##################
Dense <- function(p, diag, offdiag) {
  Sigma_inv <- matrix(offdiag, nrow = p, ncol = p)
  diag(Sigma_inv) <- diag
  return(Sigma_inv)
}

##################  Setting 3: Random_Sparse ########## ##################
Random_Sparse <- function(p,sparsity=0.7,eta = 0.5) {
  B <- matrix(rnorm(p^2), nrow = p, ncol = p)
  B_sym <- (B + t(B)) / 2
  
  B_sparse <- B_sym
  num_zero <- floor(sparsity * (p*(p-1)/2 )) 

  upper_indices <- which(upper.tri(B_sym), arr.ind = TRUE)
  
  zero_indices <- sample(1:nrow(upper_indices), num_zero)
  rows <- upper_indices[zero_indices, 1]
  cols <- upper_indices[zero_indices, 2]

  for (i in 1:length(rows)) {
    B_sparse[rows[i], cols[i]] <- 0
    B_sparse[cols[i], rows[i]] <- 0
  }
  
  I_p <- diag(p)
  Theta <- B_sparse
  
  repeat {
    Theta <- Theta + eta * I_p
    if (min(eigen(Theta)$values) > 1) {
      break  
    }
  }
  return(Theta)
}





#################################### Criteria

########### True Rate  ##########

TPR <- function(esti,true,threshold = 1e-5) {
  num <- sum(abs(esti)>=threshold & true != 0 )/sum(true != 0)
  return(num)
}

FPR <- function(esti,true,threshold = 1e-5) {
  num <- sum(abs(esti)>=threshold & true == 0)/sum(true == 0)
  return(num)
}

Sparsity_ratio <- function(esti,p,threshold = 1e-5) {
  result <- sum(abs(esti)>=threshold)/(p^2)
  return(result)
}

Fnorm <- function(esti,true) {
  result <- norm(esti-true, type = "F")
  return(result)
}

########### BIC  ##########

BIC_function <- function(esti,S,n,threshold = 1e-5) {
  term1 <- -n*log(det(esti))
  term2 <- n*(sum(diag( S %*% esti )))
  term3 <- log(n) * sum(abs(esti)>=threshold)
  BIC <- term1 + term2 + term3
  return(BIC)
}


####






########################### Simulation  ###################################
library(foreach)
library(doParallel)



main_function_case1 <- function(n, p, diag, offdiag, lambda_grid = seq(0.1, 2, 0.05), thr = 1.0e-5) {
  
  set.seed(250)
  
  Theta_true <- Banded(p, diag, offdiag)
  data <- mvrnorm(n = n, mu = numeric(p), Sigma = solve(Theta_true))
  S <- cov(data)
  lambda_grid <- lambda_grid
  L <- length(lambda_grid)
  
  BIC <- list(glasso = numeric(L), dpglasso = numeric(L), QUIC = numeric(L))
  TPR <- list(glasso = numeric(L), dpglasso = numeric(L), QUIC = numeric(L))
  FPR <- list(glasso = numeric(L), dpglasso = numeric(L), QUIC = numeric(L))
  Sparsity_ratio <- list(glasso = numeric(L), dpglasso = numeric(L), QUIC = numeric(L))
  Fnorm <- list(glasso = numeric(L), dpglasso = numeric(L), QUIC = numeric(L))
  Time <- list(glasso = numeric(L), dpglasso = numeric(L), QUIC = numeric(L))
  
  ## Record the result under the best Lambda
  Best_lambda <- list(glasso = 0, dpglasso = 0, QUIC = 0)
  TPR_Best_lambda <- list(glasso = 0, dpglasso = 0, QUIC = 0)
  FPR_Best_lambda <- list(glasso = 0, dpglasso = 0, QUIC = 0)
  Sparsity_ratio_Best_lambda <- list(glasso = 0, dpglasso = 0, QUIC = 0)
  Fnorm_Best_lambda <- list(glasso = 0, dpglasso = 0, QUIC = 0)
  Time_Best_lambda <- list(glasso = 0, dpglasso = 0, QUIC = 0)
  
  #### Run all the Lambda and find the best Lambda ####
  total_cores <- detectCores()
  cluster_number <- ifelse(total_cores < 41, total_cores - 1, 41)
  cl <- makeCluster(cluster_number)
  
  registerDoParallel(cl)
  
  clusterExport(cl, varlist = c("dpglasso.new","Banded","Dense","Random_Sparse",
                                "TPR","FPR","Sparsity_ratio","Fnorm","BIC_function"))
  
  results <- foreach(i = 1:L, .packages = c("glasso", "dpglasso", "QUIC")) %dopar% {
    lambda <- lambda_grid[i]
    
    # Run methods and record timings
    Time_glasso <- system.time({
      solution_glasso <- glasso(S, lambda)$wi
    })
    Time_dpglasso <- system.time({
      solution_dpglasso <- dpglasso.new(S, rho = lambda)$X
    })
    Time_QUIC <- system.time({
      solution_QUIC <- QUIC(S, rho = lambda, msg = 0)$X
    })
    
    # Compute metrics
    list(
      BIC_glasso = BIC_function(solution_glasso, S = S, n = n),
      BIC_dpglasso = BIC_function(solution_dpglasso, S = S, n = n),
      BIC_QUIC = BIC_function(solution_QUIC, S = S, n = n),
      
      TPR_glasso = TPR(solution_glasso, Theta_true),
      TPR_dpglasso = TPR(solution_dpglasso, Theta_true),
      TPR_QUIC = TPR(solution_QUIC, Theta_true),
      
      FPR_glasso = FPR(solution_glasso, Theta_true),
      FPR_dpglasso = FPR(solution_dpglasso, Theta_true),
      FPR_QUIC = FPR(solution_QUIC, Theta_true),
      
      Sparsity_ratio_glasso = Sparsity_ratio(solution_glasso, p),
      Sparsity_ratio_dpglasso = Sparsity_ratio(solution_dpglasso, p),
      Sparsity_ratio_QUIC = Sparsity_ratio(solution_QUIC, p),
      
      Fnorm_glasso = Fnorm(solution_glasso, Theta_true),
      Fnorm_dpglasso = Fnorm(solution_dpglasso, Theta_true),
      Fnorm_QUIC = Fnorm(solution_QUIC, Theta_true),
      
      Time_glasso = Time_glasso[3],
      Time_dpglasso = Time_dpglasso[3],
      Time_QUIC = Time_QUIC[3]
    )
  }
  
  stopCluster(cl)
  
  # Combine results
  for (i in 1:L) {
    res <- results[[i]]
    BIC$glasso[i] <- res$BIC_glasso
    BIC$dpglasso[i] <- res$BIC_dpglasso
    BIC$QUIC[i] <- res$BIC_QUIC
    
    TPR$glasso[i] <- res$TPR_glasso
    TPR$dpglasso[i] <- res$TPR_dpglasso
    TPR$QUIC[i] <- res$TPR_QUIC
    
    FPR$glasso[i] <- res$FPR_glasso
    FPR$dpglasso[i] <- res$FPR_dpglasso
    FPR$QUIC[i] <- res$FPR_QUIC
    
    Sparsity_ratio$glasso[i] <- res$Sparsity_ratio_glasso
    Sparsity_ratio$dpglasso[i] <- res$Sparsity_ratio_dpglasso
    Sparsity_ratio$QUIC[i] <- res$Sparsity_ratio_QUIC
    
    Fnorm$glasso[i] <- res$Fnorm_glasso
    Fnorm$dpglasso[i] <- res$Fnorm_dpglasso
    Fnorm$QUIC[i] <- res$Fnorm_QUIC
    
    Time$glasso[i] <- res$Time_glasso
    Time$dpglasso[i] <- res$Time_dpglasso
    Time$QUIC[i] <- res$Time_QUIC
  }
  
  #### Also Output the result with the best lambda
  
  best.index_glasso <- which.min(BIC$glasso)
  best.index_dpglasso <- which.min(BIC$dpglasso)
  best.index_QUIC <- which.min(BIC$QUIC)
  
  Best_lambda$glasso <- lambda_grid[best.index_glasso]
  Best_lambda$dpglasso <- lambda_grid[best.index_dpglasso]
  Best_lambda$QUIC <- lambda_grid[best.index_QUIC]
  
  TPR_Best_lambda$glasso <- TPR$glasso[best.index_glasso]
  TPR_Best_lambda$dpglasso <- TPR$dpglasso[best.index_dpglasso]
  TPR_Best_lambda$QUIC <- TPR$QUIC[best.index_QUIC]
  
  FPR_Best_lambda$glasso <- FPR$glasso[best.index_glasso]
  FPR_Best_lambda$dpglasso <- FPR$dpglasso[best.index_dpglasso]
  FPR_Best_lambda$QUIC <- FPR$QUIC[best.index_QUIC]
  
  Sparsity_ratio_Best_lambda$glasso <- Sparsity_ratio$glasso[best.index_glasso]
  Sparsity_ratio_Best_lambda$dpglasso <- Sparsity_ratio$dpglasso[best.index_dpglasso]
  Sparsity_ratio_Best_lambda$QUIC <- Sparsity_ratio$QUIC[best.index_QUIC]
  
  Time_Best_lambda$glasso <- Time$glasso[best.index_glasso]
  Time_Best_lambda$dpglasso <- Time$dpglasso[best.index_dpglasso]
  Time_Best_lambda$QUIC <- Time$QUIC[best.index_QUIC]
  
  Fnorm_Best_lambda$glasso <- Fnorm$glasso[best.index_glasso]
  Fnorm_Best_lambda$dpglasso <- Fnorm$dpglasso[best.index_dpglasso]
  Fnorm_Best_lambda$QUIC <- Fnorm$QUIC[best.index_QUIC]
  
  return(list(
    BIC = BIC,
    TPR = TPR,
    FPR = FPR,
    Sparsity_ratio = Sparsity_ratio,
    Fnorm = Fnorm,
    Time = Time,
    Best_lambda = Best_lambda,
    TPR_Best_lambda = TPR_Best_lambda,
    FPR_Best_lambda = FPR_Best_lambda,
    Sparsity_ratio_Best_lambda = Sparsity_ratio_Best_lambda,
    Fnorm_Best_lambda = Fnorm_Best_lambda,
    Time_Best_lambda = Time_Best_lambda
  ))
}

main_function_case2 <- function(n, p, diag, offdiag, lambda_grid = seq(0.1, 2, 0.05), thr = 1.0e-5) {
  
  set.seed(250)
  
  Theta_true <- Dense(p, diag, offdiag)
  data <- mvrnorm(n = n, mu = numeric(p), Sigma = solve(Theta_true))
  S <- cov(data)
  lambda_grid <- lambda_grid
  L <- length(lambda_grid)
  
  BIC <- list(glasso = numeric(L), dpglasso = numeric(L), QUIC = numeric(L))
  TPR <- list(glasso = numeric(L), dpglasso = numeric(L), QUIC = numeric(L))
  FPR <- list(glasso = numeric(L), dpglasso = numeric(L), QUIC = numeric(L))
  Sparsity_ratio <- list(glasso = numeric(L), dpglasso = numeric(L), QUIC = numeric(L))
  Fnorm <- list(glasso = numeric(L), dpglasso = numeric(L), QUIC = numeric(L))
  Time <- list(glasso = numeric(L), dpglasso = numeric(L), QUIC = numeric(L))
  
  ## Record the result under the best Lambda
  Best_lambda <- list(glasso = 0, dpglasso = 0, QUIC = 0)
  TPR_Best_lambda <- list(glasso = 0, dpglasso = 0, QUIC = 0)
  FPR_Best_lambda <- list(glasso = 0, dpglasso = 0, QUIC = 0)
  Sparsity_ratio_Best_lambda <- list(glasso = 0, dpglasso = 0, QUIC = 0)
  Fnorm_Best_lambda <- list(glasso = 0, dpglasso = 0, QUIC = 0)
  Time_Best_lambda <- list(glasso = 0, dpglasso = 0, QUIC = 0)
  
  #### Run all the Lambda and find the best Lambda ####
  total_cores <- detectCores()
  cluster_number <- ifelse(total_cores < 41, total_cores - 1, 41)
  cl <- makeCluster(cluster_number)
  
  registerDoParallel(cl)
  
  clusterExport(cl, varlist = c("dpglasso.new","Banded","Dense","Random_Sparse",
                                "TPR","FPR","Sparsity_ratio","Fnorm","BIC_function"))
  
  results <- foreach(i = 1:L, .packages = c("glasso", "dpglasso", "QUIC")) %dopar% {
    lambda <- lambda_grid[i]
    
    # Run methods and record timings
    Time_glasso <- system.time({
      solution_glasso <- glasso(S, lambda)$wi
    })
    Time_dpglasso <- system.time({
      solution_dpglasso <- dpglasso.new(S, rho = lambda)$X
    })
    Time_QUIC <- system.time({
      solution_QUIC <- QUIC(S, rho = lambda, msg = 0)$X
    })
    
    # Compute metrics
    list(
      BIC_glasso = BIC_function(solution_glasso, S = S, n = n),
      BIC_dpglasso = BIC_function(solution_dpglasso, S = S, n = n),
      BIC_QUIC = BIC_function(solution_QUIC, S = S, n = n),
      
      TPR_glasso = TPR(solution_glasso, Theta_true),
      TPR_dpglasso = TPR(solution_dpglasso, Theta_true),
      TPR_QUIC = TPR(solution_QUIC, Theta_true),
      
      FPR_glasso = FPR(solution_glasso, Theta_true),
      FPR_dpglasso = FPR(solution_dpglasso, Theta_true),
      FPR_QUIC = FPR(solution_QUIC, Theta_true),
      
      Sparsity_ratio_glasso = Sparsity_ratio(solution_glasso, p),
      Sparsity_ratio_dpglasso = Sparsity_ratio(solution_dpglasso, p),
      Sparsity_ratio_QUIC = Sparsity_ratio(solution_QUIC, p),
      
      Fnorm_glasso = Fnorm(solution_glasso, Theta_true),
      Fnorm_dpglasso = Fnorm(solution_dpglasso, Theta_true),
      Fnorm_QUIC = Fnorm(solution_QUIC, Theta_true),
      
      Time_glasso = Time_glasso[3],
      Time_dpglasso = Time_dpglasso[3],
      Time_QUIC = Time_QUIC[3]
    )
  }
  
  stopCluster(cl)
  
  # Combine results
  for (i in 1:L) {
    res <- results[[i]]
    BIC$glasso[i] <- res$BIC_glasso
    BIC$dpglasso[i] <- res$BIC_dpglasso
    BIC$QUIC[i] <- res$BIC_QUIC
    
    TPR$glasso[i] <- res$TPR_glasso
    TPR$dpglasso[i] <- res$TPR_dpglasso
    TPR$QUIC[i] <- res$TPR_QUIC
    
    FPR$glasso[i] <- res$FPR_glasso
    FPR$dpglasso[i] <- res$FPR_dpglasso
    FPR$QUIC[i] <- res$FPR_QUIC
    
    Sparsity_ratio$glasso[i] <- res$Sparsity_ratio_glasso
    Sparsity_ratio$dpglasso[i] <- res$Sparsity_ratio_dpglasso
    Sparsity_ratio$QUIC[i] <- res$Sparsity_ratio_QUIC
    
    Fnorm$glasso[i] <- res$Fnorm_glasso
    Fnorm$dpglasso[i] <- res$Fnorm_dpglasso
    Fnorm$QUIC[i] <- res$Fnorm_QUIC
    
    Time$glasso[i] <- res$Time_glasso
    Time$dpglasso[i] <- res$Time_dpglasso
    Time$QUIC[i] <- res$Time_QUIC
  }
  
  #### Also Output the result with the best lambda
  
  best.index_glasso <- which.min(BIC$glasso)
  best.index_dpglasso <- which.min(BIC$dpglasso)
  best.index_QUIC <- which.min(BIC$QUIC)
  
  Best_lambda$glasso <- lambda_grid[best.index_glasso]
  Best_lambda$dpglasso <- lambda_grid[best.index_dpglasso]
  Best_lambda$QUIC <- lambda_grid[best.index_QUIC]
  
  TPR_Best_lambda$glasso <- TPR$glasso[best.index_glasso]
  TPR_Best_lambda$dpglasso <- TPR$dpglasso[best.index_dpglasso]
  TPR_Best_lambda$QUIC <- TPR$QUIC[best.index_QUIC]
  
  FPR_Best_lambda$glasso <- FPR$glasso[best.index_glasso]
  FPR_Best_lambda$dpglasso <- FPR$dpglasso[best.index_dpglasso]
  FPR_Best_lambda$QUIC <- FPR$QUIC[best.index_QUIC]
  
  Sparsity_ratio_Best_lambda$glasso <- Sparsity_ratio$glasso[best.index_glasso]
  Sparsity_ratio_Best_lambda$dpglasso <- Sparsity_ratio$dpglasso[best.index_dpglasso]
  Sparsity_ratio_Best_lambda$QUIC <- Sparsity_ratio$QUIC[best.index_QUIC]
  
  Time_Best_lambda$glasso <- Time$glasso[best.index_glasso]
  Time_Best_lambda$dpglasso <- Time$dpglasso[best.index_dpglasso]
  Time_Best_lambda$QUIC <- Time$QUIC[best.index_QUIC]
  
  Fnorm_Best_lambda$glasso <- Fnorm$glasso[best.index_glasso]
  Fnorm_Best_lambda$dpglasso <- Fnorm$dpglasso[best.index_dpglasso]
  Fnorm_Best_lambda$QUIC <- Fnorm$QUIC[best.index_QUIC]
  
  return(list(
    BIC = BIC,
    TPR = TPR,
    FPR = FPR,
    Sparsity_ratio = Sparsity_ratio,
    Fnorm = Fnorm,
    Time = Time,
    Best_lambda = Best_lambda,
    TPR_Best_lambda = TPR_Best_lambda,
    FPR_Best_lambda = FPR_Best_lambda,
    Sparsity_ratio_Best_lambda = Sparsity_ratio_Best_lambda,
    Fnorm_Best_lambda = Fnorm_Best_lambda,
    Time_Best_lambda = Time_Best_lambda
  ))
}

main_function_case3 <- function(n,p,sparsity,lambda_grid=seq(0.1, 2, 0.05),thr=1.0e-5) {
  
  set.seed(250)
  
  Theta_true <- Random_Sparse(p,sparsity,eta = 0.5)
  
  data <- mvrnorm(n = n, mu = numeric(p), Sigma = solve(Theta_true))
  S <- cov(data)
  lambda_grid <- lambda_grid
  L <- length(lambda_grid)
  
  BIC <- list(glasso = numeric(L), dpglasso = numeric(L), QUIC = numeric(L))
  TPR <- list(glasso = numeric(L), dpglasso = numeric(L), QUIC = numeric(L))
  FPR <- list(glasso = numeric(L), dpglasso = numeric(L), QUIC = numeric(L))
  Sparsity_ratio <- list(glasso = numeric(L), dpglasso = numeric(L), QUIC = numeric(L))
  Fnorm <- list(glasso = numeric(L), dpglasso = numeric(L), QUIC = numeric(L))
  Time <- list(glasso = numeric(L), dpglasso = numeric(L), QUIC = numeric(L))
  
  ## Record the result under the best Lambda
  Best_lambda <- list(glasso = 0, dpglasso = 0, QUIC = 0)
  TPR_Best_lambda <- list(glasso = 0, dpglasso = 0, QUIC = 0)
  FPR_Best_lambda <- list(glasso = 0, dpglasso = 0, QUIC = 0)
  Sparsity_ratio_Best_lambda <- list(glasso = 0, dpglasso = 0, QUIC = 0)
  Fnorm_Best_lambda <- list(glasso = 0, dpglasso = 0, QUIC = 0)
  Time_Best_lambda <- list(glasso = 0, dpglasso = 0, QUIC = 0)
  
  #### Run all the Lambda and find the best Lambda ####
  total_cores <- detectCores()
  cluster_number <- ifelse(total_cores < 41, total_cores - 1, 41)
  cl <- makeCluster(cluster_number)
  
  registerDoParallel(cl)
  
  clusterExport(cl, varlist = c("dpglasso.new","Banded","Dense","Random_Sparse",
                                "TPR","FPR","Sparsity_ratio","Fnorm","BIC_function"))
  
  results <- foreach(i = 1:L, .packages = c("glasso", "dpglasso", "QUIC")) %dopar% {
    lambda <- lambda_grid[i]
    
    # Run methods and record timings
    Time_glasso <- system.time({
      solution_glasso <- glasso(S, lambda)$wi
    })
    Time_dpglasso <- system.time({
      solution_dpglasso <- dpglasso.new(S, rho = lambda)$X
    })
    Time_QUIC <- system.time({
      solution_QUIC <- QUIC(S, rho = lambda, msg = 0)$X
    })
    
    # Compute metrics
    list(
      BIC_glasso = BIC_function(solution_glasso, S = S, n = n),
      BIC_dpglasso = BIC_function(solution_dpglasso, S = S, n = n),
      BIC_QUIC = BIC_function(solution_QUIC, S = S, n = n),
      
      TPR_glasso = TPR(solution_glasso, Theta_true),
      TPR_dpglasso = TPR(solution_dpglasso, Theta_true),
      TPR_QUIC = TPR(solution_QUIC, Theta_true),
      
      FPR_glasso = FPR(solution_glasso, Theta_true),
      FPR_dpglasso = FPR(solution_dpglasso, Theta_true),
      FPR_QUIC = FPR(solution_QUIC, Theta_true),
      
      Sparsity_ratio_glasso = Sparsity_ratio(solution_glasso, p),
      Sparsity_ratio_dpglasso = Sparsity_ratio(solution_dpglasso, p),
      Sparsity_ratio_QUIC = Sparsity_ratio(solution_QUIC, p),
      
      Fnorm_glasso = Fnorm(solution_glasso, Theta_true),
      Fnorm_dpglasso = Fnorm(solution_dpglasso, Theta_true),
      Fnorm_QUIC = Fnorm(solution_QUIC, Theta_true),
      
      Time_glasso = Time_glasso[3],
      Time_dpglasso = Time_dpglasso[3],
      Time_QUIC = Time_QUIC[3]
    )
  }
  
  stopCluster(cl)
  
  # Combine results
  for (i in 1:L) {
    res <- results[[i]]
    BIC$glasso[i] <- res$BIC_glasso
    BIC$dpglasso[i] <- res$BIC_dpglasso
    BIC$QUIC[i] <- res$BIC_QUIC
    
    TPR$glasso[i] <- res$TPR_glasso
    TPR$dpglasso[i] <- res$TPR_dpglasso
    TPR$QUIC[i] <- res$TPR_QUIC
    
    FPR$glasso[i] <- res$FPR_glasso
    FPR$dpglasso[i] <- res$FPR_dpglasso
    FPR$QUIC[i] <- res$FPR_QUIC
    
    Sparsity_ratio$glasso[i] <- res$Sparsity_ratio_glasso
    Sparsity_ratio$dpglasso[i] <- res$Sparsity_ratio_dpglasso
    Sparsity_ratio$QUIC[i] <- res$Sparsity_ratio_QUIC
    
    Fnorm$glasso[i] <- res$Fnorm_glasso
    Fnorm$dpglasso[i] <- res$Fnorm_dpglasso
    Fnorm$QUIC[i] <- res$Fnorm_QUIC
    
    Time$glasso[i] <- res$Time_glasso
    Time$dpglasso[i] <- res$Time_dpglasso
    Time$QUIC[i] <- res$Time_QUIC
  }
  
  #### Also Output the result with the best lambda
  
  best.index_glasso <- which.min(BIC$glasso)
  best.index_dpglasso <- which.min(BIC$dpglasso)
  best.index_QUIC <- which.min(BIC$QUIC)
  
  Best_lambda$glasso <- lambda_grid[best.index_glasso]
  Best_lambda$dpglasso <- lambda_grid[best.index_dpglasso]
  Best_lambda$QUIC <- lambda_grid[best.index_QUIC]
  
  TPR_Best_lambda$glasso <- TPR$glasso[best.index_glasso]
  TPR_Best_lambda$dpglasso <- TPR$dpglasso[best.index_dpglasso]
  TPR_Best_lambda$QUIC <- TPR$QUIC[best.index_QUIC]
  
  FPR_Best_lambda$glasso <- FPR$glasso[best.index_glasso]
  FPR_Best_lambda$dpglasso <- FPR$dpglasso[best.index_dpglasso]
  FPR_Best_lambda$QUIC <- FPR$QUIC[best.index_QUIC]
  
  Sparsity_ratio_Best_lambda$glasso <- Sparsity_ratio$glasso[best.index_glasso]
  Sparsity_ratio_Best_lambda$dpglasso <- Sparsity_ratio$dpglasso[best.index_dpglasso]
  Sparsity_ratio_Best_lambda$QUIC <- Sparsity_ratio$QUIC[best.index_QUIC]
  
  Time_Best_lambda$glasso <- Time$glasso[best.index_glasso]
  Time_Best_lambda$dpglasso <- Time$dpglasso[best.index_dpglasso]
  Time_Best_lambda$QUIC <- Time$QUIC[best.index_QUIC]
  
  Fnorm_Best_lambda$glasso <- Fnorm$glasso[best.index_glasso]
  Fnorm_Best_lambda$dpglasso <- Fnorm$dpglasso[best.index_dpglasso]
  Fnorm_Best_lambda$QUIC <- Fnorm$QUIC[best.index_QUIC]
  
  return(list(
    BIC = BIC,
    TPR = TPR,
    FPR = FPR,
    Sparsity_ratio = Sparsity_ratio,
    Fnorm = Fnorm,
    Time = Time,
    Best_lambda = Best_lambda,
    TPR_Best_lambda = TPR_Best_lambda,
    FPR_Best_lambda = FPR_Best_lambda,
    Sparsity_ratio_Best_lambda = Sparsity_ratio_Best_lambda,
    Fnorm_Best_lambda = Fnorm_Best_lambda,
    Time_Best_lambda = Time_Best_lambda
  ))
}





########################### Test Run  ##############################

cat("Test Run:","\n")

main_function_case1(n=100,p=10,diag=1,offdiag=0.5)
main_function_case2(n=100,p=10,diag=1,offdiag=0.5)
main_function_case3(n=100,p=10,sparsity=0.5)


###########################  Run  ##############################

p_range <- c(100,200,500)





cat("Case 3_2, sp=0.99:","\n")
for (k in p_range) {
  cat("n = 100,p =",k,"diag=1,offdiag=0.99","\n")
  name <- paste("Case3_p=", k, "_Sp=0.99", sep = "")
  assign(name, main_function_case3(n=100,p=k,sparsity=0.99))  
}
save.image(file = "Results3_2.RData")



cat("Case 3_1, sp=0.5:","\n")

for (k in p_range) {
  cat("n = 100,p =",k,"diag=1,offdiag=0.5","\n")
  name <- paste("Case3_p=", k, "_Sp=0.5", sep = "")
  assign(name, main_function_case3(n=100,p=k,sparsity=0.5))  
}
save.image(file = "Results3_1.RData")
