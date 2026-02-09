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


##################  Setting 4: Hub ########## ##################

Hub <- function(p, prob1=0.02, prob2=0.7, hubset) {
  A <- matrix(0, nrow = p, ncol = p)
  
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      A[i, j] <- rbinom(1, 1, prob1)
    }
  }

  A <- A + t(A)
  
  hubs <- hubset
  
  for (hub in hubs) {
    A[hub, ] <- rbinom(p, 1, prob2)
    A[, hub] <- A[hub, ]
  }
  
  repeat {
    A <- A + 0.5 * diag(p)
    if (min(eigen(A)$values) > 1) {
      break  
    }
  }
  return(A)
}




#################################### Criteria


########### Hub_accuracy_edge ###########

Hub_accuracy_edge <- function(estimator, true, hubset) {
  
  p <- nrow(estimator)
  
  numerator <- 0
  denominator <- 0
  
  for (j in hubset) {
    for (j_prime in 1:p) {
      if (j != j_prime) {
        
        if (abs(estimator[j, j_prime]) > 1e-6 && abs(true[j, j_prime]) != 0) {
          numerator <- numerator + 1
        }
        
        if (abs(true[j, j_prime]) != 0) {
          denominator <- denominator + 1
        }
      }
    }
  }
  
  accuracy <- numerator / denominator
  return(accuracy)
}

########### Hub_accuracy_node ###########

Hub_accuracy_node <- function(estimator, true_hubset, threshold = 0.2) {
  
  p <- nrow(estimator)
  
  estimated_hubset <- c()
  
  for (i in 1:p) {
    non_zero_ratio <- sum(abs(estimator[i, ]) > 1e-5) / p
    if (non_zero_ratio > threshold) {
      estimated_hubset <- c(estimated_hubset, i)
    }
  }
  
  correctly_estimated_hubset <- intersect(estimated_hubset, true_hubset)
  
  accuracy <- length(correctly_estimated_hubset) / length(true_hubset)
  
  return(list(
    accuracy_node = accuracy,
    card_estimated_hubset = length(estimated_hubset))
  )
}

############# BIC #############

BIC_function <- function(esti,S,n,threshold = 1e-5) {
  term1 <- -n*log(det(esti))
  term2 <- n*(sum(diag( S %*% esti )))
  term3 <- log(n) * sum(abs(esti)>=threshold)
  BIC <- term1 + term2 + term3
  return(BIC)
}


########### Hub_BIC  ##########

Hub_BIC_function <- function(esti,Z,V,S,n,c=0.2,v,threshold = 1e-5) {
  term1 <- -n*log(det(esti))
  term2 <- n*(sum(diag( S %*% esti )))
  term3 <- log(n) * sum(abs(esti)>=threshold)
  term4 <- log(n) * ( v + c * (sum(abs(V) > threshold)-v) )
  
  Hub_BIC <- term1 + term2 + term3 +term4
  return(Hub_BIC)
}



########################### Simulation  ###################################

main_function_case4 <- function(n,p,hubset,lambda_grid=seq(0.1, 2, 0.05),thr=1.0e-5) {
  
  set.seed(250)
  
  Theta_true <- Hub(p, prob1=0.02, prob2=0.7, hubset)
  data <- mvrnorm(n = n, mu = numeric(p), Sigma = solve(Theta_true))
  S <- cov(data)
  lambda_grid <- lambda_grid
  L <- length(lambda_grid)
  
  
  BIC <- list(glasso = numeric(L), dpglasso = numeric(L), QUIC = numeric(L))
  Time <- list(glasso = numeric(L), dpglasso = numeric(L), QUIC = numeric(L))
  Hub_accuracy_edge <- list(glasso = numeric(L), dpglasso = numeric(L), QUIC = numeric(L))
  Hub_accuracy_node <- list(glasso = numeric(L), dpglasso = numeric(L), QUIC = numeric(L))
  
  
  ## Record the result under the best Lambda
  Best_lambda <- list(glasso = 0, dpglasso = 0, QUIC = 0)
  Time_Best_lambda <- list(glasso = 0, dpglasso = 0, QUIC = 0)
  Hub_accuracy_edge_Best_lambda <- list(glasso = 0, dpglasso = 0, QUIC = 0)
  Hub_accuracy_node_Best_lambda <- list(glasso = 0, dpglasso = 0, QUIC = 0)
  
  #### Run all the Lambda and find the best Lambda ####
  total_cores <- detectCores()
  cluster_number <- ifelse(total_cores < 65, total_cores - 1, 65)
  cl <- makeCluster(cluster_number)
  
  registerDoParallel(cl)
  
  clusterExport(cl, varlist = c("dpglasso.new","Hub","BIC_function",
                                "Hub_accuracy_edge","Hub_accuracy_node"))
  
  results <- foreach(i = 1:L, .packages = c("glasso", "dpglasso", "QUIC", "hglasso")) %dopar% {
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
      
      Hub_accuracy_edge_glasso =  Hub_accuracy_edge(solution_glasso, Theta_true, hubset),
      Hub_accuracy_edge_dpglasso =  Hub_accuracy_edge(solution_dpglasso, Theta_true, hubset),
      Hub_accuracy_edge_QUIC =  Hub_accuracy_edge(solution_QUIC, Theta_true, hubset),
      
      Hub_accuracy_node_glasso =  Hub_accuracy_node(solution_glasso, hubset)$accuracy_node,
      Hub_accuracy_node_dpglasso =  Hub_accuracy_node(solution_dpglasso, hubset)$accuracy_node,
      Hub_accuracy_node_QUIC =  Hub_accuracy_node(solution_QUIC, hubset)$accuracy_node,
      
      
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
    
    Hub_accuracy_edge$glasso[i] <- res$Hub_accuracy_edge_glasso
    Hub_accuracy_edge$dpglasso[i] <- res$Hub_accuracy_edge_dpglasso
    Hub_accuracy_edge$QUIC[i] <- res$Hub_accuracy_edge_QUIC
    
    Hub_accuracy_node$glasso[i] <- res$Hub_accuracy_node_glasso
    Hub_accuracy_node$dpglasso[i] <- res$Hub_accuracy_node_dpglasso
    Hub_accuracy_node$QUIC[i] <- res$Hub_accuracy_node_QUIC
    
    BIC$glasso[i] <- res$BIC_glasso
    BIC$dpglasso[i] <- res$BIC_dpglasso
    BIC$QUIC[i] <- res$BIC_QUIC
    
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
  
  Hub_accuracy_edge_Best_lambda$glasso <- Hub_accuracy_edge$glasso[best.index_glasso]
  Hub_accuracy_edge_Best_lambda$dpglasso <- Hub_accuracy_edge$dpglasso[best.index_dpglasso]
  Hub_accuracy_edge_Best_lambda$QUIC <- Hub_accuracy_edge$QUIC[best.index_QUIC]
  
  Hub_accuracy_node_Best_lambda$glasso <- Hub_accuracy_node$glasso[best.index_glasso]
  Hub_accuracy_node_Best_lambda$dpglasso <- Hub_accuracy_node$dpglasso[best.index_dpglasso]
  Hub_accuracy_node_Best_lambda$QUIC <- Hub_accuracy_node$QUIC[best.index_QUIC]
  
  
  Time_Best_lambda$glasso <- Time$glasso[best.index_glasso]
  Time_Best_lambda$dpglasso <- Time$dpglasso[best.index_dpglasso]
  Time_Best_lambda$QUIC <- Time$QUIC[best.index_QUIC]
  
  return(list(
    BIC = BIC,
    Time = Time,
    Hub_accuracy_edge=Hub_accuracy_edge,
    Hub_accuracy_node=Hub_accuracy_node,
    Best_lambda = Best_lambda,
    Hub_accuracy_edge_Best_lambda=Hub_accuracy_edge_Best_lambda,
    Hub_accuracy_node_Best_lambda=Hub_accuracy_node_Best_lambda,
    Time_Best_lambda = Time_Best_lambda
  ))
}




main_function_case_hglasso <- function(n,p,hubset,
                                lambda1_grid,
                                lambda2_grid,
                                lambda3_grid,
                                thr=1.0e-5) {
  
  set.seed(250)
  
  Theta_true <- Hub(p, prob1=0.02, prob2=0.7, hubset)
  
  data <- mvrnorm(n = n, mu = numeric(p), Sigma = solve(Theta_true))
  S <- cov(data)

  L1 <- length(lambda1_grid)
  L2 <- length(lambda2_grid)
  L3 <- length(lambda3_grid)
  
  BIC <- numeric(L1*L2*L3)
  Time <- numeric(L1*L2*L3)
  Hub_accuracy_edge <- numeric(L1*L2*L3)
  Hub_accuracy_node <- numeric(L1*L2*L3)
  
  np.para.matrix <- as.matrix(expand.grid(lambda1_grid, 
                                          lambda2_grid,
                                          lambda3_grid))
  
  
  ## Record the result under the best Lambda
  Best_lambda <- 0
  Time_Best_lambda <- 0
  Hub_accuracy_edge_Best_lambda <- 0
  Hub_accuracy_node_Best_lambda <- 0
  
  #### Run all the Lambda and find the best Lambda ####
  total_cores <- detectCores()
  cluster_number <- ifelse(total_cores < 65, total_cores - 1, 65)
  cl <- makeCluster(cluster_number)
  
  registerDoParallel(cl)
  
  clusterExport(cl, varlist = c("dpglasso.new","Hub","Hub_BIC_function",
                                "Hub_accuracy_edge","Hub_accuracy_node"
                                ))

 
  
   
  results <- foreach(i = 1:nrow(np.para.matrix), .packages = c("glasso", "dpglasso", "QUIC", "hglasso")) %dopar% {
    
    
    # Run methods and record timings
    Time_hglasso <- system.time({
      solution_hglasso_Theta <- hglasso(S, lambda1=np.para.matrix[i,1], 
                                        lambda2=np.para.matrix[i,2], 
                                        lambda3=np.para.matrix[i,3])$Theta
    })
    
    solution_hglasso_V <- hglasso(S, lambda1=np.para.matrix[i,1], 
                                  lambda2=np.para.matrix[i,2], 
                                  lambda3=np.para.matrix[i,3])$V
    solution_hglasso_Z <- hglasso(S, lambda1=np.para.matrix[i,1], 
                                  lambda2=np.para.matrix[i,2], 
                                  lambda3=np.para.matrix[i,3])$Z
    
    # Compute metrics
    list(
      BIC_hub = Hub_BIC_function(esti=solution_hglasso_Theta,
                                 Z=solution_hglasso_Z,
                                 V=solution_hglasso_V,S,n,c=0.2,
                                 v=Hub_accuracy_node(solution_hglasso_Theta, hubset)$card_estimated_hubset,
                                 threshold = 1e-5),
      Hub_accuracy_edge_hglasso = Hub_accuracy_edge(solution_hglasso_Theta, Theta_true, hubset),
     
      Hub_accuracy_node_hglasso =  Hub_accuracy_node(solution_hglasso_Theta, hubset)$accuracy_node,
   
      Time_hglasso = Time_hglasso[3]

    )
  }
  
  stopCluster(cl)
  
  # Combine results
  for (i in 1:nrow(np.para.matrix)) {
    res <- results[[i]]
    BIC[i] <- res$BIC_hub
    Hub_accuracy_edge[i] <- res$Hub_accuracy_edge_hglasso
    Hub_accuracy_node[i] <- res$Hub_accuracy_node_hglasso
    Time[i] <- res$Time_hglasso
  }
  
  #### Also Output the result with the best lambda
  
  best.index <- which.min(BIC)
  
  
  Best_lambda <- np.para.matrix[best.index,]
  
  
  Time_Best <- Time[best.index]
  BIC_Best <- BIC[best.index]
  Hub_accuracy_edge_Best_lambda <- Hub_accuracy_edge[best.index]
  Hub_accuracy_node_Best_lambda <- Hub_accuracy_node[best.index]
  return(list(
    BIC = BIC,
    Time = Time,
    Hub_accuracy_edge=Hub_accuracy_edge,
    Hub_accuracy_node=Hub_accuracy_node,
    
    Best_lambda = Best_lambda,
    
    BIC_Best = BIC_Best,
    Hub_accuracy_edge_Best_lambda = Hub_accuracy_edge_Best_lambda,
    Hub_accuracy_node_Best_lambda = Hub_accuracy_node_Best_lambda,
    Time_Best = Time_Best
  ))
}




########################### Test Run  ##############################

cat("Test Run:","\n")


n <- 100
p <- 10
set <- sample(1:p,2)

lambda1_grid=c(0.1,0.5)
lambda2_grid=c(0.1,0.5)
lambda3_grid=c(0.1,0.5)
np.para.matrix <- as.matrix(expand.grid(lambda1_grid, 
                                        lambda2_grid,
                                        lambda3_grid))

main_function_case4(n,p,hubset=set,lambda_grid=seq(0.1, 0.3, 0.05),thr=1.0e-5)
main_function_case_hglasso(n,p,hubset=set,
                           lambda1_grid,
                           lambda2_grid,
                           lambda3_grid,
                           thr=1.0e-5)


###########################  Run  ##############################


n <- 100
p_range <- c(10, 20, 50, 100)

lambda1_grid <- seq(0.02, 0.3, 0.02)
lambda2_grid <- seq(0.02, 0.3, 0.02)
lambda3_grid <- seq(0.02, 0.3, 0.02)

set.seed(250)

results_glasso <- list()
results_hglasso <- list()

cat("Running simulations for Case 4_glasso and Case 4_hglasso:\n")
for (k in p_range) {

  set.seed(250)
  
  set <- sample(1:k, ceiling(0.2 * k))
  cat("n = 100, p =", k, "\n")
  
  results_glasso[[paste0("Case4_p=", k)]] <- main_function_case4(n, k, hubset=set, 
                                                                 lambda_grid=seq(0.1, 2, 0.05), 
                                                                 thr=1.0e-5)
  
  results_hglasso[[paste0("Case4_p=", k)]] <- main_function_case_hglasso(n, k, hubset=set,
                                                                         lambda1_grid,
                                                                         lambda2_grid,
                                                                         lambda3_grid,
                                                                         thr=1.0e-5)
}

save(results_glasso, results_hglasso, file = "Results_4.RData")


  

  




