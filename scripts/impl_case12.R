# scripts/impl_case12.R
# Case 1 (banded) + Case 2 (dense) simulation and benchmarking

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(QUIC)
  library(dpglasso)
  library(glasso)
  library(MASS)
  library(foreach)
  library(doParallel)
})

dpglasso.new <- function (Sigma, X = NULL, invX = NULL, rho, outer.Maxiter = 100,
                          obj.seq = FALSE, outer.tol = 10^-5)
{
  output <- list()
  if (!is.matrix(Sigma))
    stop("Sigma should be a covariance matrix")
  p <- nrow(Sigma)
  if (ncol(Sigma) != p)
    stop("Sigma should be a covariance matrix")
  if (any((rho < 0)))
    stop("Penalty parameters should be non-negative")
  if (length(rho) == 1) {
    rho <- matrix(rho, p, p)
    diag(rho) <- 0
  }
  if (!isSymmetric(rho))
    stop("Penalty matrix should be symmetric")
  if (nrow(rho) != p | ncol(rho) != p)
    stop("Size of penalty matrix does not match dimension of Sigma")
  if (is.null(X) && is.null(invX)) {
    X <- diag(p)
    invX <- diag(p)
  }
  if (!is.null(X) && is.null(invX))
    invX <- solve(X)
  if (!is.null(invX) && is.null(X))
    X <- solve(invX)
  if (!isSymmetric(X)) {
    warning("Initial inverse covariance matrix is not symmetric. dpglasso uses upper triangular part only")
    X <- X * upper.tri(X) + t(X * upper.tri(X))
  }
  if (!isSymmetric(invX)) {
    warning("Initial covariance matrix is not symmetric. dpglasso uses upper triangular part only")
    invX <- invX * upper.tri(invX) + t(invX * upper.tri(invX))
  }
  iter <- 0
  dgap <- 0
  obj <- 0
  while (iter < outer.Maxiter) {
    iter <- iter + 1
    start.dual <- proc.time()
    fit.dual <- QUIC::QUIC(Sigma = Sigma, rho = rho, msg = 0, maxIter = 100,
                           tol = 1e-5)
    X <- fit.dual$X
    invX <- fit.dual$W
    dual.time <- (proc.time() - start.dual)[3]
    obj <- fit.dual$opt[1]
    dgap <- fit.dual$opt[2]
    if (obj.seq)
      output$obj[iter] <- obj
    if (dgap < outer.tol)
      break
  }
  output$X <- X
  output$W <- invX
  output$dgap <- dgap
  output$obj <- obj
  output$iter <- iter
  output
}

main_function <- function(p, case_id, lambda_grid = seq(0.1, 2, 0.05), n = 200) {
  if (case_id == 1) {
    Theta <- diag(1, p)
    for (i in 2:p) {
      Theta[i, i - 1] <- 0.5
      Theta[i - 1, i] <- 0.5
    }
  } else if (case_id == 2) {
    Theta <- matrix(0.5, p, p)
    diag(Theta) <- 1
  } else {
    stop("case_id must be 1 or 2")
  }
  
  Sigma <- solve(Theta)
  X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  S <- cov(X)
  
  BIC_vec <- TPR_vec <- FPR_vec <- Sparsity_ratio_vec <- Fnorm_vec <- Time_vec <- numeric(length(lambda_grid))
  
  for (k in seq_along(lambda_grid)) {
    lambda <- lambda_grid[k]
    
    start_t <- proc.time()
    fit_glasso <- glasso::glasso(s = S, rho = lambda)
    t_glasso <- (proc.time() - start_t)[3]
    Theta_hat_glasso <- fit_glasso$wi
    
    start_t <- proc.time()
    fit_dp <- dpglasso.new(Sigma = S, rho = lambda)
    t_dp <- (proc.time() - start_t)[3]
    Theta_hat_dp <- fit_dp$X
    
    start_t <- proc.time()
    fit_quic <- QUIC::QUIC(Sigma = S, rho = lambda, msg = 0)
    t_quic <- (proc.time() - start_t)[3]
    Theta_hat_quic <- fit_quic$X
    
    eps0 <- 1e-5
    
    TPR <- function(Theta_hat) {
      idx <- which(Theta != 0 & upper.tri(Theta), arr.ind = TRUE)
      mean(abs(Theta_hat[idx]) > eps0)
    }
    FPR <- function(Theta_hat) {
      idx0 <- which(Theta == 0 & upper.tri(Theta), arr.ind = TRUE)
      if (nrow(idx0) == 0) return(NA_real_)
      mean(abs(Theta_hat[idx0]) > eps0)
    }
    sparsity_ratio <- function(Theta_hat) {
      mean(abs(Theta_hat) <= eps0)
    }
    fnorm <- function(Theta_hat) {
      norm(Theta_hat - Theta, type = "F")
    }
    bic_score <- function(Theta_hat) {
      df <- sum(abs(Theta_hat[upper.tri(Theta_hat, diag = TRUE)]) > eps0)
      val <- -n * determinant(Theta_hat, logarithm = TRUE)$modulus + n * sum(S * Theta_hat) + log(n) * df
      as.numeric(val)
    }
    
    BIC_vec[k] <- bic_score(Theta_hat_glasso)
    TPR_vec[k] <- TPR(Theta_hat_glasso)
    FPR_vec[k] <- FPR(Theta_hat_glasso)
    Sparsity_ratio_vec[k] <- sparsity_ratio(Theta_hat_glasso)
    Fnorm_vec[k] <- fnorm(Theta_hat_glasso)
    Time_vec[k] <- t_glasso
    
    # Store full method results for plotting compatibility
    if (k == 1) {
      out <- list(BIC = list(), TPR = list(), FPR = list(), Sparsity_ratio = list(), Fnorm = list(), Time = list())
    }
    out$BIC$glasso[k] <- bic_score(Theta_hat_glasso)
    out$BIC$dpglasso[k] <- bic_score(Theta_hat_dp)
    out$BIC$QUIC[k] <- bic_score(Theta_hat_quic)
    
    out$TPR$glasso[k] <- TPR(Theta_hat_glasso)
    out$TPR$dpglasso[k] <- TPR(Theta_hat_dp)
    out$TPR$QUIC[k] <- TPR(Theta_hat_quic)
    
    out$FPR$glasso[k] <- FPR(Theta_hat_glasso)
    out$FPR$dpglasso[k] <- FPR(Theta_hat_dp)
    out$FPR$QUIC[k] <- FPR(Theta_hat_quic)
    
    out$Sparsity_ratio$glasso[k] <- sparsity_ratio(Theta_hat_glasso)
    out$Sparsity_ratio$dpglasso[k] <- sparsity_ratio(Theta_hat_dp)
    out$Sparsity_ratio$QUIC[k] <- sparsity_ratio(Theta_hat_quic)
    
    out$Fnorm$glasso[k] <- fnorm(Theta_hat_glasso)
    out$Fnorm$dpglasso[k] <- fnorm(Theta_hat_dp)
    out$Fnorm$QUIC[k] <- fnorm(Theta_hat_quic)
    
    out$Time$glasso[k] <- t_glasso
    out$Time$dpglasso[k] <- t_dp
    out$Time$QUIC[k] <- t_quic
  }
  
  out
}

Case1_p100 <- main_function(p = 100, case_id = 1)
Case1_p200 <- main_function(p = 200, case_id = 1)
Case1_p500 <- main_function(p = 500, case_id = 1)

Case2_p100 <- main_function(p = 100, case_id = 2)
Case2_p200 <- main_function(p = 200, case_id = 2)
Case2_p500 <- main_function(p = 500, case_id = 2)

save.image(file = file.path("results", "Results1and2.RData"))
