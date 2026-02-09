# scripts/impl_case12.R
# Case 1 (banded) + Case 2 (dense): simulation and benchmarking

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(QUIC)
  library(glasso)
  library(MASS)
})

source(file.path("R", "dpglasso_custom.R"))

make_theta_case12 <- function(p, case_id) {
  if (case_id == 1) {
    Theta <- diag(1, p)
    for (i in 2:p) {
      Theta[i, i - 1] <- 0.5
      Theta[i - 1, i] <- 0.5
    }
    return(Theta)
  }
  
  if (case_id == 2) {
    Theta <- matrix(0.5, p, p)
    diag(Theta) <- 1
    return(Theta)
  }
  
  stop("case_id must be 1 or 2")
}

metrics_factory <- function(Theta_true, S, n, eps0 = 1e-5) {
  tpr <- function(Theta_hat) {
    idx <- which(Theta_true != 0 & upper.tri(Theta_true), arr.ind = TRUE)
    mean(abs(Theta_hat[idx]) > eps0)
  }
  
  fpr <- function(Theta_hat) {
    idx0 <- which(Theta_true == 0 & upper.tri(Theta_true), arr.ind = TRUE)
    if (nrow(idx0) == 0) return(NA_real_)
    mean(abs(Theta_hat[idx0]) > eps0)
  }
  
  sparsity_ratio <- function(Theta_hat) {
    mean(abs(Theta_hat) <= eps0)
  }
  
  fnorm <- function(Theta_hat) {
    norm(Theta_hat - Theta_true, type = "F")
  }
  
  bic_score <- function(Theta_hat) {
    df <- sum(abs(Theta_hat[upper.tri(Theta_hat, diag = TRUE)]) > eps0)
    val <- -n * determinant(Theta_hat, logarithm = TRUE)$modulus +
      n * sum(S * Theta_hat) +
      log(n) * df
    as.numeric(val)
  }
  
  list(tpr = tpr, fpr = fpr, sparsity_ratio = sparsity_ratio, fnorm = fnorm, bic = bic_score)
}

run_case12_one <- function(p, case_id, lambda_grid, n = 100, seed = 250,
                           dp_outer_maxiter = 50, dp_outer_tol = 1e-5,
                           quic_maxiter = 50, quic_tol = 1e-4) {
  set.seed(seed + 1000L * case_id + p)
  
  Theta <- make_theta_case12(p = p, case_id = case_id)
  Sigma <- solve(Theta)
  
  X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  S <- cov(X)
  
  stopifnot(is.matrix(S))
  stopifnot(isTRUE(all.equal(S, t(S), tolerance = 1e-10)))
  
  m <- metrics_factory(Theta_true = Theta, S = S, n = n)
  
  out <- list(
    lambda = lambda_grid,
    BIC = list(glasso = numeric(length(lambda_grid)),
               dpglasso = numeric(length(lambda_grid)),
               QUIC = numeric(length(lambda_grid))),
    TPR = list(glasso = numeric(length(lambda_grid)),
               dpglasso = numeric(length(lambda_grid)),
               QUIC = numeric(length(lambda_grid))),
    FPR = list(glasso = numeric(length(lambda_grid)),
               dpglasso = numeric(length(lambda_grid)),
               QUIC = numeric(length(lambda_grid))),
    Sparsity_ratio = list(glasso = numeric(length(lambda_grid)),
                          dpglasso = numeric(length(lambda_grid)),
                          QUIC = numeric(length(lambda_grid))),
    Fnorm = list(glasso = numeric(length(lambda_grid)),
                 dpglasso = numeric(length(lambda_grid)),
                 QUIC = numeric(length(lambda_grid))),
    Time = list(glasso = numeric(length(lambda_grid)),
                dpglasso = numeric(length(lambda_grid)),
                QUIC = numeric(length(lambda_grid)))
  )
  
  for (k in seq_along(lambda_grid)) {
    lambda <- lambda_grid[k]
    
    start_t <- proc.time()
    fit_glasso <- glasso::glasso(s = S, rho = lambda)
    t_glasso <- (proc.time() - start_t)[3]
    Theta_hat_glasso <- fit_glasso$wi
    
    start_t <- proc.time()
    fit_dp <- dpglasso.new(Sigma = S, rho = lambda, outer.Maxiter = dp_outer_maxiter, outer.tol = dp_outer_tol)
    t_dp <- (proc.time() - start_t)[3]
    Theta_hat_dp <- fit_dp$X
    
    start_t <- proc.time()
    fit_quic <- QUIC::QUIC(S, rho = lambda, msg = 0, maxIter = quic_maxiter, tol = quic_tol)
    t_quic <- (proc.time() - start_t)[3]
    Theta_hat_quic <- fit_quic$X
    
    out$BIC$glasso[k] <- m$bic(Theta_hat_glasso)
    out$BIC$dpglasso[k] <- m$bic(Theta_hat_dp)
    out$BIC$QUIC[k] <- m$bic(Theta_hat_quic)
    
    out$TPR$glasso[k] <- m$tpr(Theta_hat_glasso)
    out$TPR$dpglasso[k] <- m$tpr(Theta_hat_dp)
    out$TPR$QUIC[k] <- m$tpr(Theta_hat_quic)
    
    out$FPR$glasso[k] <- m$fpr(Theta_hat_glasso)
    out$FPR$dpglasso[k] <- m$fpr(Theta_hat_dp)
    out$FPR$QUIC[k] <- m$fpr(Theta_hat_quic)
    
    out$Sparsity_ratio$glasso[k] <- m$sparsity_ratio(Theta_hat_glasso)
    out$Sparsity_ratio$dpglasso[k] <- m$sparsity_ratio(Theta_hat_dp)
    out$Sparsity_ratio$QUIC[k] <- m$sparsity_ratio(Theta_hat_quic)
    
    out$Fnorm$glasso[k] <- m$fnorm(Theta_hat_glasso)
    out$Fnorm$dpglasso[k] <- m$fnorm(Theta_hat_dp)
    out$Fnorm$QUIC[k] <- m$fnorm(Theta_hat_quic)
    
    out$Time$glasso[k] <- t_glasso
    out$Time$dpglasso[k] <- t_dp
    out$Time$QUIC[k] <- t_quic
  }
  
  out
}

# Quick config
lambda_grid_quick <- c(0.2, 1.0)
n_quick <- 100

Case1_p100 <- run_case12_one(p = 100, case_id = 1, lambda_grid = lambda_grid_quick, n = n_quick)
#Case1_p200 <- run_case12_one(p = 200, case_id = 1, lambda_grid = lambda_grid_quick, n = n_quick)
# Case1_p500 <- run_case12_one(p = 500, case_id = 1, lambda_grid = lambda_grid_quick, n = n_quick)

Case2_p100 <- run_case12_one(p = 100, case_id = 2, lambda_grid = lambda_grid_quick, n = n_quick)
#Case2_p200 <- run_case12_one(p = 200, case_id = 2, lambda_grid = lambda_grid_quick, n = n_quick)
# Case2_p500 <- run_case12_one(p = 500, case_id = 2, lambda_grid = lambda_grid_quick, n = n_quick)

save.image(file = file.path("results", "Results1and2.RData"))
