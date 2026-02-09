# scripts/impl_case3.R
# Case 3 (random sparse) simulation and benchmarking
# Outputs objects:
#   Case3_p=100_Sp=0.5, Case3_p=200_Sp=0.5, Case3_p=500_Sp=0.5 (optional)
#   Case3_p=100_Sp=0.99, Case3_p=200_Sp=0.99, Case3_p=500_Sp=0.99 (optional)
# Writes:
#   results/Results3_0.5.RData
#   results/Results3_0.99.RData

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(QUIC)
  library(glasso)
  library(MASS)
})

dpglasso_quic_once <- function(Sigma, rho) {
  if (!is.matrix(Sigma)) stop("Sigma must be a matrix")
  p <- nrow(Sigma)
  if (ncol(Sigma) != p) stop("Sigma must be square")
  
  if (length(rho) == 1) {
    R <- matrix(rho, p, p)
    diag(R) <- 0
  } else {
    R <- rho
  }
  
  fit <- QUIC::QUIC(Sigma, rho = R, msg = 0, maxIter = 50, tol = 1e-4)
  list(X = fit$X, W = fit$W)
}

make_random_sparse_theta <- function(p, sp, offdiag_scale = 0.15, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # sp: proportion of zeros in the upper triangle (excluding diagonal)
  # edge_prob: proportion of non-zeros
  edge_prob <- max(min(1 - sp, 1), 0)
  
  A <- matrix(0, p, p)
  ut <- which(upper.tri(A), arr.ind = TRUE)
  
  mask <- runif(nrow(ut)) < edge_prob
  vals <- runif(sum(mask), min = 0.05, max = offdiag_scale)
  sgn <- sample(c(-1, 1), size = sum(mask), replace = TRUE)
  
  idx <- ut[mask, , drop = FALSE]
  A[idx] <- vals * sgn
  A <- A + t(A)
  
  # Make diagonally dominant -> SPD
  diag(A) <- 0
  Theta <- A
  diag(Theta) <- rowSums(abs(A)) + 0.5
  
  Theta
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
  list(
    tpr = tpr,
    fpr = fpr,
    sparsity_ratio = sparsity_ratio,
    fnorm = fnorm,
    bic = bic_score
  )
}

run_case3_one <- function(p, sp, lambda_grid, n = 100, seed = 250) {
  set.seed(seed + p + round(sp * 1000))
  
  Theta <- make_random_sparse_theta(p = p, sp = sp, seed = seed + p)
  Sigma <- solve(Theta)
  X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  S <- cov(X)
  
  stopifnot(is.matrix(S))
  stopifnot(isTRUE(all.equal(S, t(S), tolerance = 1e-10)))
  
  m <- metrics_factory(Theta_true = Theta, S = S, n = n)
  
  out <- list(
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
    fit_dp <- dpglasso_quic_once(Sigma = S, rho = lambda)
    t_dp <- (proc.time() - start_t)[3]
    Theta_hat_dp <- fit_dp$X
    
    start_t <- proc.time()
    fit_quic <- QUIC::QUIC(S, rho = lambda, msg = 0, maxIter = 50, tol = 1e-4)
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

lambda_grid_quick <- seq(0.2, 1.5, length.out = 10)

Case3_p100_Sp0.5 <- run_case3_one(p = 100, sp = 0.5, lambda_grid = lambda_grid_quick, n = 100)
Case3_p200_Sp0.5 <- run_case3_one(p = 200, sp = 0.5, lambda_grid = lambda_grid_quick, n = 100)
Case3_p500_Sp0.5 <- run_case3_one(p = 500, sp = 0.5, lambda_grid = lambda_grid_quick, n = 100)

save.image(file = file.path("results", "Results3_0.5.RData"))

Case3_p100_Sp0.99 <- run_case3_one(p = 100, sp = 0.99, lambda_grid = lambda_grid_quick, n = 100)
Case3_p200_Sp0.99 <- run_case3_one(p = 200, sp = 0.99, lambda_grid = lambda_grid_quick, n = 100)
Case3_p500_Sp0.99 <- run_case3_one(p = 500, sp = 0.99, lambda_grid = lambda_grid_quick, n = 100)

save.image(file = file.path("results", "Results3_0.99.RData"))
