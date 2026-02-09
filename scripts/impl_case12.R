# scripts/impl_case12.R
# Case 1 (banded) + Case 2 (dense) simulation and benchmarking

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
  
  fit <- QUIC::QUIC(Sigma, rho = R, msg = 0)
  list(X = fit$X, W = fit$W)
}

main_function <- function(p, case_id, lambda_grid = seq(0.1, 1.5, length.out = 2), n = 100) {
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
  
  stopifnot(is.matrix(S))
  stopifnot(nrow(S) == ncol(S))
  stopifnot(is.numeric(S))
  stopifnot(isTRUE(all.equal(S, t(S), tolerance = 1e-10)))
  
  eps0 <- 1e-5
  
  tpr <- function(Theta_hat) {
    idx <- which(Theta != 0 & upper.tri(Theta), arr.ind = TRUE)
    mean(abs(Theta_hat[idx]) > eps0)
  }
  
  fpr <- function(Theta_hat) {
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
    val <- -n * determinant(Theta_hat, logarithm = TRUE)$modulus +
      n * sum(S * Theta_hat) +
      log(n) * df
    as.numeric(val)
  }
  
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
    fit_quic <- QUIC::QUIC(S, rho = lambda, msg = 0)
    t_quic <- (proc.time() - start_t)[3]
    Theta_hat_quic <- fit_quic$X
    
    out$BIC$glasso[k] <- bic_score(Theta_hat_glasso)
    out$BIC$dpglasso[k] <- bic_score(Theta_hat_dp)
    out$BIC$QUIC[k] <- bic_score(Theta_hat_quic)
    
    out$TPR$glasso[k] <- tpr(Theta_hat_glasso)
    out$TPR$dpglasso[k] <- tpr(Theta_hat_dp)
    out$TPR$QUIC[k] <- tpr(Theta_hat_quic)
    
    out$FPR$glasso[k] <- fpr(Theta_hat_glasso)
    out$FPR$dpglasso[k] <- fpr(Theta_hat_dp)
    out$FPR$QUIC[k] <- fpr(Theta_hat_quic)
    
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
#Case1_p200 <- main_function(p = 200, case_id = 1)
#Case1_p500 <- main_function(p = 500, case_id = 1)

Case2_p100 <- main_function(p = 100, case_id = 2)
#Case2_p200 <- main_function(p = 200, case_id = 2)
#Case2_p500 <- main_function(p = 500, case_id = 2)

save.image(file = file.path("results", "Results1and2.RData"))
