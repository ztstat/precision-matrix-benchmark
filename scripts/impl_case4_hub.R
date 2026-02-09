# scripts/impl_case4_hub.R
# Case 4 (hub network) simulation and benchmarking
# Writes:
#   results/Results_4_hub=0.2.RData

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

make_hub_theta <- function(p, hub_ratio = 0.2, hub_degree_ratio = 0.6,
                           offdiag_scale = 0.12, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  k_hub <- max(1, ceiling(hub_ratio * p))
  hubs <- sample(seq_len(p), size = k_hub, replace = FALSE)
  
  A <- matrix(0, p, p)
  
  # Connect hubs to many nodes
  for (h in hubs) {
    candidates <- setdiff(seq_len(p), h)
    deg <- max(1, ceiling(hub_degree_ratio * length(candidates)))
    nbrs <- sample(candidates, size = deg, replace = FALSE)
    vals <- runif(length(nbrs), min = 0.05, max = offdiag_scale) * sample(c(-1, 1), length(nbrs), TRUE)
    A[h, nbrs] <- vals
    A[nbrs, h] <- vals
  }
  
  # Add a small number of non-hub edges
  ut <- which(upper.tri(A), arr.ind = TRUE)
  extra_prob <- 0.02
  mask <- runif(nrow(ut)) < extra_prob
  idx <- ut[mask, , drop = FALSE]
  if (nrow(idx) > 0) {
    vals <- runif(nrow(idx), min = 0.03, max = 0.08) * sample(c(-1, 1), nrow(idx), TRUE)
    A[idx] <- vals
    A <- A + t(A)
  }
  
  diag(A) <- 0
  
  Theta <- A
  diag(Theta) <- rowSums(abs(A)) + 0.5
  
  list(Theta = Theta, hubs = hubs)
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

make_hub_penalty <- function(p, hubs, lambda, hub_weight = 0.3) {
  W <- matrix(1, p, p)
  W[hubs, ] <- hub_weight
  W[, hubs] <- hub_weight
  diag(W) <- 0
  lambda * W
}

run_case4_one <- function(p, hub_ratio, lambda_grid, n = 100, seed = 250) {
  set.seed(seed + p + round(hub_ratio * 1000))
  
  g <- make_hub_theta(p = p, hub_ratio = hub_ratio, seed = seed + p)
  Theta <- g$Theta
  hubs <- g$hubs
  
  Sigma <- solve(Theta)
  X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  S <- cov(X)
  
  stopifnot(is.matrix(S))
  stopifnot(isTRUE(all.equal(S, t(S), tolerance = 1e-10)))
  
  m <- metrics_factory(Theta_true = Theta, S = S, n = n)
  
  out <- list(
    hubs = hubs,
    BIC = list(glasso = numeric(length(lambda_grid)),
               dpglasso = numeric(length(lambda_grid)),
               QUIC = numeric(length(lambda_grid)),
               hubQUIC = numeric(length(lambda_grid))),
    TPR = list(glasso = numeric(length(lambda_grid)),
               dpglasso = numeric(length(lambda_grid)),
               QUIC = numeric(length(lambda_grid)),
               hubQUIC = numeric(length(lambda_grid))),
    FPR = list(glasso = numeric(length(lambda_grid)),
               dpglasso = numeric(length(lambda_grid)),
               QUIC = numeric(length(lambda_grid)),
               hubQUIC = numeric(length(lambda_grid))),
    Sparsity_ratio = list(glasso = numeric(length(lambda_grid)),
                          dpglasso = numeric(length(lambda_grid)),
                          QUIC = numeric(length(lambda_grid)),
                          hubQUIC = numeric(length(lambda_grid))),
    Fnorm = list(glasso = numeric(length(lambda_grid)),
                 dpglasso = numeric(length(lambda_grid)),
                 QUIC = numeric(length(lambda_grid)),
                 hubQUIC = numeric(length(lambda_grid))),
    Time = list(glasso = numeric(length(lambda_grid)),
                dpglasso = numeric(length(lambda_grid)),
                QUIC = numeric(length(lambda_grid)),
                hubQUIC = numeric(length(lambda_grid)))
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
    
    start_t <- proc.time()
    R_hub <- make_hub_penalty(p = p, hubs = hubs, lambda = lambda, hub_weight = 0.3)
    fit_hub <- QUIC::QUIC(S, rho = R_hub, msg = 0, maxIter = 50, tol = 1e-4)
    t_hub <- (proc.time() - start_t)[3]
    Theta_hat_hub <- fit_hub$X
    
    out$BIC$glasso[k] <- m$bic(Theta_hat_glasso)
    out$BIC$dpglasso[k] <- m$bic(Theta_hat_dp)
    out$BIC$QUIC[k] <- m$bic(Theta_hat_quic)
    out$BIC$hubQUIC[k] <- m$bic(Theta_hat_hub)
    
    out$TPR$glasso[k] <- m$tpr(Theta_hat_glasso)
    out$TPR$dpglasso[k] <- m$tpr(Theta_hat_dp)
    out$TPR$QUIC[k] <- m$tpr(Theta_hat_quic)
    out$TPR$hubQUIC[k] <- m$tpr(Theta_hat_hub)
    
    out$FPR$glasso[k] <- m$fpr(Theta_hat_glasso)
    out$FPR$dpglasso[k] <- m$fpr(Theta_hat_dp)
    out$FPR$QUIC[k] <- m$fpr(Theta_hat_quic)
    out$FPR$hubQUIC[k] <- m$fpr(Theta_hat_hub)
    
    out$Sparsity_ratio$glasso[k] <- m$sparsity_ratio(Theta_hat_glasso)
    out$Sparsity_ratio$dpglasso[k] <- m$sparsity_ratio(Theta_hat_dp)
    out$Sparsity_ratio$QUIC[k] <- m$sparsity_ratio(Theta_hat_quic)
    out$Sparsity_ratio$hubQUIC[k] <- m$sparsity_ratio(Theta_hat_hub)
    
    out$Fnorm$glasso[k] <- m$fnorm(Theta_hat_glasso)
    out$Fnorm$dpglasso[k] <- m$fnorm(Theta_hat_dp)
    out$Fnorm$QUIC[k] <- m$fnorm(Theta_hat_quic)
    out$Fnorm$hubQUIC[k] <- m$fnorm(Theta_hat_hub)
    
    out$Time$glasso[k] <- t_glasso
    out$Time$dpglasso[k] <- t_dp
    out$Time$QUIC[k] <- t_quic
    out$Time$hubQUIC[k] <- t_hub
  }
  
  out
}

lambda_grid_quick <- seq(0.2, 1.5, length.out = 10)

Case4_p100_hub0.2 <- run_case4_one(p = 100, hub_ratio = 0.2, lambda_grid = lambda_grid_quick, n = 100)
Case4_p200_hub0.2 <- run_case4_one(p = 200, hub_ratio = 0.2, lambda_grid = lambda_grid_quick, n = 100)
Case4_p500_hub0.2 <- run_case4_one(p = 500, hub_ratio = 0.2, lambda_grid = lambda_grid_quick, n = 100)

save.image(file = file.path("results", "Results_4_hub=0.2.RData"))
