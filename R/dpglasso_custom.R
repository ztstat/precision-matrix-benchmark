# R/dpglasso_custom.R
# Custom dpglasso implementation with defensive numerical guards

options(stringsAsFactors = FALSE)

box_qp_f <- function(Q, u, b, rho, Maxiter = 1000, tol = 1e-7) {
  if (!is.matrix(Q)) stop("Q must be a matrix")
  m <- nrow(Q)
  if (ncol(Q) != m) stop("Q must be square")
  if (length(b) != m) stop("b dimension mismatch")
  if (length(u) != m) stop("u dimension mismatch")
  if (!is.numeric(rho) || length(rho) != 1 || rho < 0) stop("rho must be a non-negative scalar")
  
  if (any(!is.finite(Q)) || any(!is.finite(u)) || any(!is.finite(b))) {
    stop("Non-finite inputs in Q/u/b")
  }
  
  Qs <- 0.5 * (Q + t(Q))
  
  ridge <- 1e-8 * (mean(diag(Qs)) + 1)
  if (!is.finite(ridge) || ridge <= 0) ridge <- 1e-6
  
  tries <- 0
  repeat {
    ok <- TRUE
    tryCatch(chol(Qs + diag(ridge, m)), error = function(e) { ok <<- FALSE })
    if (ok) break
    ridge <- ridge * 10
    tries <- tries + 1
    if (tries > 6) break
  }
  Qs <- Qs + diag(ridge, m)
  
  u0 <- pmax(pmin(u, rho), -rho)
  
  fn <- function(x) {
    v <- as.numeric(0.5 * crossprod(x, Qs %*% x) + crossprod(b, x))
    if (!is.finite(v)) 1e100 else v
  }
  
  gr <- function(x) {
    g <- as.numeric(Qs %*% x + b)
    g[!is.finite(g)] <- 0
    g
  }
  
  res <- optim(
    par = u0,
    fn = fn,
    gr = gr,
    method = "L-BFGS-B",
    lower = rep(-rho, m),
    upper = rep(rho, m),
    control = list(maxit = Maxiter, factr = tol / .Machine$double.eps)
  )
  
  u_hat <- res$par
  grad_vec <- as.numeric(Qs %*% u_hat + b)
  
  list(u = u_hat, grad_vec = grad_vec, value = res$value, conv = res$convergence)
}

dpglasso.new <- function(Sigma, X = NULL, invX = NULL, rho,
                         outer.Maxiter = 100, obj.seq = FALSE, outer.tol = 1e-5) {
  if (is.null(Sigma)) stop("Sigma is required")
  if (!is.matrix(Sigma)) stop("Sigma must be a matrix")
  
  p <- nrow(Sigma)
  if (ncol(Sigma) != p) stop("Sigma must be square")
  if (!isTRUE(all.equal(Sigma, t(Sigma), tolerance = 1e-10))) stop("Sigma must be symmetric")
  if (!is.numeric(rho) || length(rho) != 1 || rho < 0) stop("rho must be a non-negative scalar")
  
  Sigma <- 0.5 * (Sigma + t(Sigma)) + diag(1e-8, p)
  
  if (is.null(X)) {
    d <- diag(Sigma)
    d[!is.finite(d)] <- 1
    X <- diag(1 / pmax(d + rho, 1e-8), p)
  }
  if (is.null(invX)) {
    U.mat <- diag(rep(rho, p))
    invX <- Sigma + U.mat
  } else {
    if (!is.matrix(invX) || nrow(invX) != p || ncol(invX) != p) stop("invX dimension mismatch")
    U.mat <- invX - Sigma
  }
  
  if (!is.matrix(X) || nrow(X) != p || ncol(X) != p) stop("X dimension mismatch")
  
  if (any(!is.finite(X)) || any(!is.finite(invX)) || any(!is.finite(U.mat))) {
    X <- diag(1, p)
    invX <- Sigma + diag(rho, p)
    U.mat <- invX - Sigma
  }
  
  ids <- rep(seq_len(p), outer.Maxiter)
  rel.err <- numeric(outer.Maxiter)
  obj.vals <- numeric(outer.Maxiter)
  sparse.nos <- numeric(outer.Maxiter)
  
  ii <- 0
  kk <- 0
  tol_now <- Inf
  X_prev <- X
  
  eps0 <- 1e-8
  
  for (cd.iter in ids) {
    ii <- ii + 1
    
    Q <- X[-cd.iter, -cd.iter, drop = FALSE]
    u <- U.mat[-cd.iter, cd.iter]
    b <- Sigma[-cd.iter, cd.iter]
    
    if (any(!is.finite(Q)) || any(!is.finite(u)) || any(!is.finite(b))) {
      X[cd.iter, ] <- 0
      X[, cd.iter] <- 0
      X[cd.iter, cd.iter] <- 1 / pmax(Sigma[cd.iter, cd.iter] + rho, eps0)
      U.mat[-cd.iter, cd.iter] <- 0
      U.mat[cd.iter, -cd.iter] <- 0
      U.mat[cd.iter, cd.iter] <- rho
      next
    }
    
    obj <- box_qp_f(
      Q = Q,
      u = u,
      b = b,
      rho = rho,
      Maxiter = 1000,
      tol = 1e-7
    )
    
    denom <- pmax(Sigma[cd.iter, cd.iter] + rho, eps0)
    theta.hat <- -as.numeric(obj$grad_vec) * 0.5 / denom
    
    theta.hat[!is.finite(theta.hat)] <- 0
    theta.hat[abs(obj$u) < rho * 0.99999] <- 0
    theta.hat[abs(theta.hat) < 1e-6] <- 0
    
    X[cd.iter, -cd.iter] <- theta.hat
    X[-cd.iter, cd.iter] <- theta.hat
    
    diag_term <- sum((as.numeric(obj$u) + Sigma[cd.iter, -cd.iter]) * X[cd.iter, -cd.iter])
    diag_term <- if (is.finite(diag_term)) diag_term else 0
    Xii <- (1 - diag_term) / denom
    if (!is.finite(Xii)) Xii <- 1 / denom
    X[cd.iter, cd.iter] <- pmax(Xii, eps0)
    
    U.mat[-cd.iter, cd.iter] <- as.numeric(obj$u)
    U.mat[cd.iter, -cd.iter] <- U.mat[-cd.iter, cd.iter]
    U.mat[cd.iter, cd.iter] <- rho
    
    if (any(!is.finite(X))) {
      X[!is.finite(X)] <- 0
      diag(X) <- pmax(diag(X), eps0)
    }
    if (any(!is.finite(U.mat))) {
      U.mat[!is.finite(U.mat)] <- 0
      diag(U.mat) <- rho
    }
    
    if (cd.iter == p) {
      kk <- kk + 1
      num <- max(abs(X - X_prev))
      den <- pmax(max(abs(X_prev)), eps0)
      tol_now <- num / den
      
      rel.err[kk] <- tol_now
      sparse.nos[kk] <- sum(abs(X) <= 1e-9)
      
      if (isTRUE(obj.seq)) {
        v <- -sum(log(abs(diag(qr(X)$qr)))) + sum(Sigma * X) + rho * sum(abs(X))
        obj.vals[kk] <- as.numeric(v)
      }
      
      X_prev <- X
      
      if ((tol_now < outer.tol) && (kk > 1)) break
    }
  }
  
  invX_out <- Sigma + U.mat
  invX_out <- 0.5 * (invX_out + t(invX_out))
  
  if (isTRUE(obj.seq)) {
    list(
      rel.err = rel.err[seq_len(max(kk, 1))],
      obj.vals = obj.vals[seq_len(max(kk, 1))],
      X = X,
      invX = invX_out,
      sparse.nos = sparse.nos[seq_len(max(kk, 1))],
      iter = kk
    )
  } else {
    list(
      rel.err = rel.err[seq_len(max(kk, 1))],
      X = X,
      invX = invX_out,
      sparse.nos = sparse.nos[seq_len(max(kk, 1))],
      iter = kk
    )
  }
}
