# R/dpglasso_custom.R
# Custom dpglasso implementation (standalone)
# - box_qp_f solves a box-constrained quadratic program via L-BFGS-B
# - dpglasso.new performs coordinate descent updates on the precision matrix

options(stringsAsFactors = FALSE)

box_qp_f <- function(Q, u, b, rho, Maxiter = 1000, tol = 1e-7) {
  if (!is.matrix(Q)) stop("Q must be a matrix")
  m <- nrow(Q)
  if (ncol(Q) != m) stop("Q must be square")
  if (length(b) != m) stop("b dimension mismatch")
  if (length(u) != m) stop("u dimension mismatch")
  if (length(rho) != 1 || rho < 0) stop("rho must be a non-negative scalar")
  
  Qs <- 0.5 * (Q + t(Q))
  
  fn <- function(x) {
    as.numeric(0.5 * crossprod(x, Qs %*% x) + crossprod(b, x))
  }
  
  gr <- function(x) {
    as.numeric(Qs %*% x + b)
  }
  
  res <- optim(
    par = u,
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
  if (!isTRUE(all.equal(Sigma, t(Sigma), tolerance = 1e-12))) {
    stop("Sigma must be symmetric")
  }
  if (length(rho) != 1 || rho < 0) stop("rho must be a non-negative scalar")
  
  if (is.null(X)) {
    X <- diag(1 / (rep(rho, p) + diag(Sigma)))
  }
  if (is.null(invX)) {
    U.mat <- diag(rep(rho, p))
    invX <- Sigma + U.mat
  } else {
    U.mat <- invX - Sigma
  }
  
  if (nrow(X) != p || ncol(X) != p) stop("X dimension mismatch")
  if (nrow(invX) != p || ncol(invX) != p) stop("invX dimension mismatch")
  
  ids <- rep(seq_len(p), outer.Maxiter)
  time.counter.QP <- array(0, dim = c(length(ids), 3))
  rel.err <- numeric(outer.Maxiter)
  obj.vals <- numeric(outer.Maxiter)
  sparse.nos <- numeric(outer.Maxiter)
  
  ii <- 0
  kk <- 0
  tol <- 10
  vec.diagsold <- X
  
  for (cd.iter in ids) {
    ii <- ii + 1
    
    t0 <- proc.time()
    obj <- box_qp_f(
      Q = X[-cd.iter, -cd.iter, drop = FALSE],
      u = U.mat[-cd.iter, cd.iter],
      b = Sigma[-cd.iter, cd.iter],
      rho = rho,
      Maxiter = 1000,
      tol = 1e-7
    )
    dt <- proc.time() - t0
    time.counter.QP[ii, ] <- as.numeric(c(dt[1], dt[2], dt[3]))
    
    theta.hat <- -obj$grad_vec * 0.5 / (Sigma[cd.iter, cd.iter] + rho)
    
    theta.hat[abs(obj$u) < rho * 0.99999] <- 0
    theta.hat[abs(theta.hat) < 1e-5] <- 0
    
    X[cd.iter, -cd.iter] <- theta.hat
    X[-cd.iter, cd.iter] <- X[cd.iter, -cd.iter]
    
    X[cd.iter, cd.iter] <- sum((obj$u + Sigma[cd.iter, -cd.iter]) * X[cd.iter, -cd.iter])
    X[cd.iter, cd.iter] <- (1 - X[cd.iter, cd.iter]) / (Sigma[cd.iter, cd.iter] + rho)
    
    U.mat[-cd.iter, cd.iter] <- as.numeric(obj$u)
    U.mat[cd.iter, -cd.iter] <- U.mat[-cd.iter, cd.iter]
    U.mat[cd.iter, cd.iter] <- rho
    
    if (cd.iter == p) {
      kk <- kk + 1
      vec.diagsnew <- X
      tol <- max(abs(vec.diagsnew - vec.diagsold)) / max(abs(vec.diagsold))
      rel.err[kk] <- tol
      vec.diagsold <- vec.diagsnew
      sparse.nos[kk] <- sum(abs(X) <= 1e-9)
      
      if (isTRUE(obj.seq)) {
        obj.vals[kk] <- -sum(log(abs(diag(qr(X)$qr)))) + sum(Sigma * X) + rho * sum(abs(X))
      }
      
      if ((tol < outer.tol) && (kk > 1)) break
    }
  }
  
  time.counter.QP <- colSums(time.counter.QP[1:ii, , drop = FALSE])
  
  if (isTRUE(obj.seq)) {
    list(
      time.counter.QP = time.counter.QP,
      rel.err = rel.err[1:kk],
      obj.vals = obj.vals[1:kk],
      X = X,
      invX = Sigma + U.mat,
      sparse.nos = sparse.nos[1:kk]
    )
  } else {
    list(
      time.counter.QP = time.counter.QP,
      rel.err = rel.err[1:kk],
      X = X,
      invX = Sigma + U.mat,
      sparse.nos = sparse.nos[1:kk]
    )
  }
}
