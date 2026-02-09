# Wrapper functions that call each method with a consistent API.
# Each method should return:
# - theta_hat (estimated precision matrix)
# - runtime_sec
# - selected_lambda (or tuning params)
# - diagnostics (optional)

fit_glasso <- function(S, lambda) {
  stop("TODO")
}

fit_quic <- function(S, lambda) {
  stop("TODO")
}

fit_dpglasso <- function(S, lambda) {
  stop("TODO")
}

fit_hubglasso <- function(S, lambda1, lambda2, lambda3) {
  stop("TODO")
}
