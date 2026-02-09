# Data-generating mechanisms (DGMs) for precision matrices and data

#' Generate banded precision matrix (Case 1)
generate_theta_banded <- function(p, off_diag = 0.5, diag = 1) {
  stop("TODO")
}

#' Generate dense precision matrix (Case 2)
generate_theta_dense <- function(p, off_diag = 0.5, diag = 1) {
  stop("TODO")
}

#' Generate random sparse precision matrix (Case 3)
generate_theta_random_sparse <- function(p, p0, eta = 0.5, seed = NULL) {
  stop("TODO")
}

#' Generate hub-structured precision matrix (Case 4)
generate_theta_hub <- function(p, base_edge_prob = 0.02, hub_edge_prob = 0.7, hub_ratio = 0.2, eta = 0.5, seed = NULL) {
  stop("TODO")
}

#' Sample data X ~ N(0, Sigma) given Theta = Sigma^{-1}
sample_gaussian_from_theta <- function(theta, n, seed = NULL) {
  stop("TODO")
}
