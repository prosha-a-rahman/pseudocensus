#' Solve a linear system
#'
#' Solves the linear system `A %*% x = b`, or failing that, identifies the
#' linear least squares solution by computing the Moore–Penrose inverse of `A`
#' through Singular Value Decomposition.
#'
#' @param A A numeric matrix.
#' @param b A numeric vector of length `nrow(A)`.
#' @param tol Singular values at or below `tol` are treated as zero when forming
#'   the pseudoinverse. Defaults to
#'   `max(dim(A)) * max(singular values) * .Machine$double.eps`.
#'
#'
#' @returns A numeric vector of length `ncol(A)`.
#'
#' @noRd
solve_linear_system <- function(A, b, tol = NULL) {
  s <- svd(A)
  d <- s$d
  if (is.null(tol)) tol <- max(dim(A)) * max(d, 0) * .Machine$double.eps
  dinv <- ifelse(d > tol, 1 / d, 0)
  drop(s$v %*% (dinv * crossprod(s$u, b)))
}




#' Transpose a list of matrices by swapping list and matrix row indexes
#'
#' @param L a list of `K`-by-`p` numeric matrices.
#'
#' @returns A length `K` list of `length(L)`-by-`p` numeric matrices.
#'
#' @noRd
transpose_list <- function(L) {
  K <- nrow(L[[1]])
  p <- ncol(L[[1]])
  r <- length(L)

  # Map L to 3-d tensor where A[k, j, i] = L[[i]][k, j]
  A <- array(unlist(L, use.names = FALSE), dim = c(K, p, r))
  # Transpose dimensions of 3-d tensor where B[i, j, k] = L[[i]][k, j]
  B <- aperm(A, c(3, 2, 1))

  # Write tensor to list
  lapply(seq_len(K), function(k) matrix(B[, , k], nrow = r, ncol = p))
}



#' Prepend a design matrix with a column of ones if required
#'
#' Add a column of ones to a numeric matrix so that it is a valid design matrix
#' for ordinary least squares computation
#'
#' @param x A numeric matrix.
#'
#' @returns
#' If the first column of `x` are all ones, then retun `x`, else prepend a
#' column of ones to `x`.
#'
#' @noRd
prepend_ones <- function(x) {
  has_ones <- ncol(x) >= 1 && all(x[, 1] == 1)

  if (has_ones) x else cbind(1, x)
}

