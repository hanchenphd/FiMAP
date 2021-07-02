#' Make a positive definite IBD matrix
#'
#' @param ibd The raw IBD matrix
#' @param tol Relative tolerance added to the diagonal entries compared to the maximum value of each column (default = 1e-5)
#' @return The positive definite IBD matrix
#' @export
makepd_add_diag <- function(ibd, tol = 1e-5) {
	diag(ibd) <- sapply(1:nrow(ibd), function(x) max(ibd[,x])*(1+tol))
	ibd
}
