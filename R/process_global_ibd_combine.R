#' Assemble pairwise global IBD proportions into a sparse symmetric matrix
#'
#' @param ids The full ID vector
#' @param infile The input file name
#' @param cut.off Threshold to zero out small IBD values (default = 0)
#' @return A dsCMatrix
#' @export
process_global_ibd_combine <- function(ids, infile, cut.off = 0) {
	N <- length(ids)
	res <- read.table(infile, as.is = T)
	ibd <- sparseMatrix(i = res$V1, j = res$V2, x = res$V3, dims = c(N, N), symmetric = TRUE) + 2 * Diagonal(N)
	colnames(ibd) <- rownames(ibd) <- ids
	as(drop0(ibd, tol = cut.off), "symmetricMatrix")
}
