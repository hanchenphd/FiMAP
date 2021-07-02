#' Get clusters from a dsCMatrix
#'
#' @param x The dsCMatrix
#' @return An index vector showing which observations are in the same cluster
#' @export
get_clusters <- function(x) {
        if(class(x)!="dsCMatrix") stop("Error: the input must be a dsCMatrix.")
	n <- nrow(x)
	v <- seq_len(n)
	v.list <- lapply(v, "[")
	names(v.list) <- v
	cols <- which(diff(x@p) > 0)
	curr.idx <- 0
	for(ii in cols) {
		rows <- x@i[(curr.idx+1):x@p[ii+1]] + 1
		curr.idx <- x@p[ii+1]
		min.v <- min(v[c(rows, ii)])
		check.list <- unique(v[c(rows, ii)])
		check.list <- check.list[check.list > min.v]
		for(jj in check.list) {
			v.list[[as.character(min.v)]] <- unique(c(v.list[[as.character(min.v)]], v.list[[as.character(jj)]]))
			v.list[[as.character(jj)]] <- NULL
		}
		v[v.list[[as.character(min.v)]]] <- min.v
	}
	v.uniq <- unique(v)
	v.new <- seq_along(v.uniq)
	v.new[match(v,v.uniq)]
}
