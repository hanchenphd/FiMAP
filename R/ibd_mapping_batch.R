#' Conduct IBD mapping for multiple local IBD matrices in batches
#'
#' @param glmmkin.randomvec.obj The glmmkin.randomvec object
#' @param Sigma An N by N sparse covariance matrix for the fitted glmmkin (default = NULL)
#' @param offset An N.randomvec by N.randomvec dense matrix computed as the cross-product of Sigma with glmmkin.randomvec matrix of random vectors on both sides (default = NULL)
#' @param ibd.ids The ID vector for row and column names of the local IBD matrices
#' @param ibd.file.list A character vector of input file names for local IBD matrices
#' @param ncores A positive integer indicating the number of cores to be used in parallel computing (default = 1).
#' @return A data.frame
#' \item{chr}{The chromosome}
#' \item{start}{The start position}
#' \item{end}{The end position}
#' \item{n.nonzero}{Number of non-zero entries in the upper triangle of the local IBD matrix}
#' \item{p.value}{Finite-sample p-value}
#' \item{p.value.asymptotic}{Asymptotic p-value}
#' @export
ibd_mapping_batch <- function(glmmkin.randomvec.obj, Sigma = NULL, offset = NULL, ibd.ids, ibd.file.list, ncores = 1) {
	if(is.null(Sigma) && is.null(offset)) stop("Error: at least one of Sigma and offset must be specified.")
	ids <- glmmkin.randomvec.obj$id_include
	if(any(!ids %in% ibd.ids)) stop("Error: id_include from glmmkin.randomvec.obj should be a subset of ibd.ids.")
	N <- length(ids)
	ncovars <- glmmkin.randomvec.obj$ncovar
	match.idx <- match(ids, ibd.ids)
	N.randomvec <- ncol(glmmkin.randomvec.obj$random.vectors)
	if(is.null(offset)) {
		offset <- crossprod(glmmkin.randomvec.obj$random.vectors, crossprod(Sigma, glmmkin.randomvec.obj$random.vectors))
		rm(Sigma)
	}
	Q.offset <- 2*sum(glmmkin.randomvec.obj$scaled.residuals^2)
	offset0 <- 2*crossprod(glmmkin.randomvec.obj$random.vectors)
	offset1 <- offset0-Q.offset/(N-ncovars)*offset
	n.p.all <- length(ibd.file.list)
	ncores <- min(c(ncores, n.p.all, parallel::detectCores(logical = TRUE)))
	if(ncores > 1) {
		doMC::registerDoMC(cores = ncores)
		n.p.percore <- (n.p.all-1) %/% ncores + 1
    		n.p.percore_1 <- n.p.percore * ncores - n.p.all
		b <- NULL
		resout <- foreach(b = 1:ncores, .combine = c, .multicombine = TRUE, .inorder=TRUE, .options.multicore = list(preschedule = FALSE, set.seed = FALSE)) %dopar% {
    		        p.idx <- if(b <= n.p.percore_1) ((b-1)*(n.p.percore-1)+1):(b*(n.p.percore-1)) else (n.p.percore_1*(n.p.percore-1)+(b-n.p.percore_1-1)*n.p.percore+1):(n.p.percore_1*(n.p.percore-1)+(b-n.p.percore_1)*n.p.percore)
			lapply(p.idx, function(ii) {
				ibd.list <- get(load(ibd.file.list[ii]))
				ibd <- ibd.list$matrix
				if(is.null(ibd)) return(data.frame(chr = ibd.list$chr, start = ibd.list$start, end = ibd.list$end, n.nonzero = 0, p.value = NA))
				ibd <- ibd[match.idx, match.idx]
				Q <- as.numeric(crossprod(glmmkin.randomvec.obj$scaled.residuals, crossprod(ibd, glmmkin.randomvec.obj$scaled.residuals)))
				eigenmat <- crossprod(glmmkin.randomvec.obj$random.vectors, crossprod(ibd, glmmkin.randomvec.obj$random.vectors))
				lambda <- zapsmall(eigen(eigenmat+offset1-Q/(N-ncovars)*offset, only.values = TRUE, symmetric = TRUE)$values/N.randomvec)
				p.value <- try(.Q_pval(0, lambda))
				if(class(p.value) == "try-error") p.value <- NA
				if(all(lambda < 0)) p.value <- 0
				lambda1 <- zapsmall(eigen(eigenmat+offset0, only.values = TRUE, symmetric = TRUE)$values/N.randomvec)
				p.value1 <- try(.Q_pval(Q+Q.offset, lambda1))
				if(class(p.value1) == "try-error") p.value1 <- NA
				return(data.frame(chr = ibd.list$chr, start = ibd.list$start, end = ibd.list$end, n.nonzero = length(ibd@i), p.value = p.value, p.value.asymptotic = p.value1))
			})
		}
	} else {
		resout <- lapply(1:n.p.all, function(ii) {
			ibd.list <- get(load(ibd.file.list[ii]))
			ibd <- ibd.list$matrix
			if(is.null(ibd)) return(data.frame(chr = ibd.list$chr, start = ibd.list$start, end = ibd.list$end, n.nonzero = 0, p.value = NA))
			ibd <- ibd[match.idx, match.idx]
			Q <- as.numeric(crossprod(glmmkin.randomvec.obj$scaled.residuals, crossprod(ibd, glmmkin.randomvec.obj$scaled.residuals)))
			eigenmat <- crossprod(glmmkin.randomvec.obj$random.vectors, crossprod(ibd, glmmkin.randomvec.obj$random.vectors))
			lambda <- zapsmall(eigen(eigenmat+offset1-Q/(N-ncovars)*offset, only.values = TRUE, symmetric = TRUE)$values/N.randomvec)
			p.value <- try(.Q_pval(0, lambda))
			if(class(p.value) == "try-error") p.value <- NA
			if(all(lambda < 0)) p.value <- 0
			lambda1 <- zapsmall(eigen(eigenmat+offset0, only.values = TRUE, symmetric = TRUE)$values/N.randomvec)
			p.value1 <- try(.Q_pval(Q+Q.offset, lambda1))
			if(class(p.value1) == "try-error") p.value1 <- NA
			return(data.frame(chr = ibd.list$chr, start = ibd.list$start, end = ibd.list$end, n.nonzero = length(ibd@i), p.value = p.value, p.value.asymptotic = p.value1))
		})
	}
	do.call(rbind, resout)
}
