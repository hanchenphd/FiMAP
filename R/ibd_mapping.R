#' @keywords internal
.Q_pval <- function(Q, lambda, method = "davies") {
    if(method == "davies") {
        tmp <- suppressWarnings(CompQuadForm::davies(q = Q, lambda = lambda, acc = 1e-6))
        pval <- tmp$Qq
	if((tmp$ifault > 0) | (pval <= 1e-5) | (pval >= 1)) method <- "kuonen"
    }
    if(method == "kuonen") {
    	pval <- .pKuonen(x = Q, lambda = lambda)
	if(is.na(pval)) method <- "liu"
    }
    if(method == "liu") pval <- CompQuadForm::liu(q = Q, lambda = lambda)
    return(pval)
}

#' @keywords internal
.pKuonen <- function (x, lambda, delta = rep(0, length(lambda)), df = rep(1, length(lambda)))
{
    delta <- delta[lambda != 0]
    df <- df[lambda != 0]
    lambda <- lambda[lambda != 0]
    if(length(lambda) != length(delta)) stop("Error: inconsistent length in lambda and delta!")
    if(length(lambda) != length(df)) stop("Error: inconsistent length in lambda and df!")
    if (length(lambda) == 1) {
        pchisq(x/lambda, df = df, ncp = delta, lower.tail = FALSE)
    }
    d <- max(lambda)
    lambda <- lambda/d
    x <- x/d
    k0 <- function(zeta) {
        -sum(df * log(1 - 2 * zeta * lambda))/2 + sum((delta * lambda *
            zeta)/(1 - 2 * zeta * lambda))
    }
    kprime0 <- function(zeta) {
        sapply(zeta, function(zz) {
            sum(((delta + df) * lambda)/(1 - 2 * zz * lambda) + 2 * (delta *
                zz * lambda^2)/(1 - 2 * zz * lambda)^2)
        })
    }
    kpprime0 <- function(zeta) {
        sum((2 * (2 * delta + df) * lambda^2)/(1 - 2 * zeta * lambda)^2 + 8 *
            delta * zeta * lambda^3/(1 - 2 * zeta * lambda)^3)
    }
    if (any(lambda < 0)) {
        lmin <- max(1/(2 * lambda[lambda < 0])) * 0.99999
    }
    else if (x > sum((df+delta)*lambda)) {
        lmin <- -0.01
    }
    else {
        lmin <- -length(lambda)*max(df+delta)/(2 * x)
    }
    lmax <- min(1/(2 * lambda[lambda > 0])) * 0.99999
    hatzeta <- uniroot(function(zeta) kprime0(zeta) - x, lower = lmin,
        upper = lmax, tol = 1e-08)$root
    w <- sign(hatzeta) * sqrt(2 * (hatzeta * x - k0(hatzeta)))
    v <- hatzeta * sqrt(kpprime0(hatzeta))
    if (abs(hatzeta) < 1e-04)
        NA
    else pnorm(w + log(v/w)/w, lower.tail = FALSE)
}

#' Conduct IBD mapping for a single local IBD matrix
#'
#' @param glmmkin.randomvec.obj The glmmkin.randomvec object
#' @param Sigma An N by N sparse covariance matrix for the fitted glmmkin (default = NULL)
#' @param offset An N.randomvec by N.randomvec dense matrix computed as the cross-product of Sigma with glmmkin.randomvec matrix of random vectors on both sides (default = NULL)
#' @param ibd.ids The ID vector for row and column names of the local IBD matrix
#' @param ibd.file The input file names for the local IBD matrix
#' @return A data.frame
#' \item{chr}{The chromosome}
#' \item{start}{The start position}
#' \item{end}{The end position}
#' \item{n.nonzero}{Number of non-zero entries in the upper triangle of the local IBD matrix}
#' \item{p.value}{Finite-sample p-value}
#' @export
ibd_mapping <- function(glmmkin.randomvec.obj, Sigma = NULL, offset = NULL, ibd.ids, ibd.file) {
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
	offset1 <- 2*crossprod(glmmkin.randomvec.obj$random.vectors)-Q.offset/(N-ncovars)*offset
	ibd.list <- get(load(ibd.file))
	ibd <- ibd.list$matrix
	if(is.null(ibd)) return(data.frame(chr = ibd.list$chr, start = ibd.list$start, end = ibd.list$end, n.nonzero = 0, p.value = NA))
	ibd <- ibd[match.idx, match.idx]
	Q <- as.numeric(crossprod(glmmkin.randomvec.obj$scaled.residuals, crossprod(ibd, glmmkin.randomvec.obj$scaled.residuals)))
	lambda <- zapsmall(eigen(crossprod(glmmkin.randomvec.obj$random.vectors, crossprod(ibd, glmmkin.randomvec.obj$random.vectors))+offset1-Q/(N-ncovars)*offset, only.values = TRUE, symmetric = TRUE)$values/N.randomvec)
	p.value <- try(.Q_pval(0, lambda))
	if(class(p.value) == "try-error") p.value <- NA
	if(all(lambda < 0)) p.value <- 0
	return(data.frame(chr = ibd.list$chr, start = ibd.list$start, end = ibd.list$end, n.nonzero = length(ibd@i), p.value = p.value))
}
