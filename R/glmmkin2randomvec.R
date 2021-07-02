#' Create random vectors for a glmmkin object
#'
#' @param obj The glmmkin object
#' @param Z A list of design matrices for the random effects. The length must match the number of variance components
#' @param N.randomvec The number of random vectors to generate
#' @param group.idx A length N index vector showing which observation belongs to which variance group, for heteroscedastic linear mixed models (default = NULL for homoscedastic linear mixed models)
#' @param cluster.idx A length N index vector showing which observation belongs to which cluster (default = NULL for no clusters)
#' @param robust A logical switch: whether robust variance should be used (default = FALSE)
#' @return A list of class glmmkin.randomvec
#' \item{theta}{Variance estimates, inherited from the glmmkin object}
#' \item{scaled.residuals}{Scaled residuals, inherited from the glmmkin object}
#' \item{random.vectors}{An N by N.randomvec matrix of the random vectors generated}
#' \item{ncovar}{Number of fixed-effects covariates (including the intercept) from the glmmkin object}
#' \item{id_include}{The ID vector of included samples, inherited from the glmmkin object}
#' @export
glmmkin2randomvec <- function(obj, Z = NULL, N.randomvec = 1000, group.idx=NULL,cluster.idx=NULL,robust = FALSE) {
        if(class(obj) != "glmmkin") stop("Error: \"obj\" must be a class glmmkin object.")
	N <- length(obj$id_include)
       	random.vectors <- matrix(rnorm(N*N.randomvec),nrow=N,ncol=N.randomvec)
	if(!is.null(obj$P) && !robust) {
        	eig <- eigen(obj$P, symmetric = TRUE)
        	random.vectors <- tcrossprod(eig$vectors, t(random.vectors * sqrt(pmax(eig$values, 0))))
        	rm(eig)
	} else {
	        if(obj$n.groups != 1 && (is.null(group.idx) || !all.equal(seq_len(obj$n.groups), sort(unique(group.idx))))) stop("Error: heteroscedastic linear mixed models should include a valid group.idx argument.")
		if(is.null(group.idx)) group.idx <- rep(1, N)
		if(!robust) random.vectors <- sqrt(obj$theta[group.idx]) * random.vectors
		else {
			res <- as.numeric(obj$Y - tcrossprod(obj$X, t(obj$coefficient)))
			if(is.null(cluster.idx)) random.vectors <- random.vectors * res
			else random.vectors <- random.vectors[match(cluster.idx,unique(cluster.idx)),] * res
	        }
		if(!is.null(Z)) {
			if(class(Z) != "list") stop("Error: \"Z\" must be a list of matrices.")
			if(length(Z) != length(obj$theta) - obj$n.groups) stop("Error: number of matrices in \"Z\" does not match the number of variance components in \"obj\".")
			for(i in 1:length(Z)) {
			      	if(nrow(Z[[i]]) != N) stop("Error: \"Z\" matrix ", i, " is not compatible in sample size with \"obj\".")
				p <- ncol(Z[[i]])
				if(obj$theta[i+obj$n.groups] < 0) stop("Error: negative variance component estimates are not allowed.")
				if(obj$theta[i+obj$n.groups] == 0) next
				random.vectors2 <- matrix(rnorm(p*N.randomvec), nrow=N.randomvec, ncol=p)
				random.vectors <- random.vectors + sqrt(obj$theta[i+obj$n.groups]) * tcrossprod(Z[[i]], random.vectors2)
			}
		}
		if(!is.null(obj$P)) random.vectors <- crossprod(obj$P, random.vectors)
		else random.vectors <- crossprod(obj$Sigma_i, random.vectors) - tcrossprod(obj$Sigma_iX, tcrossprod(crossprod(random.vectors, obj$Sigma_iX), obj$cov))
	}
       	out <- list(theta = obj$theta, scaled.residuals = obj$scaled.residuals, random.vectors = as.matrix(random.vectors), ncovar = ncol(obj$X), id_include = obj$id_include)
       	class(out) <- "glmmkin.randomvec"
       	return(out)
}
