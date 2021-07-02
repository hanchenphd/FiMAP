#' Calculate local IBD proportions from IBD segment summary files
#'
#' @param ids The full ID vector
#' @param infiles A character vector of input file names
#' @param chr.length Total length of the chromosome (can be in cM or bp)
#' @param window.size Chunking window size (can be in cM or bp)
#' @param id1.col Column number of individual 1 ID in \code{infiles}
#' @param id2.col Column number of individual 2 ID in \code{infiles}
#' @param start.col Column number of start positions for each IBD segment in \code{infiles}
#' @param end.col Column number of end positions for each IBD segment in \code{infiles}
#' @param is.integer A logical switch: whether the start and end positions of each IBD segment are integers (default = FALSE for cM; use \code{is.integer = TRUE} for bp)
#' @param cut.off Threshold to zero out small IBD values (default = 0)
#' @param ncores A positive integer indicating the number of cores to be used in parallel computing (default = 1).
#' @return A list of length equal to the number of windows. Each window is also a list with the following three components:
#' \item{start}{The start position}
#' \item{end}{The end position}
#' \item{matrix}{A dsCMatrix}
#' @export
get_local_ibd <- function(ids, infiles, chr.length, window.size, id1.col = 2, id2.col = 3, start.col, end.col, is.integer = FALSE, cut.off = 0, ncores = 1) {
	N <- length(ids)
	n.p.all <- length(infiles)
	n.q.all <- ceiling(chr.length/window.size)
	ncores <- min(c(ncores, n.p.all, parallel::detectCores(logical = TRUE)))
	if(ncores > 1) {
		doMC::registerDoMC(cores = ncores)
		n.p.percore <- (n.p.all-1) %/% ncores + 1
    		n.p.percore_1 <- n.p.percore * ncores - n.p.all
		b <- NULL
		RaPID.list <- foreach(b = 1:ncores, .combine = c, .multicombine = TRUE, .inorder=TRUE, .options.multicore = list(preschedule = FALSE, set.seed = FALSE)) %dopar% {
    		        p.idx <- if(b <= n.p.percore_1) ((b-1)*(n.p.percore-1)+1):(b*(n.p.percore-1)) else (n.p.percore_1*(n.p.percore-1)+(b-n.p.percore_1-1)*n.p.percore+1):(n.p.percore_1*(n.p.percore-1)+(b-n.p.percore_1)*n.p.percore)
			lapply(p.idx, function(ii) {
	                	res <- read.table(infiles[ii], as.is = T)
				res <- res[res[,id1.col] %in% ids & res[,id2.col] %in% ids,]
				names(res)[c(start.col, end.col)] <- c("start.pos", "end.pos")
				idx1 <- match(res[,id1.col], ids)
				idx2 <- match(res[,id2.col], ids)
				res$idxmin <- pmin(idx1, idx2)
				res$idxmax <- pmax(idx1, idx2)
				rm(idx1, idx2)
				res <- res[, c("idxmin", "idxmax", "start.pos", "end.pos")]
				lapply(1:n.q.all, function(jj) {
					if(is.integer) loc <- (jj-1)*window.size+1
					else loc <- (jj-1)*window.size
					end <- min(chr.length, jj*window.size)
		                	if(is.integer) overlap_len <- pmin(res$end.pos-res$start.pos+1,res$end.pos-loc+1,end-res$start.pos+1,end-loc+1)
		                	else overlap_len <- pmin(res$end.pos-res$start.pos,res$end.pos-loc,end-res$start.pos,end-loc+1)
                			overlap_len[overlap_len < 0] <- 0
                			if(all(overlap_len==0)) return(NULL)
                			resout <- res[overlap_len > 0,]
                			overlap_len <- overlap_len[overlap_len > 0]
					if(is.integer) resout$len <- overlap_len/(end-loc+1)
					else resout$len <- overlap_len/(end-loc)
					return(resout[, c("idxmin", "idxmax", "len")])
				})
			})
		}
		ncores <- min(ncores, n.q.all)
		doMC::registerDoMC(cores = ncores)
		n.p.percore <- (n.q.all-1) %/% ncores + 1
    		n.p.percore_1 <- n.p.percore * ncores - n.q.all
		ibd.list <- foreach(b = 1:ncores, .combine = c, .multicombine = TRUE, .inorder=TRUE, .options.multicore = list(preschedule = FALSE, set.seed = FALSE)) %dopar% {
    		        p.idx <- if(b <= n.p.percore_1) ((b-1)*(n.p.percore-1)+1):(b*(n.p.percore-1)) else (n.p.percore_1*(n.p.percore-1)+(b-n.p.percore_1-1)*n.p.percore+1):(n.p.percore_1*(n.p.percore-1)+(b-n.p.percore_1)*n.p.percore)
			lapply(p.idx, function(jj) {
				if(is.integer) loc <- (jj-1)*window.size+1
				else loc <- (jj-1)*window.size
				end <- min(chr.length, jj*window.size)
				RaPID <- do.call(rbind, lapply(RaPID.list, "[[", jj))
				if(is.null(RaPID)) return(list(start = as.integer(loc), end = as.integer(end), matrix = NULL))
				ibd <- sparseMatrix(RaPID$idxmin, RaPID$idxmax, x = RaPID$len, dims = c(N, N), symmetric = TRUE)
				colnames(ibd) <- rownames(ibd) <- ids
				ibd <- as(drop0(ibd, tol = cut.off), "symmetricMatrix")
				if(is.integer) return(list(start = as.integer(loc), end = as.integer(end), matrix = ibd))
				else return(list(start = loc, end = end, matrix = ibd))
			})
		}
	} else {
		RaPID.list <- lapply(1:n.p.all, function(ii) {
	                res <- read.table(infiles[ii], as.is = T)
			res <- res[res[,id1.col] %in% ids & res[,id2.col] %in% ids,]
			names(res)[c(start.col, end.col)] <- c("start.pos", "end.pos")
			idx1 <- match(res[,id1.col], ids)
			idx2 <- match(res[,id2.col], ids)
			res$idxmin <- pmin(idx1, idx2)
			res$idxmax <- pmax(idx1, idx2)
			rm(idx1, idx2)
			res <- res[, c("idxmin", "idxmax", "start.pos", "end.pos")]
			lapply(1:n.q.all, function(jj) {
				if(is.integer) loc <- (jj-1)*window.size+1
				else loc <- (jj-1)*window.size
				end <- min(chr.length, jj*window.size)
		                if(is.integer) overlap_len <- pmin(res$end.pos-res$start.pos+1,res$end.pos-loc+1,end-res$start.pos+1,end-loc+1)
		                else overlap_len <- pmin(res$end.pos-res$start.pos,res$end.pos-loc,end-res$start.pos,end-loc+1)
                		overlap_len[overlap_len < 0] <- 0
                		if(all(overlap_len==0)) return(NULL)
                		resout <- res[overlap_len > 0,]
                		overlap_len <- overlap_len[overlap_len > 0]
				if(is.integer) resout$len <- overlap_len/(end-loc+1)
				else resout$len <- overlap_len/(end-loc)
				return(resout[, c("idxmin", "idxmax", "len")])
			})
		})
		ibd.list <- lapply(1:n.q.all, function(jj) {
			if(is.integer) loc <- (jj-1)*window.size+1
			else loc <- (jj-1)*window.size
			end <- min(chr.length, jj*window.size)
			RaPID <- do.call(rbind, lapply(RaPID.list, "[[", jj))
			if(is.null(RaPID)) return(list(start = as.integer(loc), end = as.integer(end), matrix = NULL))
			ibd <- sparseMatrix(RaPID$idxmin, RaPID$idxmax, x = RaPID$len, dims = c(N, N), symmetric = TRUE)
			colnames(ibd) <- rownames(ibd) <- ids
			ibd <- as(drop0(ibd, tol = cut.off), "symmetricMatrix")
			if(is.integer) return(list(start = as.integer(loc), end = as.integer(end), matrix = ibd))
			else return(list(start = loc, end = end, matrix = ibd))
		})
	}
	ibd.list
}
