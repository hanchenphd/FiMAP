#' Calculate pairwise global IBD proportions from IBD segment summary files
#'
#' @param ids The full ID vector
#' @param infiles A character vector of input file names
#' @param total.length Total length of the genome (can be in cM or bp; default = 3391.3543292972 cM)
#' @param outfile.prefix Prefix of output file names. The output files are indexed from 1 to \code{length(infiles)} and gzipped.
#' @param id1.col Column number of individual 1 ID in \code{infiles}
#' @param id2.col Column number of individual 2 ID in \code{infiles}
#' @param start.col Column number of start positions for each IBD segment in \code{infiles}
#' @param end.col Column number of end positions for each IBD segment in \code{infiles}
#' @param is.integer A logical switch: whether the start and end positions of each IBD segment are integers (default = FALSE for cM; use \code{is.integer = TRUE} for bp)
#' @param ncores A positive integer indicating the number of cores to be used in parallel computing (default = 1).
#' @return NULL
#' @export
process_global_ibd <- function(ids, infiles, total.length = 3391.3543292972, outfile.prefix, id1.col = 2, id2.col = 3, start.col, end.col, is.integer = FALSE, ncores = 1) {
	N <- length(ids)
	n.p.all <- length(infiles)
	ncores <- min(c(ncores, n.p.all, parallel::detectCores(logical = TRUE)))
	if(ncores > 1) {
		doMC::registerDoMC(cores = ncores)
		n.p.percore <- (n.p.all-1) %/% ncores + 1
    		n.p.percore_1 <- n.p.percore * ncores - n.p.all
		b <- NULL
		ibd <- foreach(b = 1:ncores, .combine = c, .multicombine = TRUE, .inorder=FALSE, .options.multicore = list(preschedule = FALSE, set.seed = FALSE)) %dopar% {
    		        p.idx <- if(b <= n.p.percore_1) ((b-1)*(n.p.percore-1)+1):(b*(n.p.percore-1)) else (n.p.percore_1*(n.p.percore-1)+(b-n.p.percore_1-1)*n.p.percore+1):(n.p.percore_1*(n.p.percore-1)+(b-n.p.percore_1)*n.p.percore)
			lapply(p.idx, function(ii) {
	                	res <- read.table(infiles[ii], as.is = T)
				res <- res[res[,id1.col] %in% ids & res[,id2.col] %in% ids,]
				idx1 <- match(res[,id1.col], ids)
				idx2 <- match(res[,id2.col], ids)
				idxmin <- pmin(idx1, idx2)
				idxmax <- pmax(idx1, idx2)
				if(is.integer) len <- (res[,end.col] - res[,start.col] + 1)/total.length
				else len <- (res[,end.col] - res[,start.col])/total.length
				rm(res, idx1, idx2)
				idx <- order(idxmin, idxmax)
				idxmin <- idxmin[idx]
				idxmax <- idxmax[idx]
				len <- len[idx]
				res.id <- paste(idxmin, idxmax, sep = ":")
				idx <- duplicated(res.id)
				unique.id <- res.id[!idx]
				idxmin <- idxmin[!idx]
				idxmax <- idxmax[!idx]
				idx.end <- findInterval(1:length(unique.id), match(res.id, unique.id))
				rm(idx, res.id)
                                idx.start <- c(1, idx.end[-length(idx.end)] + 1)
                                len <- sapply(1:length(unique.id), function(xx) sum(len[idx.start[xx]:idx.end[xx]]))
                                rm(idx.end, idx.start)
				outfile <- gzfile(paste0(outfile.prefix, ii, ".gz"), "w")
				write.table(cbind(unique.id, idxmin, idxmax, len), outfile, quote=F, row.names=F, col.names=F, sep="\t")
				close(outfile)
				invisible(NULL)
			})
		}
	} else {
		ibd <- lapply(1:n.p.all, function(ii) {
	                res <- read.table(infiles[ii], as.is = T)
			res <- res[res[,id1.col] %in% ids & res[,id2.col] %in% ids,]
			idx1 <- match(res[,id1.col], ids)
			idx2 <- match(res[,id2.col], ids)
			idxmin <- pmin(idx1, idx2)
			idxmax <- pmax(idx1, idx2)
			if(is.integer) len <- (res[,end.col] - res[,start.col] + 1)/total.length
			else len <- (res[,end.col] - res[,start.col])/total.length
			rm(res, idx1, idx2)
			idx <- order(idxmin, idxmax)
			idxmin <- idxmin[idx]
			idxmax <- idxmax[idx]
			len <- len[idx]
			res.id <- paste(idxmin, idxmax, sep = ":")
			idx <- duplicated(res.id)
			unique.id <- res.id[!idx]
			idxmin <- idxmin[!idx]
			idxmax <- idxmax[!idx]
			idx.end <- findInterval(1:length(unique.id), match(res.id, unique.id))
			rm(idx, res.id)
                        idx.start <- c(1, idx.end[-length(idx.end)] + 1)
                        len <- sapply(1:length(unique.id), function(xx) sum(len[idx.start[xx]:idx.end[xx]]))
                        rm(idx.end, idx.start)
			outfile <- gzfile(paste0(outfile.prefix, ii, ".gz"), "w")
			write.table(cbind(unique.id, idxmin, idxmax, len), outfile, quote=F, row.names=F, col.names=F, sep="\t")
			close(outfile)
			invisible(NULL)
		})
	}
	invisible(NULL)
}

