########################################

###  Code for: Contingency and determinism in rapid adaptation shape evolutionary rescue and extinction in Drosophila
###  Authors: René S. Shahmohamadloo*, Joshua J. Faber-Hammond, Amir R. Gabidulin, Catherine I. Clare, Seth M. Rudman*

###  In this file: liftover_dgrp2_frq.R — converts the DGRP2 founder-line allele frequencies (dgrp2_founder_list.frq,
###  FlyBase Release 5 / dm3 coordinates) to Release 6 / dm6 coordinates, the coordinates of the HAF-pipe table (BDGP6.46).
###  Uses the UCSC dm3ToDm6 chain (get_liftover_chain.sh). Positions are mapped through the highest-scoring chain that
###  contains them; positions that map to another chromosome arm, to the minus strand, or to the same position as
###  another SNP are dropped. Reads the chain from, and writes dgrp2_founder_list.dm6.frq to, PATHS$liftover (input to 11_hafpipe_validation.R).

########################################

f_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
RF_HOME <- if (length(f_arg)) normalizePath(dirname(sub("^--file=", "", f_arg))) else getwd()
source(file.path(RF_HOME, "00_config.R"))

## ---- Chain file: blocks of (dm3 start, dm6 start, length) per chain, chains in descending score order
ln <- readLines(gzfile(file.path(PATHS$liftover, "dm3ToDm6.over.chain.gz")))
h <- grep("^chain", ln); ends <- c(h[-1] - 1, length(ln))
chains <- lapply(seq_along(h), function(k) {
	f <- strsplit(ln[h[k]], "[ \t]+")[[1]]
	b <- ln[(h[k] + 1):ends[k]]; b <- b[nzchar(b)]
	v <- lapply(strsplit(b, "[ \t]+"), as.numeric)
	size <- sapply(v, `[`, 1); dt <- sapply(v, function(x) if (length(x) > 1) x[2] else 0); dq <- sapply(v, function(x) if (length(x) > 1) x[3] else 0)
	list(score = as.numeric(f[2]), t = f[3], q = f[8], qstrand = f[10],
		tstart = as.numeric(f[6]) + c(0, cumsum(size + dt)[-length(size)]),
		qstart = as.numeric(f[11]) + c(0, cumsum(size + dq)[-length(size)]), size = size)
})
chains <- chains[order(-sapply(chains, `[[`, "score"))]
msg("Chain file: %d chains", length(chains))

## ---- Lift the .frq positions -----------------------------------------------------------
dg <- read.delim(file.path(PATHS$raw_input, "dgrp2_founder_list.frq"), header = TRUE, sep = "\t", colClasses = c(CHROM = "character"))
dg <- dg[order(dg$CHROM, dg$POS), ]
new <- rep(NA_real_, nrow(dg)); p0 <- dg$POS - 1    # chain coordinates are 0-based
by_arm <- split(seq_len(nrow(dg)), paste0("chr", dg$CHROM))   # row indices per arm, in position order
for (ch in chains) {
	if (ch$qstrand != "+" || ch$q != ch$t || is.null(by_arm[[ch$t]])) next
	rows <- by_arm[[ch$t]]
	r <- findInterval(c(ch$tstart[1], max(ch$tstart + ch$size)) - c(0.5, 0.5), p0[rows])   # rows inside the chain's span
	if (r[2] <= r[1]) next
	i <- rows[(r[1] + 1):r[2]]; i <- i[is.na(new[i])]; if (!length(i)) next
	j <- findInterval(p0[i], ch$tstart)
	ok <- j > 0; ok[ok] <- p0[i][ok] < ch$tstart[j[ok]] + ch$size[j[ok]]
	new[i[ok]] <- ch$qstart[j[ok]] + (p0[i][ok] - ch$tstart[j[ok]]) + 1
}
out <- dg; out$POS <- as.integer(new); out <- out[!is.na(out$POS), ]
key <- paste(out$CHROM, out$POS); out <- out[!(duplicated(key) | duplicated(key, fromLast = TRUE)), ]
out <- out[order(out$CHROM, out$POS), ]
rep <- do.call(rbind, lapply(sort(unique(dg$CHROM)), function(c) data.frame(arm = c, dm3_SNPs = sum(dg$CHROM == c),
	lifted = sum(out$CHROM == c), pct = round(100 * sum(out$CHROM == c) / sum(dg$CHROM == c), 2))))
cat("\nSNPs lifted to dm6 per arm:\n"); print(rep, row.names = FALSE)
write.table(out, file.path(PATHS$liftover, "dgrp2_founder_list.dm6.frq"), sep = "\t", quote = FALSE, row.names = FALSE)
msg("Wrote %s (%s SNPs)", file.path(PATHS$liftover, "dgrp2_founder_list.dm6.frq"), format(nrow(out), big.mark = ","))
