########################################

###  Code for: Contingency and determinism in rapid adaptation shape evolutionary rescue and extinction in Drosophila
###  Authors: René S. Shahmohamadloo*, Joshua J. Faber-Hammond, Amir R. Gabidulin, Catherine I. Clare, Seth M. Rudman*
###  *Corresponding authors:
###  René S. Shahmohamadloo, School of Biological Sciences, Washington State University, 14204 NE Salmon Creek Ave, Vancouver, WA 98686, USA. Email: rene.shahmohamadloo@wsu.edu
###  Seth M. Rudman, School of Biological Sciences, Washington State University, 14204 NE Salmon Creek Ave, Vancouver, WA 98686, USA. Email: seth.rudman@wsu.edu

###  In this file: 07_concordance_LOO.R.

###  Produces:
###  — parallel (concordant) allele-frequency change in exposed populations: SNP clusters, leave-one-out concordance t-tests, Fig. S9.             

###  Outputs:
###  Tables S3, S4, S6   concordance t-tests (extinct Aug; rescued Aug; rescued Nov)
###  SI 2.4              cluster counts (all-sample GLMs) and leave-one-out ranges
###  Fig. S9             Euclidean distances among August populations at the 441 November convergent candidate SNPs (from module 05)

########################################

f_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
RF_HOME <- if (length(f_arg)) normalizePath(dirname(sub("^--file=", "", f_arg))) else getwd()
source(file.path(RF_HOME, "00_config.R"))
suppressMessages({ library(dplyr); library(ggplot2) })
OUT <- module_dir("07_concordance_LOO"); setwd(OUT)
sink("07_concordance_LOO_report.txt", split = TRUE)
dir.create("clusters", FALSE); dir.create("ttests", FALSE)
FIX <- toupper(Sys.getenv("RF_FIX_SCORE_PAIRING", "FALSE")) == "TRUE"
if (FIX) cat("*** RF_FIX_SCORE_PAIRING=TRUE: scores use matching columns; results will NOT match the paper ***\n")
rd <- function(f) read.table(file.path(PATHS$glm_tables, f), header = TRUE)
CHRS <- c("2L", "2R", "3L", "3R", "X")
N_ITER <- 100

## The matched-SNP draw below uses sample.int() directly; confirm it gives the same draw as
## dplyr::sample_n() with this dplyr version.
same_draw <- all(vapply(c(1, 2, 7, 5000), function(n) {
	set.seed(7); a <- dplyr::sample_n(data.frame(i = seq_len(n)), 1)$i; r1 <- runif(1)
	set.seed(7); b <- sample.int(n, 1); r2 <- runif(1); a == b && r1 == r2 }, logical(1)))
cat("sample.int() reproduces dplyr::sample_n():", same_draw, "\n")
draw1 <- if (same_draw) function(ii) ii[sample.int(length(ii), 1)] else
	function(ii) dplyr::sample_n(data.frame(i = ii), 1)$i

##########################################################################################
## 1. Allele frequencies and allele-frequency differences                              ##
##########################################################################################
haf <- load_haf(exclude = TRUE)
sites_filt <- field_samples(haf)$sites          # haf.sites.filt (all time points, variance filter)
## Founder means: sample columns 3-5 (E) and 6-8 (PA) of the HAF-pipe table
cat("Founder columns used (E; PA):", paste(names(haf$freq)[3:5], collapse = ","), ";", paste(names(haf$freq)[6:8], collapse = ","), "\n")
Fmean <- cbind(haf$freq[, 1:2], E = rowMeans(haf$freq[, 3:5]), PA = rowMeans(haf$freq[, 6:8]))
is_sorted <- function(s) identical(order(s$CHROM, s$POS), seq_len(nrow(s)))

## August: four groups; cage order of the leave-one-out loops = metadata sorted by sample
mt <- merge(haf$meta, COND_T1, by = c("cage", "tpt"))
mt <- unique(mt[mt$batch == "a" & mt$treat %in% c("S", "PA", "E"), ])
mt <- mt[mt$samp %in% KEEP_SAMP, ]
mt <- mt[order(mt$samp), ]
SE_loo <- as.character(mt$cage[mt$condition == "SE"])
SP_loo <- as.character(mt$cage[mt$condition == "SP"])
cat("Leave-one-out cage order (sample order): extinct", paste(SE_loo, collapse = ","), "| rescued", paste(SP_loo, collapse = ","), "\n")
ft <- haf$freq[, names(haf$freq) %in% mt$samp]
cat("August samples: frequency columns in metadata order:", identical(as.character(mt$samp), names(ft)), "\n")
m1 <- mt[match(names(ft), as.character(mt$samp)), ]           # metadata in frequency-column order
v1 <- rowVars(as.matrix(ft)) > 0.001
sites1 <- haf$freq[v1, 1:2]; ft <- ft[v1, ]
cat("August SNPs (variance filter):", format(nrow(sites1), big.mark = ","), "; already sorted by CHROM/POS:", is_sorted(sites1), "\n")

## November: E, PA and S samples at time point 4
m4 <- haf$meta[haf$meta$tpt == "4" & haf$meta$batch == "a" & haf$meta$treat.fix %in% c("PA", "S", "E"), ]
f4 <- haf$freq[, names(haf$freq) %in% m4$samp]
m4 <- m4[m4$samp %in% names(f4), ]
cat("November samples: frequency columns in metadata order:", identical(as.character(m4$samp), names(f4)), "\n")
m4 <- m4[match(names(f4), as.character(m4$samp)), ]
v4 <- rowVars(as.matrix(f4)) > 0.001
sites4 <- haf$freq[v4, 1:2]; f4 <- f4[v4, ]
cat("November SNPs (variance filter):", format(nrow(sites4), big.mark = ","), "; already sorted:", is_sorted(sites4), "\n")
cat("November samples per exposed cage:", paste(names(table(m4$cage[m4$treat.fix == "S"])), table(m4$cage[m4$treat.fix == "S"]), sep = ":", collapse = " "), "\n")

rm_ <- function(x, sel) rowMeans(x[, sel, drop = FALSE])
one <- function(x, sel, what) {     # a single left-in population (the first column if there are two)
	y <- x[, sel, drop = FALSE]
	if (ncol(y) != 1) cat("  note:", what, "has", ncol(y), "samples; the first (", names(y)[1], ") is used\n")
	y[, 1]
}
to_bed <- function(sites, d) { b <- cbind(sites, STOP = sites$POS + 1, d); b[order(b[, 1], b[, 2]), ] }
c1 <- m1$condition; k1 <- as.character(m1$cage)
E1 <- rm_(ft, c1 == "E"); PA1 <- rm_(ft, c1 == "PA"); SE1 <- rm_(ft, c1 == "SE"); SP1 <- rm_(ft, c1 == "SP")
fd_t1 <- cbind(EvSE.T1.diff = E1 - SE1, EvSP.T1.diff = E1 - SP1, EvPA.T1.diff = E1 - PA1,
	SEvPA.T1.diff = SE1 - PA1, SPvPA.T1.diff = SP1 - PA1, SEvSP.T1.diff = SE1 - SP1)
no1 <- only1 <- list()
for (i in SE_loo) {
	se_no <- rm_(ft, c1 == "SE" & k1 != i); se_in <- one(ft, c1 == "SE" & k1 == i, paste("August cage", i))
	no1[[i]] <- cbind(EvSE.T1.diff = E1 - se_no, EvSP.T1.diff = E1 - SP1, EvPA.T1.diff = E1 - PA1,
		SEvPA.T1.diff = se_no - PA1, SPvPA.T1.diff = SP1 - PA1, SEvSP.T1.diff = se_no - SP1)
	only1[[i]] <- cbind(EvSE.T1.diff = E1 - se_in, EvSP.T1.diff = E1 - SP1, EvPA.T1.diff = E1 - PA1,
		SEvPA.T1.diff = se_in - PA1, SPvPA.T1.diff = SP1 - PA1, SEvSP.T1.diff = se_in - SP1)
}
for (i in SP_loo) {
	sp_no <- rm_(ft, c1 == "SP" & k1 != i); sp_in <- one(ft, c1 == "SP" & k1 == i, paste("August cage", i))
	no1[[i]] <- cbind(EvSE.T1.diff = E1 - SE1, EvSP.T1.diff = E1 - sp_no, EvPA.T1.diff = E1 - PA1,
		SEvPA.T1.diff = SE1 - PA1, SPvPA.T1.diff = sp_no - PA1, SEvSP.T1.diff = SE1 - sp_no)
	only1[[i]] <- cbind(EvSE.T1.diff = E1 - SE1, EvSP.T1.diff = E1 - sp_in, EvPA.T1.diff = E1 - PA1,
		SEvPA.T1.diff = SE1 - PA1, SPvPA.T1.diff = sp_in - PA1, SEvSP.T1.diff = SE1 - sp_in)
}
no1_bed <- lapply(no1, function(d) to_bed(sites1, d)); only1_bed <- lapply(only1, function(d) to_bed(sites1, d))

c4 <- m4$treat.fix; k4 <- as.character(m4$cage)
E4 <- rm_(f4, c4 == "E"); PA4 <- rm_(f4, c4 == "PA"); S4 <- rm_(f4, c4 == "S")
fd_t4 <- cbind(EvSP.T4.diff = E4 - S4, EvPA.T4.diff = E4 - PA4, PAvSP.T4.diff = PA4 - S4)
no4 <- only4 <- list()
for (i in SP_loo) {
	s_no <- rm_(f4, c4 == "S" & k4 != i); s_in <- one(f4, c4 == "S" & k4 == i, paste("November cage", i))
	no4[[i]] <- cbind(EvSP.T4.diff = E4 - s_no, EvPA.T4.diff = E4 - PA4, SPvPA.T4.diff = s_no - PA4)
	only4[[i]] <- cbind(EvSP.T4.diff = E4 - s_in, EvPA.T4.diff = E4 - PA4, SPvPA.T4.diff = s_in - PA4)
}
no4_bed <- lapply(no4, function(d) to_bed(sites4, d)); only4_bed <- lapply(only4, function(d) to_bed(sites4, d))

## Sign-corrected allele-frequency difference of each left-out population: its difference
## from unexposed, multiplied by the sign of the difference of the remaining populations
sign_table <- function(no_bed, only_bed, cages, f) {
	st <- sites_filt
	for (i in cages) {
		s <- cbind(no_bed[[i]][, 1:2], only_bed[[i]][, f] * (abs(no_bed[[i]][, f]) / no_bed[[i]][, f]))
		names(s)[3] <- paste0("AFdiff_no", i)
		st <- merge(st, s, by = c("CHROM", "POS"))
	}
	st
}
sign_EvSET1 <- sign_table(no1_bed, only1_bed, SE_loo, 4)    # column 4 = EvSE.T1.diff
sign_EvSPT1 <- sign_table(no1_bed, only1_bed, SP_loo, 5)    # column 5 = EvSP.T1.diff
sign_EvSPT4 <- sign_table(no4_bed, only4_bed, SP_loo, 4)    # column 4 = EvSP.T4.diff
rm(ft, f4); invisible(gc())

##########################################################################################
## 2. GLM tables and VEP consequence priority (for matched SNPs)                       ##
##########################################################################################
T1 <- rd("rudflies_2023_PAvSvSEvE.wLoci.GLMcontrast.txt")
names(T1) <- c("CHROM", "POS", "EvPA.T1", "EvSE.T1", "EvSP.T1", "PAvSE.T1", "PAvSP.T1", "SEvSP.T1")
T1_all <- na.omit(T1)                                        # used for the leave-one-out table
T1 <- T1[T1$CHROM != "4", ]
PAvE <- merge(rd("rudflies_2023_PAvE_treat.all.table.GLMcontrast.txt")[, -3],
	rd("rudflies_2023_PAvE_treat.table.GLMcontrast.founders.txt")[, -3], by = c("CHROM", "POS"))
PAvE <- merge(PAvE, rd("rudflies_2023_PAvE_treat.table.GLMcontrast2.txt")[, -3], by = c("CHROM", "POS"))
names(PAvE) <- c("CHROM", "POS", "PAvE", "PAvE.F", paste0("PAvE.T", 1:4))
SvE <- merge(rd("rudflies_2023_SvE_treat.all.table.GLMcontrast.txt")[, c(1, 2, 4)],
	rd("rudflies_2023_SvE_treat.GLMcontrast.table.txt"), by = c("CHROM", "POS"))
SvE <- merge(SvE, rd("rudflies_2023_SvE_int.GLMcontrast.table.txt"), by = c("CHROM", "POS"))
names(SvE)[1:7] <- c("CHROM", "POS", "SvE", paste0("SvE.T", 1:4))
PAvS_trt <- rd("rudflies_2023_PAvS_treat.GLMcontrast.txt"); PAvS_tpt <- rd("rudflies_2023_PAvS_tpt.GLMcontrast.txt"); PAvS_int <- rd("rudflies_2023_PAvS_int.GLMcontrast.txt")
stopifnot(identical(PAvS_trt[, 1:2], PAvS_tpt[, 1:2]), identical(PAvS_trt[, 1:2], PAvS_int[, 1:2]))
PAvS_all <- rd("rudflies_2023_PAvS_treat.all.table.GLMcontrast.txt")[, -3]; names(PAvS_all)[3] <- "PAvS"
PAvS <- merge(PAvS_all, cbind(PAvS_trt, PAvS_tpt[, -c(1, 2)], PAvS_int[, -c(1, 2)]), by = c("CHROM", "POS"))
PAvS <- PAvS[order(PAvS[, 1], PAvS[, 2]), ]

## August leave-one-out table (all-sample + ten leave-one-out GLMs)
loo1 <- T1_all
for (cg in sort(as.integer(c(SE_CAGES, SP_CAGES)))) {
	x <- rd(paste0("rudflies_2023_PAvSvSEvE.wLociGLMcontrast.LOO.no", cg, ".txt"))
	names(x) <- c("CHROM", "POS", paste0(c("EvPA.T1", "EvSE.T1", "EvSP.T1", "PAvSE.T1", "PAvSP.T1", "SEvSP.T1"), "_no", cg))
	loo1 <- merge(loo1, x, by = c("CHROM", "POS"))
}
loo1 <- loo1[loo1$CHROM != "4", ]
## columns: extinct cages, then rescued (cage-number order)
loo1_cols <- c(paste0("EvSE.T1_no", SE_CAGES), paste0("EvSP.T1_no", SP_CAGES))
loo1_fdr <- loo1[, 1:2]
for (cc in loo1_cols) loo1_fdr[[paste0(cc, ".fdr")]] <- p.adjust(as.numeric(loo1[[cc]]), method = "fdr")

## November table (all-sample E vs S, E vs PA, PA vs S + five leave-one-out GLMs)
tpt4 <- merge(merge(cbind(SvE[, 1:2], EvSP.T4 = SvE$SvE.T4), cbind(PAvE[, 1:2], EvPA.T4 = PAvE$PAvE.T4), by = c("CHROM", "POS")),
	cbind(PAvS[, 1:2], PAvSP.T4 = PAvS$PAvS.T4), by = c("CHROM", "POS"))
for (cg in SP_CAGES) {
	x <- rd(paste0("rudflies_2023_PAvSvE.wLociGLMcontrast.TPT4.LOO.no", cg, ".txt"))
	names(x) <- c("CHROM", "POS", paste0(c("EvPA.T4", "EvSP.T4", "PAvSP.T4"), "_no", cg))
	tpt4 <- merge(tpt4, x, by = c("CHROM", "POS"))
}
tpt4 <- tpt4[tpt4$CHROM != "4", ]
for (cc in c("EvSP.T4", "EvPA.T4", "PAvSP.T4", paste0("EvSP.T4_no", SP_CAGES))) tpt4[[paste0(cc, ".fdr")]] <- p.adjust(as.numeric(tpt4[[cc]]), method = "fdr")
msg("GLM tables: August %s SNPs | August leave-one-out %s | November %s",
	format(nrow(T1), big.mark = ","), format(nrow(loo1), big.mark = ","), format(nrow(tpt4), big.mark = ","))

## One consequence per SNP, chosen by VEP impact then mean FlyCADD score of the consequence
glm_keys <- load_rds("glm_all")[, c("CHROM", "POS")]
vep <- load_vep()
annot <- merge(glm_keys, vep, by = c("CHROM", "POS"))
vep_cons <- unique(cbind(annot[, c(1, 2)], Consequence = annot$Consequence)); rm(annot)
vep_score <- as.data.frame(cbind(Consequence = vep$Consequence, FLYCADD = vep$FLYCADD, Extra = vep$Extra))
vep_score$Extra <- sub(".*=", "", sub(";.*", "", vep_score$Extra))
vep_score$FLYCADD <- suppressWarnings(as.numeric(vep_score$FLYCADD))
vep_score_mean <- as.data.frame(vep_score %>% group_by(Consequence, Extra) %>% summarise_all(mean))
vep_score_mean$Impact <- as.numeric(c(MODIFIER = 1, LOW = 2, MODERATE = 3, HIGH = 4)[vep_score_mean$Extra])
vep_score_mean <- vep_score_mean[order(-vep_score_mean[, 4], -vep_score_mean[, 3]), ]
rm(vep, vep_score); invisible(gc())
vep_priority <- c(); vep_unique <- merge(vep_cons, Fmean, by = c("CHROM", "POS"))
for (r in seq_len(nrow(vep_score_mean))) {
	tmp <- unique(vep_unique[vep_unique$Consequence == vep_score_mean[r, 1], ])
	vep_priority <- rbind(vep_priority, tmp)
	vep_unique <- anti_join(vep_unique, tmp, by = c("CHROM", "POS"))
}
msg("SNPs with a consequence for matching: %s", format(nrow(vep_priority), big.mark = ","))

##########################################################################################
## 3. Cluster helpers (501-SNP windows, 100-SNP steps)                                 ##
##########################################################################################
## Score 0-6 from FDR and |delta p| (Rudman et al. 2022)
score_of <- function(fdr, dp, six_on_score = FALSE) {
	s <- rep(NA_real_, length(fdr))
	s[which(fdr > 0.2)] <- 0
	s[which(fdr < 0.2 | (fdr > 0.2 & dp > 0.02))] <- 1
	for (k in 2:5) s[which(fdr < c(0.05, 0.01, 0.001, 1e-4)[k - 1] & dp > 0.02)] <- k
	s[which(fdr < 1e-5 & (if (six_on_score) s else dp) > 0.02)] <- 6
	s
}
## Shuffle loci within chromosomes (one sample() per chromosome)
shuffle <- function(sc) {
	r <- c()
	for (ch in unique(sc$CHROM)) { t <- sc[sc$CHROM == ch, ]; t[, c("CHROM", "POS")] <- t[sample(1:nrow(t)), c("CHROM", "POS")]; r <- rbind(r, t) }
	r[order(r[, 1], r[, 2]), ]
}
by_chrom <- function(sc) setNames(lapply(CHRS, function(ch) { t <- unique(sc[sc$CHROM == ch, ]); t <- t[order(t$POS), ]; t$POS <- as.integer(t$POS); t }), CHRS)
## One-sided Welch t-tests: true-locus window scores > shuffled-locus window scores
window_tests <- function(tr, ra, col) {
	out <- list()
	for (ch in CHRS) {
		x <- tr[[ch]][[col]]; y <- ra[[ch]][[col]]; p <- tr[[ch]]$POS; n <- length(x)
		if (n == 0) next
		for (j in seq(1, n, by = 100)) {
			w <- max(j - 250, 1):min(j + 250, n)
			tt <- tryCatch(t.test(x[w], y[w], alternative = "greater"), error = function(e) NULL)   # skipped when both windows are constant
			if (is.null(tt)) next
			out[[length(out) + 1]] <- data.frame(CHROM = ch, POS = p[j], START = min(p[w]), STOP = max(p[w]),
				statistic = unname(tt$statistic), parameter = unname(tt$parameter), p.value = tt$p.value, stderr = tt$stderr,
				mean_score = mean(x[w]), med_score = median(x[w]), max_score = max(x[w]), rand_score = mean(y[w]))
		}
	}
	d <- na.omit(do.call(rbind, out))
	d$fdr <- p.adjust(as.numeric(vapply(d$p.value, toString, "")), method = "fdr")   # FDR computed on the toString() values
	d
}
## Merge significant windows that overlap or touch (GenomicRanges::reduce rules)
reduce_iv <- function(d) {
	res <- list()
	for (ch in CHRS) {
		s <- d[d$CHROM == ch, ]; if (!nrow(s)) next
		s <- s[order(s$START, s$STOP), ]; st <- s$START[1]; en <- s$STOP[1]
		for (k in seq_len(nrow(s))[-1]) {
			if (s$START[k] <= en + 1) en <- max(en, s$STOP[k]) else { res[[length(res) + 1]] <- c(ch, st, en); st <- s$START[k]; en <- s$STOP[k] }
		}
		res[[length(res) + 1]] <- c(ch, st, en)
	}
	if (!length(res)) return(data.frame(seqnames = character(0), start = integer(0), end = integer(0), width = integer(0), strand = character(0)))
	m <- do.call(rbind, res)
	data.frame(seqnames = m[, 1], start = as.integer(m[, 2]), end = as.integer(m[, 3]), width = as.integer(m[, 3]) - as.integer(m[, 2]) + 1L, strand = "+")
}
## SNPs in clusters: [POS, POS+1] overlaps [start, end]; cluster named CHROM_start_end
assign_clusters <- function(loci, iv) {
	if (!nrow(iv)) return(data.frame(CHROM = character(0), POS = integer(0), clust = character(0), color = character(0)))
	iv$clust <- paste(iv$seqnames, iv$start, iv$end, sep = "_"); iv$color <- rep(c("blue4", "orange3"), length.out = nrow(iv))
	out <- list()
	for (ch in unique(iv$seqnames)) {
		s <- iv[iv$seqnames == ch, ]; s <- s[order(s$start), ]
		l <- loci[loci$CHROM == ch, c("CHROM", "POS")]
		k <- findInterval(l$POS + 1, s$start)
		ok <- k > 0; ok[ok] <- s$end[k[ok]] >= l$POS[ok]
		out[[ch]] <- data.frame(CHROM = l$CHROM[ok], POS = l$POS[ok], clust = s$clust[k[ok]], color = s$color[k[ok]])
	}
	do.call(rbind, out)
}
## Full clustering for a set of score columns; returns intervals and cluster-SNP tables
make_clusters <- function(sc, cols, loci, prefix) {
	set.seed(42)
	ra <- shuffle(sc)
	tr <- by_chrom(sc); rb <- by_chrom(ra)
	res <- list()
	for (cc in cols) {
		tt <- window_tests(tr, rb, cc)
		nm <- file.path("clusters", paste0(prefix, "_", sub("\\.score$", "", cc)))
		write.table(tt, paste0(nm, "_window_ttests.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
		iv <- reduce_iv(tt[tt$fdr < 0.05 & tt$max_score >= 2, ])
		write.table(iv, paste0(nm, "_intervals.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
		cl <- merge(loci[, c("CHROM", "POS")], assign_clusters(loci, iv), by = c("CHROM", "POS"))
		write.table(cl, paste0(nm, "_cluster_SNPs.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
		res[[cc]] <- list(intervals = iv, clusters = cl)
		msg("  %-16s %4d windows tested, %4d significant, %4d clusters, %s SNPs in clusters", cc, nrow(tt),
			sum(tt$fdr < 0.05 & tt$max_score >= 2), nrow(iv), format(nrow(cl), big.mark = ","))
	}
	res
}
pairing <- list()
note_pair <- function(set, score, fdr_col, dp_col, same, note = "") {
	pairing[[length(pairing) + 1]] <<- data.frame(cluster_set = set, score = score, FDR_column = fdr_col,
		dp_column = dp_col, same_contrast_or_cage = same, note = note)
}

##########################################################################################
## 4. All-sample clusters: August (E vs extinct, E vs rescued) and November (E vs rescued)
##########################################################################################
msg("August all-sample clusters ...")
t1 <- merge(T1, PAvE[, c(1, 2, 5)], by = c("CHROM", "POS"))   # E vs PA from the time-point model (more P precision)
t1$EvPA.T1 <- t1$PAvE.T1; t1 <- t1[, -9]
t1_fdr <- t1; for (i in 3:8) t1_fdr[, i] <- p.adjust(as.numeric(t1_fdr[, i]), method = "fdr")
sc1 <- merge(t1_fdr, cbind(sites1, abs(fd_t1)), by = c("CHROM", "POS"))
sc1 <- sc1[order(sc1[, 1], sc1[, 2]), ]
## FDR column f (EvPA, EvSE, EvSP, ...) is paired with difference column f + 6 (EvSE, EvSP, EvPA, ...)
dp_for <- if (FIX) c(EvSE.T1 = "EvSE.T1.diff", EvSP.T1 = "EvSP.T1.diff") else c(EvSE.T1 = names(sc1)[4 + 6], EvSP.T1 = names(sc1)[5 + 6])
sc1s <- sc1[, 1:2]
for (cc in names(dp_for)) {
	sc1s[[cc]] <- score_of(sc1[[cc]], sc1[[dp_for[cc]]])
	note_pair("August all-sample", cc, paste0(cc, ".fdr"), dp_for[cc], dp_for[cc] == paste0(cc, ".diff"))
}
cl1 <- make_clusters(sc1s, names(dp_for), t1, "Aug")

msg("November all-sample clusters ...")
sc4 <- merge(tpt4[, c("CHROM", "POS", "EvSP.T4.fdr")], cbind(sites4, abs(fd_t4)), by = c("CHROM", "POS"))
sc4 <- sc4[order(sc4[, 1], sc4[, 2]), ]
sc4$EvSP.T4 <- score_of(sc4$EvSP.T4.fdr, sc4$EvSP.T4.diff, six_on_score = !FIX)
note_pair("November all-sample", "EvSP.T4", "EvSP.T4.fdr", "EvSP.T4.diff", TRUE,
	if (FIX) "" else "score 6 = FDR < 1e-5 without the |delta p| > 0.02 condition")
cl4 <- make_clusters(sc4[, c("CHROM", "POS", "EvSP.T4")], "EvSP.T4", tpt4, "Nov")
## make_clusters() set the seed for the November shuffle; the cross-cluster tests follow it
## directly, so the matched-SNP draws continue from that state.

##########################################################################################
## 5. Concordance t-tests                                                              ##
##########################################################################################
## fdr_tab: CHROM, POS, fdr. clusters: CHROM, POS, clust. sign: CHROM, POS, AFdiff columns.
## cols: sign columns tested (each tested against the same focal/matched SNP draws).
concordance <- function(fdr_tab, clusters, sign, cols, cages, label) {
	vep_temp <- merge(fdr_tab, vep_priority, by = c("CHROM", "POS"))
	vep_temp2 <- merge(vep_temp, clusters[, c("CHROM", "POS", "clust")], by = c("CHROM", "POS"))
	vep_sig <- c()
	for (c in unique(vep_temp2$clust)) {   # top SNP per cluster (lowest FDR; ties broken at random)
		vep_sig <- rbind(vep_sig, sample_n(unique(vep_temp2[vep_temp2$fdr == min(vep_temp2[vep_temp2$clust == c, ]$fdr) & vep_temp2$clust == c, ]), 1))
	}
	if (is.null(vep_sig)) { cat("  ", label, ": no clusters - skipped\n"); return(NULL) }
	vep_sig <- merge(sign[, c("CHROM", "POS", cols)], vep_sig, by = c("CHROM", "POS"))
	vep_sig <- unique(vep_sig[vep_sig$fdr < 0.05, ])
	vep_nonsig <- unique(vep_temp[vep_temp$fdr > 0.05, ])
	vep_nonsig <- vep_nonsig[!vep_nonsig$POS %in% (vep_temp2$POS), ]
	## matched-SNP pool of each focal SNP: same consequence and chromosome, unexposed founder
	## frequency within 25%, > 50 kb away (row order kept)
	idx <- split(seq_len(nrow(vep_nonsig)), paste(vep_nonsig$CHROM, vep_nonsig$Consequence, sep = "\r"))
	pool <- lapply(seq_len(nrow(vep_sig)), function(j) {
		ii <- idx[[paste(vep_sig$CHROM[j], vep_sig$Consequence[j], sep = "\r")]]; if (is.null(ii)) return(integer(0))
		found <- vep_sig$E[j]
		ii[which((vep_nonsig$E[ii] < found * 1.25) & (vep_nonsig$E[ii] > found * .75) & abs(vep_nonsig$POS[ii] - vep_sig$POS[j]) > 50000)]
	})
	use <- which(lengths(pool) > 4)        # focal SNPs with at least 5 possible matches
	if (!length(use)) { cat("  ", label, ": no focal SNP has 5 matches - skipped\n"); return(NULL) }
	## sign-table rows of the background SNPs (merge() order: sorted by the CHROM/POS key)
	skey <- paste(sign$CHROM, sign$POS, sep = "\r"); bkey <- paste(vep_nonsig$CHROM, vep_nonsig$POS, sep = "\r")
	srow <- match(bkey, skey)
	res <- list()
	for (z in seq_len(N_ITER)) {
		b <- vapply(use, function(j) draw1(pool[[j]]), 1L)
		b <- b[!is.na(srow[b])]; b <- b[sort.list(bkey[b])]
		for (x in seq_along(cols)) {
			bg <- sign[srow[b], cols[x]]
			tt <- tryCatch(t.test(vep_sig[[cols[x]]], bg, alternative = "greater"), error = function(e) {
				cat("  note:", label, cols[x], "iteration", z, "t-test failed:", conditionMessage(e), "\n")
				list(statistic = NA_real_, parameter = NA_real_, p.value = NA_real_, stderr = NA_real_) })
			res[[length(res) + 1]] <- data.frame(cage = paste0("no", cages[x]), clusters = label, iteration = z,
				mean_focal = mean(vep_sig[[cols[x]]]), mean_BG = mean(bg), t = unname(tt$statistic), df = unname(tt$parameter),
				p.value = tt$p.value, stderr = tt$stderr)
		}
	}
	msg("  %-22s %3d clusters, %3d focal SNPs (FDR < 0.05), %3d with >= 5 matches", label, length(unique(vep_temp2$clust)), nrow(vep_sig), length(use))
	d <- do.call(rbind, res)
	d[order(as.integer(sub("^no", "", d$cage)), d$iteration), ]   # cage-number order
}
f_T1 <- function(col) cbind(T1[, 1:2], fdr = p.adjust(T1[, col], method = "fdr"))
f_T4 <- function() cbind(tpt4[, 1:2], fdr = p.adjust(tpt4[, "EvSP.T4"], method = "fdr"))
seCols <- paste0("AFdiff_no", SE_loo); spCols <- paste0("AFdiff_no", SP_loo)
cross <- list(
	SET1clust_EvSET1 = list(f_T1(4), cl1$EvSE.T1$clusters, sign_EvSET1, seCols, SE_loo, "SET1clust"),
	SPT1clust_EvSET1 = list(f_T1(5), cl1$EvSP.T1$clusters, sign_EvSET1, seCols, SE_loo, "SPT1clust"),
	SPT4clust_EvSET1 = list(f_T4(),  cl4$EvSP.T4$clusters, sign_EvSET1, seCols, SE_loo, "SPT4clust"),
	SET1clust_EvSPT1 = list(f_T1(4), cl1$EvSE.T1$clusters, sign_EvSPT1, spCols, SP_loo, "SET1clust"),
	SPT1clust_EvSPT1 = list(f_T1(5), cl1$EvSP.T1$clusters, sign_EvSPT1, spCols, SP_loo, "SPT1clust"),
	SPT4clust_EvSPT1 = list(f_T4(),  cl4$EvSP.T4$clusters, sign_EvSPT1, spCols, SP_loo, "SPT4clust"),
	SET1clust_EvSPT4 = list(f_T1(4), cl1$EvSE.T1$clusters, sign_EvSPT4, spCols, SP_loo, "SET1clust"),
	SPT1clust_EvSPT4 = list(f_T1(5), cl1$EvSP.T1$clusters, sign_EvSPT4, spCols, SP_loo, "SPT1clust"),
	SPT4clust_EvSPT4 = list(f_T4(),  cl4$EvSP.T4$clusters, sign_EvSPT4, spCols, SP_loo, "SPT4clust"))
msg("Cross-cluster concordance t-tests (9 x 5 populations x %d iterations) ...", N_ITER)
tests <- list()
for (nm in names(cross)) {
	a <- cross[[nm]]
	tests[[nm]] <- concordance(a[[1]], a[[2]], a[[3]], a[[4]], a[[5]], a[[6]])
	if (!is.null(tests[[nm]])) write.table(tests[[nm]], paste0("ttests/ttest_", nm, "_topsig_results.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
}

##########################################################################################
## 6. Leave-one-out clusters and their concordance t-tests                             ##
##########################################################################################
## Leave-one-out FDR columns are in cage-number order; the allele-frequency difference
## columns in the sample order of the cages (SE_loo / SP_loo). They are paired by
## position; RF_FIX_SCORE_PAIRING pairs them by cage.
msg("August leave-one-out clusters ...")
dpL <- sites1[order(sites1$CHROM, sites1$POS), 1:2]   # loci of a sorted bed, values unsorted
for (i in SE_loo) dpL[[paste0("EvSE.T1_no", i)]] <- abs(no1[[i]][, 1])
for (i in SP_loo) dpL[[paste0("EvSP.T1_no", i)]] <- abs(no1[[i]][, 2])
scL <- merge(loo1_fdr, dpL, by = c("CHROM", "POS"))
scL <- scL[order(scL[, 1], scL[, 2]), ]
lcolsL <- c()
for (k in seq_along(loo1_cols)) {
	fc <- paste0(loo1_cols[k], ".fdr")
	dc <- if (FIX) loo1_cols[k] else names(dpL)[k + 2]
	nm <- paste0(loo1_cols[k], ".score"); scL[[nm]] <- score_of(scL[[fc]], scL[[dc]]); lcolsL <- c(lcolsL, nm)
	note_pair("August leave-one-out", nm, fc, paste0(dc, " (|delta p| without that cage)"), sub("\\.fdr$", "", fc) == dc)
}
clL1 <- make_clusters(scL[, c("CHROM", "POS", lcolsL)], lcolsL, loo1, "AugLOO")

msg("November leave-one-out clusters ...")
dpL4 <- sites4[order(sites4$CHROM, sites4$POS), 1:2]
for (i in SP_loo) dpL4[[paste0("EvSP.T4_no", i)]] <- abs(no4[[i]][, 1])
scL4 <- merge(tpt4[, c("CHROM", "POS", paste0("EvSP.T4_no", SP_CAGES, ".fdr"))], dpL4, by = c("CHROM", "POS"))
scL4 <- scL4[order(scL4[, 1], scL4[, 2]), ]
lcols4 <- c()
for (k in seq_along(SP_CAGES)) {
	fc <- paste0("EvSP.T4_no", SP_CAGES[k], ".fdr")
	dc <- if (FIX) paste0("EvSP.T4_no", SP_CAGES[k]) else names(dpL4)[k + 2]
	nm <- paste0("EvSP.T4_no", SP_CAGES[k], ".score"); scL4[[nm]] <- score_of(scL4[[fc]], scL4[[dc]]); lcols4 <- c(lcols4, nm)
	note_pair("November leave-one-out", nm, fc, paste0(dc, " (|delta p| without that cage)"), sub("\\.fdr$", "", fc) == dc)
}
clL4 <- make_clusters(scL4[, c("CHROM", "POS", lcols4)], lcols4, tpt4, "NovLOO")
## the leave-one-out t-tests follow the November leave-one-out shuffle directly

msg("Leave-one-out concordance t-tests (3 x 5 populations x %d iterations) ...", N_ITER)
loo_test <- function(fdr_cols, fdr_src, clus, cages, sign, label) {
	out <- list()
	for (i in seq_along(cages)) {
		## FDR column i of the leave-one-out table; clusters of cages[i]; sign column i
		fc <- if (FIX) grep(paste0("_no", cages[i], "\\.fdr$"), fdr_cols, value = TRUE) else fdr_cols[i]
		cc <- grep(paste0("_no", cages[i], "\\.score$"), names(clus), value = TRUE)
		if (fc != sub("\\.score$", ".fdr", cc)) cat("  note:", label, "left-out cage", cages[i], "uses FDR column", fc, "with clusters", cc, "(paired by position)\n")
		out[[i]] <- concordance(cbind(fdr_src[, 1:2], fdr = fdr_src[[fc]]), clus[[cc]]$clusters,
			sign[, c("CHROM", "POS", paste0("AFdiff_no", cages[i]))], paste0("AFdiff_no", cages[i]), cages[i], label)
	}
	d <- do.call(rbind, out)
	if (!is.null(d)) d[order(as.integer(sub("^no", "", d$cage)), d$iteration), ] else NULL
}
tests$LOO_SET1clust_EvSET1 <- loo_test(paste0("EvSE.T1_no", SE_CAGES, ".fdr"), loo1_fdr, clL1, SE_loo, sign_EvSET1, "SET1clust")
tests$LOO_SPT1clust_EvSPT1 <- loo_test(paste0("EvSP.T1_no", SP_CAGES, ".fdr"), loo1_fdr, clL1, SP_loo, sign_EvSPT1, "SPT1clust")
tests$LOO_SPT4clust_EvSPT4 <- loo_test(paste0("EvSP.T4_no", SP_CAGES, ".fdr"), tpt4, clL4, SP_loo, sign_EvSPT4, "SPT4clust")
for (nm in grep("^LOO_", names(tests), value = TRUE)) if (!is.null(tests[[nm]])) write.table(tests[[nm]], paste0("ttests/ttest_", nm, "_topsig_results.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

pairing <- do.call(rbind, pairing)
write.csv(pairing, "score_pairing.csv", row.names = FALSE)
cat("\nCluster scores: which FDR column was paired with which |delta p| column (FALSE = different contrast or cage)\n")
for (i in seq_len(nrow(pairing))) cat(sprintf("  %-24s %-20s FDR %-18s |delta p| %-14s %-5s %s\n", pairing$cluster_set[i], pairing$score[i],
	pairing$FDR_column[i], sub(" .*", "", pairing$dp_column[i]), pairing$same_contrast_or_cage[i], pairing$note[i]))

##########################################################################################
## 7. Tables S3, S4, S6 and cluster counts                                             ##
##########################################################################################
fa <- function(x) sprintf("%.4f", x)
ft2 <- function(x) sprintf("%.2f", x)
fp <- function(x) ifelse(x < 0.001, sub("e-0?", "E-", sprintf("%.2e", x)), formatC(signif(x, 3), format = "fg", digits = 3, flag = "#"))
fs <- function(x) ifelse(abs(x - round(x)) < 1e-9, sprintf("%d", as.integer(round(x))), sprintf("%.1f", x))
trow <- function(lab, cl, d) data.frame(`Left-out cage` = lab, clusters = cl,
	`median focal AF diff` = fa(median(d$mean_focal)), `median background AF diff` = fa(median(d$mean_BG, na.rm = TRUE)),
	`median t-test statistic` = ft2(median(d$t)), `median t-test P-value` = fp(median(d$p.value)),
	`perc. sig` = fs(100 * mean(d$p.value < 0.05)), check.names = FALSE)
build <- function(test, cages, sets) {
	out <- NULL
	for (cl in names(sets)) {
		d <- tests[[sets[[cl]]]]; if (is.null(d)) next
		out <- rbind(out, trow(paste0("Across all ", test, " cages"), cl, d))
		for (cg in cages) out <- rbind(out, trow(paste0("Cage ", cg, " (", test, ")"), cl, d[d$cage == paste0("no", cg), ]))
	}
	out
}
tabs <- list(
	S3 = build("SE T1", SE_CAGES, list(`SE T1` = "LOO_SET1clust_EvSET1", `SP T1` = "SPT1clust_EvSET1", `SP T4` = "SPT4clust_EvSET1")),
	S4 = build("SP T1", SP_CAGES, list(`SE T1` = "SET1clust_EvSPT1", `SP T1` = "LOO_SPT1clust_EvSPT1", `SP T4` = "SPT4clust_EvSPT1")),
	S6 = build("SP T4", SP_CAGES, list(`SE T1` = "SET1clust_EvSPT4", `SP T1` = "SPT1clust_EvSPT4", `SP T4` = "LOO_SPT4clust_EvSPT4")))
obs <- c()
cat("\n===== Tables S3, S4, S6 compared with data/published_tables =====\n")
for (n in names(tabs)) {
	write.csv(tabs[[n]], paste0("Table", n, ".csv"), row.names = FALSE)
	pf <- file.path(RF_HOME, "data", "published_tables", paste0("Table", n, ".csv"))
	if (file.exists(pf)) {
		pub <- read.csv(pf, check.names = FALSE, colClasses = "character")
		same <- identical(dim(pub), dim(tabs[[n]])) && all(as.matrix(pub) == as.matrix(tabs[[n]]))
		cat(sprintf("  Table %s: %s\n", n, if (same) "identical to the published table" else "DIFFERS from the published table"))
		if (!same && identical(dim(pub), dim(tabs[[n]]))) {
			w <- which(as.matrix(pub) != as.matrix(tabs[[n]]), arr.ind = TRUE)
			for (k in seq_len(min(nrow(w), 20))) cat(sprintf("    row %d (%s, %s) %s: published %s, now %s\n", w[k, 1], pub[w[k, 1], 1], pub[w[k, 1], 2],
				names(pub)[w[k, 2]], pub[w[k, 1], w[k, 2]], tabs[[n]][w[k, 1], w[k, 2]]))
		}
		obs[paste0("Table", n, "_matches_published")] <- as.numeric(same)
	}
	for (cl in c("SE T1", "SP T1", "SP T4")) {
		r <- tabs[[n]][grepl("^Across", tabs[[n]][, 1]) & tabs[[n]]$clusters == cl, ]
		if (!nrow(r)) next
		k <- paste0(n, "_", gsub(" ", "", cl), "clust")
		obs[paste0(k, "_median_t")] <- as.numeric(r$`median t-test statistic`)
		obs[paste0(k, "_median_P")] <- as.numeric(sub("E", "e", r$`median t-test P-value`))
	}
}
counts <- data.frame(cluster_set = c("August extinct (all samples)", "August rescued (all samples)", "November rescued (all samples)",
	paste("August extinct, without cage", SE_CAGES), paste("August rescued, without cage", SP_CAGES), paste("November rescued, without cage", SP_CAGES)),
	clusters = c(nrow(cl1$EvSE.T1$intervals), nrow(cl1$EvSP.T1$intervals), nrow(cl4$EvSP.T4$intervals),
		vapply(paste0("EvSE.T1_no", SE_CAGES, ".score"), function(k) nrow(clL1[[k]]$intervals), 1L),
		vapply(paste0("EvSP.T1_no", SP_CAGES, ".score"), function(k) nrow(clL1[[k]]$intervals), 1L),
		vapply(paste0("EvSP.T4_no", SP_CAGES, ".score"), function(k) nrow(clL4[[k]]$intervals), 1L)))
write.csv(counts, "cluster_counts.csv", row.names = FALSE)
cat("\nCluster counts (SI 2.4)\n"); print(counts, row.names = FALSE)
obs["clusters_EvSE.T1"] <- counts$clusters[1]; obs["clusters_EvSP.T1"] <- counts$clusters[2]; obs["clusters_EvSP.T4"] <- counts$clusters[3]
for (g in list(c("LOO_EvSE.T1", "August extinct, without"), c("LOO_EvSP.T1", "August rescued, without"), c("LOO_EvSP.T4", "November rescued, without"))) {
	v <- counts$clusters[startsWith(counts$cluster_set, g[2])]
	obs[paste0(g[1], "_min")] <- min(v); obs[paste0(g[1], "_max")] <- max(v)
}

##########################################################################################
## 8. Fig. S9 — Euclidean distances among August populations at the November           ##
##    convergent candidate SNPs                                                        ##
##########################################################################################
cf <- file.path(PATHS$out_root, "05_convergent_candidates", "convergent_candidate_SNPs_TPT4.txt")
if (file.exists(cf)) {
	cand <- read.table(cf, header = TRUE, colClasses = c(CHROM = "character"))
	af <- merge(unique(cand[, 1:2]), haf$freq, by = c("CHROM", "POS"))
	obs["FigS9_n_snps"] <- nrow(af)
	x <- unique(subset(af, select = as.character(mt$samp)))       # August samples, in sample order
	D <- as.matrix(dist(t(x), method = "euclidean", diag = TRUE, upper = FALSE)); D[D == 0] <- NA
	w <- which(!is.na(D), arr.ind = TRUE)
	dd <- data.frame(row = rownames(D)[w[, 1]], col = colnames(D)[w[, 2]], value = D[w])
	sl <- data.frame(samp = as.character(mt$samp), condition = as.character(mt$condition))
	dd <- merge(merge(dd, sl, by.x = "row", by.y = "samp"), sl, by.x = "col", by.y = "samp")
	pr <- paste(pmin(dd$condition.x, dd$condition.y), pmax(dd$condition.x, dd$condition.y))
	lab <- c("E E" = "EvE", "SE SP" = "SEvSP", "SE SE" = "SEvSE", "SP SP" = "SPvSP", "PA PA" = "PAvPA",
		"E SE" = "EvSE", "E SP" = "EvSP", "E PA" = "EvPA", "PA SE" = "PAvSE", "PA SP" = "PAvSP")
	dd$contrast <- unname(lab[pr])
	s9 <- unique(dd[dd$contrast %in% c("EvE", "SEvSP"), c("contrast", "value")])
	write.csv(s9, "FigS9_distances.csv", row.names = FALSE)
	summ <- do.call(rbind, lapply(split(s9$value, s9$contrast), function(v) data.frame(n_pairs = length(v), median = median(v), min = min(v), max = max(v))))
	write.csv(cbind(contrast = rownames(summ), summ), "FigS9_summary.csv", row.names = FALSE)
	cat("\nFig. S9 (", nrow(af), " SNPs): distances among August populations\n", sep = ""); print(summ)
	set.seed(SEED)   # reproducible jitter
	pdf("FigS9_euclidean_distances_jitter.pdf")
	print(ggplot(s9, aes(x = contrast, y = value)) + geom_jitter(aes(colour = contrast)) +
		ggtitle("Pairwise Euclidean distances: T4 candidate loci (AF) among T1 populations") + theme_classic())
	invisible(dev.off())
} else cat("\nFig. S9 skipped: run 05_convergent_candidates.R first (", cf, " not found)\n", sep = "")

ok <- check_expected("07_concordance_LOO", obs)
sink()
if (!ok) quit(status = 1)
