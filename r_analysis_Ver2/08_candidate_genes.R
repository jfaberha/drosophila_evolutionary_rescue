########################################

###  Code for: Contingency and determinism in rapid adaptation shape evolutionary rescue and extinction in Drosophila
###  Authors: René S. Shahmohamadloo*, Joshua J. Faber-Hammond, Amir R. Gabidulin, Catherine I. Clare, Seth M. Rudman*
###  *Corresponding authors:
###  René S. Shahmohamadloo, School of Biological Sciences, Washington State University, 14204 NE Salmon Creek Ave, Vancouver, WA 98686, USA. Email: rene.shahmohamadloo@wsu.edu
###  Seth M. Rudman, School of Biological Sciences, Washington State University, 14204 NE Salmon Creek Ave, Vancouver, WA 98686, USA. Email: seth.rudman@wsu.edu

###  In this file: 08_candidate_genes.R — are the a priori spinosad candidate genes (Table S2) enriched among GLM outliers?
###  Outputs:
###  SI 1.8.6   SNP counts in candidate-gene regions
###  SI 2.1.4   hypergeometric overlap tests (all time points; all / missense SNPs; FDR < 0.05 / 0.01) and Wilcoxon rank-sum tests vs matched genes      ###
###  SI 2.1.4   the same for the August contrasts, and unsmoothed counts (Fig. S4) (the 21-SNP-window ['smoothed'] and unsmoothed versions of the hypergeometric tests; Wilcoxon on the smoothed values, as in the text)

###  Each gene's minimum FDR (and P) across its SNPs is taken numerically and compared with the threshold.

########################################

f_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
RF_HOME <- if (length(f_arg)) normalizePath(dirname(sub("^--file=", "", f_arg))) else getwd()
source(file.path(RF_HOME, "00_config.R"))
OUT <- module_dir("08_candidate_genes"); setwd(OUT)
sink("08_candidate_genes_report.txt", split = TRUE)
N_ITER <- 1000

glm.all <- load_rds("glm_all"); sm <- load_rds("glm_all_rolwin21")
vep <- load_vep(); spino <- load_spino()
gff <- read.delim(file.path(PATHS$raw_input, FILES$genic_bed), sep = "\t", header = TRUE)
msg("Candidate genes listed: %d | genic-feature table: %s rows", nrow(spino), format(nrow(gff), big.mark = ","))
gs <- unique(vep[, c("Gene", "CHROM", "POS", "Consequence")]); rm(vep); invisible(gc())

## Gene-SNP links for one GLM table (SNPs present in the table)
links <- function(tab) {
	L <- merge(gs, tab[, c("CHROM", "POS")], by = c("CHROM", "POS"))
	L$idx <- match(paste(L$CHROM, L$POS), paste(tab$CHROM, tab$POS))
	L
}
## Minimum value per gene across its SNPs
gene_min <- function(L, tab, cols) {
	g <- sort(unique(L$Gene))
	out <- data.frame(Gene = g, stringsAsFactors = FALSE)
	for (cc in cols) out[[cc]] <- unname(tapply(tab[[cc]][L$idx], L$Gene, min)[g])
	out
}

## Hypergeometric overlap (candidate genes vs genes with an outlier SNP), every contrast
hyper <- function(cand, bg, thr, bh) {
	cols <- names(cand)[-1]
	r <- do.call(rbind, lapply(cols, function(cc) {
		q <- length(unique(cand$Gene[which(cand[[cc]] < thr)]))
		m <- length(unique(cand$Gene)); n <- length(unique(bg$Gene)) - m
		k <- length(unique(bg$Gene[which(bg[[cc]] < thr)]))
		data.frame(contrast = cc, q = q, m = m, k = k, n = n, P = phyper(q - 1, m, n, k, lower.tail = FALSE))
	}))
	r$P_BH <- if (bh) p.adjust(r$P, method = "fdr") else NA   # BH across contrasts (smoothed tables)
	r
}

##########################################################################################
## 1. Hypergeometric tests: smoothed (21-SNP windows) and unsmoothed                  ##
##########################################################################################
hyp <- list(); sets <- list()
for (smooth in c(TRUE, FALSE)) {
	tab <- if (smooth) sm else glm.all
	lab <- if (smooth) "smoothed" else "unsmoothed"
	fcols <- grep("fdr", names(tab), value = TRUE)
	L <- links(tab); Ls <- L[L$Gene %in% spino$Gene, ]
	sets[[lab]] <- list(L = L, Ls = Ls, tab = tab)
	msg("%s: %s gene-SNP links; %d candidate genes with SNPs", lab, format(nrow(L), big.mark = ","), length(unique(Ls$Gene)))
	cand <- gene_min(Ls, tab, fcols); bg <- gene_min(L, tab, fcols)
	Lm <- L[L$Consequence == "missense_variant", ]; Lsm <- Ls[Ls$Consequence == "missense_variant", ]
	candm <- gene_min(Lsm, tab, fcols); bgm <- gene_min(Lm, tab, fcols)
	write.table(cand, paste0("candidate_gene_minFDR_", lab, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
	for (thr in c(0.05, 0.01)) {
		hyp[[length(hyp) + 1]] <- cbind(data = lab, SNPs = "all", FDR_threshold = thr, hyper(cand, bg, thr, smooth))
		hyp[[length(hyp) + 1]] <- cbind(data = lab, SNPs = "missense", FDR_threshold = thr, hyper(candm, bgm, thr, smooth))
	}
}
hyp <- do.call(rbind, hyp)
hyp$contrast <- sub("\\.fdr(\\.rolwin21)?$", "", hyp$contrast)
write.csv(hyp, "hypergeometric_candidate_genes.csv", row.names = FALSE)

##########################################################################################
## 2. Counts (SI 1.8.6; SI 2.1.4 unsmoothed paragraph)                                 ##
##########################################################################################
Ls <- sets$smoothed$Ls
flank <- c("upstream_gene_variant", "downstream_gene_variant")
counts <- c(SI186_candidate_SNPs = nrow(unique(Ls[, c("CHROM", "POS")])),
	SI186_candidate_SNPs_noflank = nrow(unique(Ls[!Ls$Consequence %in% flank, c("CHROM", "POS")])),
	SI186_candidate_missense_SNPs = nrow(unique(Ls[Ls$Consequence == "missense_variant", c("CHROM", "POS")])),
	candidate_genes_in_data = length(unique(Ls$Gene)))
Lu <- sets$unsmoothed$Ls; tu <- sets$unsmoothed$tab
for (cc in c("EvSE.T1", "EvSP.T1", "EvPA.T1", "PAvE.T1")) {
	s <- Lu[which(tu[[paste0(cc, ".fdr")]][Lu$idx] < 0.05), ]
	counts[paste0("SI214_unsmoothed_", cc, "_SNPs")] <- nrow(unique(s[, c("CHROM", "POS")]))
	counts[paste0("SI214_unsmoothed_", cc, "_genes")] <- length(unique(s$Gene))
}
counts["SI214_unsmoothed_candidate_genes_in_data"] <- length(unique(Lu$Gene))
cat("\nCounts\n"); print(data.frame(value = counts))

##########################################################################################
## 3. Wilcoxon rank-sum tests: candidate genes vs matched genes (smoothed P values)    ##
##    Matched gene: same type and chromosome, length and SNP count within 25%.          ##
##########################################################################################
tab <- sm; L <- sets$smoothed$L; Ls <- sets$smoothed$Ls
pcols <- grep("rolwin21$", names(tab), value = TRUE); pcols <- pcols[!grepl("fdr|logp", pcols)]
cmin <- gene_min(Ls, tab, pcols); bmin <- gene_min(L, tab, pcols)
cg <- merge(gff, cmin, by.x = "ID", by.y = "Gene"); bg <- merge(gff, bmin, by.x = "ID", by.y = "Gene")
bg <- bg[!bg$ID %in% cg$ID, ]
pool <- lapply(seq_len(nrow(cg)), function(i) which(bg$TYPE == cg$TYPE[i] & bg$CHROM == cg$CHROM[i] &
	bg$LENGTH < cg$LENGTH[i] * 1.25 & bg$LENGTH > cg$LENGTH[i] * 0.75 &
	bg$SNP_COUNT < cg$SNP_COUNT[i] * 1.25 & bg$SNP_COUNT > cg$SNP_COUNT[i] * 0.75))
if (any(lengths(pool) == 0)) cat("  note:", sum(lengths(pool) == 0), "candidate gene(s) have no matched gene and are skipped\n")
msg("Wilcoxon: %d candidate genes with feature information, %s background genes, %d iterations", nrow(cg), format(nrow(bg), big.mark = ","), N_ITER)
set.seed(SEED)
res <- matrix(NA_real_, N_ITER, length(pcols), dimnames = list(NULL, pcols))
for (it in seq_len(N_ITER)) {
	b <- unlist(lapply(pool, function(ii) if (length(ii)) ii[sample.int(length(ii), 1)]))   # one matched gene per candidate
	for (j in seq_along(pcols)) {
		cand <- as.numeric(cg[[pcols[j]]]); bs <- as.numeric(bg[[pcols[j]]][b])
		if (length(cand) > length(bs)) cand <- sample(cand, length(bs)) else bs <- sample(bs, length(cand))
		res[it, j] <- suppressWarnings(wilcox.test(cand, bs, alternative = "less")$p.value)
	}
}
write.table(res, "wilcoxon_P_by_iteration.txt", sep = "\t", quote = FALSE, row.names = FALSE)
wil <- data.frame(contrast = sub("\\.rolwin21$", "", pcols), mean_P = colMeans(res),
	median_P = apply(res, 2, median), percent_P_below_0.05 = 100 * colMeans(res < 0.05, na.rm = TRUE))
rownames(wil) <- NULL
write.csv(wil, "wilcoxon_candidate_genes.csv", row.names = FALSE)

##########################################################################################
## 4. Values cited in SI 1.8.6, 2.1.3 and 2.1.4                                         ##
##########################################################################################
hv <- function(d, cc, snps, thr, what) { x <- hyp[hyp$data == d & hyp$contrast == cc & hyp$SNPs == snps & hyp$FDR_threshold == thr, ]; if (nrow(x)) x[[what]] else NA }
wv <- function(cc, what) { x <- wil[wil$contrast == cc, ]; if (nrow(x)) x[[what]] else NA }
low_aug <- function(d) min(hyp$P[hyp$data == d & hyp$SNPs == "all" & hyp$contrast %in% c("EvSE.T1", "EvSP.T1", "EvPA.T1", "PAvE.T1")])
items <- list(
	list("SI213_PAvE_FDR01_genes",     "2.1.3 pre-adapted vs unexposed, FDR < 0.01: candidate genes", hv("smoothed", "PAvE", "all", 0.01, "q")),
	list("SI213_PAvE_FDR01_P",         "2.1.3 pre-adapted vs unexposed, FDR < 0.01: uncorrected P", hv("smoothed", "PAvE", "all", 0.01, "P")),
	list("SI213_SvE_missense_FDR05_genes", "2.1.3 exposed vs unexposed missense, FDR < 0.05: candidate genes", hv("smoothed", "SvE", "missense", 0.05, "q")),
	list("SI213_SvE_missense_FDR05_P", "2.1.3 exposed vs unexposed missense, FDR < 0.05: uncorrected P", hv("smoothed", "SvE", "missense", 0.05, "P")),
	list("SI213_SvE_FDR01_genes",      "2.1.3 exposed vs unexposed (all time points), FDR < 0.01: candidate genes", hv("smoothed", "SvE", "all", 0.01, "q")),
	list("SI213_SvE_FDR05_genes",      "2.1.3 exposed vs unexposed (all time points), FDR < 0.05: candidate genes", hv("smoothed", "SvE", "all", 0.05, "q")),
	list("SI213_SvE_T1toT3_FDR05_genes", "2.1.3 exposed vs unexposed Aug-Oct, FDR < 0.05: candidate genes (sum)", sum(sapply(c("SvE.T1", "SvE.T2", "SvE.T3"), function(cc) hv("smoothed", cc, "all", 0.05, "q")))),
	list("SI213_SvE_T4_FDR05_genes",   "2.1.3 exposed vs unexposed November, FDR < 0.05: candidate genes", hv("smoothed", "SvE.T4", "all", 0.05, "q")),
	list("SI213_wilcox_mean_PAvE",     "2.1.3 Wilcoxon mean P, pre-adapted vs unexposed", wv("PAvE", "mean_P")),
	list("SI213_wilcox_mean_SvE",      "2.1.3 Wilcoxon mean P, exposed vs unexposed", wv("SvE", "mean_P")),
	list("SI213_wilcox_pctsig_PAvE",   "2.1.3 Wilcoxon % iterations P < 0.05, pre-adapted vs unexposed", wv("PAvE", "percent_P_below_0.05")),
	list("SI213_wilcox_pctsig_SvE",    "2.1.3 Wilcoxon % iterations P < 0.05, exposed vs unexposed", wv("SvE", "percent_P_below_0.05")),
	list("SI213_unsmoothed_SvE_FDR01_genes", "2.1.3 unsmoothed exposed vs unexposed, FDR < 0.01: candidate genes", hv("unsmoothed", "SvE", "all", 0.01, "q")),
	list("SI213_unsmoothed_SvE_FDR01_P",     "2.1.3 unsmoothed exposed vs unexposed, FDR < 0.01: P", hv("unsmoothed", "SvE", "all", 0.01, "P")),
	list("SI213_unsmoothed_PAvE_FDR01_genes", "2.1.3 unsmoothed pre-adapted vs unexposed, FDR < 0.01: candidate genes", hv("unsmoothed", "PAvE", "all", 0.01, "q")),
	list("SI213_unsmoothed_PAvE_FDR01_P",    "2.1.3 unsmoothed pre-adapted vs unexposed, FDR < 0.01: P", hv("unsmoothed", "PAvE", "all", 0.01, "P")),
	list("SI214_lowest_P_August",      "2.1.4 lowest uncorrected hypergeometric P, August contrasts", low_aug("smoothed")),
	list("SI214_wilcox_mean_EvSE.T1",  "2.1.4 Wilcoxon mean P, extinct vs unexposed (Aug)", wv("EvSE.T1", "mean_P")),
	list("SI214_wilcox_mean_EvSP.T1",  "2.1.4 Wilcoxon mean P, rescued vs unexposed (Aug)", wv("EvSP.T1", "mean_P")),
	list("SI214_wilcox_mean_EvPA.T1",  "2.1.4 Wilcoxon mean P, pre-adapted vs unexposed (Aug; four-group model)", wv("EvPA.T1", "mean_P")),
	list("SI214_wilcox_max_pctsig_August", "2.1.4 Wilcoxon % iterations P < 0.05, largest of the August contrasts", max(sapply(c("EvSE.T1", "EvSP.T1", "PAvE.T1", "EvPA.T1"), function(cc) wv(cc, "percent_P_below_0.05")))),
	list("SI214_unsmoothed_lowest_P_August", "2.1.4 unsmoothed lowest uncorrected P, August contrasts", low_aug("unsmoothed")))
summ <- do.call(rbind, lapply(items, function(it) data.frame(key = it[[1]], item = it[[2]], value = it[[3]])))
summ <- rbind(summ, data.frame(key = names(counts), item = "count", value = unname(counts)))
write.csv(summ, "SI_values.csv", row.names = FALSE)
cat("\nValues cited in the SI\n")
for (i in seq_len(nrow(summ))) cat(sprintf("  %-38s %12s   %s\n", summ$key[i], format(summ$value[i], digits = 3, big.mark = ","), summ$item[i]))

ok <- check_expected("08_candidate_genes", setNames(summ$value, summ$key))
sink()
if (!ok) quit(status = 1)
