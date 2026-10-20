########################################

###  Code for: Contingency and determinism in rapid adaptation shape evolutionary rescue and extinction in Drosophila
###  Authors: René S. Shahmohamadloo*, Joshua J. Faber-Hammond, Amir R. Gabidulin, Catherine I. Clare, Seth M. Rudman*
###  *Corresponding authors:
###  René S. Shahmohamadloo, School of Biological Sciences, Washington State University, 14204 NE Salmon Creek Ave, Vancouver, WA 98686, USA. Email: rene.shahmohamadloo@wsu.edu
###  Seth M. Rudman, School of Biological Sciences, Washington State University, 14204 NE Salmon Creek Ave, Vancouver, WA 98686, USA. Email: seth.rudman@wsu.edu

###  In this file: 03_genome_scans.R — genome-wide divergence scans and allele-frequency change at outlier and candidate SNPs.

###  Produces:
###  Fig. 2D-E  smoothed Manhattans, all time points (exposed; pre-adapted)
###  Fig. S3    the same without smoothing
###  Fig. 2F-G  |Δp| at spinosad candidate SNPs and 'exposed vs unexposed' outliers vs matched background SNPs (Welch tests: main text; SI 2.1.3)
###  Fig. 3D-F  smoothed August Manhattans (extinct; rescued; pre-adapted)
###  Fig. S8    the same without smoothing
###  Fig. S4    unsmoothed August Manhattans highlighting candidate-gene SNPs
###  Fig. 3G-I  August |Δp| at outliers vs matched SNPs (main text; SI 2.1.4)
###  SI 2.1.2   hypergeometric overlap of exposed and pre-adapted outliers

########################################

f_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
RF_HOME <- if (length(f_arg)) normalizePath(dirname(sub("^--file=", "", f_arg))) else getwd()
source(file.path(RF_HOME, "00_config.R"))
suppressMessages({ library(dplyr); library(ggplot2) })
OUT <- module_dir("03_genome_scans"); setwd(OUT)
sink("03_genome_scans_report.txt", split = TRUE)

glm.all <- load_rds("glm_all")
glm.all.rolwin21 <- load_rds("glm_all_rolwin21")
freq_means <- load_rds("freq_means")
freq_means_t1 <- load_rds("freq_means_t1")
vep <- load_vep(); spino.cand <- load_spino()
msg("Loaded: %s SNPs (unsmoothed), %s (smoothed); VEP rows %s; candidate genes %d",
	format(nrow(glm.all), big.mark = ","), format(nrow(glm.all.rolwin21), big.mark = ","), format(nrow(vep), big.mark = ","), nrow(spino.cand))

## SNPs in and around (±5 kb, as assigned by VEP) spinosad candidate genes: every VEP row
vep_spino <- vep[vep$Gene %in% spino.cand$Gene, ]
glm.all.annot.spino <- merge(spino.cand, merge(glm.all, vep_spino, by = c("CHROM", "POS")), by = "Gene")
glm.all.rolwin21.annot.spino <- merge(spino.cand, merge(glm.all.rolwin21, vep_spino, by = c("CHROM", "POS")), by = "Gene")
## One VEP row per SNP (the first)
vep_nr <- vep[!duplicated(vep[, c("CHROM", "POS")]), ]
rm(vep); invisible(gc())

## ======================================================================================
## Manhattan plots
## ======================================================================================
## d: table; y: -log10(P) column; layers: list(list(rows, colour, na_omit)); hline: FDR
## column for a dashed line at the smallest -log10(P) among SNPs below `thr`
manh <- function(d, y, title, layers, hline = NULL, thr = 0.05, ymax_add = 0, base_na_omit = FALSE) {
	g <- ggplot(if (base_na_omit) na.omit(d) else d, aes(POS, .data[[y]])) + geom_line(alpha = 1, colour = "#CCCCCC")
	if (!is.null(hline)) g <- g + geom_hline(yintercept = min(d[[y]][which(d[[hline]] < thr)], na.rm = TRUE), color = "black", linetype = "dashed", linewidth = .25)
	for (L in layers) {
		pd <- L$rows; if (isTRUE(L$na)) pd <- na.omit(pd)
		g <- g + geom_point(data = pd, aes(POS, .data[[y]]), color = L$col, size = 1)
	}
	g + facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
		scale_y_continuous(limits = c(0, max(na.omit(d[[y]])) + ymax_add)) +
		scale_x_continuous(breaks = c(0, 5e6, 1e7, 1.5e7, 2e7, 2.5e7, 3e7), guide = guide_axis(angle = 45)) +
		labs(col = "candidate\ngene\n-log10(p)") + xlab("chromosome position") + ylab("-log10(p)") +
		theme_classic() + theme(legend.position = "none") + ggtitle(title)
}
lay <- function(rows, col, na = TRUE) list(rows = rows, col = col, na = na)
w <- function(d, cond) d[which(cond), ]          # rows where cond is TRUE
RED <- "#BD2D23"; BROWN <- "#3F0100"; LBLUE <- "#809BCA"; DBLUE <- "#025197"
RED3 <- "#BF2D23"; PINK <- "#DE9691"; BROWN3 <- "#3F0000"

## ---- Fig. 2D-E: all time points, 21-SNP windows (outliers FDR < 0.01) ----
r <- glm.all.rolwin21; rs <- glm.all.rolwin21.annot.spino
p2d <- manh(r, "SvE.logp.rolwin21", "E vs S", base_na_omit = TRUE, layers = list(
	lay(w(r, r$SvE.fdr.rolwin21 < 0.01), RED), lay(rs, LBLUE), lay(w(rs, rs$SvE.fdr.rolwin21 < 0.01), DBLUE)))
p2e <- manh(r, "PAvE.logp.rolwin21", "E vs PA", base_na_omit = TRUE, layers = list(
	lay(w(r, r$PAvE.fdr.rolwin21 < 0.01), BROWN), lay(w(r, r$SvE.fdr.rolwin21 < 0.01), RED),
	lay(rs, LBLUE, na = FALSE), lay(w(rs, rs$PAvE.fdr.rolwin21 < 0.01), DBLUE)))
save_panels("Fig2DE_manhattan_alltpts_rolwin21.pdf", list(p2d, p2e), width = 7.5, height = 4)

## ---- Fig. S3: the same without smoothing ----
u <- glm.all; us <- glm.all.annot.spino
ps3a <- manh(u, "SvE.logp", "E vs S", base_na_omit = TRUE, hline = "SvE.fdr", thr = 0.01, layers = list(
	lay(w(u, u$SvE.fdr < 0.01), RED), lay(us, LBLUE), lay(w(us, us$SvE.fdr < 0.01), DBLUE)))
ps3b <- manh(u, "PAvE.logp", "E vs PA", base_na_omit = TRUE, hline = "PAvE.fdr", thr = 0.01, layers = list(
	lay(w(u, u$PAvE.fdr < 0.01), BROWN), lay(w(u, u$SvE.fdr < 0.01), RED),
	lay(us, LBLUE, na = FALSE), lay(w(us, us$PAvE.fdr < 0.01), DBLUE)))
save_panels("FigS3_manhattan_alltpts_unsmoothed.pdf", list(ps3a, ps3b), width = 7.5, height = 4)

## ---- Fig. S8: August, unsmoothed (outliers FDR < 0.05) ----
both <- function(d) w(d, d$EvSE.T1.fdr < 0.05 & d$EvSP.T1.fdr < 0.05)
ps8a <- manh(u, "EvSE.T1.logp", "E vs SE (TPT1)", hline = "EvSE.T1.fdr", layers = list(
	lay(w(u, u$EvSE.T1.fdr < 0.05), RED3), lay(w(u, u$EvSP.T1.fdr < 0.05), PINK), lay(both(u), "yellow")))
ps8b <- manh(u, "EvSP.T1.logp", "E vs SP (TPT1)", hline = "EvSP.T1.fdr", ymax_add = 0.5, layers = list(
	lay(w(u, u$EvSP.T1.fdr < 0.05), PINK), lay(w(u, u$EvSE.T1.fdr < 0.05), RED3), lay(both(u), "yellow")))
ps8c <- manh(u, "PAvE.T1.logp", "E vs PA (TPT1)", hline = "PAvE.T1.fdr", layers = list(
	lay(w(u, u$PAvE.T1.fdr < 0.05), BROWN3), lay(w(u, u$EvSE.T1.fdr < 0.05 & u$PAvE.T1.fdr < 0.05), RED3),
	lay(w(u, u$EvSP.T1.fdr < 0.05), PINK), lay(both(u), "yellow")))
save_panels("FigS8_manhattan_August_unsmoothed.pdf", list(ps8a, ps8b, ps8c), width = 7.5, height = 6)

## ---- Fig. S4: August, unsmoothed, candidate-gene SNPs highlighted ----
s4 <- function(y, f, title, add = 0) manh(u, y, title, hline = f, ymax_add = add, layers = list(
	lay(us, LBLUE, na = FALSE), lay(w(us, us[[f]] < 0.05), DBLUE, na = FALSE)))
save_panels("FigS4_manhattan_August_candidate_genes.pdf", list(
	s4("EvSE.T1.logp", "EvSE.T1.fdr", "E vs SE (TPT1)"), s4("EvSP.T1.logp", "EvSP.T1.fdr", "E vs SP (TPT1)", 0.5),
	s4("PAvE.T1.logp", "PAvE.T1.fdr", "E vs PA (TPT1)")), width = 7.5, height = 6)

## ---- Fig. 3D-F: August, 21-SNP windows ----
r3a <- manh(r, "EvSE.T1.logp.rolwin21", "E vs SE (TPT1): 21-SNP sliding window", layers = list(
	lay(w(r, r$EvSE.T1.fdr.rolwin21 < 0.05), RED3)))
r3b <- manh(r, "EvSP.T1.logp.rolwin21", "E vs SP (TPT1): 21-SNP sliding window", ymax_add = 0.5, layers = list(
	lay(w(r, r$EvSE.T1.fdr.rolwin21 < 0.05 & r$EvSP.T1.fdr.rolwin21 > 0.05), PINK),
	lay(w(r, r$EvSE.T1.fdr.rolwin21 < 0.05 & r$EvSP.T1.fdr.rolwin21 < 0.05), RED3)))
r3c <- manh(r, "PAvE.T1.logp.rolwin21", "E vs PA (TPT1): 21-SNP sliding window", layers = list(
	lay(w(r, r$PAvE.T1.fdr.rolwin21 < 0.05), BROWN3),
	lay(w(r, r$EvSE.T1.fdr.rolwin21 < 0.05 & r$PAvE.T1.fdr.rolwin21 > 0.05), PINK),
	lay(w(r, r$EvSE.T1.fdr.rolwin21 < 0.05 & r$PAvE.T1.fdr.rolwin21 < 0.05), RED3)))
save_panels("Fig3DF_manhattan_August_rolwin21.pdf", list(r3a, r3b, r3c), width = 7.5, height = 6)
msg("Manhattan plots written")

## ======================================================================================
## SI 2.1.2: overlap of 'exposed vs unexposed' and 'pre-adapted vs unexposed' outliers
## ======================================================================================
hyp <- function(d, a, b, thr, label) {
	d <- d[!is.na(d[[a]]) & !is.na(d[[b]]), ]
	A <- d[[a]] < thr; B <- d[[b]] < thr; N <- nrow(d); m <- sum(A); k <- sum(B); q <- sum(A & B)
	e <- as.numeric(m) * k / N
	lp <- phyper(q - 1, m, N - m, k, lower.tail = FALSE, log.p = TRUE) / log(10)
	data.frame(test = label, N = N, outliers_A = m, outliers_B = k, overlap = q, pct_of_A = round(100 * q / m, 1),
		expected = round(e, 1), fold = round(q / e, 2), log10_P = round(lp, 1))
}
ov <- rbind(hyp(glm.all.rolwin21, "SvE.fdr.rolwin21", "PAvE.fdr.rolwin21", 0.01, "exposed (A) vs pre-adapted (B), all time points, 21-SNP windows, FDR<0.01"),
	hyp(glm.all, "SvE.fdr", "PAvE.fdr", 0.01, "exposed (A) vs pre-adapted (B), all time points, unsmoothed, FDR<0.01"))
cat("\n===== SI 2.1.2 hypergeometric overlap =====\n"); print(ov, row.names = FALSE)
write.csv(ov, "SI2.1.2_overlap_hypergeometric.csv", row.names = FALSE)

## ======================================================================================
## Allele-frequency change at outlier and candidate SNPs vs matched background SNPs
## ======================================================================================
## For each focal SNP: background SNPs on the same chromosome with the same VEP consequence
## and unexposed mean frequency within ±25% (of min(p, 1-p)); up to 100 sampled; mean and
## variance of |Δp| for each contrast. One row per focal SNP, in focal-SNP order.
## (A focal SNP with no matches gives NaN.)
matched_bg <- function(focal, bg, contrasts) {
	idx <- split(seq_len(nrow(bg)), paste(bg$CHROM, bg$Consequence, sep = "\r"))
	out <- matrix(NA_real_, nrow(focal), 4 * length(contrasts))
	for (g in seq_len(nrow(focal))) {
		e <- focal$E.af[g]; tol <- 0.25 * min(e, 1 - e)
		ii <- idx[[paste(focal$CHROM[g], focal$Consequence[g], sep = "\r")]]
		if (is.null(ii)) ii <- integer(0)
		ii <- ii[bg$E.af[ii] > (e - tol) & bg$E.af[ii] < (e + tol)]
		tb <- bg[ii[sample(length(ii), min(c(length(ii), 100)))], ]
		af <- do.call(cbind, c(lapply(contrasts, function(cc) abs(tb$E.af - tb[[cc]])),
			lapply(contrasts, function(cc) abs((tb$E.af - tb[[cc]]) / tb$E.af))))
		out[g, ] <- c(colMeans(af), colVars(af))
	}
	nm <- c(paste0("bg.E2", sub("\\.af$", "", contrasts)), paste0("bg.E2", sub("\\.af$", "", contrasts), ".perc"))
	colnames(out) <- c(paste0(nm, ".mean"), paste0(nm, ".var"))
	as.data.frame(out)
}
## Welch test: candidate (focal) |Δp| greater than the mean |Δp| of their matched SNPs
welch <- function(panel, focal_set, cand, matched) {
	tt <- t.test(cand, matched, alternative = "greater")   # same as t.test(af_diff ~ set, alternative = 'greater')
	data.frame(panel = panel, focal_SNPs = focal_set, n = length(cand), mean_focal = mean(cand), mean_matched = mean(matched, na.rm = TRUE),
		t = unname(tt$statistic), df = unname(tt$parameter), P = tt$p.value)
}
boxp <- function(cand_list, matched_list, title) {
	d <- rbind(do.call(rbind, lapply(names(cand_list), function(k) data.frame(set = "candidate", category = k, af_diff = cand_list[[k]]))),
		do.call(rbind, lapply(names(matched_list), function(k) data.frame(set = "matched", category = k, af_diff = matched_list[[k]]))))
	ggplot(d, aes(x = set, y = af_diff)) + geom_boxplot(aes(colour = category)) +
		scale_color_manual(values = c("#39568A", "#C45C5C")) + ggtitle(title) + theme_classic()
}
set.seed(SEED)   # one seed before the seven matching loops, which run in this order

## ---- All time points (Fig. 2F-G) ----
fm_nr <- merge(freq_means, vep_nr, by = c("CHROM", "POS"))
fm_spino <- fm_nr[fm_nr$Gene %in% spino.cand$Gene, ]
fm_bg <- fm_nr[!fm_nr$Gene %in% spino.cand$Gene, ]
## all-sample GLM rows with an annotation, in locus order; top N outliers (N = candidate SNPs)
glm_nr <- merge(glm.all, vep_nr[, c("CHROM", "POS", "Uploaded_variation")], by = c("CHROM", "POS"))
topN <- function(fdr_col, N) head(glm_nr[order(glm_nr[[fdr_col]]), ], N)
N2 <- nrow(fm_spino)
top_EvS <- topN("SvE.fdr", N2); top_EvPA <- topN("PAvE.fdr", N2)
fm_EvS <- merge(top_EvS[, 1:2], fm_nr, by = c("CHROM", "POS"))
fm_EvS_bg <- fm_nr[!fm_nr$Uploaded_variation %in% top_EvS$Uploaded_variation, ]
fm_EvPA <- merge(top_EvPA[, 1:2], fm_nr, by = c("CHROM", "POS"))
fm_EvPA_bg <- fm_nr[!fm_nr$Uploaded_variation %in% top_EvPA$Uploaded_variation, ]
msg("All time points: %d candidate-gene SNPs; matching background ...", N2)
c2 <- c("S.af", "PA.af")
res_spino <- matched_bg(fm_spino, fm_bg, c2);     names(res_spino) <- paste0(names(res_spino), ".spino")
res_EvS   <- matched_bg(fm_EvS, fm_EvS_bg, c2);   names(res_EvS)   <- paste0(names(res_EvS), ".EvS")
res_EvPA  <- matched_bg(fm_EvPA, fm_EvPA_bg, c2); names(res_EvPA)  <- paste0(names(res_EvPA), ".EvPA")  # not plotted; run to keep the random draws in order
write.table(res_spino, "freq_means_results.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(res_EvS, "freq_means_results_EvS.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(res_EvPA, "freq_means_results_EvPA.txt", sep = "\t", quote = FALSE, row.names = FALSE)

E2S  <- list(spino = abs(fm_spino$E.af - fm_spino$S.af),  EvS = abs(fm_EvS$E.af - fm_EvS$S.af))
E2PA <- list(spino = abs(fm_spino$E.af - fm_spino$PA.af), EvS = abs(fm_EvS$E.af - fm_EvS$PA.af))
mS   <- list(spino = res_spino$bg.E2S.mean.spino,  EvS = res_EvS$bg.E2S.mean.EvS)
mPA  <- list(spino = res_spino$bg.E2PA.mean.spino, EvS = res_EvS$bg.E2PA.mean.EvS)
save_panels("Fig2FG_dp_boxplots.pdf", list(boxp(E2S, mS, "AF differences: E vs S"), boxp(E2PA, mPA, "AF differences: E vs PA")), width = 3, height = 6)
stats <- rbind(
	welch("Fig. 2F (unexposed vs exposed)", "candidate genes", E2S$spino, mS$spino),
	welch("Fig. 2F (unexposed vs exposed)", "exposed outliers", E2S$EvS, mS$EvS),
	welch("Fig. 2G (unexposed vs pre-adapted)", "candidate genes", E2PA$spino, mPA$spino),
	welch("Fig. 2G (unexposed vs pre-adapted)", "exposed outliers", E2PA$EvS, mPA$EvS))

## ---- August (Fig. 3G-I) ----
f1_nr <- merge(freq_means_t1, vep_nr, by = c("CHROM", "POS"))
f1_spino <- f1_nr[f1_nr$Gene %in% spino.cand$Gene, ]
f1_bg <- f1_nr[!f1_nr$Gene %in% spino.cand$Gene, ]   # background for all August loops
N1 <- nrow(f1_spino)
f1_EvSE <- merge(topN("EvSE.T1.fdr", N1)[, 1:2], f1_nr, by = c("CHROM", "POS"))
f1_EvSP <- merge(topN("EvSP.T1.fdr", N1)[, 1:2], f1_nr, by = c("CHROM", "POS"))
f1_EvPA <- merge(topN("PAvE.T1.fdr", N1)[, 1:2], f1_nr, by = c("CHROM", "POS"))
msg("August: %d candidate-gene SNPs; matching background ...", N1)
c1 <- c("SE.af", "SP.af", "PA.af")
r1_spino <- matched_bg(f1_spino, f1_bg, c1); names(r1_spino) <- paste0(names(r1_spino), ".spino")
r1_EvSE  <- matched_bg(f1_EvSE, f1_bg, c1);  names(r1_EvSE)  <- paste0(names(r1_EvSE), ".EvSE.T1")
r1_EvSP  <- matched_bg(f1_EvSP, f1_bg, c1);  names(r1_EvSP)  <- paste0(names(r1_EvSP), ".EvSP.T1")  # not plotted; keeps the draw order
r1_EvPA  <- matched_bg(f1_EvPA, f1_bg, c1);  names(r1_EvPA)  <- paste0(names(r1_EvPA), ".EvPA.T1")
write.table(r1_spino, "freq_means_t1_results.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(r1_EvSE, "freq_means_t1_results_EvSE.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(r1_EvSP, "freq_means_t1_results_EvSP.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(r1_EvPA, "freq_means_t1_results_EvPA.txt", sep = "\t", quote = FALSE, row.names = FALSE)

E2SE <- list(spino = abs(f1_spino$E.af - f1_spino$SE.af), EvSE = abs(f1_EvSE$E.af - f1_EvSE$SE.af))
E2SP <- list(spino = abs(f1_spino$E.af - f1_spino$SP.af), EvSE = abs(f1_EvSE$E.af - f1_EvSE$SP.af))
E2PAt1 <- list(spino = abs(f1_spino$E.af - f1_spino$PA.af), EvPA = abs(f1_EvPA$E.af - f1_EvPA$PA.af))
mSE <- list(spino = r1_spino$bg.E2SE.mean.spino, EvSE = r1_EvSE$bg.E2SE.mean.EvSE.T1)
mSP <- list(spino = r1_spino$bg.E2SP.mean.spino, EvSE = r1_EvSE$bg.E2SP.mean.EvSE.T1)
mPAt1 <- list(spino = r1_spino$bg.E2PA.mean.spino, EvPA = r1_EvPA$bg.E2PA.mean.EvPA.T1)
save_panels("Fig3GI_dp_boxplots_August.pdf", list(boxp(E2SE, mSE, "AF differences: E vs SE"),
	boxp(E2SP, mSP, "AF differences: E vs SP"), boxp(E2PAt1, mPAt1, "AF differences: E vs PA")), width = 3, height = 6)
stats <- rbind(stats,
	welch("Fig. 3G (unexposed vs extinct)", "extinct outliers", E2SE$EvSE, mSE$EvSE),
	welch("Fig. 3G (unexposed vs extinct)", "candidate genes", E2SE$spino, mSE$spino),
	welch("Fig. 3H (unexposed vs rescued)", "extinct outliers", E2SP$EvSE, mSP$EvSE),
	welch("Fig. 3H (unexposed vs rescued)", "candidate genes", E2SP$spino, mSP$spino),
	welch("Fig. 3I (unexposed vs pre-adapted)", "pre-adapted outliers", E2PAt1$EvPA, mPAt1$EvPA),
	welch("Fig. 3I (unexposed vs pre-adapted)", "candidate genes", E2PAt1$spino, mPAt1$spino))
cat("\n===== Δp at focal vs matched SNPs (one-sided Welch t-tests) =====\n")
print(transform(stats, mean_focal = round(mean_focal, 4), mean_matched = round(mean_matched, 4), t = round(t, 2), df = round(df, 1), P = signif(P, 3)), row.names = FALSE)
write.csv(stats, "dp_matched_welch_tests.csv", row.names = FALSE)

## ---- Checks against the paper ----
g <- function(panel, set) stats[stats$panel == panel & stats$focal_SNPs == set, ]
obs <- list(
	Fig2F_P_candidate_genes = g("Fig. 2F (unexposed vs exposed)", "candidate genes")$P,
	Fig2G_P_candidate_genes = g("Fig. 2G (unexposed vs pre-adapted)", "candidate genes")$P,
	Fig3G_mean_dp_extinct_outliers = round(g("Fig. 3G (unexposed vs extinct)", "extinct outliers")$mean_focal, 2),
	Fig3G_mean_dp_matched = round(g("Fig. 3G (unexposed vs extinct)", "extinct outliers")$mean_matched, 2),
	"Fig3G_P_below_2.2e-16" = as.numeric(g("Fig. 3G (unexposed vs extinct)", "extinct outliers")$P < 2.2e-16),
	Fig3H_mean_dp_extinct_outliers = round(g("Fig. 3H (unexposed vs rescued)", "extinct outliers")$mean_focal, 2),
	Fig3H_mean_dp_matched = round(g("Fig. 3H (unexposed vs rescued)", "extinct outliers")$mean_matched, 2),
	Fig3H_P = g("Fig. 3H (unexposed vs rescued)", "extinct outliers")$P,
	SI214_P_candidates_extinct = g("Fig. 3G (unexposed vs extinct)", "candidate genes")$P,
	SI214_P_candidates_rescued = g("Fig. 3H (unexposed vs rescued)", "candidate genes")$P,
	SI214_P_candidates_preadapted = g("Fig. 3I (unexposed vs pre-adapted)", "candidate genes")$P,
	SI212_overlap = ov$overlap[1], SI212_log10P = ov$log10_P[1])
pass <- check_expected("03_genome_scans", obs)
msg("Finished")
sink()
if (!pass) quit(status = 1)
