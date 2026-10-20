########################################

###  Code for: Contingency and determinism in rapid adaptation shape evolutionary rescue and extinction in Drosophila
###  Authors: René S. Shahmohamadloo*, Joshua J. Faber-Hammond, Amir R. Gabidulin, Catherine I. Clare, Seth M. Rudman*
###  *Corresponding authors:
###  René S. Shahmohamadloo, School of Biological Sciences, Washington State University, 14204 NE Salmon Creek Ave, Vancouver, WA 98686, USA. Email: rene.shahmohamadloo@wsu.edu
###  Seth M. Rudman, School of Biological Sciences, Washington State University, 14204 NE Salmon Creek Ave, Vancouver, WA 98686, USA. Email: seth.rudman@wsu.edu

###  In this file: 02_build_glm_master.R — builds the genome-wide GLM master table.

###  Usage:
###   rds/glm_all.rds          per-SNP P, FDR and -log10(P) for all 58 contrasts
###   rds/glm_all_rolwin21.rds the same, averaged over 21-SNP windows
###   rds/freq_means.rds       mean allele frequency per treatment (x time point)
###   rds/freq_means_t1.rds    August means (E, S, SE, SP, PA), August variance filter
###   rds/freq_t1_allmean.rds  August mean over all samples (RRHO frequency filter)
###   rds/meta_T1.rds          August sample table with SE/SP/PA/E condition

### Steps:
###   1. merge every GLM table by CHROM/POS (SNPs must pass the variance filter in every model), drop chromosome 4, sort;
###   2. Benjamini-Hochberg FDR within each contrast over all merged SNPs;
###   3. -log10(P), with P = 0 set to half the smallest non-zero P (fixes Inf);
###   4. 21-SNP centred rolling means of P, FDR and -log10(P) within each chromosome arm (Inf replaced by max finite + 100 first, as in the original).
### Then checks the SNP total and the counts reported in the paper (expected_values).

########################################

f_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
RF_HOME <- if (length(f_arg)) normalizePath(dirname(sub("^--file=", "", f_arg))) else getwd()
source(file.path(RF_HOME, "00_config.R"))
suppressMessages(library(dplyr))
OUT <- module_dir("02_build_glm_master"); dir.create(PATHS$rds, FALSE, TRUE)
sink(file.path(OUT, "02_build_glm_master_report.txt"), split = TRUE)
msg("GLM tables from: %s", PATHS$glm_tables)
rd <- function(f) read.table(file.path(PATHS$glm_tables, f), header = TRUE)

## ---- 1. Load and merge the GLM tables ----
## PA vs S: all time points; by time point (additive model); time-point pairs; interaction
PAvS_all <- rd("rudflies_2023_PAvS_treat.all.table.GLMcontrast.txt")[, -3]; names(PAvS_all)[3] <- "PAvS"
PAvS_trt <- rd("rudflies_2023_PAvS_treat.GLMcontrast.txt")
PAvS_tpt <- rd("rudflies_2023_PAvS_tpt.GLMcontrast.txt")
PAvS_int <- rd("rudflies_2023_PAvS_int.GLMcontrast.txt")
## these three come from one model run and are combined by position:
## confirm their loci line up row for row before doing so
stopifnot(identical(PAvS_trt[, 1:2], PAvS_tpt[, 1:2]), identical(PAvS_trt[, 1:2], PAvS_int[, 1:2]))
PAvS <- cbind(PAvS_trt, PAvS_tpt[, -c(1, 2)], PAvS_int[, -c(1, 2)])
PAvS <- merge(PAvS_all, PAvS, by = c("CHROM", "POS"))
PAvS <- PAvS[order(PAvS[, 1], PAvS[, 2]), ]

## Within-treatment time-point contrasts (S, PA, E)
tS  <- rd("rudflies_2023_S_tpt.table.GLMcontrast.txt")
tPA <- rd("rudflies_2023_PA_tpt.table.GLMcontrast.txt")
tE  <- rd("rudflies_2023_E_tpt.table.GLMcontrast.txt")
names(tE) <- c("CHROM", "POS", "E.af.mean", paste0("E.", c("T1vT2", "T1vT3", "T1vT4", "T2vT3", "T2vT4", "T3vT4")))
tpt <- merge(merge(tS, tPA, by = c("CHROM", "POS")), tE, by = c("CHROM", "POS"))[, -c(3, 10, 17)]   # drop AF means

## August four-group model (E, PA, SE, SP)
T1 <- rd("rudflies_2023_PAvSvSEvE.wLoci.GLMcontrast.txt")

## PA vs E: all time points; founders; by time point
PAvE <- merge(rd("rudflies_2023_PAvE_treat.all.table.GLMcontrast.txt")[, -3],
	rd("rudflies_2023_PAvE_treat.table.GLMcontrast.founders.txt")[, -3], by = c("CHROM", "POS"))
PAvE <- merge(PAvE, rd("rudflies_2023_PAvE_treat.table.GLMcontrast2.txt")[, -3], by = c("CHROM", "POS"))
names(PAvE) <- c("CHROM", "POS", "PAvE", "PAvE.F", paste0("PAvE.T", 1:4))

## S vs E: all time points; by time point; interaction
SvE <- merge(rd("rudflies_2023_SvE_treat.all.table.GLMcontrast.txt")[, c(1, 2, 4)],
	rd("rudflies_2023_SvE_treat.GLMcontrast.table.txt"), by = c("CHROM", "POS"))
SvE <- merge(SvE, rd("rudflies_2023_SvE_int.GLMcontrast.table.txt"), by = c("CHROM", "POS"))
names(SvE) <- c("CHROM", "POS", "SvE", paste0("SvE.T", 1:4), paste0("int.SvE.", c("T1vT2", "T1vT3", "T1vT4", "T2vT3", "T2vT4", "T3vT4")))

msg("SNPs per table: PAvS %s | S/PA/E time points %s | August 4-group %s | PAvE %s | SvE %s",
	format(nrow(PAvS), big.mark = ","), format(nrow(tpt), big.mark = ","), format(nrow(T1), big.mark = ","),
	format(nrow(PAvE), big.mark = ","), format(nrow(SvE), big.mark = ","))

glm.all <- merge(PAvS, tpt, by = c("CHROM", "POS"))
glm.all <- merge(glm.all, T1, by = c("CHROM", "POS"))
glm.all <- merge(glm.all, PAvE, by = c("CHROM", "POS"))
glm.all <- merge(glm.all, SvE, by = c("CHROM", "POS"))
glm.all <- glm.all[glm.all$CHROM != "4", ]
glm.all <- glm.all %>% arrange(CHROM, POS)
x <- ncol(glm.all)
msg("Merged master table: %s SNPs x %d contrasts", format(nrow(glm.all), big.mark = ","), x - 2)
cat("Contrasts:", paste(names(glm.all)[-(1:2)], collapse = " "), "\n")

## Diagnostics: failed fits (NA) per contrast
na <- colSums(is.na(glm.all[, 3:x])); if (any(na > 0)) { cat("NA P-values per contrast:\n"); print(na[na > 0]) } else cat("No NA P-values.\n")

## ---- 2. FDR and 3. -log10(P) ----------------------------------------------------------
for (i in 3:x) glm.all[[paste0(names(glm.all)[i], ".fdr")]] <- p.adjust(glm.all[[i]], method = "fdr")
for (i in 3:x) glm.all[[paste0(names(glm.all)[i], ".logp")]] <- neglog10(glm.all[[i]])

## ---- 4. 21-SNP rolling windows, per chromosome arm ------------------------------------
roll_arm <- function(arm) {
	g <- unique(glm.all[glm.all$CHROM == arm, ]) %>% arrange(POS)
	g$POS <- as.integer(g$POS)
	r <- g[, 1:2]
	for (i in 3:ncol(g)) {
		v <- g[[i]]
		v[is.infinite(v)] <- max(v[is.finite(v)]) + 100   # Inf -> max finite + 100
		r[[paste0(names(g)[i], ".rolwin21")]] <- slide_mean(v, before = 10, after = 10)
	}
	r
}
glm.all.rolwin21 <- do.call(rbind, lapply(c("2L", "2R", "3L", "3R", "X"), roll_arm))
glm.all.rolwin21 <- na.omit(glm.all.rolwin21)
msg("Rolling-window table: %s SNPs", format(nrow(glm.all.rolwin21), big.mark = ","))

## ---- Shared allele-frequency summaries (used by modules 03-10) --------------------------
haf <- load_haf(exclude = TRUE)
## All time points: PA, S and E samples of the spinosad experiment
fs <- field_samples(haf); m <- fs$meta; sites <- fs$sites; fr <- fs$freq
## samples are paired with frequency columns by name; report whether the file order already agreed
cat("All-time-point samples: frequency columns in metadata order:", fs$positional_ok, "\n")
stopifnot(identical(as.character(m$samp), names(fr)))
mn <- function(sel) rowMeans(fr[, sel, drop = FALSE])
freq_means <- cbind(sites, E.af = mn(m$treat.fix == "E"), S.af = mn(m$treat.fix == "S"), PA.af = mn(m$treat.fix == "PA"))
for (tr in c("E", "S", "PA")) for (t in 1:4) freq_means[[paste0(tr, ".T", t, ".af")]] <- mn(m$treat.fix == tr & m$tpt == t)

## August: four groups, August variance filter
mt <- merge(haf$meta, COND_T1, by = c("cage", "tpt"))
mt <- unique(mt[mt$batch == "a" & mt$treat %in% c("S", "PA", "E"), ])
mt <- mt[mt$samp %in% KEEP_SAMP, ]
mt <- mt[order(mt$samp), ]
ft <- haf$freq[, names(haf$freq) %in% mt$samp]
cat("August samples: frequency columns in metadata order:", identical(as.character(mt$samp), names(ft)), "\n")
mt <- mt[match(names(ft), as.character(mt$samp)), ]
stopifnot(identical(as.character(mt$samp), names(ft)))
v1 <- rowVars(as.matrix(ft)) > 0.001
sites1 <- haf$freq[v1, 1:2]; ft <- ft[v1, ]
mn1 <- function(sel) rowMeans(ft[, sel, drop = FALSE])
freq_means_t1 <- cbind(sites1, E.af = mn1(mt$treat.fix == "E"), S.af = mn1(mt$treat.fix == "S"),
	SE.af = mn1(mt$condition == "SE"), SP.af = mn1(mt$condition == "SP"), PA.af = mn1(mt$treat.fix == "PA"))
## mean over all August samples (used by module 04 to drop SNPs with mean frequency outside 0.15-0.85)
freq_t1_allmean <- cbind(sites1, af_mean = rowMeans(ft))
cat("August samples per group:\n"); print(table(mt$condition))
rm(haf, fr, ft); invisible(gc())

## ---- Save ------------------------------------------------------------------------------
saveRDS(glm.all, file.path(PATHS$rds, "glm_all.rds"))
saveRDS(glm.all.rolwin21, file.path(PATHS$rds, "glm_all_rolwin21.rds"))
saveRDS(freq_means, file.path(PATHS$rds, "freq_means.rds"))
saveRDS(freq_means_t1, file.path(PATHS$rds, "freq_means_t1.rds"))
saveRDS(freq_t1_allmean, file.path(PATHS$rds, "freq_t1_allmean.rds"))
saveRDS(mt, file.path(PATHS$rds, "meta_T1.rds"))
msg("Saved to %s", PATHS$rds)

## ---- Checks against the numbers reported in the paper -------------------------------------
sm <- glm.all.rolwin21
on3 <- function(chr) round(100 * mean(chr %in% c("3L", "3R")), 1)
obs <- c(n_snps = nrow(glm.all))
for (ct in c("PAvE", "SvE")) {
	pu <- pr <- c()
	for (thr in c(0.05, 0.01)) {
		u <- glm.all[[paste0(ct, ".fdr")]] < thr
		s <- sm[[paste0(ct, ".fdr.rolwin21")]] < thr
		obs[sprintf("%s_unsmoothed_fdr%02d", ct, round(thr * 100))] <- sum(u, na.rm = TRUE)
		obs[sprintf("%s_smoothed_fdr%02d", ct, round(thr * 100))] <- sum(s, na.rm = TRUE)
		pu <- c(pu, on3(glm.all$CHROM[which(u)])); pr <- c(pr, on3(sm$CHROM[which(s)]))
	}
	obs[paste0("pct3_", ct, "_unsmoothed_min")] <- min(pu); obs[paste0("pct3_", ct, "_unsmoothed_max")] <- max(pu)
	obs[paste0("pct3_", ct, "_smoothed_min")] <- min(pr); obs[paste0("pct3_", ct, "_smoothed_max")] <- max(pr)
}
obs["EvSE.T1_smoothed_fdr05"] <- sum(sm$EvSE.T1.fdr.rolwin21 < 0.05)
obs["EvSP.T1_smoothed_fdr05"] <- sum(sm$EvSP.T1.fdr.rolwin21 < 0.05)
obs["EvSE.T1_unsmoothed_fdr05"] <- sum(glm.all$EvSE.T1.fdr < 0.05, na.rm = TRUE)
obs["EvSP.T1_unsmoothed_fdr05"] <- sum(glm.all$EvSP.T1.fdr < 0.05, na.rm = TRUE)
obs["EvSE.T1_EvSP.T1_unsmoothed_overlap"] <- sum(glm.all$EvSE.T1.fdr < 0.05 & glm.all$EvSP.T1.fdr < 0.05, na.rm = TRUE)
obs["SvE.T4_smoothed_fdr05"] <- sum(sm$SvE.T4.fdr.rolwin21 < 0.05)
pass <- check_expected("02_build_glm_master", as.list(obs))
msg("Finished")
sink()
if (!pass) quit(status = 1)
