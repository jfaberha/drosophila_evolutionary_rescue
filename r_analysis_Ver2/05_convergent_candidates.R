########################################

###  Code for: Contingency and determinism in rapid adaptation shape evolutionary rescue and extinction in Drosophila
###  Authors: René S. Shahmohamadloo*, Joshua J. Faber-Hammond, Amir R. Gabidulin, Catherine I. Clare, Seth M. Rudman*
###  *Corresponding authors:
###  René S. Shahmohamadloo, School of Biological Sciences, Washington State University, 14204 NE Salmon Creek Ave, Vancouver, WA 98686, USA. Email: rene.shahmohamadloo@wsu.edu
###  Seth M. Rudman, School of Biological Sciences, Washington State University, 14204 NE Salmon Creek Ave, Vancouver, WA 98686, USA. Email: seth.rudman@wsu.edu

###  In this file: 05_convergent_candidates.R — loci converging between exposed and pre-adapted populations.

###  Produces:
###  Candidate sets (SI 1.8.10), all on 21-SNP-window FDR:
###  November (rescued): PAvE Aug & Nov FDR < 0.01, exposed vs unexposed Nov < 0.05
###  August (extinct):   PAvE Aug & Nov FDR < 0.01, extinct vs unexposed Aug < 0.05
###  Relaxed three-way:  all four contrasts FDR < 0.1
###  each with unexposed August frequency 0.1-0.9 and |PA - E| in August > 0.1

###  Outputs:
###  SNP lists of the three sets (used by modules 07 and 09)
###  gene lists (used by 09 GO enrichment and the tables below)
###  Tables S7 (November), S10 (August), S11 (relaxed): features within 5 kb
###  Fig. S11  residual Manhattan plots
###  Fig. 4D   allele-frequency trajectories of the November set counts reported in the main text and SI 2.3.1-2.3.3

########################################

f_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
RF_HOME <- if (length(f_arg)) normalizePath(dirname(sub("^--file=", "", f_arg))) else getwd()
source(file.path(RF_HOME, "00_config.R"))
suppressMessages({ library(dplyr); library(ggplot2) })
OUT <- module_dir("05_convergent_candidates"); setwd(OUT)
sink("05_convergent_candidates_report.txt", split = TRUE)

glm.all <- load_rds("glm_all"); r <- load_rds("glm_all_rolwin21")
freq_means <- load_rds("freq_means"); fm1 <- load_rds("freq_means_t1")
vep <- load_vep()
key <- function(d) paste(d$CHROM, d$POS)
vkeys <- unique(key(vep))
fly_max <- vep %>% group_by(CHROM, POS) %>% summarise(FLYCADD = suppressWarnings(max(FLYCADD, na.rm = TRUE)), .groups = "drop") %>% as.data.frame()
on3 <- function(d) round(100 * mean(d$CHROM %in% c("3L", "3R")), 1)

## ---- Candidate sets ------------------------------------------------------------------
af_ok <- fm1[which(fm1$E.af < 0.9 & fm1$E.af > 0.1 & abs(fm1$E.af - fm1$PA.af) > 0.1), c("CHROM", "POS")]
ra <- r[key(r) %in% vkeys, ]                       # windowed SNPs with a VEP annotation (as glm.all.rolwin21.annot)
pa_both <- ra$PAvE.T1.fdr.rolwin21 < 0.01 & ra$PAvE.T4.fdr.rolwin21 < 0.01
sets <- list(
	TPT4 = merge(unique(ra[which(pa_both & ra$SvE.T4.fdr.rolwin21 < 0.05), c("CHROM", "POS")]), af_ok, by = c("CHROM", "POS")),
	TPT1 = merge(unique(ra[which(pa_both & ra$EvSE.T1.fdr.rolwin21 < 0.05), c("CHROM", "POS")]), af_ok, by = c("CHROM", "POS")),
	relaxed_fdr10 = merge(unique(ra[which(pmax(ra$PAvE.T1.fdr.rolwin21, ra$PAvE.T4.fdr.rolwin21, ra$SvE.T4.fdr.rolwin21, ra$EvSE.T1.fdr.rolwin21) < 0.1), c("CHROM", "POS")]), af_ok, by = c("CHROM", "POS")))
for (s in names(sets)) {
	x <- merge(sets[[s]], fly_max, by = c("CHROM", "POS"))
	write.table(x, paste0("convergent_candidate_SNPs_", s, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
	msg("%s: %d SNPs, %d with FlyCADD > 0.6, %.1f%% on 3L/3R", s, nrow(x), sum(x$FLYCADD > 0.6, na.rm = TRUE), on3(x))
}
fly_n <- function(s) sum(merge(sets[[s]], fly_max, by = c("CHROM", "POS"))$FLYCADD > 0.6, na.rm = TRUE)

## ---- Gene lists (features within 5 kb, via VEP) ----------------------------------------
genes_of <- function(snps, v = vep) { x <- merge(snps, v, by = c("CHROM", "POS")); unique(x$Gene) }
noflank <- vep[!vep$Consequence %in% c("upstream_gene_variant", "downstream_gene_variant"), ]
fly_rows <- vep[which(vep$FLYCADD > 0.6), ]
gl <- list()
for (s in names(sets)) gl[[s]] <- list(all = genes_of(sets[[s]]), contains = genes_of(sets[[s]], noflank), flycadd = genes_of(sets[[s]], fly_rows))
## gene-list files read by 09_GO_enrichment.R
wl <- function(v, f) write.table(v, f, sep = "\t", quote = FALSE, row.names = FALSE)
wl(gl$TPT1$all, "glm.all.rolwin21.TPT1.convergent_candidate_genes.txt")
wl(gl$TPT4$all, "glm.all.rolwin21.TPT4.convergent_candidate_genes.txt")
wl(gl$relaxed_fdr10$all, "glm.all.rolwin21.adaptive.fdr10.txt")

## ---- Tables S7, S10, S11 ---------------------------------------------------------------
## Annotation text (symbol, type, coordinates, names) comes from the published tables in
## data/published_tables (FlyBase r6.69); asterisk = feature contains a candidate SNP (not only
## up/downstream); FLYCADD = feature has a candidate SNP with FlyCADD > 0.6.
pub <- lapply(list.files(file.path(RF_HOME, "data", "published_tables"), pattern = "^TableS(7|10|11)\\.csv$", full.names = TRUE), read.csv, check.names = FALSE, colClasses = "character")
look <- do.call(rbind, pub); look$ID0 <- sub("\\*$", "", look$ID); look <- look[!duplicated(look$ID0), ]
fb_ids <- function(v) { v <- sub("_df_nrg$", "", trimws(v)); unique(v[grepl("^FB(gn|ti|tr)[0-9]", v)]) }
build_table <- function(g, tab) {
	ids <- sort(fb_ids(g$all)); star <- fb_ids(g$contains); yes <- fb_ids(g$flycadd)
	a <- look[match(ids, look$ID0), c("SYMBOL", "TYPE", "CHROM", "START", "STOP", "NAME", "GENE SYMBOL")]
	miss <- is.na(a$SYMBOL)
	if (any(miss)) { a[miss, ] <- "-"; cat(sprintf("  %s: %d features lack annotation text in data/published_tables: %s\n", tab, sum(miss), paste(ids[miss], collapse = " "))) }
	d <- data.frame(ID = paste0(ids, ifelse(ids %in% star, "*", "")), a, `FLYCADD > 0.60` = ifelse(ids %in% yes, "YES", "NO"), check.names = FALSE)
	write.csv(d, paste0("Table", tab, ".csv"), row.names = FALSE)
	ref <- file.path(RF_HOME, "data", "published_tables", paste0("Table", tab, ".csv"))
	same <- if (file.exists(ref)) { p <- read.csv(ref, check.names = FALSE, colClasses = "character"); identical(unname(as.matrix(p)), unname(as.matrix(d))) } else NA
	cat(sprintf("  Table %s: %d features (%d with *, %d FlyCADD YES); identical to published table: %s\n", tab, nrow(d), sum(endsWith(d$ID, "*")), sum(d[[9]] == "YES"), same))
	nrow(d)
}
cat("\n===== SI tables =====\n")
nS7 <- build_table(gl$TPT4, "S7"); nS10 <- build_table(gl$TPT1, "S10"); nS11 <- build_table(gl$relaxed_fdr10, "S11")

## ---- Counts reported in the text ----------------------------------------------------------
## 'near features' = SNP has a VEP gene assignment (genes within 5 kb); FlyCADD counted among all SNPs
ann <- function(keys) {
	v <- merge(unique(keys[, c("CHROM", "POS")]), vep, by = c("CHROM", "POS"))
	near <- unique(v[!is.na(v$Gene) & !v$Gene %in% c("-", ""), c("CHROM", "POS")])
	f <- merge(unique(keys[, c("CHROM", "POS")]), fly_max, by = c("CHROM", "POS"))
	c(n = nrow(unique(keys[, c("CHROM", "POS")])), near = nrow(near), fly = sum(f$FLYCADD > 0.6, na.rm = TRUE))
}
a1 <- ann(r[which(r$EvSE.T1.fdr.rolwin21 < 0.05), ])
a2 <- ann(glm.all[which(glm.all$EvSP.T1.fdr < 0.05), ])
a3 <- ann(r[which(r$SvE.T4.fdr.rolwin21 < 0.05), ])
## November rescued outliers that diverge in the same direction as pre-adapted (main text);
## three possible definitions are reported
nv <- merge(r[which(r$SvE.T4.fdr.rolwin21 < 0.05), c("CHROM", "POS", "PAvE.T4.fdr.rolwin21")],
	freq_means[, c("CHROM", "POS", "E.T4.af", "S.T4.af", "PA.T4.af")], by = c("CHROM", "POS"))
same <- sign(nv$S.T4.af - nv$E.T4.af) == sign(nv$PA.T4.af - nv$E.T4.af)
sd3 <- c(any_PA = sum(same), PA_fdr05 = sum(same & nv$PAvE.T4.fdr.rolwin21 < 0.05), PA_fdr01 = sum(same & nv$PAvE.T4.fdr.rolwin21 < 0.01))
cat("\nNovember rescued outliers in the same direction as pre-adapted (main text 2,456):\n"); print(sd3)
shared <- nrow(merge(sets$TPT1, sets$TPT4, by = c("CHROM", "POS")))
counts <- data.frame(item = c("Aug extinct outliers (21-SNP FDR<0.05)", "  near features", "  FlyCADD > 0.6",
	"Aug rescued outliers (unsmoothed FDR<0.05)", "  near features", "  FlyCADD > 0.6",
	"Nov rescued outliers (21-SNP FDR<0.05)", "  near features", "  FlyCADD > 0.6",
	"Nov same direction as PA (any / PA FDR<0.05 / PA FDR<0.01)",
	"Nov convergent set", "  FlyCADD > 0.6", "  on 3L/3R", "Aug convergent set", "  FlyCADD > 0.6", "  % on 3L/3R",
	"Relaxed three-way set", "  FlyCADD > 0.6", "SNPs in both Aug and Nov sets", "Table S7 / S10 / S11 features"),
	value = c(a1, a2, a3, paste(sd3, collapse = " / "), nrow(sets$TPT4), fly_n("TPT4"), sum(sets$TPT4$CHROM %in% c("3L", "3R")),
		nrow(sets$TPT1), fly_n("TPT1"), on3(sets$TPT1), nrow(sets$relaxed_fdr10), fly_n("relaxed_fdr10"), shared, paste(nS7, nS10, nS11, sep = " / ")))
cat("\n===== Counts =====\n"); print(counts, row.names = FALSE)
write.csv(counts, "convergent_counts.csv", row.names = FALSE)

## ---- Fig. S11: residual Manhattan plots ------------------------------------------------
## residual = A - B / max(B) * max(A) on window -log10(P) (A: vs unexposed; B: pre-adapted vs
## exposed), floored at 0. Points: November set (green; red if FlyCADD > 0.6), August set
## (cyan; yellow if FlyCADD > 0.6), relaxed three-way set (dark).
fly_keys <- unique(key(fly_rows))
relaxed_plot <- merge(na.omit(r[, c("CHROM", "POS", "PAvE.T1.fdr.rolwin21", "PAvE.T4.fdr.rolwin21", "SvE.T4.fdr.rolwin21", "EvSE.T1.fdr.rolwin21")]), af_ok, by = c("CHROM", "POS"))
relaxed_plot <- relaxed_plot[pmax(relaxed_plot[[3]], relaxed_plot[[4]], relaxed_plot[[5]], relaxed_plot[[6]]) < 0.1, c("CHROM", "POS")]
hl <- list(list(sets$TPT4, "#97AA71"), list(sets$TPT4[key(sets$TPT4) %in% fly_keys, ], "#C95B5B"),
	list(sets$TPT1, "#78B0BF"), list(sets$TPT1[key(sets$TPT1) %in% fly_keys, ], "#E8B53E"), list(relaxed_plot, "#4C4C4C"))
xs <- function(g) g + facet_grid(~ CHROM, scales = "free_x", space = "free_x") +
	scale_x_continuous(breaks = c(0, 5e6, 1e7, 1.5e7, 2e7, 2.5e7, 3e7), guide = guide_axis(angle = 45)) +
	labs(col = "candidate\ngene\n-log10(p)") + xlab("chromosome position") + ylab("-log10(p)") + theme_classic() + theme(legend.position = "none")
resid_plot <- function(a, b, title) {
	A <- r[[paste0(a, ".logp.rolwin21")]]; B <- r[[paste0(b, ".logp.rolwin21")]]
	d <- na.omit(data.frame(CHROM = r$CHROM, POS = r$POS, res = A - (B / max(na.omit(B)) * max(na.omit(A))),
		r[[paste0(a, ".fdr.rolwin21")]], A, r[[paste0(b, ".fdr.rolwin21")]], B))
	d$res[d$res < 0] <- 0
	g <- ggplot(d, aes(POS, res)) + geom_line(alpha = 0.5, colour = "grey")
	for (h in hl) g <- g + geom_point(data = merge(h[[1]], d, by = c("CHROM", "POS")), aes(POS, res), color = h[[2]], size = 1)
	xs(g + scale_y_continuous(limits = c(0, round(max(d$res)))) + ggtitle(title))
}
ref_plot <- function(a, title) {
	y <- paste0(a, ".logp.rolwin21")
	xs(ggplot(na.omit(r), aes(POS, .data[[y]])) + geom_line(alpha = 1, colour = "#CCCCCC") +
		geom_point(data = na.omit(r[which(r[[paste0(a, ".fdr.rolwin21")]] < 0.01), ]), aes(POS, .data[[y]]), color = "darkmagenta", size = 1) +
		scale_y_continuous(limits = c(0, max(na.omit(r[[y]])))) + ggtitle(title))
}
save_panels("FigS11_residual_manhattans.pdf", list(
	ref_plot("PAvE", "E vs PA"), ref_plot("SvE", "E vs S"),
	resid_plot("PAvE.T1", "PAvSE.T1", "E vs PA minus PA vs SE (TPT1): 21-SNP sliding window"),
	resid_plot("EvSE.T1", "PAvSE.T1", "E vs SE minus PA vs SE (TPT1): 21-SNP sliding window"),
	resid_plot("PAvE.T4", "PAvS.T4", "PA vs E minus PA vs S (TPT4): 21-SNP sliding window"),
	resid_plot("SvE.T4", "PAvS.T4", "S vs E minus PA vs S (TPT4): 21-SNP sliding window"),
	resid_plot("PAvE", "PAvS", "PA vs E minus PA vs S: 21-SNP sliding window"),
	resid_plot("SvE", "PAvS", "S vs E minus PA vs S: 21-SNP sliding window")), ncol = 2, nrow = 4, width = 11.25, height = 8)
msg("Fig. S11 written")

## ---- Fig. 4D: allele-frequency trajectories of the November set ------------------------
rm(glm.all); invisible(gc())
haf <- load_haf(exclude = TRUE)
mf <- haf$meta[(haf$meta$batch == "a" | haf$meta$batch == "F") & haf$meta$experiment == "spino" & haf$meta$treat.fix %in% c("PA", "S", "E"), ]
ff <- haf$freq[, names(haf$freq) %in% mf$samp]; mf <- mf[mf$samp %in% names(ff), ]
keep <- rowVars(as.matrix(ff)) > 0.001
sites <- haf$freq[keep, 1:2]; ff <- ff[keep, ]; ok <- complete.cases(ff); sites <- sites[ok, ]; ff <- ff[ok, ]
ff <- ff[, match(as.character(mf$samp), names(ff))]                    # pair by name
idx <- match(key(sets$TPT4), key(sites)); af <- as.matrix(ff[idx[!is.na(idx)], ])
msg("Fig. 4D: %d of %d November candidates in the founder + field frequency table", nrow(af), nrow(sets$TPT4))
if (nrow(af) < 2) stop("Fewer than two November candidates found in the frequency table; cannot draw Fig. 4D.")
grp <- as.character(mf$treat.fix)
grp[mf$tpt == 1 & mf$cage %in% as.numeric(SE_CAGES)] <- "SE"; grp[grp == "S"] <- "SP"
gk <- paste(grp, mf$tpt)
means <- sapply(split(seq_along(gk), gk), function(i) rowMeans(af[, i, drop = FALSE]))
tpts <- sort(unique(as.numeric(mf$tpt))); f0 <- min(tpts); E0 <- paste("E", f0); PA0 <- paste("PA", f0)
means <- cbind(means, `SP F` = means[, E0])                              # exposed founders = unexposed founders
flip <- means[, PA0] < means[, E0]; means[flip, ] <- 1 - means[flip, ]   # oriented to the direction of pre-adapted allelic change
rel <- means - means[, E0]
conv <- rel[, paste("SP", max(tpts))] - rel[, paste("E", max(tpts))] >= 0
msg("Fig. 4D: %d SNPs reoriented; %d of %d shifted toward the pre-adapted allele (rescued vs unexposed, November)", sum(flip), sum(conv), length(conv))
long <- do.call(rbind, lapply(c("E", "PA", "SP"), function(g) {
	cl <- c(if (g == "SP") "SP F" else paste(g, f0), paste(g, tpts[tpts > f0]))
	data.frame(treat.fix = g, tpt = factor(rep(c(f0, tpts[tpts > f0]), each = nrow(rel))), value = as.vector(rel[, cl]), snp = rep(seq_len(nrow(rel)), length(cl)))
}))
p4d <- ggplot(long, aes(x = tpt, y = value, group = paste(treat.fix, snp), color = treat.fix)) + geom_path(alpha = 0.8, linewidth = 0.2) +
	theme_classic() + ggtitle("convergent candidate loci") + scale_color_manual(values = unname(SAMPLE_COLS[c("E", "PA", "SP")]))
ggsave("Fig4D_convergent_candidate_trajectories.pdf", p4d, width = 5, height = 4)

## ---- Checks ----------------------------------------------------------------------------
pass <- check_expected("05_convergent_candidates", list(
	Aug_extinct_outliers = a1[["n"]], Aug_extinct_near = a1[["near"]], Aug_extinct_flycadd = a1[["fly"]],
	Aug_rescued_outliers_unsmoothed = a2[["n"]], Aug_rescued_near = a2[["near"]], Aug_rescued_flycadd = a2[["fly"]],
	Nov_rescued_outliers = a3[["n"]], Nov_rescued_near = a3[["near"]], Nov_rescued_flycadd = a3[["fly"]],
	Nov_same_direction_2456_found = as.numeric(2456 %in% sd3),
	Nov_convergent = nrow(sets$TPT4), Nov_convergent_flycadd = fly_n("TPT4"), Nov_convergent_on_3L3R = sum(sets$TPT4$CHROM %in% c("3L", "3R")),
	Aug_convergent = nrow(sets$TPT1), Aug_convergent_flycadd = fly_n("TPT1"), Aug_convergent_pct_3L3R = on3(sets$TPT1),
	Relaxed_convergent = nrow(sets$relaxed_fdr10), Relaxed_convergent_flycadd = fly_n("relaxed_fdr10"), Shared_Aug_Nov = shared,
	Fig4D_snps_toward_preadapted = sum(conv), TableS7_rows = nS7, TableS10_rows = nS10, TableS11_rows = nS11))
msg("Finished")
sink()
if (!pass) quit(status = 1)
